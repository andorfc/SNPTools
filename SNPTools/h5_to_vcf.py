"""
h5_to_vcf.py  --  region + accession list  ->  VCF, pulled from the HDF5 store.

    python h5_to_vcf.py <db.h5> <out.vcf> <start> <end> <accessions>

<accessions> is either a path to a JSON file holding an array of accession
ids (how processForm.php passes it, to dodge Windows shell-quoting of '"'),
or a raw JSON array string (any direct caller).

The heavy lifting is vectorized: genotype codes (stored as int8) are mapped
to VCF strings through a NumPy lookup table over the whole variant x
accession block at once, and rows are assembled and written in chunks. The
per-cell Python loop it replaced was O(variants x accessions) with a
function call + dict lookup each, which is what made wide / full-panel
queries take minutes. POS is decoded with a single astype() rather than a
Python loop over every position in the file.

The INFO column is stored verbatim in the HDF5 and passed straight through,
so every annotation key -- GENEMODEL, TYPE, EFFECT, SUB, MQ, CVP, MAXR2,
MAF, plantcad1/2_score, ESM1/2/3_score -- survives unchanged.

Memory is bounded by a block of rows, not by the request: the fixed columns
and the genotypes are read one block at a time (blocks aligned to the
store's 65,536-row chunks). Reading the whole slice at once held every
column for the whole interval -- a 6.6 GB peak footprint for chromosome 2
with only 5 lines, the INFO column alone padded to its longest entry over
5.18 M rows; now ~0.3 GB for any interval. (On macOS the resident size reads
higher: freed memory stays mapped as reusable until the OS reclaims it.)

A request larger than SNPTOOLS_MAX_CELLS (variants x accessions, default
2e9, about 8 GB of genotype text) is refused before anything is written:
the script prints "TOO_LARGE: ..." and exits 3, and processForm.php passes
the message on.
"""
import gzip
import h5py
import os
import sys
import json
import datetime
import numpy as np

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
hdf5_file_path = sys.argv[1]
output_vcf_path = sys.argv[2]
lower_bound = int(sys.argv[3])
upper_bound = int(sys.argv[4])

arg5 = sys.argv[5]
if os.path.isfile(arg5):
    with open(arg5, 'r') as _f:
        json_string = _f.read()
else:
    json_string = arg5
genome_list = json.loads(json_string)
if not isinstance(genome_list, list):
    genome_list = []
# ids are column names: anything that is not a string becomes a missing column's name
genome_list = [g if isinstance(g, str) else str(g) for g in genome_list]

# HDF5 genotype code -> VCF genotype string. Anything outside 0..3 is clamped
# to 3 ('./.'), matching the old dict-with-default behaviour.
GT_LUT = np.array([b'0/0', b'1/0', b'1/1', b'./.'], dtype='S3')
FIXED_COLS = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO']
ROW_CHUNK = 20000                 # rows per write (string assembly), fewer for wide requests
GT_BLOCK_CELLS = 50_000_000       # genotype cells read at once (~50 MB of int8)
MAX_BLOCK_CHUNKS = 4              # and at most 4 store chunks (262,144 rows) per read
try:
    MAX_CELLS = int(float(os.environ.get('SNPTOOLS_MAX_CELLS') or 2e9))
except ValueError:
    MAX_CELLS = int(2e9)
current_date = datetime.date.today().strftime('%Y%m%d')


def as_bytes_col(arr):
    """h5 column (object-of-bytes, 'S', or 'U') -> contiguous fixed-width 'S' array."""
    if arr.dtype.kind == 'U':
        return np.char.encode(arr)
    return arr.astype('S')


# ---------------------------------------------------------------------------
# Locate the requested slice in the HDF5 store
# ---------------------------------------------------------------------------
# No chunk cache: blocks are cut on chunk boundaries, so each chunk is read once, and the
# default 1 MB cache per dataset doubled the footprint for a full panel's 933 genotype
# columns (chr2:100-130 Mb: 711 MB with it, 356 MB without).
hdf5_file = h5py.File(hdf5_file_path, 'r', rdcc_nbytes=0)
if 'POS' not in hdf5_file:
    print("No 'POS' dataset found in the file.")
    sys.exit(1)

pos_raw = hdf5_file['POS'][:]
try:
    pos_data = pos_raw.astype(np.int64)
except (ValueError, TypeError):
    pos_data = np.array(
        [int(p.decode('utf-8')) if isinstance(p, bytes) else int(p) for p in pos_raw],
        dtype=np.int64,
    )
del pos_raw

lower_index = int(np.searchsorted(pos_data, lower_bound, side='left'))
upper_index = int(np.searchsorted(pos_data, upper_bound, side='right'))
n_rows = upper_index - lower_index
del pos_data

if n_rows == 0:
    print("No data found in the specified position range.")
    sys.exit(0)

n_acc = len(genome_list)
if n_rows * max(1, n_acc) > MAX_CELLS:
    print(f"TOO_LARGE: {n_rows:,} variants x {n_acc:,} accessions is more than this server builds "
          f"in one request ({MAX_CELLS:,} genotype cells). Choose a smaller interval or fewer accessions.")
    sys.exit(3)


def genotype_dataset(g):
    """The store's genotype column for accession id g, or None. The fixed columns (POS,
    INFO, ...) and anything that is not a top-level dataset are not accessions."""
    if g in FIXED_COLS or '/' in g or g not in hdf5_file:
        return None
    d = hdf5_file[g]
    return d if isinstance(d, h5py.Dataset) else None


# A stray/foreign accession id in the request is filled with './.' rather than aborting.
gt_dsets = [genotype_dataset(g) for g in genome_list]
missing = [g for g, d in zip(genome_list, gt_dsets) if d is None]
if missing:
    print(f"Note: {len(missing)} requested accession(s) not in {hdf5_file_path}; "
          f"filling with ./.  e.g. {missing[:5]}")
fixed_dsets = {c: hdf5_file[c] for c in FIXED_COLS}

# Rows per read block: whole store chunks (65,536 rows in the vcf_to_h5.py stores), as many
# as GT_BLOCK_CELLS genotype cells allow, from 1 to MAX_BLOCK_CHUNKS of them.
STORE_CHUNK = (hdf5_file['POS'].chunks or (65536,))[0]
BLOCK = STORE_CHUNK * min(MAX_BLOCK_CHUNKS, max(1, GT_BLOCK_CELLS // (STORE_CHUNK * max(1, n_acc))))


def blocks():
    """(a, b) store-index ranges covering [lower_index, upper_index), cut at multiples of
    BLOCK so each store chunk is decompressed once."""
    a = lower_index
    while a < upper_index:
        b = min(upper_index, (a // BLOCK + 1) * BLOCK)
        yield a, b
        a = b


def read_block(a, b):
    """Fixed columns as stored and the genotype code matrix (rows, accessions) in request
    order. The fixed columns become fixed-width bytes per write chunk, not per block: the
    conversion pads every value to the longest INFO in the slice it is given."""
    cols = {c: d[a:b] for c, d in fixed_dsets.items()}
    gt = np.empty((b - a, n_acc), dtype=np.int8)
    for k, d in enumerate(gt_dsets):
        gt[:, k] = d[a:b] if d is not None else 3
    np.clip(gt, 0, 3, out=gt)
    return cols, gt


# A write chunk's strings take about rows x (4 bytes per accession + the INFO); keep the
# genotype text of one chunk near 40 MB, so a full panel writes ~10,000 rows at a time.
ROW_CHUNK = max(1000, min(ROW_CHUNK, 40_000_000 // (4 * max(1, n_acc) + 512)))
ID_COL = np.full(ROW_CHUNK, b'.', dtype='S1')
FILTER_COL = np.full(ROW_CHUNK, b'.', dtype='S1')
FORMAT_COL = np.full(ROW_CHUNK, b'GT', dtype='S2')


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
def common_info_header():
    return ["##fileformat=VCFv4.2", "##fileDate=" + current_date]


def info_defs():
    return [
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">',
        '##INFO=<ID=CVC,Number=1,Type=Integer,Description="The number of accessions that have genotype data for a particular variant">',
        '##INFO=<ID=CVP,Number=1,Type=Float,Description="The percent of accessions that have genotype data for a particular variant.">',
        '##INFO=<ID=TYPE,Number=.,Type=String,Description="The type of effect using Sequence Ontology terms">',
        '##INFO=<ID=EFFECT,Number=.,Type=String,Description="An estimation of putative impact/deleteriousness">',
        '##INFO=<ID=GENEMODEL,Number=.,Type=String,Description="The name of the gene model affected by the variant">',
        '##INFO=<ID=SUB,Number=.,Type=String,Description="The amino acid substitution for missense and non-synonymous variants">',
        '##INFO=<ID=MAXR2,Number=1,Type=Float,Description="The maximum R2 for a given loci">',
        '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + genome_list),
    ]


def header_text():
    lines = []
    if "maizegdb" in hdf5_file_path:
        lines += common_info_header()
        lines.append("##source=MaizeGDB")
        lines.append("##reference=Andorf CM, Ross-Ibarra J, Seetharam AS, Hufford MB, Woodhouse MR. (2024) A unified VCF data set from nearly 1,500 diverse maize accessions and resources to explore the genomic landscape of maize. G3 Genes|Genomes|Genetics.")
        lines.append("##doi=https://doi.org/10.1101/2024.04.30.591904")
        lines += info_defs()
    if "nam" in hdf5_file_path:
        lines += common_info_header()
        lines.append("##source=MaizeGDB+NAM")
        lines.append("##reference=Hufford MB, Seetharam AS, Woodhouse MR, et al. De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes. Science. 2021;373(6555):655-662.")
        lines.append("##doi=https://doi.org/10.1126/science.abg5289")
        lines += info_defs()
    if "schnable" in hdf5_file_path:
        lines += common_info_header()
        lines.append("##source=Schnable2023")
        lines.append("##reference=Grzybowski MW, Mural RV, Xu G, Turkus J, Yang J, Schnable JC. A common resequencing-based genetic marker data set for global maize diversity. Plant J. 2023;113(6):1109-1121.")
        lines.append("##doi=https://doi.org/10.1111/tpj.16123")
        lines += info_defs()
    if "zmgrin" in hdf5_file_path:
        lines += common_info_header()
        lines.append("##source=MaizeGDB GRIN-linked 2026 (release v1.4; Grzybowski et al. 2023 sites, Beagle-imputed + direct-call companion lines)")
        lines.append("##reference=Grzybowski MW, Mural RV, Xu G, Turkus J, Yang J, Schnable JC. A common resequencing-based genetic marker data set for global maize diversity. Plant J. 2023;113(6):1109-1121.")
        lines.append("##doi=https://doi.org/10.1111/tpj.16123")
        lines += info_defs()
    if not lines:
        # Unknown family: still emit a minimal header. Without the #CHROM line the
        # browser's parseVcf() cannot map sample columns and reads every call as missing.
        lines += common_info_header()
        lines += info_defs()
    return ("\n".join(lines) + "\n").encode()


# ---------------------------------------------------------------------------
# Write the VCF, one row-chunk at a time. A '.gz' output path is gzip-
# compressed transparently. compresslevel=1: the genotype text is ~90% '0/0',
# so LZ77 alone gets ~10x and level 1 is ~30x faster than the default level 9
# for a barely-larger file (a 250 MB VCF: L1 = 25 MB in ~1s, L9 = 11 MB in ~38s).
# h5_to_vcf.py called directly with a plain '.vcf' path still writes plain text.
# ---------------------------------------------------------------------------
def _open(path, mode):
    return gzip.open(path, mode, compresslevel=1) if path.endswith('.gz') else open(path, mode)
with _open(output_vcf_path, 'wb') as vcf_file:
    vcf_file.write(header_text())

    for blo, bhi in blocks():
        fixed_cols, gt_codes = read_block(blo, bhi)
        for a in range(0, bhi - blo, ROW_CHUNK):
            b = min(a + ROW_CHUNK, bhi - blo)
            m = b - a

            # 9 fixed columns, tab-joined: only 9 vectorized adds.
            fc = {c: as_bytes_col(v[a:b]) for c, v in fixed_cols.items()}
            line = fc['CHROM']
            for col in (fc['POS'], ID_COL[:m], fc['REF'], fc['ALT'], fc['QUAL'], FILTER_COL[:m],
                        fc['INFO'], FORMAT_COL[:m]):
                line = np.char.add(np.char.add(line, b'\t'), col)

            # Genotype block: one vectorized LUT map over the whole (m, n_acc)
            # slice, prefix each cell with a tab, then reinterpret each row's
            # contiguous 4-byte cells ('\t' + 'x/y') as a single string.
            if n_acc:
                cells = np.ascontiguousarray(np.char.add(b'\t', GT_LUT[gt_codes[a:b]]))
                gt_lines = cells.view(f'S{4 * n_acc}').reshape(m)
                line = np.char.add(line, gt_lines)

            vcf_file.write(b'\n'.join(line.tolist()))
            vcf_file.write(b'\n')
        del fixed_cols, gt_codes

hdf5_file.close()
print(f"variants: {n_rows}")
print(f"VCF data has been saved to {output_vcf_path}")
