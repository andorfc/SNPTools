#!/usr/bin/env python3
"""
h5_to_vcf.py — region + accession list -> VCF, for the Fusarium SNPTools store.

Differences from the maize build:
  * Fusarium isolates are HAPLOID, so genotype codes map to single alleles
    {0:"0", 1:"1", 2:"2", 3:"."}  (not diploid "0/0" etc.).
  * Three reference genomes are supported; the ##source / ##reference header
    block is chosen from the HDF5 path (fusarium2025 = F. graminearum,
    vert7600 = F. verticillioides 7600, vertMRC826 = F. verticillioides MRC826).
  * A requested accession that is not a column in this .h5 is filled with "."
    instead of aborting the whole query (robustness borrowed from the maize build).

Usage:  h5_to_vcf.py <db.h5> <out.vcf> <start> <end> <genotypesJson>
"""
import h5py
import sys
import json
import datetime
import gzip
import os
import shutil
import tempfile
import atexit
import numpy as np

hdf5_file_path  = sys.argv[1]
output_vcf_path = sys.argv[2]
lower_bound     = int(sys.argv[3])
upper_bound     = int(sys.argv[4])
json_string     = sys.argv[5]

# The Fusarium test stores ship gzip-compressed (.h5.gz). h5py cannot read a gzip
# stream directly (that is whole-file gzip, not HDF5's internal filters), so if we
# were handed a .gz, transparently decompress the single requested chromosome file
# to a temp .h5 and read that. Only the one chromosome needed for this query is
# expanded, and the temp file is removed on exit.
if hdf5_file_path.endswith('.gz'):
    _tmp_fd, _tmp_path = tempfile.mkstemp(suffix='.h5')
    os.close(_tmp_fd)
    atexit.register(lambda p=_tmp_path: os.path.exists(p) and os.remove(p))
    with gzip.open(hdf5_file_path, 'rb') as _gz, open(_tmp_path, 'wb') as _out:
        shutil.copyfileobj(_gz, _out)
    hdf5_file_path = _tmp_path

var_list = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO']
genome_list = json.loads(json_string)

# Haploid: one allele index per call. 3 (or anything unmapped) -> missing.
reverse_genotype_mapping = {0: "0", 1: "1", 2: "2", 3: "."}
def gt_string(val):
    try:
        return reverse_genotype_mapping.get(int(val), ".")
    except (ValueError, TypeError):
        return "."

current_date = datetime.date.today().strftime('%Y%m%d')
ID = "."
FILTER = "."
FORMAT = "GT"

with h5py.File(hdf5_file_path, 'r') as hdf5_file:
    if 'POS' not in hdf5_file:
        print("No 'POS' dataset found in the file.")
        sys.exit(1)

    pos_data = hdf5_file['POS'][:]
    try:
        pos_data = np.array([int(pos.decode('utf-8')) if isinstance(pos, bytes) else int(pos)
                             for pos in pos_data])
    except ValueError as e:
        print(f"Error converting 'POS' data to integers: {e}")
        sys.exit(1)

    lower_index = np.searchsorted(pos_data, lower_bound, side='left')
    upper_index = np.searchsorted(pos_data, upper_bound, side='right')

    variant_data = {var: hdf5_file[var][lower_index:upper_index] for var in var_list}

    # Read only columns that exist; stray/foreign ids are filled with "." below,
    # so one bad accession id never aborts the whole request.
    present = [g for g in genome_list if g in hdf5_file]
    missing = [g for g in genome_list if g not in hdf5_file]
    if missing:
        print(f"Note: {len(missing)} requested accession(s) not in {hdf5_file_path}; "
              f"filling with .  e.g. {missing[:5]}")
    genome_data = {g: hdf5_file[g][lower_index:upper_index] for g in present}

if len(variant_data['POS']) == 0:
    print("No data found in the specified position range.")
else:
    path = sys.argv[1]   # ORIGINAL path (hdf5_file_path may be a temp file after gz decompression)
    # Genome-key prefixes match the hdf5/fusarium/ filenames (2026 rebuild):
    #   Fvert_7600_* , Fvert_mrc_* , fusarium2026_* (graminearum).
    # Older aliases kept so pre-2026 files still resolve.
    if ("Fvert_7600" in path) or ("vert7600" in path):
        source = "Fusarium_verticillioides_7600_2026"
        reference = "Fusarium verticillioides 7600 (FungiDB / FVEG gene models)"
    elif ("Fvert_mrc" in path) or ("vertMRC826" in path):
        source = "Fusarium_verticillioides_MRC826_2026"
        reference = "Fusarium verticillioides MRC826 (FVERT4 gene models)"
    else:
        source = "Fusarium_graminearum_2026"
        reference = "Fusarium graminearum PH-1 (FGSG gene models)"

    with open(output_vcf_path, 'w') as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write("##fileDate=" + current_date + "\n")
        vcf_file.write("##source=" + source + "\n")
        vcf_file.write("##reference=" + reference + "\n")
        vcf_file.write("##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS mapping quality.\">\n")
        vcf_file.write("##INFO=<ID=CVC,Number=1,Type=Integer,Description=\"The number of accessions that have genotype data for a particular variant.\">\n")
        vcf_file.write("##INFO=<ID=CVP,Number=1,Type=Float,Description=\"The percent of accessions that have genotype data for a particular variant.\">\n")
        vcf_file.write("##INFO=<ID=TYPE,Number=.,Type=String,Description=\"The type of effect using Sequence Ontology terms.\">\n")
        vcf_file.write("##INFO=<ID=EFFECT,Number=.,Type=String,Description=\"An estimation of putative impact/deleteriousness.\">\n")
        vcf_file.write("##INFO=<ID=GENEMODEL,Number=.,Type=String,Description=\"The name of the gene model affected by the variant.\">\n")
        vcf_file.write("##INFO=<ID=SUB,Number=.,Type=String,Description=\"The amino acid substitution for missense and non-synonymous variants.\">\n")
        vcf_file.write("##INFO=<ID=MAXR2,Number=1,Type=Float,Description=\"The maximum R2 for a given loci.\">\n")
        # 2026 language-model scores (graminearum). INFO is passed through verbatim
        # from the HDF5, so these headers are declarative; the values come from the store.
        vcf_file.write("##INFO=<ID=FUNDLM_SCORE,Number=1,Type=Float,Description=\"FunDLM DNA language-model variant-effect score; more negative = more conserved / potentially deleterious.\">\n")
        vcf_file.write("##INFO=<ID=EVO2_SCORE,Number=1,Type=Float,Description=\"Evo2 DNA language-model variant-effect score.\">\n")
        vcf_file.write("##INFO=<ID=ESM1_SCORE,Number=1,Type=Float,Description=\"ESM1 protein language-model variant-effect score; lower = higher likelihood of functional impact.\">\n")
        vcf_file.write("##INFO=<ID=ESM2_SCORE,Number=1,Type=Float,Description=\"ESM2 protein language-model variant-effect score.\">\n")
        vcf_file.write("##INFO=<ID=ESM3_SCORE,Number=1,Type=Float,Description=\"ESM3 protein language-model variant-effect score.\">\n")
        vcf_file.write("##INFO=<ID=ESMC_SCORE,Number=1,Type=Float,Description=\"ESM-C protein language-model variant-effect score.\">\n")
        # 2025 legacy scores (verticillioides, until its 2026 rebuild).
        vcf_file.write("##INFO=<ID=DNA_SCORE,Number=1,Type=Float,Description=\"Legacy DNA language-model (DNABERT-2) mutation-effect score.\">\n")
        vcf_file.write("##INFO=<ID=AA_SCORE,Number=1,Type=Float,Description=\"Legacy ESM1b protein language-model variant-effect score.\">\n")
        vcf_file.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor Allele Frequency.\">\n")
        vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        vcf_file.write("#" + "\t".join(["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"] + genome_list) + "\n")

        def dec(v):
            return v.decode('utf-8') if isinstance(v, bytes) else str(v)

        for i in range(len(variant_data['POS'])):
            row = [
                dec(variant_data['CHROM'][i]),
                dec(variant_data['POS'][i]),
                ID,
                dec(variant_data['REF'][i]),
                dec(variant_data['ALT'][i]),
                dec(variant_data['QUAL'][i]),
                FILTER,
                dec(variant_data['INFO'][i]),
                FORMAT,
            ]
            genotypes = [gt_string(genome_data[g][i]) if g in genome_data else "." for g in genome_list]
            vcf_file.write("\t".join(row + genotypes) + "\n")

    print(f"VCF data has been saved to {output_vcf_path}")
