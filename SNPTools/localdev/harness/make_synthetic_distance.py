#!/usr/bin/env python3
"""make_synthetic_distance.py <outdir> <site_root> -- small synthetic genome-wide IBS files
for the SNPCompare / SNPTree tests (the real matrices are computed on Ceres).

<outdir>/zmgrin2026/   all 933 release sample_ids (data/zmgrin2026_samples.tsv), in a fixed
                       shuffled order, with similarity.csv, missing_pct.csv (fraction 0..1),
                       the *_snp.csv pair, ids.txt, tree_nj.nwk, tree_upgma.nwk, and the files
                       the endpoint does not serve (allele_distance.csv, co_called_sites.csv).
<outdir>/              the unchanged MaizeGDB 2026 layout (maizegdb_allchr_final_*.csv + ids.txt)
                       for the first 60 mgdb2026 catalogue ids, to check that path still works.
Values are deterministic functions of the matrix indices (see sim()/mis()), so a checker can
recompute any row. Standard library only."""
import csv, json, os, random, sys


def sim(i, j, snp=False):
    if i == j:
        return 1.0
    v = 0.80 + 0.19 * (((i * j) + i + j) % 1000) / 1000
    return round(v - (0.0015 if snp else 0), 4)


def mis(i, j, snp=False):
    return 0.0 if i == j else round(((i + j) % 50) / 1000 + (0.002 if snp else 0), 4)


def write_matrix(path, n, f):
    with open(path, 'w') as fh:
        for i in range(n):
            fh.write(','.join(f'{f(i, j):.4f}' for j in range(n)) + '\n')


def main(out, root):
    z = os.path.join(out, 'zmgrin2026'); os.makedirs(z, exist_ok=True)
    ids = [r['sample_id'] for r in csv.DictReader(open(f'{root}/data/zmgrin2026_samples.tsv'), delimiter='\t')]
    random.Random(20260930).shuffle(ids)
    n = len(ids)
    open(f'{z}/ids.txt', 'w').write('\n'.join(ids) + '\n')
    write_matrix(f'{z}/similarity.csv', n, sim)
    write_matrix(f'{z}/missing_pct.csv', n, mis)
    write_matrix(f'{z}/similarity_snp.csv', n, lambda i, j: sim(i, j, True))
    write_matrix(f'{z}/missing_pct_snp.csv', n, lambda i, j: mis(i, j, True))
    write_matrix(f'{z}/allele_distance.csv', n, lambda i, j: 1 - sim(i, j))
    write_matrix(f'{z}/co_called_sites.csv', n, lambda i, j: 1000)
    open(f'{z}/tree_nj.nwk', 'w').write('(' + ','.join(f'{x}:0.01' for x in ids[:5]) + ');\n')
    open(f'{z}/tree_upgma.nwk', 'w').write('((' + ids[0] + ':0.01,' + ids[1] + ':0.01):0.02,' + ids[2] + ':0.03);\n')
    src = open(f'{root}/js/accessions.catalog.js', encoding='utf-8').read()
    cat = json.loads(src[src.index('{'):src.rindex('}') + 1])
    mids = [a['id'] for p in cat['families']['mgdb2026']['projects'] for g in p['groups'] for a in g['accessions']][:60]
    open(f'{out}/ids.txt', 'w').write('\n'.join(mids) + '\n')
    write_matrix(f'{out}/maizegdb_allchr_final_similarity.csv', len(mids), sim)
    write_matrix(f'{out}/maizegdb_allchr_final_missing_pct.csv', len(mids), mis)
    print(json.dumps({'zmgrin2026': n, 'mgdb2026': len(mids)}))


if __name__ == '__main__':
    main(*sys.argv[1:3])
