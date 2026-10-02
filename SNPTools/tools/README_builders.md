# Maize SNPTrait / SNPGeo data builders

All builders are Python 3.8+ standard library unless noted. Outputs mirror `SNPTools/`.

| Script | Output | Inputs |
|---|---|---|
| `build_maize_samples.py` | `data/zmgrin2026_samples.tsv` | T0 population members, `grz2023_grin_crosswalk.tsv`, optional `--release-manifest` (vcf_name -> sample_id as written to the release VCFs) and `--extra-samples` (separately called lines) |
| `build_maize_trait_catalog.py` | `js/zmgrin.catalog.js` (`SNP_CATALOG.families.zmgrin2026`, `SNPTRAIT_SCHEMA.zmgrin2026`), `data/traits/zmgrin2026.traits.json` | samples TSV, T5 passport/geo, T4b trait summary, T4c trait dictionary, crosswalk |
| `build_maize_snpgeo_data.py` | `js/snpgeo.regions.js` (and `data/geo/countries.geo.json`) | Natural Earth admin-0, samples TSV, T5, crosswalk |
| `build_geo_layers.py` | `data/geo/countries.geo.json`, `data/geo/admin1_na.geo.json` | Natural Earth v5.1.2 admin-0 110m, admin-1 10m (see `data/geo/PROVENANCE.md`) |
| `vcf_to_h5.py` (h5py, numpy) | `hdf5/version3/zmgrin2026_<chr>_impute.h5` | one release VCF per chromosome |

The files committed in this branch were built from the release v1.3 manifest
(932 samples = 926 imputed + 6 companion lines). Run order and pinned sources are in
`../Makefile` (`make help`). The prototype validators from the design phase
(`validate_maize_port.py`, `smoke_test_snptools.js`) ran the unmodified fusarium modules
with proposed patches; they are superseded by `localdev/harness/` now that the changes are
applied. `validate_maize_port.py` is kept for its structural checks of the catalogue shape.

Natural Earth: https://github.com/nvkelso/natural-earth-vector tag v5.1.2
(f1890d9f152c896d250a77557a5751a93d494776), public domain.
