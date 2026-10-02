SHELL := /bin/bash
# SNPTools maize (zmgrin2026) data build + local instance.   make help
# Default to the conda env built for this port when it exists; override with PHP_BIN=/PYTHON_PATH=.
LOCAL_ENV ?= $(HOME)/.claude-science/conda/envs/snptools-local
PY        ?= $(if $(wildcard $(LOCAL_ENV)/bin/python),$(LOCAL_ENV)/bin/python,python3)
PHP_BIN   ?= $(if $(wildcard $(LOCAL_ENV)/bin/php),$(LOCAL_ENV)/bin/php,php)
PYTHON_PATH ?= $(PY)
NE_TAG    ?= f1890d9f152c896d250a77557a5751a93d494776   # natural-earth-vector v5.1.2
NE_RAW    := https://raw.githubusercontent.com/nvkelso/natural-earth-vector/$(strip $(NE_TAG))/geojson
NE_DIR    ?= build/ne
# project tables (set on the command line)
T0 ?=
CROSSWALK ?=
T4B ?=
T4C ?=
T5 ?=
RELEASE_MANIFEST ?=
EXTRA_SAMPLES ?=
# release VCFs, one per chromosome: $(RELEASE_VCF_DIR)/grz2023_release_v1.3_all_chr<N>.vcf.gz
RELEASE_VCF_DIR ?=
CHRS := 1 2 3 4 5 6 7 8 9 10

.PHONY: help samples catalog regions geo-layers store local-store start stop test
help:
	@echo "samples     data/zmgrin2026_samples.tsv       (T0= CROSSWALK= [RELEASE_MANIFEST=] [EXTRA_SAMPLES=])"
	@echo "catalog     js/zmgrin.catalog.js + data/traits (T5= T4B= T4C= CROSSWALK=)"
	@echo "regions     js/snpgeo.regions.js               (T5= CROSSWALK=)"
	@echo "geo-layers  data/geo/*.geo.json from Natural Earth v5.1.2"
	@echo "store       hdf5/version3/zmgrin2026_chr*_impute.h5 from RELEASE_VCF_DIR (full release; large)"
	@echo "local-store small test store from localdev/fixtures"
	@echo "start/stop  local PHP server on 127.0.0.1:8877;  test  headless end-to-end checks"

samples:
	$(PY) tools/build_maize_samples.py --t0 $(T0) --crosswalk $(CROSSWALK) \
	  $(if $(RELEASE_MANIFEST),--release-manifest $(RELEASE_MANIFEST)) $(if $(EXTRA_SAMPLES),--extra-samples $(EXTRA_SAMPLES)) \
	  --out data/zmgrin2026_samples.tsv
catalog:
	$(PY) tools/build_maize_trait_catalog.py --samples data/zmgrin2026_samples.tsv --t5 $(T5) --t4b $(T4B) --t4c $(T4C) \
	  --crosswalk $(CROSSWALK) --outdir .
$(NE_DIR)/%.geojson:
	mkdir -p $(NE_DIR) && curl -fsSL -o $@ $(NE_RAW)/$*.geojson
regions: $(NE_DIR)/ne_110m_admin_0_countries.geojson
	$(PY) tools/build_maize_snpgeo_data.py --ne $< --samples data/zmgrin2026_samples.tsv --t5 $(T5) \
	  --crosswalk $(CROSSWALK) --outdir . --catalog js/zmgrin.catalog.js
geo-layers: $(NE_DIR)/ne_110m_admin_0_countries.geojson $(NE_DIR)/ne_10m_admin_1_states_provinces.geojson
	$(PY) tools/build_geo_layers.py --ne0 $(word 1,$^) --out data/geo/countries.geo.json \
	  --ne1 $(word 2,$^) --out1 data/geo/admin1_na.geo.json
store:
	@test -n "$(RELEASE_VCF_DIR)" || { echo "set RELEASE_VCF_DIR"; exit 1; }
	for c in $(CHRS); do $(PYTHON_PATH) tools/vcf_to_h5.py $(RELEASE_VCF_DIR)/grz2023_release_v1.3_all_chr$$c.vcf.gz \
	  hdf5/version3/zmgrin2026_chr$${c}_impute.h5 --samples <(cut -f1 data/zmgrin2026_samples.tsv | tail -n +2) || exit 1; done
local-store:
	PYTHON_PATH=$(PYTHON_PATH) localdev/build_test_store.sh
start:
	PHP_BIN=$(PHP_BIN) PYTHON_PATH=$(PYTHON_PATH) localdev/start.sh
stop:
	localdev/stop.sh
test:
	cd localdev/harness && npm install --no-audit --no-fund && PHP_BIN=$(PHP_BIN) PYTHON_PATH=$(PYTHON_PATH) npm test
