# data/geo — provenance and licence

| File | Source | Processing |
|---|---|---|
| `countries.geo.json` | Natural Earth 1:110m Cultural Vectors, Admin 0 – Countries (`geojson/ne_110m_admin_0_countries.geojson`, 838,726 bytes, sha256 `6866c877d39cba9c357620878839b336d569f8c662d3cfab4cb1dbe2d39c977f`) | `tools/build_geo_layers.py`: 2-dp coordinates, Antarctica dropped, properties `{iso3, name}`; 176 features |

Source repository: https://github.com/nvkelso/natural-earth-vector, tag `v5.1.2`
(commit `f1890d9f152c896d250a77557a5751a93d494776`), fetched 2026-09-28 from
raw.githubusercontent.com.

Licence: Natural Earth data are in the **public domain** ("All versions of Natural
Earth raster + vector map data found on this website are in the public domain",
https://www.naturalearthdata.com/about/terms-of-use/). No permission or attribution
is required; attribution "Made with Natural Earth" is appreciated.

Boundaries follow Natural Earth's de facto policy and do not imply any position by
USDA-ARS or MaizeGDB on disputed territories.

Regenerate:

    python3 tools/build_geo_layers.py --ne0 ne_110m_admin_0_countries.geojson \
        --out data/geo/countries.geo.json
