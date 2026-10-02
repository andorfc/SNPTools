# data/geo — provenance and licence

| File | Source | Processing |
|---|---|---|
| `countries.geo.json` | Natural Earth 1:110m Cultural Vectors, Admin 0 – Countries (`geojson/ne_110m_admin_0_countries.geojson`, 838,726 bytes, sha256 `6866c877d39cba9c357620878839b336d569f8c662d3cfab4cb1dbe2d39c977f`) | `tools/build_geo_layers.py`: 2-dp coordinates, Antarctica dropped, properties `{iso3, name}`; 176 features |
| `admin1_na.geo.json` | Natural Earth 1:10m Cultural Vectors, Admin 1 – States, Provinces (`geojson/ne_10m_admin_1_states_provinces.geojson`, 40,726,851 bytes, sha256 `22d0e3ad85eb3e27f17cabf8ba2d50e554fbc27a87796ff891d958185da62fb5`) | `tools/build_geo_layers.py --ne1`: USA, CAN, MEX only; Douglas-Peucker 0.04 deg per ring; rings with bounding box < 0.15 deg dropped; U.S. rings east of 0 deg longitude (Aleutians across the antimeridian) dropped; 2-dp coordinates; properties `{iso3, code (ISO 3166-2), name}`; 97 features (51 + 13 + 33), 274,757 bytes |

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
        --out data/geo/countries.geo.json \
        --ne1 ne_10m_admin_1_states_provinces.geojson --out1 data/geo/admin1_na.geo.json

The 110 m admin-1 file covers only the USA and the 50 m file lacks Mexico, which is
why the 10 m file is simplified here. All 39 ISO 3166-2 codes carried by
`SNPGEO_REGIONS[...].admin1Code` for USA/CAN/MEX samples in the v1.3 set have a
polygon in `admin1_na.geo.json`.
