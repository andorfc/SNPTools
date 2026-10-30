# The SNPMaize MCP server

`mcp.php` makes the SNPMaize instance a [Model Context Protocol](https://modelcontextprotocol.io) server, so an
LLM client (Claude, Claude Code, Cursor, any MCP client) can ask it about variants, lines, passport and
trait data, curated alleles, identity-by-state and GWAS, and get answers that agree with the pages and
carry permalinks back to them. It is built to work alongside the MaizeGDB MCP server
(https://claude.maizegdb.org/mcp), which holds the gene records, names, sequences, expression, GO terms,
QTL and literature: resolve a gene name there, query its variants here.

## Connecting

Endpoint: `https://<host>/mcp.php` (Streamable HTTP, stateless, no authentication; JSON-RPC 2.0 over
POST; `GET` answers a short plain-text page; `OPTIONS` answers CORS).

    claude mcp add --transport http snpmaize https://snpmaize.maizegdb.org/mcp.php
    claude mcp add --transport http maizegdb https://claude.maizegdb.org/mcp

    curl -s -X POST -H 'Content-Type: application/json' https://snpmaize.maizegdb.org/mcp.php \
         -d '{"jsonrpc":"2.0","id":1,"method":"tools/list"}'

Protocol versions 2025-06-18, 2025-03-26 and 2024-11-05 are accepted; `initialize`, `ping`,
`tools/list`, `tools/call`, `prompts/list` and `resources/list` (both empty) are served; notifications get
202 and no body; batches of up to 20 are answered as batches; bodies are capped at 256 KB.

## Tools

All tools are read-only. Every answer is JSON in `structuredContent` (and the same text in `content[0].text`)
with `ok`, `tool`, the payload, `links` (permalinks that open the same view on the site: `#snpcurate/SC0001`,
`#snpfunction/<gene>`, `#snpversity?ds=&chr=&start=&end=` ...) and `meta` (release, seconds). An expected
failure is an `isError` result with `error.code` one of `invalid_params`, `not_found`, `too_large`,
`unavailable`, `timeout`, `busy`, and a message the model can act on.

| tool | what it answers | main arguments |
|---|---|---|
| `snpmaize_datasets` | the release, the two variant sets, site-QC codes and rule, line-QC classes, the scores and the composite/tier calibration, registry counts, trait descriptor keys, GWAS studies, caps, conventions | none |
| `snpmaize_lines` | lines matching passport facets, a name query and/or a trait descriptor value or range; facet counts; the descriptor's distribution | `query`, `subpop`, `country`, `state`, `panel`, `kernel_type`, `improvement`, `era`, `species`, `trait`, `trait_value`, `trait_min`, `trait_max`, `limit` |
| `snpmaize_line` | one line: passport, every trait decoded, line QC, NAM founder, curated alleles carried | `line` |
| `snpmaize_variants` | the sites of a gene or region with consequence, scores, composite, percentile, tier, MAF, max R², carrier counts, site QC; filters; optional restriction to sites carried by given lines | `gene` or `region`, `flank`, `dataset`, `consequence`, `tier`, `qc`, `min_maf`, `lines`, `gene_only`, `sort`, `limit` |
| `snpmaize_gene` | one gene: coordinates, counts by QC/class/tier, TOP/HIGH variants, registry alleles, GWAS within 10 kb, protein sites, links | `gene`, `dataset` |
| `snpmaize_carriers` | one site: annotation, genotype counts, carriers with passport, by subpopulation and country, B73, carrier line QC, registry entry, trait by genotype | `site` or `curated`, `ref`, `alt`, `dataset`, `trait`, `limit` |
| `snpmaize_curated` | the SNPCurate registry: one record, or the list by gene, text, mark, role | `id`, `gene`, `query`, `mark`, `role`, `limit` |
| `snpmaize_alleles` | genotypes of lines at sites, with favorable-allele carry calls for curated sites | `lines` or `founders`, `sites`, `dataset` |
| `snpmaize_similarity` | identity-by-state: genome-wide nearest lines from the precomputed matrix, or computed in a gene/region (focal vs panel, or all pairs) | `line`, `lines`, `gene`, `region`, `mode`, `sites`, `top`, `dataset` |
| `snpmaize_gwas` | the GWAS studies; hits in a gene/region, by study or trait; a study's top hits | `gene`, `region`, `flank`, `study`, `trait`, `significant_only`, `limit` |

Conventions: coordinates are 1-based on B73 v5 (chr1..chr10); genes are `Zm00001eb` ids (a transcript id or a
registry symbol such as `tga1` is accepted); lines resolve by id (`ZmG_B73`), name (`B73`), GRIN accession
(`PI 550473`) or genotype id (`ZmG:B73`); genotypes are `ref`, `het`, `alt`, `missing`; `dataset` is
`zmgrin2026_imp` (default) or `zmgrin2026_hq`.

Caps (per call): regions up to 250,000 bp (genes up to 1 Mb for the summary), flank up to 50,000 bp, 500
rows, 100 lines, 50 sites, 200 carrier names, 100 similar lines. Counts in an answer cover the whole window
even when rows are cut (`truncated: true`).

## How it works

`mcp.php` validates the JSON-RPC request and the tool arguments (types, enums, caps), then runs
`tools/snpq.py` once per call (`PYTHON_PATH`, as processForm.php runs h5_to_vcf.py) with the request on
stdin and relays its JSON. The library reads what the pages read, the variant stores and their site-QC
sidecars, `js/zmgrin.catalog.js`, `data/traits`, `js/zmgrin.lineqc.js`, `js/snpcurate.data.js`,
`data/qc/zmgrin2026.scorecal.json`, `gff/genes_index.txt`, `data/gwas` and `distance/zmgrin2026`, through
the page mirrors in `tools/build_snpcurate.py` (site-QC rule, consequence classes, composite z and tiers,
annotation parsing), so an answer agrees with the same question on the page; `localdev/harness/check_mcp.py`
checks that against `js/snpcurate.data.js` on every `make test`.

Environment (set where processForm.php's `PYTHON_PATH` is set):

| variable | default | meaning |
|---|---|---|
| `PYTHON_PATH` | `python3` | the interpreter with h5py and numpy |
| `SNPMAIZE_SITE` | the request's scheme, host and folder | public URL used in permalinks (set it behind a proxy) |
| `SNPMAIZE_MCP_MAX_CONCURRENT` | 4 | tool calls running at once; more answer `busy` with `retry_after_s` |
| `SNPMAIZE_MCP_TIMEOUT` | 60 | seconds per call; the library stops itself and answers `timeout` |
| `SNPTOOLS_DEBUG` | unset | `1` adds the library's stderr to internal-error results |

Every call is one line in the PHP error log (`snpmaize-mcp {"event":"tools/call","tool":...,"ms":...}`),
arguments shortened, no server paths in replies.

## Example

> Which NAM founders carry the teosinte allele of tga1, and what else in the gene looks damaging?

1. `snpmaize_gene {"gene":"tga1"}`: coordinates, 88 sites, the registry entry SC0001 (gold), the TOP/HIGH
   variants (S376F, E68K).
2. `snpmaize_alleles {"founders":true,"sites":["SC0001"]}`: all 26 founders are `ref` (the maize
   asparagine), the alternate is confined to the 7 parviglumis lines.
3. `snpmaize_carriers {"site":"chr4:46652277:C:T"}`: who carries S376F, by subpopulation and country.
4. The `links` open SNPCurate, SNPFunction and SNPGeo on the same gene; the MaizeGDB MCP's
   `maizegdb_record` gives the gene's function and literature.
