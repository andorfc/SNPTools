#!/usr/bin/env python3
"""check_mcp.py -- the SNPMaize MCP server (mcp.php + tools/snpq.py) over HTTP.

    python3 check_mcp.py <root>

Starts PHP's built-in server on the tree (PHP_BIN, PYTHON_PATH as the other checks), then checks
the transport (GET page, OPTIONS, 405, bad JSON, notifications, batches, initialize, tools/list,
unknown method/tool, argument validation, the wall-time limit) and the tools against the data the
pages are built from: carriers, composite, tier and site QC of registry entries must equal
js/snpcurate.data.js; permalinks must point at this server; caps must refuse; rows must be sorted
as documented. Tool cases that need a variant store skip when that chromosome's store is not
installed. Exit 1 on a failure."""
import json, os, re, socket, subprocess, sys, time, urllib.error, urllib.request

root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', '..'))
PHP = os.environ.get('PHP_BIN') or 'php'
PY = os.environ.get('PYTHON_PATH') or sys.executable
bad = 0
skipped = 0


def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1


def skip(msg):
    global skipped
    skipped += 1
    print('skip ' + msg)


def free_port():
    s = socket.socket(); s.bind(('127.0.0.1', 0)); p = s.getsockname()[1]; s.close(); return p


class Server:
    def __init__(self, env=None):
        self.port = free_port()
        self.url = f'http://127.0.0.1:{self.port}/mcp.php'
        e = dict(os.environ, PYTHON_PATH=PY, **(env or {}))
        self.log = open(os.path.join(os.path.dirname(__file__), '..', 'out', f'mcp_server_{self.port}.log'), 'w')
        self.p = subprocess.Popen([PHP, '-S', f'127.0.0.1:{self.port}', '-t', root], cwd=root, env=e, stdout=self.log, stderr=subprocess.STDOUT)
        for _ in range(100):
            try:
                urllib.request.urlopen(self.url, timeout=1).read(); return
            except urllib.error.HTTPError:
                return
            except Exception:
                time.sleep(0.1)
        raise SystemExit('FAIL the PHP server did not start')

    def http(self, method='POST', body=None, headers=None):
        data = body.encode() if isinstance(body, str) else body
        req = urllib.request.Request(self.url, data=data, method=method, headers={'Content-Type': 'application/json', **(headers or {})})
        try:
            r = urllib.request.urlopen(req, timeout=180)
            return r.status, dict(r.headers), r.read().decode()
        except urllib.error.HTTPError as e:
            return e.code, dict(e.headers), e.read().decode()

    def rpc(self, method, params=None, id=1):
        st, h, body = self.http('POST', json.dumps({'jsonrpc': '2.0', 'id': id, 'method': method, 'params': params or {}}))
        return st, (json.loads(body) if body else None)

    def call(self, tool, args=None):
        st, r = self.rpc('tools/call', {'name': tool, 'arguments': args or {}})
        res = (r or {}).get('result') or {}
        return st, res, res.get('structuredContent') or {}

    def stop(self):
        self.p.terminate()
        try: self.p.wait(5)
        except subprocess.TimeoutExpired: self.p.kill()
        self.log.close()


def js_object(path, assign_re):
    s = open(os.path.join(root, path), encoding='utf-8').read()
    m = re.search(assign_re, s)
    return json.JSONDecoder().raw_decode(s, m.end())[0]


def has_store(chrom, suffix='_impute'):
    return os.path.isfile(os.path.join(root, 'hdf5', 'grin2026', f'zmgrin2026_{chrom}{suffix}.h5'))


os.makedirs(os.path.join(os.path.dirname(__file__), '..', 'out'), exist_ok=True)
for f in ('mcp.php', 'tools/snpq.py'):
    check(os.path.isfile(os.path.join(root, f)), f'{f} present')
S = Server()
try:
    # ---------------- transport ----------------
    st, h, body = S.http('GET')
    check(st == 200 and h.get('Content-Type', '').startswith('text/plain') and 'SNPMaize MCP server' in body and 'claude.maizegdb.org/mcp' in body,
          f'GET: {st} {h.get("Content-Type")}, a plain-text page naming the companion server')
    st, h, body = S.http('GET', headers={'Accept': 'text/event-stream'})
    check(st == 405 and json.loads(body)['error']['code'] == -32600, f'GET for an event stream: {st} -32600 (no SSE stream)')
    st, h, body = S.http('GET', headers={'Accept': 'text/html,application/xhtml+xml,*/*;q=0.8'})
    check(st == 200 and 'SNPMaize MCP server' in body, f'GET from a browser: {st}, the page')
    st, h, body = S.http('OPTIONS')
    check(st == 204 and h.get('Access-Control-Allow-Origin') == '*' and 'POST' in h.get('Access-Control-Allow-Methods', ''), f'OPTIONS: {st}, CORS headers')
    st, h, body = S.http('PUT', '{}')
    check(st == 405 and json.loads(body)['error']['code'] == -32600, f'PUT: {st} -32600')
    st, h, body = S.http('POST', '{not json')
    check(st == 400 and json.loads(body)['error']['code'] == -32700, f'bad JSON: {st} -32700')
    st, h, body = S.http('POST', json.dumps({'jsonrpc': '2.0', 'method': 'notifications/initialized'}))
    check(st == 202 and body == '', f'a notification alone: {st}, empty body')
    st, h, body = S.http('POST', json.dumps({'jsonrpc': '2.0', 'id': 1, 'method': 'nope'}))
    check(st == 200 and json.loads(body)['error']['code'] == -32601 and h.get('MCP-Protocol-Version') == '2025-06-18', f'unknown method: -32601, MCP-Protocol-Version header')
    st, r = S.rpc('initialize', {'protocolVersion': '2025-03-26', 'clientInfo': {'name': 'check_mcp', 'version': '1'}})
    res = r['result']
    check(res['protocolVersion'] == '2025-03-26' and res['serverInfo']['name'] == 'snpmaize' and 'claude.maizegdb.org/mcp' in res['instructions']
          and res['capabilities']['tools'] == {'listChanged': False},
          f"initialize: protocol {res['protocolVersion']}, server {res['serverInfo']['name']} {res['serverInfo']['version']}, instructions name the MaizeGDB MCP")
    st, r = S.rpc('initialize', {'protocolVersion': '1999-01-01'})
    check(r['result']['protocolVersion'] == '2025-06-18', 'initialize with an unknown protocol version answers the latest')
    st, r = S.rpc('ping')
    check(r['result'] == {}, 'ping -> {}')
    st, r = S.rpc('tools/list')
    tools = r['result']['tools']
    names = [t['name'] for t in tools]
    expect = ['snpmaize_datasets', 'snpmaize_lines', 'snpmaize_line', 'snpmaize_variants', 'snpmaize_gene', 'snpmaize_carriers',
              'snpmaize_curated', 'snpmaize_alleles', 'snpmaize_similarity', 'snpmaize_gwas']
    check(names == expect, f'tools/list: {len(tools)} tools, the ten in order')
    check(all(t['inputSchema'].get('type') == 'object' and t['description'] and t['title'] and t['annotations'].get('readOnlyHint') is True
              and t['annotations'].get('destructiveHint') is False for t in tools), 'every tool: object schema, title, description, read-only annotations')
    listed = subprocess.run([PY, '-I', os.path.join(root, 'tools', 'snpq.py'), '--list'], capture_output=True, text=True).stdout.split()
    check(listed == expect, 'tools/snpq.py --list names the same ten tools')
    st, h, body = S.http('POST', json.dumps([{'jsonrpc': '2.0', 'id': 1, 'method': 'ping'}, {'jsonrpc': '2.0', 'id': 2, 'method': 'tools/list'},
                                             {'jsonrpc': '2.0', 'method': 'notifications/x'}]))
    b = json.loads(body)
    check(st == 200 and isinstance(b, list) and [x['id'] for x in b] == [1, 2], f'a batch of two requests and a notification answers two: {st}')
    st, r = S.rpc('tools/call', {'name': 'snpmaize_nothing', 'arguments': {}})
    check(r.get('error', {}).get('code') == -32602, 'tools/call of an unknown tool: -32602')

    # ---------------- validation (PHP side, no python run) ----------------
    st, res, sc = S.call('snpmaize_variants', {'gene': 'tga1', 'bogus': 1})
    check(res.get('isError') is True and sc['error']['code'] == 'invalid_params' and 'bogus' in sc['error']['message'], f"unknown argument: isError, {sc['error']['message'][:60]}")
    st, res, sc = S.call('snpmaize_variants', {})
    check(res.get('isError') is True and 'gene, region' in sc['error']['message'], 'variants without gene or region: isError naming both')
    st, res, sc = S.call('snpmaize_variants', {'region': 'chr1:1-1000', 'limit': 5000})
    check(res.get('isError') is True and sc['error']['code'] == 'too_large', 'limit over the cap: too_large')
    st, res, sc = S.call('snpmaize_variants', {'region': 'chr1:1-1000', 'tier': 'HUGE'})
    check(res.get('isError') is True and sc['error']['code'] == 'invalid_params', 'tier outside the enum: invalid_params')
    st, res, sc = S.call('snpmaize_line', {})
    check(res.get('isError') is True and 'line' in sc['error']['message'], 'line without a name: isError')
    st, res, sc = S.call('snpmaize_variants', {'region': 'chr1:1-300000'})
    check(res.get('isError') is True and sc['error']['code'] == 'too_large' and '250,000' in sc['error']['message'], 'a 300 kb region: too_large (the library)')
    check(json.loads(res['content'][0]['text']) == sc, 'content[0].text is the structuredContent as JSON')

    # ---------------- static-data tools ----------------
    st, res, sc = S.call('snpmaize_datasets')
    check(res.get('isError') is False and sc['ok'] is True and sc['tool'] == 'snpmaize_datasets' and sc['release']['lines'] == 933
          and sc['site_qc']['codes'] == ['NO_CARRIER', 'HET_ONLY', 'HET_EXCESS', 'HET_ELEVATED', 'PASS'] and sc['scores']['composite']['tiers']['top'] == -1.5353
          and sc['registry']['entries'] == 76 and sc['traits']['descriptors'] == 78 and sc['gwas']['studies'] == 38 and sc['caps']['max_span_bp'] == 250000,
          f"datasets: {sc['release']['lines']} lines, {sc['registry']['entries']} registry entries, {sc['traits']['descriptors']} descriptors, {sc['gwas']['studies']} studies, meta {sc['meta']['seconds']} s")
    check(sc['links']['site'] == f'http://127.0.0.1:{S.port}/' and sc['links']['snpcurate'] == f'http://127.0.0.1:{S.port}/#snpcurate',
          f"permalinks point at this server: {sc['links']['snpcurate']}")
    st, res, sc = S.call('snpmaize_line', {'line': 'PI 550473'})
    check(sc.get('ok') and sc['line']['id'] == 'ZmG_B73' and sc['passport']['grin'] == 'PI 550473' and sc['is_nam_founder'] and sc['n_traits'] >= 50
          and sc['line']['line_qc']['class'] == 'I' and any(t['descriptor'] == 'KERNEL-TYPE' and t['label'] == 'Dent' for t in sc['traits']),
          f"line by GRIN accession: {sc.get('line', {}).get('label')}, {sc.get('n_traits')} traits decoded, kernel type Dent")
    st, res, sc = S.call('snpmaize_line', {'line': 'no such line'})
    check(res.get('isError') and sc['error']['code'] == 'not_found', 'an unknown line: not_found')
    st, res, sc = S.call('snpmaize_lines', {'subpop': 'Tropical', 'trait': 'KERNEL-WEIGHT-1000', 'trait_min': 300, 'limit': 5})
    check(sc.get('ok') and sc['matched'] >= 5 and len(sc['lines']) == 5 and all(l['subpop'] == 'Tropical' and l['trait']['mean'] >= 300 for l in sc['lines'])
          and sc['descriptor']['scale'] == 'numeric' and sc['descriptor']['lines_measured'] > 800,
          f"lines: Tropical with 1000-kernel weight >= 300 g: {sc.get('matched')} matched, descriptor measured on {sc.get('descriptor', {}).get('lines_measured')} lines")
    st, res, sc = S.call('snpmaize_lines', {'trait': 'Kernel Type', 'trait_value': 'Flint', 'panel': 'Ames282', 'limit': 200})
    check(sc.get('ok') and sc['matched'] > 0 and all(l['trait']['label'] == 'Flint' and 'Ames282' in l['panels'] for l in sc['lines']),
          f"lines: coded descriptor by label + panel: {sc.get('matched')} Flint lines in MAP")
    st, res, sc = S.call('snpmaize_lines', {'query': 'Mo17'})
    check(sc.get('ok') and sc['lines'] and sc['lines'][0]['id'] == 'ZmG_MO17', f"lines: query Mo17 puts Mo17 first ({sc.get('matched')} matched)")
    st, res, sc = S.call('snpmaize_curated', {'id': 'SC0001'})
    check(sc.get('ok') and sc['gene'] == 'Zm00001eb175150' and sc['mark'] == 'gold' and sc['links']['snpcurate'].endswith('/#snpcurate/SC0001')
          and sc['links']['maizegdb_gene'].endswith('/Zm00001eb175150') and sc['site']['nHom'] == 7, 'curated SC0001: the record with its permalink and MaizeGDB link')
    st, res, sc = S.call('snpmaize_curated', {'mark': 'gold'})
    cu = js_object('js/snpcurate.data.js', r'window\.SNP_CURATE\s*=\s*')
    gold = [e for e in cu['entries'] if e.get('mark') == 'gold']
    check(sc.get('ok') and sc['matched'] == len(gold) and all(e['mark'] == 'gold' for e in sc['entries']), f"curated mark=gold: {sc.get('matched')} = {len(gold)} in snpcurate.data.js")
    st, res, sc = S.call('snpmaize_curated', {'gene': 'tga1'})
    check(sc.get('ok') and sc['matched'] >= 1 and all(e['gene'] == 'Zm00001eb175150' for e in sc['entries']), f"curated by symbol: {sc.get('matched')} tga1 entries")
    st, res, sc = S.call('snpmaize_gwas')
    check(sc.get('ok') and len(sc['studies']) == 38 and sc['studies'][0]['link'].endswith('#snpgwas/' + sc['studies'][0]['id']), 'gwas with no arguments: the 38 studies with permalinks')
    st, res, sc = S.call('snpmaize_gwas', {'study': 'tibbscortes2024_nam_dta_intcp', 'limit': 3})
    ps = [h['p'] for h in sc.get('hits', [])]
    check(sc.get('ok') and len(ps) == 3 and ps == sorted(ps) and all(h['significant'] for h in sc['hits']) and sc['truncated'] and sc['matched'] > 3,
          f"gwas by study: {sc.get('matched')} significant hits, top 3 by p")
    st, res, sc = S.call('snpmaize_similarity', {'line': 'B73', 'top': 3})
    sims = [r['ibs_similarity'] for r in sc.get('most_similar', [])]
    check(sc.get('ok') and len(sims) == 3 and sims == sorted(sims, reverse=True) and 0.9 < sims[0] <= 1 and sc['sites_compared'] == 46054265
          and sc['focal']['id'] == 'ZmG_B73' and all(r['line']['id'] != 'ZmG_B73' for r in sc['most_similar']),
          f"similarity genome-wide: B73's nearest {[r['line']['id'] for r in sc.get('most_similar', [])]} {sims}")

    # ---------------- store-backed tools ----------------
    if not has_store('chr4'):
        skip('chr4 store not installed: gene, variants, carriers, alleles, region similarity not checked')
    else:
        st, res, sc = S.call('snpmaize_gene', {'gene': 'tga1'})
        check(sc.get('ok') and sc['gene'] == 'Zm00001eb175150' and sc['symbol'] == 'tga1' and sc['location']['chr'] == 'chr4' and sc['sites_in_gene'] > 0
              and any(c['id'] == 'SC0001' for c in sc['curated_alleles']) and sc['links']['snpfunction'].endswith('#snpfunction/Zm00001eb175150')
              and sc['links']['maizegdb_gene'] == 'https://www.maizegdb.org/gene_center/gene/Zm00001eb175150' and isinstance(sc['gwas_within_10kb']['significant_hits'], int),
              f"gene tga1: {sc.get('sites_in_gene')} sites, QC {sc.get('by_site_qc')}, tiers {sc.get('by_tier_coding_usable')}, {len(sc.get('top_variants', []))} TOP/HIGH")
        gene_sites = sc['sites_in_gene']
        st, res, sc = S.call('snpmaize_variants', {'gene': 'Zm00001eb175150_T001', 'qc': 'all', 'limit': 500, 'sort': 'position'})
        rows = sc.get('variants', [])
        check(sc.get('ok') and sc['matched'] == gene_sites and len(rows) == gene_sites and [r['pos'] for r in rows] == sorted(r['pos'] for r in rows)
              and all(r['gene'] == 'Zm00001eb175150' for r in rows) and sc['links']['snpversity'].endswith('#snpversity?ds=zmgrin2026_imp&chr=chr4&start=46647932&end=46652896'),
              f"variants of the gene (qc=all, by position): {sc.get('matched')} rows = the gene summary's {gene_sites}; SNPVersity permalink")
        reg = next(r for r in rows if r['pos'] == 46648374 and r['alt'] == 'G')
        e = next(x for x in cu['entries'] if x['id'] == 'SC0001')
        check(reg['composite'] == e['site']['combined'] and reg['tier'] == e['site']['priority'] and reg['site_qc'] == e['site']['qc'] and reg['sub'] == e['site']['sub']
              and reg['n_hom'] == e['site']['nHom'] and reg['n_het'] == e['site']['nHet'] and reg['b73'] == e['site']['b73'],
              f"the tga1 N6K row agrees with snpcurate.data.js: composite {reg['composite']}, tier {reg['tier']}, QC {reg['site_qc']}, {reg['n_hom']} hom, B73 {reg['b73']}")
        st, res, sc = S.call('snpmaize_variants', {'gene': 'tga1', 'tier': 'HIGH', 'limit': 50})
        rows = sc.get('variants', [])
        comps = [r['composite'] for r in rows]
        check(sc.get('ok') and rows and all(r['tier'] in ('TOP', 'HIGH') and r['site_qc'] in ('PASS', 'HET_ELEVATED') for r in rows) and comps == sorted(comps)
              and sc['sites_in_window'] == gene_sites and sum(sc['by_site_qc_in_window'].values()) == gene_sites,
              f"variants tier>=HIGH: {len(rows)} rows, composite ascending {comps}, window counts cover all {gene_sites} sites")
        st, res, sc = S.call('snpmaize_variants', {'gene': 'tga1', 'qc': 'all', 'lines': ['TIL25', 'B73'], 'limit': 500})
        rows = sc.get('variants', [])
        check(sc.get('ok') and sc['filters']['lines'] == ['ZmG_C_Ames28398', 'ZmG_B73'] and any(r['pos'] == 46648374 for r in rows)
              and all(any(g in ('het', 'alt') for g in r['lines']['genotypes'].values()) for r in rows)
              and next(r for r in rows if r['pos'] == 46648374)['lines']['genotypes'] == {'ZmG_C_Ames28398': 'alt', 'ZmG_B73': 'ref'},
              f"variants restricted to carriers among TIL25 and B73: {len(rows)} sites, N6K alt in TIL25 and ref in B73")
        # carriers of every registry site on an installed chromosome agree with the built data file (first 6)
        done = 0
        for e in cu['entries']:
            if e.get('mark') == 'white' or not e.get('pos') or not has_store(e['chr']) or done >= 6: continue
            st, res, sc = S.call('snpmaize_carriers', {'curated': e['id'], 'limit': 200})
            s = e['site']
            ok = (sc.get('ok') and sc['genotypes']['alt_hom'] == s['nHom'] and sc['genotypes']['het'] == s['nHet'] and sc['genotypes']['missing'] == s['nMiss']
                  and sc['site']['composite'] == s['combined'] and sc['site']['tier'] == s['priority'] and sc['site']['site_qc'] == s['qc'] and sc['b73'] == s['b73']
                  and sc['carriers']['alt_hom']['n'] == len(s['carriers']['hom']) and set(l['id'] for l in sc['carriers']['alt_hom']['lines']) <= set(s['carriers']['hom'])
                  and sc['carriers']['alt_hom']['truncated'] == (len(s['carriers']['hom']) > 200)
                  and sc['carriers']['het']['n'] == len(s['carriers']['het']) and set(l['id'] for l in sc['carriers']['het']['lines']) <= set(s['carriers']['het'])
                  and sc['curated']['id'] == e['id'] and sc['links']['snpcurate'].endswith('#snpcurate/' + e['id']))
            check(ok, f"carriers {e['id']} {e.get('symbol')} {e['chr']}:{e['pos']}: {sc.get('genotypes')} = registry {s['nHom']}/{s['nHet']}/{s['nMiss']}, "
                      f"composite {sc.get('site', {}).get('composite')} tier {sc.get('site', {}).get('tier')} QC {sc.get('site', {}).get('site_qc')}, carrier ids from the registry")
            done += 1
        st, res, sc = S.call('snpmaize_carriers', {'site': 'chr4:46648374', 'trait': 'KERNEL-TYPE', 'limit': 0})
        check(sc.get('ok') and sc['carriers']['alt_hom']['n'] == 7 and sc['carriers']['alt_hom']['lines'] == [] and sc['by_subpop'][0]['subpop'] == 'Z.parviglumis'
              and sc['by_subpop'][0]['alt_allele_freq'] == 1.0 and sc['trait_by_genotype']['descriptor'] == 'KERNEL-TYPE' and 'codes' in sc['trait_by_genotype'],
              f"carriers by chr:pos with limit 0 and a coded trait: {sc.get('carriers', {}).get('alt_hom', {}).get('n')} hom, top subpop {sc.get('by_subpop', [{}])[0]}")
        st, res, sc = S.call('snpmaize_carriers', {'site': 'chr4:46648375'})
        check(res.get('isError') and sc['error']['code'] == 'not_found' and 'Nearby' in sc['error']['message'], 'carriers at a position with no site: not_found with nearby sites')
        st, res, sc = S.call('snpmaize_alleles', {'founders': True, 'sites': ['SC0001', 'chr4:46652277:C:T', 'SC9999']})
        col = sc.get('sites', [{}])[0]
        check(sc.get('ok') and len(sc['lines']) == 26 and col.get('curated', {}).get('id') == 'SC0001' and col['genotypes']['ZmG_B73'] == 'ref'
              and col['favorable_carry']['ZmG_B73'] == 'yes' and col['link'].endswith('#snpcurate/SC0001?lines=' + ','.join(l['id'] for l in sc['lines']))
              and sc['sites'][1]['curated'] is None and sc['sites'][1]['favorable_carry'] is None and 'error' in sc['sites'][2]
              and sc['favorable_sites_per_line']['ZmG_B73'] == 1,
              f"alleles of the 26 founders at SC0001 + a plain site + a bad id: B73 {col.get('genotypes', {}).get('ZmG_B73')}/{col.get('favorable_carry', {}).get('ZmG_B73')}, "
              f"bad id reported inline")
        st, res, sc = S.call('snpmaize_similarity', {'region': 'chr4:46600000-46700000', 'lines': ['B73', 'Mo17', 'CML103', 'W22']})
        check(sc.get('ok') and sc['mode'] == 'region' and sc['n_pairs'] == 6 and len(sc['pairs']) == 6 and all(0 <= p['ibs_similarity'] <= 1 and p['sites_co_called'] > 0 for p in sc['pairs'])
              and sc['sites_in_region'] > 100, f"similarity in a region, all pairs of 4 lines: 6 pairs over {sc.get('sites_in_region')} sites")
        st, res, sc = S.call('snpmaize_similarity', {'gene': 'tga1', 'line': 'B73', 'top': 5})
        check(sc.get('ok') and sc['lines_compared'] == 932 and len(sc['most_similar']) == 5 and sc['most_similar'][0]['ibs_similarity'] >= sc['most_similar'][-1]['ibs_similarity'],
              f"similarity in a gene, focal B73 against the panel: {sc.get('lines_compared')} lines, top {[r['line']['id'] for r in sc.get('most_similar', [])]}")
        st, res, sc = S.call('snpmaize_gwas', {'gene': 'tga1', 'flank': 50000, 'significant_only': False, 'limit': 5})
        check(sc.get('ok') and sc['region']['gene'] == 'Zm00001eb175150' and sc['region']['flank'] == 50000 and isinstance(sc['matched'], int),
              f"gwas around tga1 ±50 kb, all rows: {sc.get('matched')} rows in {sum(1 for _ in sc.get('studies', []))} studies")
        if has_store('chr4', '_HQ'):
            st, res, sc = S.call('snpmaize_variants', {'gene': 'tga1', 'dataset': 'zmgrin2026_hq', 'limit': 1})
            check(sc.get('ok') and sc['dataset'] == 'zmgrin2026_hq' and sc['links']['snpversity'].endswith('ds=zmgrin2026_hq&chr=chr4&start=46647932&end=46652896'),
                  f"variants in the High Quality set: {sc.get('sites_in_window')} sites in the window")
        else:
            skip('chr4 HQ store not installed')
finally:
    S.stop()

# ---------------- wall-time limit (a server of its own with a 0.01 s budget) ----------------
T = Server({'SNPMAIZE_MCP_TIMEOUT': '0.01'})
try:
    st, res, sc = T.call('snpmaize_datasets')
    check(res.get('isError') is True and sc['error']['code'] == 'timeout', f"a call over the time budget is stopped: {sc.get('error', {}).get('code')}")
finally:
    T.stop()

print(f'{"FAIL" if bad else "ok"}: check_mcp {bad} failures, {skipped} skipped')
sys.exit(1 if bad else 0)
