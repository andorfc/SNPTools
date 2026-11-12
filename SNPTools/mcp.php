<?php
/* =====================================================================
 *  mcp.php -- the SNPMaize MCP server (Model Context Protocol, Streamable HTTP, stateless).
 *
 *  POST <site>/mcp.php   JSON-RPC 2.0: initialize, ping, tools/list, tools/call (prompts/list
 *                        and resources/list answer empty lists; notifications get 202 and no body).
 *  GET  <site>/mcp.php   a short plain-text page: what the server is and how to connect.
 *  OPTIONS               CORS preflight (204).
 *
 *  Same transport as the MaizeGDB MCP (https://claude.maizegdb.org/mcp) so a client configured for
 *  one connects to the other unchanged; the two are meant to be used together: MaizeGDB for gene
 *  records, names, sequences, expression, GO and literature, SNPMaize for the variants, lines,
 *  passport, traits and curated alleles of the GRIN-linked 2026 release.
 *
 *  Every tool runs tools/snpq.py (one process per call, as processForm.php runs h5_to_vcf.py) and
 *  relays its JSON. This file validates the request (types, enums, caps), limits concurrency and
 *  wall time, and shapes the JSON-RPC answer; the data logic lives in snpq.py, which reuses the
 *  page mirrors of tools/build_snpcurate.py so an answer here agrees with the page.
 *
 *  Configuration (environment, as processForm.php): PYTHON_PATH (default python3 on PATH);
 *  SNPMAIZE_MCP_MAX_CONCURRENT (default 4) simultaneous tool calls, the rest answer 429-style
 *  isError results; SNPMAIZE_MCP_TIMEOUT seconds per call (default 60); SNPMAIZE_SITE the public
 *  URL used in permalinks (default: derived from the request). SNPTOOLS_DEBUG=1 adds the script's
 *  stderr to error results.
 * ===================================================================== */

final class SnpMcp {
  const PROTOCOL_LATEST = '2025-06-18';
  const SERVER_NAME = 'snpmaize';
  const SERVER_TITLE = 'SNPMaize MCP server';
  const VERSION = '1.0.0';
  const E_PARSE = -32700;
  const E_INVALID_REQUEST = -32600;
  const E_METHOD = -32601;
  const E_PARAMS = -32602;
  const E_INTERNAL = -32603;
  const MAX_BODY = 262144;
  public static $protocols = array('2025-06-18', '2025-03-26', '2024-11-05');

  /* ---- transport ------------------------------------------------------ */
  public static function serve() {
    $method = isset($_SERVER['REQUEST_METHOD']) ? $_SERVER['REQUEST_METHOD'] : 'GET';
    ob_start();                                   // nothing but the answer reaches the client (notices go to the log)
    ini_set('display_errors', '0');
    header('Access-Control-Allow-Origin: *');
    header('Access-Control-Allow-Methods: POST, GET, OPTIONS');
    header('Access-Control-Allow-Headers: Content-Type, Accept, Authorization, Mcp-Session-Id, MCP-Protocol-Version, Last-Event-ID');
    header('Access-Control-Expose-Headers: MCP-Protocol-Version');
    header('Access-Control-Max-Age: 86400');
    header('X-Content-Type-Options: nosniff');
    header('Cache-Control: no-store');
    header('Vary: Accept, Origin');
    if ($method === 'OPTIONS') { ob_end_clean(); http_response_code(204); return; }
    if ($method === 'GET' || $method === 'HEAD') {
      ob_end_clean();
      // A GET asking for text/event-stream (and not a page) is a client opening the server-sent event
      // stream. This server is stateless and has none, which Streamable HTTP says to answer with 405.
      $accept = isset($_SERVER['HTTP_ACCEPT']) ? $_SERVER['HTTP_ACCEPT'] : '';
      if (stripos($accept, 'text/event-stream') !== false && stripos($accept, 'text/html') === false) {
        http_response_code(405);
        header('Allow: POST, GET, OPTIONS');
        header('Content-Type: application/json; charset=utf-8');
        echo json_encode(self::error(null, self::E_INVALID_REQUEST, 'This server opens no event stream; send JSON-RPC requests with POST.'));
        return;
      }
      header('Content-Type: text/plain; charset=utf-8');
      echo self::page();
      return;
    }
    if ($method !== 'POST') {
      ob_end_clean();
      http_response_code(405);
      header('Allow: POST, GET, OPTIONS');
      header('Content-Type: application/json; charset=utf-8');
      echo json_encode(self::error(null, self::E_INVALID_REQUEST, 'Use POST with a JSON-RPC 2.0 body.'));
      return;
    }
    $raw = file_get_contents('php://input', false, null, 0, self::MAX_BODY + 1);
    if (strlen($raw) > self::MAX_BODY) {
      ob_end_clean();
      http_response_code(413);
      header('Content-Type: application/json; charset=utf-8');
      echo json_encode(self::error(null, self::E_INVALID_REQUEST, 'Request body larger than 256 KB.'));
      return;
    }
    $r = self::handleRaw($raw);
    ob_end_clean();
    http_response_code($r['status']);
    header('MCP-Protocol-Version: ' . self::PROTOCOL_LATEST);
    if ($r['body'] !== null) {
      header('Content-Type: application/json; charset=utf-8');
      header('Content-Length: ' . strlen($r['body']));
      echo $r['body'];
    }
  }

  /* A raw body -> {status, body|null}. Batches are answered as batches; notifications alone give 202. */
  public static function handleRaw($raw) {
    $msg = json_decode($raw, true);
    if ($raw === '' || $msg === null && json_last_error() !== JSON_ERROR_NONE) {
      return array('status' => 400, 'body' => json_encode(self::error(null, self::E_PARSE, 'Body is not JSON.')));
    }
    $batch = is_array($msg) && array_keys($msg) === range(0, count($msg) - 1) && count($msg) > 0;
    $items = $batch ? $msg : array($msg);
    if (count($items) > 20) {
      return array('status' => 400, 'body' => json_encode(self::error(null, self::E_INVALID_REQUEST, 'At most 20 requests per batch.')));
    }
    $out = array();
    foreach ($items as $m) {
      $res = self::handle($m);
      if ($res !== null) $out[] = $res;
    }
    if (!count($out)) return array('status' => 202, 'body' => null);
    return array('status' => 200, 'body' => json_encode($batch ? $out : $out[0], JSON_UNESCAPED_SLASHES | JSON_UNESCAPED_UNICODE));
  }

  /* One message -> one response array, or null for a notification. */
  public static function handle($m) {
    if (!is_array($m) || !isset($m['jsonrpc']) || $m['jsonrpc'] !== '2.0' || !isset($m['method']) || !is_string($m['method'])) {
      return self::error(is_array($m) && isset($m['id']) ? $m['id'] : null, self::E_INVALID_REQUEST, 'Not a JSON-RPC 2.0 request.');
    }
    $id = array_key_exists('id', $m) ? $m['id'] : null;
    $notification = !array_key_exists('id', $m);
    $params = isset($m['params']) && is_array($m['params']) ? $m['params'] : array();
    if (strpos($m['method'], 'notifications/') === 0) return null;
    switch ($m['method']) {
      case 'initialize':
        $asked = isset($params['protocolVersion']) ? (string) $params['protocolVersion'] : '';
        $version = in_array($asked, self::$protocols, true) ? $asked : self::PROTOCOL_LATEST;
        $result = array(
          'protocolVersion' => $version,
          'capabilities' => array('tools' => array('listChanged' => false), 'prompts' => array('listChanged' => false),
                                  'resources' => array('subscribe' => false, 'listChanged' => false)),
          'serverInfo' => array('name' => self::SERVER_NAME, 'title' => self::SERVER_TITLE, 'version' => self::VERSION),
          'instructions' => self::instructions());
        break;
      case 'ping': $result = new stdClass(); break;
      case 'tools/list': $result = array('tools' => self::toolList()); break;
      case 'tools/call': return $notification ? null : self::toolsCall($id, $params);
      case 'prompts/list': $result = array('prompts' => array()); break;
      case 'resources/list': $result = array('resources' => array()); break;
      case 'resources/templates/list': $result = array('resourceTemplates' => array()); break;
      case 'logging/setLevel': $result = new stdClass(); break;
      default:
        return $notification ? null : self::error($id, self::E_METHOD, 'Method not found: ' . $m['method']);
    }
    return $notification ? null : array('jsonrpc' => '2.0', 'id' => $id, 'result' => $result);
  }

  private static function error($id, $code, $message, $data = null) {
    $e = array('code' => $code, 'message' => $message);
    if ($data !== null) $e['data'] = $data;
    return array('jsonrpc' => '2.0', 'id' => $id, 'error' => $e);
  }

  /* ---- tools/call ------------------------------------------------------ */
  private static function toolsCall($id, $params) {
    $name = isset($params['name']) ? (string) $params['name'] : '';
    $args = isset($params['arguments']) && is_array($params['arguments']) ? $params['arguments'] : array();
    $tools = self::tools();
    if (!isset($tools[$name])) return self::error($id, self::E_PARAMS, 'Unknown tool: ' . $name);
    $t0 = microtime(true);
    $isError = false;
    try {
      $clean = self::validate($tools[$name], $args);
      $payload = self::runQuery($name, $clean);
      if (!isset($payload['ok']) || !$payload['ok']) {
        $isError = true;
        $err = isset($payload['error']) && is_array($payload['error']) ? $payload['error'] : array('code' => 'internal', 'message' => 'The query failed.');
        $payload = array('ok' => false, 'tool' => $name, 'error' => $err);
      } else {
        unset($payload['ok']);
        $payload = array_merge(array('ok' => true, 'tool' => $name), $payload);
      }
    } catch (SnpMcpToolError $e) {
      $isError = true;
      $payload = array('ok' => false, 'tool' => $name, 'error' => array('code' => $e->kind, 'message' => $e->getMessage()));
      if ($e->extra) $payload['error'] = array_merge($payload['error'], $e->extra);
    }
    $ms = (int) round((microtime(true) - $t0) * 1000);
    self::log(array('event' => 'tools/call', 'tool' => $name, 'ok' => !$isError, 'ms' => $ms, 'args' => self::logArgs($args)));
    $text = json_encode($payload, JSON_UNESCAPED_SLASHES | JSON_UNESCAPED_UNICODE | JSON_PRETTY_PRINT);
    if ($text === false) {
      $isError = true;
      $payload = array('ok' => false, 'tool' => $name, 'error' => array('code' => 'internal', 'message' => 'The answer could not be encoded.'));
      $text = json_encode($payload);
    }
    return array('jsonrpc' => '2.0', 'id' => $id, 'result' => array(
      'content' => array(array('type' => 'text', 'text' => $text)),
      'structuredContent' => $payload,
      'isError' => $isError));
  }

  /* Run tools/snpq.py for one call: concurrency gate, wall-time limit, JSON back. */
  private static function runQuery($name, $args) {
    $root = __DIR__;
    $script = $root . '/tools/snpq.py';
    if (!is_file($script)) throw new SnpMcpToolError('unavailable', 'The query library (tools/snpq.py) is not installed.');
    $python = getenv('PYTHON_PATH'); if (!$python) $python = 'python3';
    $max = (int) getenv('SNPMAIZE_MCP_MAX_CONCURRENT'); if ($max < 1) $max = 4;
    $timeout = (float) getenv('SNPMAIZE_MCP_TIMEOUT'); if ($timeout <= 0) $timeout = 60.0;
    $slot = self::acquireSlot($max);
    if ($slot === null) throw new SnpMcpToolError('busy', 'The server is answering ' . $max . ' queries already; try again in a few seconds.', array('retry_after_s' => 3));
    $req = json_encode(array('tool' => $name, 'args' => $args, 'site' => self::siteUrl()), JSON_UNESCAPED_SLASHES | JSON_UNESCAPED_UNICODE);
    $errf = tempnam(sys_get_temp_dir(), 'snpq');
    $spec = array(0 => array('pipe', 'r'), 1 => array('pipe', 'w'), 2 => array('file', $errf, 'w'));   // child reads 0, writes 1 and 2
    $env = array();
    foreach (array('PATH', 'HOME', 'LANG', 'LC_ALL', 'TMPDIR', 'HDF5_USE_FILE_LOCKING') as $k) { $v = getenv($k); if ($v !== false) $env[$k] = $v; }
    $env['PYTHONIOENCODING'] = 'utf-8';
    $env['SNPQ_TIMEOUT'] = (string) $timeout;         // the library stops itself at the budget and answers {"error": {"code": "timeout"}}
    if (!isset($env['HDF5_USE_FILE_LOCKING'])) $env['HDF5_USE_FILE_LOCKING'] = 'FALSE';
    @set_time_limit((int) $timeout + 30);
    $proc = proc_open(escapeshellarg($python) . ' -I ' . escapeshellarg($script), $spec, $pipes, $root, $env);
    if (!is_resource($proc)) { self::releaseSlot($slot); @unlink($errf); throw new SnpMcpToolError('unavailable', 'The query process could not be started.'); }
    fwrite($pipes[0], $req); fclose($pipes[0]);
    $out = stream_get_contents($pipes[1]); fclose($pipes[1]);
    $code = proc_close($proc);
    $err = (string) @file_get_contents($errf); @unlink($errf);
    self::releaseSlot($slot);
    $json = json_decode($out, true);
    if (!is_array($json)) {
      self::log(array('event' => 'snpq_failed', 'tool' => $name, 'exit' => $code, 'stderr' => substr($err, -2000)));
      $extra = getenv('SNPTOOLS_DEBUG') === '1' ? array('stderr' => substr($err, -2000), 'exit' => $code) : array();
      throw new SnpMcpToolError('internal', 'The query failed on the server (exit ' . $code . ').', $extra);
    }
    return $json;
  }

  /* Up to $max calls at once: the first free of $max lock files (non-blocking flock). */
  private static function acquireSlot($max) {
    $dir = sys_get_temp_dir();
    for ($i = 0; $i < $max; $i++) {
      $fh = @fopen($dir . '/snpmaize-mcp-slot-' . $i . '.lock', 'c');
      if (!$fh) continue;
      if (flock($fh, LOCK_EX | LOCK_NB)) return $fh;
      fclose($fh);
    }
    return null;
  }

  private static function releaseSlot($fh) { if ($fh) { flock($fh, LOCK_UN); fclose($fh); } }

  /* The public URL of the site for permalinks: SNPMAIZE_SITE, else this request's scheme, host and folder. */
  public static function siteUrl() {
    $s = getenv('SNPMAIZE_SITE');
    if ($s) return rtrim($s, '/') . '/';
    if (empty($_SERVER['HTTP_HOST'])) return '';
    $https = (!empty($_SERVER['HTTPS']) && $_SERVER['HTTPS'] !== 'off')
          || (!empty($_SERVER['HTTP_X_FORWARDED_PROTO']) && $_SERVER['HTTP_X_FORWARDED_PROTO'] === 'https');
    $host = preg_replace('/[^A-Za-z0-9.:\-\[\]]/', '', $_SERVER['HTTP_HOST']);
    $dir = isset($_SERVER['SCRIPT_NAME']) ? rtrim(str_replace('\\', '/', dirname($_SERVER['SCRIPT_NAME'])), '/') : '';
    return ($https ? 'https' : 'http') . '://' . $host . $dir . '/';
  }

  /* ---- argument validation (types, enums, caps; the library re-checks and reports) ---- */
  private static function validate($tool, $args) {
    $schema = $tool['inputSchema'];
    $props = isset($schema['properties']) ? $schema['properties'] : array();
    $clean = array();
    foreach ($args as $k => $v) {
      if (!isset($props[$k])) throw new SnpMcpToolError('invalid_params', 'Unknown argument "' . $k . '" for ' . $tool['name'] . '; it takes ' . implode(', ', array_keys($props)) . '.');
      if ($v === null) continue;
      $p = $props[$k];
      $type = isset($p['type']) ? $p['type'] : 'string';
      if ($type === 'integer') {
        if (!is_int($v) && !(is_string($v) && preg_match('/^-?\d+$/', trim($v))) && !(is_float($v) && floor($v) == $v)) throw new SnpMcpToolError('invalid_params', $k . ': an integer is required.');
        $v = (int) $v;
        if (isset($p['minimum']) && $v < $p['minimum']) throw new SnpMcpToolError('invalid_params', $k . ': at least ' . $p['minimum'] . '.');
        if (isset($p['maximum']) && $v > $p['maximum']) throw new SnpMcpToolError('too_large', $k . ': at most ' . $p['maximum'] . '.');
      } elseif ($type === 'number') {
        if (!is_numeric($v)) throw new SnpMcpToolError('invalid_params', $k . ': a number is required.');
        $v = (float) $v;
      } elseif ($type === 'boolean') {
        if (is_string($v)) $v = in_array(strtolower(trim($v)), array('1', 'true', 'yes'), true);
        $v = (bool) $v;
      } elseif ($type === 'array') {
        if (is_string($v)) $v = preg_split('/[,\n;]+/', $v, -1, PREG_SPLIT_NO_EMPTY);
        if (!is_array($v)) throw new SnpMcpToolError('invalid_params', $k . ': a list is required.');
        $v = array_values(array_filter(array_map(function ($x) { return is_scalar($x) ? trim((string) $x) : ''; }, $v), 'strlen'));
        if (isset($p['maxItems']) && count($v) > $p['maxItems']) throw new SnpMcpToolError('too_large', $k . ': at most ' . $p['maxItems'] . ' items.');
        if (isset($p['items']['enum'])) foreach ($v as $x) if (!in_array($x, $p['items']['enum'], true)) throw new SnpMcpToolError('invalid_params', $k . ': each item one of ' . implode(', ', $p['items']['enum']) . '.');
      } else {
        if (!is_scalar($v)) throw new SnpMcpToolError('invalid_params', $k . ': a string is required.');
        $v = trim((string) $v);
        if (strlen($v) > 4000) throw new SnpMcpToolError('too_large', $k . ': at most 4000 characters.');
        if ($v === '') continue;
        if (isset($p['enum']) && !in_array($v, $p['enum'], true)) throw new SnpMcpToolError('invalid_params', $k . ': one of ' . implode(', ', $p['enum']) . '.');
      }
      $clean[$k] = $v;
    }
    foreach ((isset($schema['required']) ? $schema['required'] : array()) as $k) {
      if (!isset($clean[$k])) throw new SnpMcpToolError('invalid_params', $tool['name'] . ' needs "' . $k . '".');
    }
    if (isset($tool['oneOfArgs'])) {
      $have = 0; foreach ($tool['oneOfArgs'] as $k) if (isset($clean[$k])) $have++;
      if (!$have) throw new SnpMcpToolError('invalid_params', $tool['name'] . ' needs one of ' . implode(', ', $tool['oneOfArgs']) . '.');
    }
    return $clean;
  }

  /* ---- the tools ------------------------------------------------------ */
  private static function str($d, $extra = array()) { return array_merge(array('type' => 'string', 'description' => $d), $extra); }
  private static function int($d, $min, $max, $extra = array()) { return array_merge(array('type' => 'integer', 'description' => $d, 'minimum' => $min, 'maximum' => $max), $extra); }
  private static function bool($d) { return array('type' => 'boolean', 'description' => $d); }
  private static function lst($d, $max, $extra = array()) { return array_merge(array('type' => 'array', 'items' => array('type' => 'string'), 'description' => $d, 'maxItems' => $max), $extra); }
  private static function ro($extra = array()) { return array_merge(array('readOnlyHint' => true, 'destructiveHint' => false, 'idempotentHint' => true, 'openWorldHint' => false), $extra); }

  public static function tools() {
    static $tools = null;
    if ($tools !== null) return $tools;
    $ds = self::str('Variant set: zmgrin2026_imp (Imputed, default: every line called at every site) or zmgrin2026_hq (High Quality subset of the same sites).', array('enum' => array('zmgrin2026_imp', 'zmgrin2026_hq')));
    $gene = self::str('B73 v5 gene model id (Zm00001eb175150; a transcript id like Zm00001eb175150_T001 is accepted) or a symbol from the SNPCurate registry (tga1, su1). Other names: resolve them first with the MaizeGDB MCP (maizegdb_resolve).');
    $region = self::str('Region on B73 v5 as chr:start-end, e.g. chr4:46640000-46660000 (commas allowed). At most 250,000 bp per call.');
    $flank = self::int('Base pairs added on each side of the gene (0-50,000; default 0).', 0, 50000);
    $limit = self::int('Rows to return (default 100, at most 500).', 1, 500);
    $defs = array(
      array('name' => 'snpmaize_datasets', 'title' => 'SNPMaize: release, datasets and conventions',
        'description' => 'What SNPMaize holds and how to ask: the GRIN-linked 2026 release (933 lines, 46 M B73 v5 sites), the two variant sets, the site-QC codes and rule, the line-QC classes, the per-site AI scores and the composite/tier calibration, the SNPCurate registry counts, the GRIN trait descriptors (keys), the GWAS studies, the per-call caps and the name conventions. Call once before the other tools when unsure.',
        'inputSchema' => array('type' => 'object', 'properties' => new stdClass()), 'annotations' => self::ro()),
      array('name' => 'snpmaize_lines', 'title' => 'Find lines by passport or trait',
        'description' => 'Lines of the panel that match passport facets (subpopulation, country, state, NAM/Ames282/WiDiv panel, kernel type, improvement status, era, species), a name query (id, line name, GRIN accession) and/or a GRIN trait descriptor with a value or range. Returns the matching lines with line QC, facet counts of the match and, with a trait, the descriptor\'s panel-wide distribution. A trait without a value just describes the descriptor.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'query' => self::str('Text matched against line id, name, GRIN accession (PI 550473) and genotype id; exact names come first.'),
          'pedigree' => self::bool('Also match the query against pedigrees (default false).'),
          'subpop' => self::str('Subpopulation: SS, NSS, IDT, Tropical, Mixed, Broad origin-public, Landrace, Popcorn, Sweet, Z.parviglumis ...'),
          'country' => self::str('Country of origin, e.g. Mexico.'), 'state' => self::str('US state of origin.'),
          'panel' => self::str('Membership: NAM, Ames282 or WiDiv.', array('enum' => array('NAM', 'Ames282', 'WiDiv'))),
          'kernel_type' => self::str('Dent, Semident, Semiflint, Flint, Pop, Sweet, Floury, Waxy ...'),
          'improvement' => self::str('Improvement status, e.g. Breeding material, Landrace, Cultivar.'),
          'era' => self::str('PublicBreeding, Landrace, ExPVP ...'), 'species' => self::str('Z.mays, Z.parviglumis ...'),
          'trait' => self::str('GRIN descriptor key or name (KERNEL-WEIGHT-1000, Plant Height, Kernel Type); snpmaize_datasets lists the keys.'),
          'trait_value' => self::str('For coded descriptors: the code or its label (Flint, Red).'),
          'trait_min' => array('type' => 'number', 'description' => 'For numeric descriptors: lowest accepted mean.'),
          'trait_max' => array('type' => 'number', 'description' => 'For numeric descriptors: highest accepted mean.'),
          'limit' => self::int('Lines to return (default 50, at most 200).', 1, 200))),
        'annotations' => self::ro()),
      array('name' => 'snpmaize_line', 'title' => 'One line: passport, traits, QC',
        'description' => 'Everything SNPMaize knows about one line: passport (GRIN accession, origin, pedigree, subpopulation, panels), every GRIN trait observation decoded, its line QC, whether it is a NAM founder, and the curated (SNPCurate) alleles it carries, with permalinks.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'line' => self::str('Line id (ZmG_B73), name (B73, Mo17, CML 103), GRIN accession (PI 550473) or genotype id (ZmG:B73).')), 'required' => array('line')),
        'annotations' => self::ro()),
      array('name' => 'snpmaize_variants', 'title' => 'Variants in a gene or region',
        'description' => 'The variant sites of a gene or a region with consequence, amino-acid change, SnpEff impact, the per-site AI scores (PlantCAD1/2, Evo2, ESM1/2/3/C), the composite z-score, percentile and priority tier, MAF, LD max R², carrier counts and site QC. Filter by consequence class, minimum tier, QC set and MAF; optionally restrict to sites where given lines carry the alternate allele (and see their genotypes). Sorted most damaging first by default; counts cover the whole window even when rows are capped.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'gene' => $gene, 'region' => $region, 'flank' => $flank, 'dataset' => $ds,
          'consequence' => self::lst('Keep these classes only: missense, lof (stop gained, frameshift, start/stop lost), splice, indel (in-frame), syn, other.', 6, array('items' => array('type' => 'string', 'enum' => array('missense', 'lof', 'splice', 'indel', 'syn', 'other')))),
          'tier' => self::str('Keep this priority tier and above: TOP (top 5 % of the composite, or a severe consequence), HIGH (top 25 %), MODERATE (top 80 %), LOW (all).', array('enum' => array('TOP', 'HIGH', 'MODERATE', 'LOW'))),
          'qc' => self::str('Site-QC set: usable (PASS or HET_ELEVATED, default), pass (PASS only) or all.', array('enum' => array('usable', 'pass', 'all'))),
          'min_maf' => array('type' => 'number', 'description' => 'Lowest minor-allele frequency kept.'),
          'lines' => self::lst('Keep only sites where at least one of these lines carries the alternate allele; their genotypes are reported per site (at most 100 names).', 100),
          'gene_only' => self::bool('With a gene: keep only sites annotated to that gene (default true).'),
          'sort' => self::str('composite (most damaging first, default), position or maf.', array('enum' => array('composite', 'position', 'maf'))),
          'limit' => $limit)),
        'oneOfArgs' => array('gene', 'region'), 'annotations' => self::ro()),
      array('name' => 'snpmaize_gene', 'title' => 'Gene summary',
        'description' => 'One gene at a glance: B73 v5 coordinates, how many sites it has by QC and consequence class, the coding-site tier counts, the TOP/HIGH variants, its SNPCurate registry alleles, significant GWAS hits within 10 kb, whether protein site annotations exist, and permalinks to SNPFunction, SNPFold, SNPGeo, PanEffect and the MaizeGDB gene record.',
        'inputSchema' => array('type' => 'object', 'properties' => array('gene' => $gene, 'dataset' => $ds), 'required' => array('gene')),
        'annotations' => self::ro()),
      array('name' => 'snpmaize_carriers', 'title' => 'Who carries a variant',
        'description' => 'One site: its annotation and scores, genotype counts, the lines homozygous and heterozygous for the alternate allele (with passport), the breakdown by subpopulation and country, B73\'s genotype, the carriers\' line QC, the registry entry if it is a curated allele, and optionally the mean of a GRIN trait by genotype class.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'site' => self::str('chr:pos or chr:pos:ref:alt on B73 v5, e.g. chr4:46648374:C:G.'),
          'curated' => self::str('A SNPCurate registry id (SC0001) instead of a site.'),
          'ref' => self::str('Reference allele, when chr:pos alone is ambiguous.'), 'alt' => self::str('Alternate allele.'),
          'dataset' => $ds,
          'trait' => self::str('GRIN descriptor to summarise by genotype class (numeric: means; coded: counts).'),
          'limit' => self::int('Carrier names per class (default 50, at most 200; 0 for counts only).', 0, 200))),
        'oneOfArgs' => array('site', 'curated'), 'annotations' => self::ro()),
      array('name' => 'snpmaize_curated', 'title' => 'SNPCurate registry',
        'description' => 'Published, curated maize alleles (SNPCurate): one record by id with its sources, note, expected and observed carriers, favorable allele and scores; or the list filtered by gene/symbol, text (trait, label, note), mark (gold = published causal change at a site of this release; silver = published variant or marker at a site; white = published variant not in this marker set) and role.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'id' => self::str('Registry id, e.g. SC0001.'), 'gene' => self::str('Gene id or symbol.'),
          'query' => self::str('Text matched against symbol, gene, label, trait, note and role.'),
          'mark' => self::str('gold, silver or white.', array('enum' => array('gold', 'silver', 'white'))),
          'role' => self::str('causal, proposed, allele, associated, lead_snp or marker.'),
          'limit' => self::int('Entries to return (default 50, at most 200).', 1, 200))),
        'annotations' => self::ro()),
      array('name' => 'snpmaize_alleles', 'title' => 'Allele passport of lines',
        'description' => 'Genotypes of chosen lines at chosen sites: registry ids (SC0001) and/or chr:pos[:ref:alt] sites, for named lines or the 26 NAM founders. For curated sites with a favorable allele, says per line whether it carries it (yes/het/no). Mirrors the "Your lines" panel of SNPCurate; the permalinks open it.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'lines' => self::lst('Line names, ids or GRIN accessions (at most 100).', 100),
          'founders' => self::bool('Use the 26 NAM founders as the lines.'),
          'sites' => self::lst('Registry ids and/or chr:pos:ref:alt sites (at most 50).', 50), 'dataset' => $ds),
        'required' => array('sites')), 'oneOfArgs' => array('lines', 'founders'), 'annotations' => self::ro()),
      array('name' => 'snpmaize_similarity', 'title' => 'Identity-by-state between lines',
        'description' => 'Genome-wide: the lines most and least similar to a focal line from the precomputed IBS matrix (all sites or SNPs only), optionally within a named subset. In a gene or region: IBS computed from the genotypes for a focal line against the panel or a subset, or all pairs among up to 100 named lines.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'line' => self::str('Focal line.'), 'lines' => self::lst('Lines to compare (subset, or all pairs when no focal line is given; at most 100).', 100),
          'gene' => $gene, 'region' => $region, 'flank' => $flank, 'dataset' => $ds,
          'mode' => self::str('genome (default without a gene/region) or region.', array('enum' => array('genome', 'region'))),
          'sites' => self::str('Genome-wide matrix: all (default) or snp.', array('enum' => array('all', 'snp'))),
          'top' => self::int('Most-similar lines to return (default 20, at most 100).', 1, 100))),
        'oneOfArgs' => array('line', 'lines'), 'annotations' => self::ro()),
      array('name' => 'snpmaize_gwas', 'title' => 'GWAS hits',
        'description' => 'Published GWAS results carried by SNPGWAS (Tibbs-Cortes et al. 2024 NAM flowering, yield-component and morphology traits and the other studies in the manifest). With no arguments: the study list. With a gene (±10 kb by default) or region, and optionally a study or trait: the hits there, significant by the study\'s own threshold unless significant_only is false. With a study alone: its top hits genome-wide. Positions are B73 v5; the GWAS panels are not the SNPMaize lines.',
        'inputSchema' => array('type' => 'object', 'properties' => array(
          'gene' => $gene, 'region' => $region, 'flank' => self::int('Base pairs around the gene (default 10,000; at most 50,000).', 0, 50000),
          'study' => self::str('Study id (tibbscortes2024_nam_dta_intcp) or a distinctive part of it.'),
          'trait' => self::str('Trait or category text, e.g. Days to Anthesis, Flowering Time.'),
          'significant_only' => self::bool('Only hits at or below the study threshold (default true).'), 'limit' => $limit)),
        'annotations' => self::ro()),
    );
    $tools = array();
    foreach ($defs as $d) $tools[$d['name']] = $d;
    return $tools;
  }

  public static function toolList() {
    $out = array();
    foreach (self::tools() as $t) {
      $out[] = array('name' => $t['name'], 'title' => $t['title'], 'description' => $t['description'],
                     'inputSchema' => $t['inputSchema'], 'annotations' => $t['annotations']);
    }
    return $out;
  }

  public static function instructions() {
    return "SNPMaize is MaizeGDB's variant browser for the GRIN-linked 2026 release: 933 maize lines, each a USDA GRIN accession, "
         . "called at 46 million B73 v5 sites, with per-site AI damage scores (PlantCAD, Evo2, ESM), site QC, GRIN passport and trait "
         . "descriptors, a registry of published alleles (SNPCurate), genome-wide identity-by-state and published GWAS. "
         . "Coordinates are 1-based on B73 v5 (chr1..chr10). Start with snpmaize_datasets when unsure what a field means. "
         . "Genes are B73 v5 ids (Zm00001eb...); turn a gene name into an id with the MaizeGDB MCP (maizegdb_resolve, maizegdb_record) "
         . "at " . 'https://claude.maizegdb.org/mcp' . ", then query here. Lines resolve by id, name or GRIN accession. "
         . "Answers are capped (regions up to 250 kb, 500 rows, 100 lines); counts in an answer cover the whole window even when rows are cut. "
         . "Every answer carries `links`: permalinks that open the same view on the site; give them to the user. "
         . "Damage scores and tiers are predictions, not phenotypes; SNPCurate gold marks are the published causal changes. "
         . "Trait means by genotype describe this panel and are not association tests.";
  }

  private static function page() {
    $site = self::siteUrl();
    $names = array_keys(self::tools());
    return self::SERVER_TITLE . " " . self::VERSION . "\n\n"
      . "A Model Context Protocol server (JSON-RPC 2.0 over HTTP, stateless) for the SNPMaize instance of SNPTools.\n"
      . "Connect an MCP client to this URL (POST). Tools: " . implode(', ', $names) . ".\n\n"
      . "Site: " . ($site ? $site : '(this host)') . "\nDocumentation: " . ($site ? $site : '') . "MCP.md\nCompanion server for gene records, names, sequences, expression, GO and literature: https://claude.maizegdb.org/mcp\n\n"
      . "Try:  curl -s -X POST -H 'Content-Type: application/json' " . ($site ? $site . 'mcp.php' : 'mcp.php')
      . " -d '{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"tools/list\"}'\n";
  }

  /* ---- log ------------------------------------------------------------ */
  private static function log($fields) {
    $fields['ts'] = gmdate('c');
    $fields['ip'] = isset($_SERVER['REMOTE_ADDR']) ? $_SERVER['REMOTE_ADDR'] : null;
    error_log('snpmaize-mcp ' . json_encode($fields, JSON_UNESCAPED_SLASHES));
  }

  private static function logArgs($args) {
    $out = array();
    foreach ($args as $k => $v) {
      if (is_array($v)) $out[$k] = count($v) . ' items';
      elseif (is_string($v)) $out[$k] = strlen($v) > 80 ? substr($v, 0, 77) . '...' : $v;
      else $out[$k] = $v;
    }
    return $out;
  }
}

final class SnpMcpToolError extends Exception {
  public $kind; public $extra;
  public function __construct($code, $message, $extra = array()) { parent::__construct($message); $this->kind = $code; $this->extra = $extra; }
}

if (!defined('SNPMAIZE_MCP_NO_SERVE')) SnpMcp::serve();
