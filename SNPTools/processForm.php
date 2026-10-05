<?php
/* =====================================================================
 *  processForm.php — region + accession list  ->  VCF (via h5_to_vcf.py)
 *
 *  Local testing (MAMP) vs production (Linux/Docker): nothing in this file
 *  changes between machines. The Python interpreter is resolved from the
 *  PYTHON_PATH environment variable, falling back to 'python3' on PATH.
 *  Python deps are declared in requirements.txt (pip-installed at image build).
 * ===================================================================== */

header('Content-Type: application/json');

/* Building a VCF for a wide region / large accession set legitimately takes
 * more than PHP's default 30 s max_execution_time (which on Windows counts
 * the blocking shell_exec below as wall-clock). Lift it so the request runs
 * to completion instead of dying with a Fatal error mid-stream. */
set_time_limit(0);

/* ---------------------------------------------------------------------
 *  CONFIG — provided by the environment; no machine paths committed to git
 * ------------------------------------------------------------------- */
// 1) Python interpreter that has h5py + numpy installed. Resolution order:
//      a. PYTHON_PATH environment variable (first hit wins). Set it per host,
//         not in code, so nothing here is machine-specific:
//           - Docker:  ENV PYTHON_PATH=python3   (or just rely on the fallback)
//           - MAMP:    SetEnv PYTHON_PATH /path/to/venv/bin/python  in an Apache
//                      conf / .htaccess, or export it in the shell that starts MAMP
//                      — keeps your local absolute path out of the repo.
//      b. 'python3' on the system PATH (the default). In the Docker image the
//         deps from requirements.txt are pip-installed globally, so python3 works
//         out of the box with no configuration.
$PYTHON_PATH = getenv('PYTHON_PATH');
if (!$PYTHON_PATH) {
    $PYTHON_PATH = 'python3';                     // resolved via PATH — portable across OSes
}

// 2) Where the .h5 files live, relative to this PHP file.
$VERSION_PATH = './hdf5/grin2026/';

// 3) Where VCFs are written (must be web-served AND writable). Matches CFG.vcfDir in data.js.
$VCF_DIR = './vcf/';

// 4) How long a built VCF stays downloadable. Every query (and every SNPGeo / SNPFunction /
//    SNPFold gene view) writes one; older ones are removed, at most once every 10 minutes.
//    SNPTOOLS_VCF_TTL_HOURS overrides the default of 24 hours; 0 turns the clean-up off.
$VCF_TTL_HOURS = getenv('SNPTOOLS_VCF_TTL_HOURS');
$VCF_TTL_HOURS = ($VCF_TTL_HOURS === false || trim($VCF_TTL_HOURS) === '') ? 24.0 : (float) $VCF_TTL_HOURS;

// 5) The largest request h5_to_vcf.py builds (variants x accessions) is SNPTOOLS_MAX_CELLS,
//    read by the script itself (default 2e9). A list longer than this many ids is refused here.
$MAX_IDS = 20000;

// 6) Replies never name server paths. The script's command line and output (which do) go to
//    the PHP error log; SNPTOOLS_DEBUG=1 also returns them in the reply, for local debugging.
$DEBUG = getenv('SNPTOOLS_DEBUG') === '1';
function reply($fields, $debug = array()) {
    global $DEBUG;
    echo json_encode($DEBUG ? array_merge($fields, $debug) : $fields);
}

/* ---------------------------------------------------------------------
 *  INPUT
 * ------------------------------------------------------------------- */
if ($_SERVER['REQUEST_METHOD'] !== 'POST') {
    echo json_encode(array('status' => 'error', 'message' => 'POST required'));
    exit;
}

$start         = isset($_POST['start'])     ? $_POST['start']     : '';
$end           = isset($_POST['end'])       ? $_POST['end']       : '';
$chr           = isset($_POST['chr'])       ? $_POST['chr']       : '';
$dataset       = isset($_POST['dataSet'])   ? $_POST['dataSet']   : '';
$genotypesJson = isset($_POST['genotypes']) ? $_POST['genotypes'] : '[]';
$outName       = isset($_POST['outName'])   ? $_POST['outName']   : '';

// integer interval (is_numeric() also passed "1e5" and "2.5", which the script cannot read)
$start = trim((string) $start); $end = trim((string) $end);
if (!preg_match('/^-?[0-9]{1,12}$/', $start) || !preg_match('/^-?[0-9]{1,12}$/', $end)) {
    echo json_encode(array('status' => 'error', 'message' => 'Invalid interval (start/end must be whole numbers).'));
    exit;
}
// chromosome token must look like chr10 (used verbatim in the .h5 filename)
if (!preg_match('/^chr[0-9]{1,2}$/', $chr)) {
    echo json_encode(array('status' => 'error', 'message' => "Invalid chromosome '$chr'."));
    exit;
}

/* ---------------------------------------------------------------------
 *  DATASET  ->  (family, quality)  ->  <family>_<chr>_<quality>.h5
 * ------------------------------------------------------------------- */
$ds_part0 = null;   // family        (every dataset id is listed below; anything else is refused)
$ds_part2 = null;   // quality tier

switch ($dataset) {
    case 'mgdb2026_hq':  $ds_part0 = 'maizegdb2026'; $ds_part2 = 'HQ';     break;
    case 'mgdb2026_hc':  $ds_part0 = 'maizegdb2026'; $ds_part2 = 'HC';     break;
    case 'mgdb2024_hq':  $ds_part0 = 'maizegdb2024'; $ds_part2 = 'HQ';     break;
    case 'mgdb2024_hc':  $ds_part0 = 'maizegdb2024'; $ds_part2 = 'HC';     break;
    case 'schnable2023': $ds_part0 = 'schnable2023'; $ds_part2 = 'impute'; break;
    // MaizeGDB GRIN-linked 2026 release (v1.4): Grzybowski et al. 2023 sites,
    // Beagle-imputed; sample columns are release sample_ids (ZmG_<genotype>).
    // Store: hdf5/grin2026/zmgrin2026_<chr>_impute.h5
    case 'zmgrin2026_imp': $ds_part0 = 'zmgrin2026'; $ds_part2 = 'impute'; break;
    case 'nam2021':      // new UI sends the bare id
    case 'nam2021_hq':   $ds_part0 = 'nam2021';      $ds_part2 = 'HQ';     break;
    case 'nam2021_hc':   $ds_part0 = 'nam2021';      $ds_part2 = 'HC';     break;
    default:
        echo json_encode(array('status' => 'error', 'message' => "Unknown dataset '$dataset'."));
        exit;
}

$db_filename = $VERSION_PATH . $ds_part0 . '_' . $chr . '_' . $ds_part2 . '.h5';
if (!is_file($db_filename)) {
    error_log("processForm.php: no store $db_filename");
    reply(array('status' => 'error',
        'message' => "No variant store for $dataset $chr on this server."), array('store' => $db_filename));
    exit;
}

/* ---------------------------------------------------------------------
 *  OUTPUT PATH — force it inside $VCF_DIR (no path traversal)
 * ------------------------------------------------------------------- */
if (!is_dir($VCF_DIR)) { @mkdir($VCF_DIR, 0775, true); }
if (!is_writable($VCF_DIR)) {
    error_log("processForm.php: VCF directory not writable: $VCF_DIR");
    reply(array('status' => 'error',
        'message' => 'The server cannot write the VCF (its output folder is not writable).'), array('vcfDir' => $VCF_DIR));
    exit;
}
prune_old_vcfs($VCF_DIR, $VCF_TTL_HOURS);

/* Output is gzip-compressed (.vcf.gz) — genotype text compresses ~15-20x.
 * h5_to_vcf.py writes gzip transparently when the path ends in .gz, and
 * data.js inflates the response with DecompressionStream.
 * The browser proposes a name (snpv_<time>_<random>_<start>_<end>.vcf.gz); it is used only
 * when it has that shape and names no existing file, so a request can never overwrite
 * another query's VCF. Otherwise the server names the file. */
$base = basename((string) $outName);
if (!preg_match('/^snpv_[A-Za-z0-9_]{1,120}\.vcf\.gz$/', $base) || file_exists(rtrim($VCF_DIR, '/') . '/' . $base)) {
    $base = 'snpv_' . time() . '_' . bin2hex(random_bytes(8)) . '.vcf.gz';
}
$vcf_path = rtrim($VCF_DIR, '/') . '/' . $base;

/* ---------------------------------------------------------------------
 *  GENOTYPES — pass the selected accession IDs through as a JSON array
 * ------------------------------------------------------------------- */
$genotypesArray = json_decode($genotypesJson);
if (!is_array($genotypesArray)) { $genotypesArray = array(); }
// column names only: strings, each once, in the order given
$genotypesArray = array_values(array_unique(array_filter($genotypesArray, 'is_string')));
if (count($genotypesArray) > $MAX_IDS) {
    echo json_encode(array('status' => 'error',
        'message' => 'Too many accessions in one request (' . count($genotypesArray) . '; at most ' . $MAX_IDS . ').'));
    exit;
}

/* The accession list is handed to Python through a sidecar JSON file, not
 * as a command-line argument. Windows' escapeshellarg() strips '"' (and
 * '%', '!') from arguments, which shreds a JSON array on the command line;
 * a file sidesteps shell quoting entirely and behaves the same on every
 * OS. h5_to_vcf.py's 5th arg is now that file's path. */
$acc_path = $vcf_path . '.acc.json';
file_put_contents($acc_path, json_encode(array_values($genotypesArray)));

/* ---------------------------------------------------------------------
 *  RUN  h5_to_vcf.py  <db> <out.vcf> <start> <end> <accessions.json>
 * ------------------------------------------------------------------- */
$command = escapeshellarg($PYTHON_PATH) . ' ' . escapeshellarg('h5_to_vcf.py') . ' '
         . escapeshellarg($db_filename) . ' '
         . escapeshellarg($vcf_path)    . ' '
         . escapeshellarg($start)       . ' '
         . escapeshellarg($end)         . ' '
         . escapeshellarg($acc_path)    . ' 2>&1';

$output = shell_exec($command);

@unlink($acc_path);

// The script writes the VCF as a side effect; success = the file now exists.
if (is_file($vcf_path)) {
    $variants = null;
    if ($output !== null && preg_match('/^variants:\s*(\d+)/m', $output, $mm)) {
        $variants = (int) $mm[1];
    }
    reply(array(
        'status'   => 'success',
        'outFile'  => $vcf_path,   // the web path the browser fetches (./vcf/...)
        'variants' => $variants,   // exact site count — lets the client size the result before parsing
        'message'  => 'VCF written',
    ), array('output' => $output));
} else if ($output !== null && preg_match('/^TOO_LARGE:\s*(.+)$/m', $output, $tl)) {
    // Over the server's build limit (h5_to_vcf.py checks variants x accessions first).
    reply(array(
        'status'  => 'error',
        'message' => 'This query is too large to build here: ' . trim($tl[1]),
        'tooLarge'=> true,
    ));
} else if ($output !== null && strpos($output, 'No data found in the specified position range') !== false) {
    // Python ran fine, the interval simply contained no variants.
    reply(array(
        'status'  => 'empty',
        'message' => 'No variants in this interval.',
    ), array('output' => $output));
} else {
    // Real failure (bad Python path, missing h5py/numpy, dataset key error, ...).
    $detail = ($output === null ? '(no output — check that $PYTHON_PATH is correct and executable)' : $output);
    error_log('processForm.php: no VCF produced. command: ' . $command . ' | output: ' . substr($detail, 0, 4000));
    reply(array(
        'status'  => 'error',
        'message' => 'No VCF produced: the variant extraction failed on the server.',
    ), array('command' => $command, 'output' => $detail));
}

/* Remove VCFs (and stray accession lists) older than $ttlHours from the output folder. Runs
 * at most once every 10 minutes, marked by the folder's .last_prune file. Only names this
 * endpoint writes are touched. */
function prune_old_vcfs($dir, $ttlHours) {
    if ($ttlHours <= 0) return;
    $dir = rtrim($dir, '/');
    $mark = $dir . '/.last_prune';
    if (is_file($mark) && filemtime($mark) > time() - 600) return;
    @touch($mark);
    $cutoff = time() - (int) round($ttlHours * 3600);
    foreach ((array) @scandir($dir) as $f) {
        if (!preg_match('/^snpv_[A-Za-z0-9_]+\.vcf(\.gz)?(\.acc\.json)?$/', $f)) continue;
        $p = $dir . '/' . $f;
        $m = @filemtime($p);
        if ($m !== false && $m < $cutoff) @unlink($p);
    }
}
?>
