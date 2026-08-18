<?php
/* =====================================================================
 *  processForm.php — region + accession list -> VCF (via h5_to_vcf.py)
 *  FUSARIUM build.  Three reference genomes; haploid genotypes are handled
 *  inside h5_to_vcf.py.  The Python interpreter is resolved from PYTHON_PATH
 *  (set per host) falling back to 'python3' on PATH.
 * ===================================================================== */

header('Content-Type: application/json');

$PYTHON_PATH = getenv('PYTHON_PATH');
if (!$PYTHON_PATH) { $PYTHON_PATH = 'python3'; }

// Where the Fusarium .h5 files live, relative to this PHP file.
$VERSION_PATH = './hdf5/fusarium/';
// Where VCFs are written (web-served AND writable). Matches CFG.vcfDir in data.js.
$VCF_DIR = './vcf/';

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

if (!is_numeric($start) || !is_numeric($end)) {
    echo json_encode(array('status' => 'error', 'message' => 'Invalid interval (start/end must be numeric).'));
    exit;
}
// chromosome token: chr1 .. chr12 (used verbatim in the .h5 filename).
// graminearum has chr1–4, verticillioides 7600 chr1–11, MRC826 chr1–12.
if (!preg_match('/^chr(1[0-2]|[1-9])$/', $chr)) {
    echo json_encode(array('status' => 'error', 'message' => "Invalid chromosome '$chr'."));
    exit;
}

/* DATASET -> (genome key, quality) -> <genome>_<chr>_<quality>.h5
 * All three references are the 2026 language-model rebuild (FunDLM/EVO2 +
 * ESM1/2/3/ESM-C). The FULL per-chromosome store is preferred; if it is not
 * present the 500kb test subset (_test500kb) is used instead (see resolution
 * below). Genome-key file prefixes: graminearum = fusarium2026, verticillioides
 * 7600 = Fvert_7600, verticillioides MRC826 = Fvert_mrc (match hdf5/fusarium/). */
$ds_part0   = 'fusarium2026';   // genome-key file prefix
$ds_part2   = 'HQ';             // quality tier

switch ($dataset) {
    case 'gram_hq':     $ds_part0 = 'fusarium2026'; $ds_part2 = 'HQ'; break;
    case 'gram_hc':     $ds_part0 = 'fusarium2026'; $ds_part2 = 'HC'; break;
    case 'vert7600_hq': $ds_part0 = 'Fvert_7600';   $ds_part2 = 'HQ'; break;
    case 'vert7600_hc': $ds_part0 = 'Fvert_7600';   $ds_part2 = 'HC'; break;
    case 'vertmrc_hq':  $ds_part0 = 'Fvert_mrc';    $ds_part2 = 'HQ'; break;
    case 'vertmrc_hc':  $ds_part0 = 'Fvert_mrc';    $ds_part2 = 'HC'; break;
    default:
        echo json_encode(array('status' => 'error', 'message' => "Unknown dataset '$dataset'."));
        exit;
}

$db_base = $VERSION_PATH . $ds_part0 . '_' . $chr . '_' . $ds_part2;
/* Resolution order: prefer the FULL store, fall back to the 500kb test subset.
 * Within each, prefer a decompressed .h5, else the gzip (.h5.gz) which
 * h5_to_vcf.py expands on the fly. First existing candidate wins. */
$db_candidates = array(
    $db_base . '.h5',              // full, decompressed
    $db_base . '.h5.gz',           // full, gzipped
    $db_base . '_test500kb.h5',    // test subset, decompressed
    $db_base . '_test500kb.h5.gz', // test subset, gzipped
);
$db_filename = null;
foreach ($db_candidates as $cand) {
    if (is_file($cand)) { $db_filename = $cand; break; }
}
if ($db_filename === null) {
    echo json_encode(array('status' => 'error',
        'message' => 'HDF5 file not found for ' . $db_base .
                     ' (.h5 / .h5.gz, full or _test500kb) — this reference/chromosome/quality ' .
                     'combination has not been generated yet.'));
    exit;
}

/* OUTPUT PATH — forced inside $VCF_DIR (no path traversal) */
if (!is_dir($VCF_DIR)) { @mkdir($VCF_DIR, 0775, true); }
if (!is_writable($VCF_DIR)) {
    echo json_encode(array('status' => 'error', 'message' => 'VCF directory is not writable: ' . $VCF_DIR));
    exit;
}
$base = basename($outName ? $outName : ('snpv_' . time() . '_' . mt_rand() . '.vcf'));
if (substr($base, -4) !== '.vcf') { $base .= '.vcf'; }
$vcf_path = rtrim($VCF_DIR, '/') . '/' . $base;

/* GENOTYPES — selected accession IDs passed through as a JSON array */
$genotypesArray = json_decode($genotypesJson);
if (!is_array($genotypesArray)) { $genotypesArray = array(); }
$jsonArray = escapeshellarg(json_encode(array_values($genotypesArray)));

/* RUN  h5_to_vcf.py <db> <out.vcf> <start> <end> <genotypesJson> */
$command = escapeshellarg($PYTHON_PATH) . ' ' . escapeshellarg('h5_to_vcf.py') . ' '
         . escapeshellarg($db_filename) . ' '
         . escapeshellarg($vcf_path)    . ' '
         . escapeshellarg($start)       . ' '
         . escapeshellarg($end)         . ' '
         . $jsonArray . ' 2>&1';

$output = shell_exec($command);

if (is_file($vcf_path)) {
    echo json_encode(array('status' => 'success', 'outFile' => $vcf_path,
        'message' => 'VCF written', 'output' => $output));
} else if ($output !== null && strpos($output, 'No data found in the specified position range') !== false) {
    echo json_encode(array('status' => 'empty', 'message' => 'No variants in this interval.', 'output' => $output));
} else {
    echo json_encode(array('status' => 'error',
        'message' => 'No VCF produced (script error). See output.',
        'command' => $command,
        'output'  => ($output === null ? '(no output — check that $PYTHON_PATH is correct and executable)' : $output)));
}
?>
