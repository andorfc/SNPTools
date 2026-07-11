<?php
/* =====================================================================
 *  processForm.php — region + accession list  ->  VCF (via h5_to_vcf.py)
 *
 *  Local testing (MAMP) vs production (Linux): the ONLY thing you should
 *  need to change between machines is $PYTHON_PATH below (or set the
 *  PYTHON_PATH environment variable). Everything else is relative.
 * ===================================================================== */

header('Content-Type: application/json');

/* ---------------------------------------------------------------------
 *  CONFIG — edit for your machine
 * ------------------------------------------------------------------- */
// 1) Python interpreter that has h5py + numpy installed.
//    MAMP example (miniconda):  /Users/<you>/opt/miniconda3/bin/python
//    Prefer setting the PYTHON_PATH env var; this is the fallback.
$PYTHON_PATH = getenv('PYTHON_PATH');
if (!$PYTHON_PATH) {
    $PYTHON_PATH = '/usr/bin/python3';           // <-- CHANGE ME for MAMP if needed
}

// 2) Where the .h5 files live, relative to this PHP file.
$VERSION_PATH = './hdf5/version3/';

// 3) Where VCFs are written (must be web-served AND writable). Matches CFG.vcfDir in data.js.
$VCF_DIR = './vcf/';

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

// numeric interval
if (!is_numeric($start) || !is_numeric($end)) {
    echo json_encode(array('status' => 'error', 'message' => 'Invalid interval (start/end must be numeric).'));
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
$ds_part0 = 'maizegdb2026';   // family
$ds_part2 = 'HQ';             // quality tier

switch ($dataset) {
    case 'mgdb2026_hq':  $ds_part0 = 'maizegdb2026'; $ds_part2 = 'HQ';     break;
    case 'mgdb2026_hc':  $ds_part0 = 'maizegdb2026'; $ds_part2 = 'HC';     break;
    case 'mgdb2024_hq':  $ds_part0 = 'maizegdb2024'; $ds_part2 = 'HQ';     break;
    case 'mgdb2024_hc':  $ds_part0 = 'maizegdb2024'; $ds_part2 = 'HC';     break;
    case 'schnable2023': $ds_part0 = 'schnable2023'; $ds_part2 = 'impute'; break;
    case 'nam2021':      // new UI sends the bare id
    case 'nam2021_hq':   $ds_part0 = 'nam2021';      $ds_part2 = 'HQ';     break;
    case 'nam2021_hc':   $ds_part0 = 'nam2021';      $ds_part2 = 'HC';     break;
    default:
        echo json_encode(array('status' => 'error', 'message' => "Unknown dataset '$dataset'."));
        exit;
}

$db_filename = $VERSION_PATH . $ds_part0 . '_' . $chr . '_' . $ds_part2 . '.h5';
if (!is_file($db_filename)) {
    echo json_encode(array('status' => 'error',
        'message' => 'HDF5 file not found: ' . $db_filename));
    exit;
}

/* ---------------------------------------------------------------------
 *  OUTPUT PATH — force it inside $VCF_DIR (no path traversal)
 * ------------------------------------------------------------------- */
if (!is_dir($VCF_DIR)) { @mkdir($VCF_DIR, 0775, true); }
if (!is_writable($VCF_DIR)) {
    echo json_encode(array('status' => 'error',
        'message' => 'VCF directory is not writable: ' . $VCF_DIR));
    exit;
}
$base = basename($outName ? $outName : ('snpv_' . time() . '_' . mt_rand() . '.vcf'));
if (substr($base, -4) !== '.vcf') { $base .= '.vcf'; }
$vcf_path = rtrim($VCF_DIR, '/') . '/' . $base;

/* ---------------------------------------------------------------------
 *  GENOTYPES — pass the selected accession IDs through as a JSON array
 * ------------------------------------------------------------------- */
$genotypesArray = json_decode($genotypesJson);
if (!is_array($genotypesArray)) { $genotypesArray = array(); }
$jsonArray = escapeshellarg(json_encode(array_values($genotypesArray)));

/* ---------------------------------------------------------------------
 *  RUN  h5_to_vcf.py  <db> <out.vcf> <start> <end> <genotypesJson>
 * ------------------------------------------------------------------- */
$command = escapeshellarg($PYTHON_PATH) . ' ' . escapeshellarg('h5_to_vcf.py') . ' '
         . escapeshellarg($db_filename) . ' '
         . escapeshellarg($vcf_path)    . ' '
         . escapeshellarg($start)       . ' '
         . escapeshellarg($end)         . ' '
         . $jsonArray . ' 2>&1';

$output = shell_exec($command);

// The script writes the VCF as a side effect; success = the file now exists.
if (is_file($vcf_path)) {
    echo json_encode(array(
        'status'  => 'success',
        'outFile' => $vcf_path,
        'message' => 'VCF written',
        'output'  => $output,
    ));
} else if ($output !== null && strpos($output, 'No data found in the specified position range') !== false) {
    // Python ran fine, the interval simply contained no variants.
    echo json_encode(array(
        'status'  => 'empty',
        'message' => 'No variants in this interval.',
        'output'  => $output,
    ));
} else {
    // Real failure (bad Python path, missing h5py/numpy, dataset key error, ...).
    echo json_encode(array(
        'status'  => 'error',
        'message' => 'No VCF produced (script error). See output.',
        'command' => $command,
        'output'  => ($output === null ? '(no output — check that $PYTHON_PATH is correct and executable)' : $output),
    ));
}
?>
