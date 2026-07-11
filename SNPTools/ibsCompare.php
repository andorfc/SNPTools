<?php
/* =====================================================================
 *  ibsCompare.php — serve the precomputed genome-wide IBS row for one
 *  focal accession from the dense matrices in ./distance/.
 *
 *  CALL:      ibsCompare.php?focal=<SNPVersity ID>
 *  RESPONSE:  { "focal":"<ID>", "rows":[ {"id":"<id>","similarity":0.9963,"missing":23.51}, ... ] }
 *
 *  Files in ./distance/ :
 *    maizegdb_allchr_final_similarity.csv   dense 2710x2710, headerless, symmetric, 0..1
 *    maizegdb_allchr_final_missing_pct.csv  dense 2710x2710, headerless, fraction 0..1
 *    ids.txt                                REQUIRED: 2710 accession IDs, ONE PER LINE,
 *                                           in the SAME ORDER as the matrix rows/columns.
 *
 *  Similarity is returned as-is (0..1); missing is returned as a PERCENT (fraction x 100)
 *  so it matches SNPCompare's local scope. Metadata (project/SRA/name) is joined in the
 *  browser from the accession catalog, so it is intentionally NOT included here.
 *
 *  A small byte-offset index (<csv>.offidx) is built once so each request seeks directly
 *  to the focal row instead of scanning the whole file.
 * ===================================================================== */

header('Content-Type: application/json');

$DIR = './distance/';
$SIM = $DIR . 'maizegdb_allchr_final_similarity.csv';
$MIS = $DIR . 'maizegdb_allchr_final_missing_pct.csv';
$IDS = $DIR . 'ids.txt';

$focal = isset($_GET['focal']) ? $_GET['focal'] : '';
if (!preg_match('/^[A-Za-z0-9_.-]+$/', $focal)) {
    echo json_encode(array('error' => 'Invalid focal id')); exit;
}
foreach (array($SIM, $MIS, $IDS) as $f) {
    if (!is_file($f)) { echo json_encode(array('error' => 'Missing file: ' . $f)); exit; }
}

$ids = file($IDS, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
$pos = array_search($focal, $ids, true);
if ($pos === false) {
    echo json_encode(array('focal' => $focal, 'rows' => array(),
        'error' => 'Focal id not found in distance/ids.txt')); exit;
}

function line_offsets($csv) {
    $idx = $csv . '.offidx';
    if (is_file($idx) && filemtime($idx) >= filemtime($csv)) {
        return unserialize(file_get_contents($idx));
    }
    $offs = array(); $fh = fopen($csv, 'rb'); $p = 0;
    while (($l = fgets($fh)) !== false) { $offs[] = $p; $p = ftell($fh); }
    fclose($fh);
    @file_put_contents($idx, serialize($offs));
    return $offs;
}
function read_row($csv, $i) {
    $offs = line_offsets($csv);
    if ($i >= count($offs)) return null;
    $fh = fopen($csv, 'rb'); fseek($fh, $offs[$i]); $line = fgets($fh); fclose($fh);
    return explode(',', rtrim($line, "\r\n"));
}

$sim = read_row($SIM, $pos);
$mis = read_row($MIS, $pos);
if ($sim === null || $mis === null) {
    echo json_encode(array('error' => 'Row not found (matrix / ids.txt length mismatch?)')); exit;
}

$n = min(count($ids), count($sim), count($mis));
$rows = array();
for ($j = 0; $j < $n; $j++) {
    $rows[] = array(
        'id'         => $ids[$j],
        'similarity' => (float)$sim[$j],
        'missing'    => ((float)$mis[$j]) * 100.0,
    );
}
echo json_encode(array('focal' => $focal, 'rows' => $rows));
?>
