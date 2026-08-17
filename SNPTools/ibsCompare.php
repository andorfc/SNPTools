<?php
/* =====================================================================
 *  ibsCompare.php — serve the precomputed genome-wide IBS row for one
 *  focal accession from the dense matrices in ./distance/.
 *
 *  Fusarium build: three reference genomes, one distance folder each.
 *
 *  CALL:      ibsCompare.php?focal=<SNPVersity ID>&ds=<dataset id>
 *  RESPONSE:  { "focal":"<ID>", "ref":"Fgram_ph1",
 *               "rows":[ {"id":"<id>","similarity":0.9963,"missing":23.51}, ... ] }
 *
 *  Reference selection (from the dataset id, mirrors h5_to_vcf.py):
 *      ds contains "vert7600" -> Fvert_7600
 *      ds contains "vertmrc"  -> Fvert_mrc
 *      otherwise              -> Fgram_ph1   (graminearum PH-1)
 *
 *  Files in ./distance/<REF>/ :
 *    <REF>_allchr_similarity.csv    dense NxN, headerless, symmetric, 0..1
 *    <REF>_allchr_missing_pct.csv   dense NxN, headerless, fraction 0..1
 *    <REF>_ids.txt                  N accession IDs, ONE PER LINE, in the
 *                                   SAME ORDER as the matrix rows/columns.
 *
 *  Similarity is returned as-is (0..1); missing is returned as a PERCENT
 *  (fraction x 100) so it matches SNPCompare's local scope. Metadata
 *  (population / species / host / country / substrate / …) is joined in the
 *  browser from the accession catalog, so it is NOT included here.
 *
 *  A small byte-offset index (<csv>.offidx) is built once (best-effort) so
 *  each request seeks directly to the focal row instead of scanning.
 * ===================================================================== */

header('Content-Type: application/json');

$DIR = './distance/';

$focal = isset($_GET['focal']) ? $_GET['focal'] : '';
if (!preg_match('/^[A-Za-z0-9_.-]+$/', $focal)) {
    echo json_encode(array('error' => 'Invalid focal id')); exit;
}

/* dataset -> reference folder */
$ds = isset($_GET['ds']) ? strtolower($_GET['ds']) : '';
if (strpos($ds, 'vert7600') !== false)      { $REF = 'Fvert_7600'; }
elseif (strpos($ds, 'vertmrc') !== false)   { $REF = 'Fvert_mrc'; }
else                                        { $REF = 'Fgram_ph1'; }

$RD  = $DIR . $REF . '/';
$SIM = $RD . $REF . '_allchr_similarity.csv';
$MIS = $RD . $REF . '_allchr_missing_pct.csv';
$IDS = $RD . $REF . '_ids.txt';

foreach (array($SIM, $MIS, $IDS) as $f) {
    if (!is_file($f)) { echo json_encode(array('error' => 'Missing file: ' . $f)); exit; }
}

$ids = file($IDS, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
$ids = array_map('trim', $ids);
$pos = array_search($focal, $ids, true);
if ($pos === false) {
    echo json_encode(array('focal' => $focal, 'ref' => $REF, 'rows' => array(),
        'error' => 'Focal id not found in ' . $REF . '_ids.txt')); exit;
}

function line_offsets($csv) {
    $idx = $csv . '.offidx';
    if (is_file($idx) && filemtime($idx) >= filemtime($csv)) {
        $cached = @unserialize(file_get_contents($idx));
        if (is_array($cached)) return $cached;
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
echo json_encode(array('focal' => $focal, 'ref' => $REF, 'rows' => $rows));
?>
