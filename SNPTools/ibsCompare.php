<?php
/* =====================================================================
 *  ibsCompare.php — serve the precomputed genome-wide IBS row for one
 *  focal accession from the dense matrices in ./distance/.
 *
 *  CALL:      ibsCompare.php?focal=<SNPVersity ID>[&dataset=<family>][&sites=snp]
 *  RESPONSE:  { "focal":"<ID>", "rows":[ {"id":"<id>","similarity":0.9963,"missing":23.51}, ... ] }
 *  Also:      ibsCompare.php?probe=1&dataset=<family>  -> {"dataset","available","n","sites":[..],"trees":[..]}
 *             ibsCompare.php?tree=nj|upgma&dataset=<family> -> the precomputed genome-wide Newick
 *
 *  dataset = SNPTools dataset FAMILY (Data.familyOf); SNPCompare and SNPTree always send it.
 *  Without it the release's family is assumed: zmgrin2026 (MaizeGDB GRIN-linked 2026).
 *  mgdb2026 (unchanged layout), files in ./distance/ :
 *    maizegdb_allchr_final_similarity.csv   dense 2710x2710, headerless, symmetric, 0..1
 *    maizegdb_allchr_final_missing_pct.csv  dense 2710x2710, headerless, fraction 0..1
 *    ids.txt                                REQUIRED: 2710 accession IDs, ONE PER LINE,
 *                                           in the SAME ORDER as the matrix rows/columns.
 *  any other family, files in ./distance/<family>/ (same formats):
 *    similarity.csv, missing_pct.csv, ids.txt          all variant sites (sites=all, default)
 *    similarity_snp.csv, missing_pct_snp.csv           optional, SNPs only (sites=snp)
 *    tree_nj.nwk, tree_upgma.nwk                       optional, genome-wide trees (tree=...)
 *  (other files there, e.g. allele_distance.csv / co_called_sites.csv, are not served.)
 *
 *  Similarity is returned as-is (0..1); missing is returned as a PERCENT (fraction x 100)
 *  so it matches SNPCompare's local scope. Metadata (project/SRA/name) is joined in the
 *  browser from the accession catalog, so it is intentionally NOT included here.
 *
 *  A small byte-offset index (<csv>.offidx) is built once so each request seeks directly
 *  to the focal row instead of scanning the whole file.
 *  The distance root can be moved with the SNPTOOLS_DISTANCE_DIR environment variable
 *  (used by the headless tests); the default is ./distance/.
 * ===================================================================== */

header('Content-Type: application/json');

$ROOT = getenv('SNPTOOLS_DISTANCE_DIR');
$ROOT = $ROOT ? rtrim($ROOT, '/') . '/' : './distance/';
$family = isset($_GET['dataset']) && $_GET['dataset'] !== '' ? (string) $_GET['dataset'] : 'zmgrin2026';
if (!preg_match('/^[a-z0-9_]+$/', $family)) {
    echo json_encode(array('error' => 'Invalid dataset')); exit;
}
$sites = (isset($_GET['sites']) && $_GET['sites'] === 'snp') ? 'snp' : 'all';
if ($family === 'mgdb2026') {
    $DIR = $ROOT;
    $SIM = $DIR . 'maizegdb_allchr_final_similarity.csv';
    $MIS = $DIR . 'maizegdb_allchr_final_missing_pct.csv';
    $VARIANTS = array('all' => array($SIM, $MIS));
    $TREES = array();
} else {
    $DIR = $ROOT . $family . '/';
    $VARIANTS = array('all' => array($DIR . 'similarity.csv', $DIR . 'missing_pct.csv'),
                      'snp' => array($DIR . 'similarity_snp.csv', $DIR . 'missing_pct_snp.csv'));
    $TREES = array('nj' => $DIR . 'tree_nj.nwk', 'upgma' => $DIR . 'tree_upgma.nwk');
    list($SIM, $MIS) = $VARIANTS[$sites];
}
$IDS = $DIR . 'ids.txt';

if (isset($_GET['probe'])) {
    $have = array();
    foreach ($VARIANTS as $k => $pair) if (is_file($pair[0]) && is_file($pair[1])) $have[] = $k;
    $trees = array();
    foreach ($TREES as $k => $f) if (is_file($f)) $trees[] = $k;
    $ok = is_file($IDS) && in_array('all', $have, true);
    $n = $ok ? count(file($IDS, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES)) : 0;
    echo json_encode(array('dataset' => $family, 'available' => $ok, 'n' => $n, 'sites' => $have, 'trees' => $trees));
    exit;
}
if (isset($_GET['tree'])) {
    $t = $_GET['tree'];
    if (!isset($TREES[$t]) || !is_file($TREES[$t])) {
        http_response_code(404); echo json_encode(array('error' => 'No precomputed tree ' . $t . ' for ' . $family)); exit;
    }
    header('Content-Type: text/plain; charset=utf-8');
    header('Content-Disposition: attachment; filename="' . $family . '_genomewide_' . $t . '.nwk"');
    readfile($TREES[$t]); exit;
}

$focal = isset($_GET['focal']) ? $_GET['focal'] : '';
if (!preg_match('/^[A-Za-z0-9_.-]+$/', $focal)) {
    echo json_encode(array('error' => 'Invalid focal id')); exit;
}
foreach (array($SIM, $MIS, $IDS) as $f) {
    if (!is_file($f)) {   // the path goes to the server log, not to the browser
        error_log("ibsCompare.php: missing $f");
        echo json_encode(array('error' => "Genome-wide matrices for $family are not installed on this server.")); exit;
    }
}

$ids = file($IDS, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
$pos = array_search($focal, $ids, true);
if ($pos === false) {
    echo json_encode(array('focal' => $focal, 'rows' => array(),
        'error' => "Focal id not found in the $family matrices.")); exit;
}

/* The index is written to a temporary file and renamed into place, so a request reading
   it while another builds it sees the old index or the whole new one, never half of one;
   an index that does not read back as an array is ignored and rebuilt. */
function line_offsets($csv) {
    $idx = $csv . '.offidx';
    if (is_file($idx) && filemtime($idx) >= filemtime($csv)) {
        $offs = @unserialize(file_get_contents($idx), array('allowed_classes' => false));
        if (is_array($offs)) return $offs;
    }
    $offs = array(); $fh = fopen($csv, 'rb'); $p = 0;
    while (($l = fgets($fh)) !== false) { $offs[] = $p; $p = ftell($fh); }
    fclose($fh);
    $tmp = $idx . '.' . getmypid() . '.' . mt_rand() . '.tmp';
    if (@file_put_contents($tmp, serialize($offs)) !== false) {
        if (!@rename($tmp, $idx)) @unlink($tmp);
    }
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
