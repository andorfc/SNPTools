<?php
/* lookupGeneModel.php — resolve a Fusarium gene-model ID to its coordinates.
 * The reference is chosen with ?database=  (graminearum | vert7600 | vertMRC826).
 *   graminearum  -> gff/genes_data.serialized            (FGSG_ gene models)
 *   vert7600     -> gff/genes_data_vert_7600.serialized  (FVEG_ gene models)
 *   vertMRC826   -> gff/genes_data_vert.serialized       (FVERT4_ gene models)
 */
$geneId   = isset($_GET['geneModelId']) ? $_GET['geneModelId'] : '';
$database = isset($_GET['database'])    ? $_GET['database']    : 'graminearum';

switch ($database) {
    case 'vert7600':
    case 'vertV':
        $file = './gff/genes_data_vert_7600.serialized'; break;
    case 'vertMRC826':
    case 'vert':
        $file = './gff/genes_data_vert.serialized'; break;
    case 'graminearum':
    case 'maizegdb':   // legacy alias
    default:
        $file = './gff/genes_data.serialized'; break;
}

$genesData = @unserialize(@file_get_contents($file));

/* Resolve a gene ID to a key in the store, tolerating zero-padding differences.
 * FungiDB writes IDs with 5-digit numbers (e.g. FVEG_03144) while these stores use
 * 6-digit numbers (FVEG_003144); FGSG uses 5-digit. We split PREFIX_ + number and
 * try the number re-padded to a few common widths (and unpadded) so a FungiDB-style
 * ID pasted from the browser resolves to the stored key. */
function resolve_gene_key($genesData, $geneId) {
    if (!is_array($genesData)) return null;
    if (isset($genesData[$geneId])) return $geneId;
    if (preg_match('/^(.*_)0*([0-9]+)$/', $geneId, $m)) {
        $prefix = $m[1];
        $num    = $m[2];
        $n      = (string) intval($num);
        $cands  = array();
        foreach (array(6, 5, 4, 3, strlen($num)) as $w) {
            $cands[$prefix . str_pad($n, $w, '0', STR_PAD_LEFT)] = true;
        }
        $cands[$prefix . $n] = true;                 // unpadded
        foreach (array_keys($cands) as $k) {
            if (isset($genesData[$k])) return $k;
        }
    }
    return null;
}

$key = resolve_gene_key($genesData, $geneId);
if ($key !== null) {
    $info = $genesData[$key];
    // Flag whether this gene sits on a chr1..chrN chromosome (queryable in the HDF5
    // variant store) or on an unplaced scaffold / mitochondrion (no variants there,
    // so gene->region autofill has nothing to load). The variant HDF5s use chr#
    // tokens only; a handful of genes remain on NW_* / contig_* scaffolds.
    $chrom = isset($info['chromosome']) ? $info['chromosome'] : '';
    $info['placed'] = (bool) preg_match('/^chr[0-9]+$/', $chrom);
    if (!$info['placed']) { $info['scaffold'] = $chrom; }
    if (!isset($info['id'])) { $info['id'] = $key; }
    $info['queryId'] = $key;                          // the actual key matched (post-normalization)
    echo json_encode($info);
} else {
    echo json_encode(array('chromosome' => 'chr1', 'start' => '0', 'end' => '0', 'id' => 'empty'));
}
?>
