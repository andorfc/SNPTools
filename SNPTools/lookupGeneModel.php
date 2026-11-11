<?php
/* lookupGeneModel.php?geneModelId=<B73 v5 gene model> -> {"chromosome","start","end","ID"}
   (all strings), or {"chromosome":"chr1","start":"0","end":"0","id":"empty"} when the id is
   unknown. The lookup is a binary search of gff/genes_index.txt (gene_index_lib.php); the
   index is rebuilt from gff/genes_data.serialized when the store changes. */
require __DIR__ . '/gene_index_lib.php';
header('Content-Type: application/json');

$geneId = isset($_GET['geneModelId']) ? (string) $_GET['geneModelId'] : '';
$geneInfo = ($geneId !== '') ? gene_lookup(__DIR__ . '/gff', $geneId) : null;

if ($geneInfo !== null) {
    echo json_encode($geneInfo);
} else {
    $data = [
        'chromosome' => 'chr1',
        'start' => '0',
        'end' => '0',
        'id' => 'empty'
    ];
    echo json_encode($data);
}
