<?php
/* php tools/build_gene_index.php [gff dir]  -- rebuild gff/genes_index.txt from
   gff/genes_data.serialized (see gene_index_lib.php). lookupGeneModel.php also rebuilds it
   by itself when the store changes and gff/ is writable; run this where it is not. */
require __DIR__ . '/../gene_index_lib.php';
$dir = isset($argv[1]) ? $argv[1] : __DIR__ . '/../gff';
list($source, $index) = gene_index_paths($dir);
$n = gene_index_build($source, $index);
if ($n === false) { fwrite(STDERR, "could not build $index from $source\n"); exit(1); }
echo "$index: $n genes, " . filesize($index) . " bytes\n";
