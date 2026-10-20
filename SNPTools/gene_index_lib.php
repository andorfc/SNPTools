<?php
/* =====================================================================
 *  gene_index_lib.php — gene model -> coordinates, without loading the gene store.
 *
 *  gff/genes_data.serialized (39,756 B73 v5 genes) is the source. Unserializing it
 *  for every lookup cost ~10 ms and ~25 MB of PHP memory to answer for one gene.
 *  gff/genes_index.txt holds the same records sorted by id in fixed-width lines,
 *  so a lookup is a binary search: ~16 reads of one line each.
 *
 *    #snptools-gene-index v1 source_bytes=<n> records=<n> id=<w> chr=<w> num=<w>
 *    <id>\t<chromosome>\t<start>\t<end>\n   one line per gene, every field padded with
 *    spaces to its column width (ids and chromosomes on the right, numbers on the left, so
 *    no line ends in a space)
 *
 *  The index answers only while its source_bytes equals the size of the source; a
 *  stale or missing index is rebuilt when gff/ is writable (written to a temp file
 *  and renamed), and otherwise the source itself is read, as before.
 *  Rebuild by hand: php tools/build_gene_index.php
 * ===================================================================== */

const GENE_INDEX_MAGIC = '#snptools-gene-index v1';

function gene_index_paths($dir) {
    $dir = rtrim($dir, '/');
    return array($dir . '/genes_data.serialized', $dir . '/genes_index.txt');
}

/* Write the index for $source to $index. Returns the number of records, or false. */
function gene_index_build($source, $index) {
    $data = @unserialize(file_get_contents($source), array('allowed_classes' => false));
    if (!is_array($data)) return false;
    ksort($data, SORT_STRING);
    $wid = 1; $wchr = 1; $wnum = 1;
    foreach ($data as $id => $g) {
        $wid = max($wid, strlen($id));
        $wchr = max($wchr, strlen((string) $g['chromosome']));
        $wnum = max($wnum, strlen((string) $g['start']), strlen((string) $g['end']));
    }
    $out = GENE_INDEX_MAGIC . ' source_bytes=' . filesize($source) . ' records=' . count($data)
         . " id=$wid chr=$wchr num=$wnum\n";
    foreach ($data as $id => $g) {
        $out .= str_pad($id, $wid) . "\t" . str_pad((string) $g['chromosome'], $wchr) . "\t"
              . str_pad((string) $g['start'], $wnum, ' ', STR_PAD_LEFT) . "\t"
              . str_pad((string) $g['end'], $wnum, ' ', STR_PAD_LEFT) . "\n";
    }
    $tmp = $index . '.' . getmypid() . '.' . mt_rand() . '.tmp';
    if (@file_put_contents($tmp, $out) === false) return false;
    if (!@rename($tmp, $index)) { @unlink($tmp); return false; }
    return count($data);
}

/* Open the index if it matches the source: {fh, off, n, len, wid, wchr, wnum} or null. */
function gene_index_open($source, $index) {
    if (!is_file($index) || !is_file($source)) return null;
    $fh = @fopen($index, 'rb');
    if (!$fh) return null;
    $head = fgets($fh);
    if ($head === false || strpos($head, GENE_INDEX_MAGIC . ' ') !== 0
        || !preg_match('/source_bytes=(\d+) records=(\d+) id=(\d+) chr=(\d+) num=(\d+)$/', rtrim($head), $m)
        || (int) $m[1] !== filesize($source)) {
        fclose($fh); return null;
    }
    $ix = array('fh' => $fh, 'off' => strlen($head), 'n' => (int) $m[2],
                'wid' => (int) $m[3], 'wchr' => (int) $m[4], 'wnum' => (int) $m[5]);
    $ix['len'] = $ix['wid'] + $ix['wchr'] + 2 * $ix['wnum'] + 4;      // 3 tabs + newline
    if ($ix['off'] + $ix['n'] * $ix['len'] !== filesize($index)) { fclose($fh); return null; }
    return $ix;
}

/* The gene's record, shaped as in the source ({chromosome, start, end, ID}, all strings), or null. */
function gene_index_find($ix, $id) {
    $lo = 0; $hi = $ix['n'] - 1;
    while ($lo <= $hi) {
        $mid = ($lo + $hi) >> 1;
        fseek($ix['fh'], $ix['off'] + $mid * $ix['len']);
        $rec = fread($ix['fh'], $ix['len']);
        $key = rtrim(substr($rec, 0, $ix['wid']), ' ');
        $c = strcmp($key, $id);
        if ($c === 0) {
            $f = explode("\t", rtrim($rec, "\n"));
            return array('chromosome' => rtrim($f[1], ' '), 'start' => ltrim($f[2], ' '),
                         'end' => ltrim($f[3], ' '), 'ID' => $key);
        }
        if ($c < 0) $lo = $mid + 1; else $hi = $mid - 1;
    }
    return null;
}

/* Look a gene up: the index when it is current (rebuilt first if stale and gff/ is
   writable), else the serialized source. Returns the record or null. */
function gene_lookup($dir, $id) {
    list($source, $index) = gene_index_paths($dir);
    $ix = gene_index_open($source, $index);
    if (!$ix && is_file($source) && is_writable(dirname($index))
        && (!is_file($index) || is_writable($index)) && gene_index_build($source, $index)) {
        $ix = gene_index_open($source, $index);
    }
    if ($ix) {
        $g = gene_index_find($ix, $id);
        fclose($ix['fh']);
        return $g;
    }
    $data = @unserialize(@file_get_contents($source), array('allowed_classes' => false));
    return (is_array($data) && isset($data[$id])) ? $data[$id] : null;
}
