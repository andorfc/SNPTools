/* =====================================================================
 *  snphelp.js — Help & FAQ page for the SNPMaize suite (the SNPTools code base).
 *
 *  Registers itself as the 'help' tool so the router can reach it, and
 *  exposes openHelp() for the MaizeGDB masthead link (which doesn't want
 *  to depend on TOOLS[] being populated). Explains every page in the
 *  suite and every dataset defined in data.js, plus a scores/annotations
 *  glossary and a short FAQ.
 *
 *  Depends on core.js (SNPTools, S, go, ICONS, renderNav) and — for the
 *  live dataset table — data.js (Data.datasets). Degrades gracefully if
 *  Data isn't present.
 * ===================================================================== */
const SNPHelp = (function () {

  /* ---- tool descriptions, in user language, kept in step with core.js ---- */
  const LIVE = 'live', SOON = 'in development';
  const PAGES = [
    { id:'snpversity', name:'SNPVersity', icon:'dna', color:'#2563eb', status:LIVE,
      tag:'Build a variant view across accessions',
      what:'The front door of the suite and the starting point for most work. Choose a dataset, type a genomic interval (or a B73 v5 gene model ID), and pick the accessions you want: a quick pick (a panel, a random sample), a pasted list (IDs, line names or GRIN accessions), or Browse & filter, which opens SNPTrait over the page to filter by panel, origin, subpopulation and GRIN traits and applies the result. SNPVersity queries the variant store and returns a color-coded genotype table plus a downloadable VCF — allele states, predicted effects, and DNA/protein language-model scores included. Each row carries its Site QC class, with a Site QC filter above the table; a curated allele (see SNPCurate) has its mark beside the position and its site and annotation cells shaded light gold, and the curated alleles of the interval are named above the table. The gene models of the region link to their MaizeGDB gene pages, and a line whose sample is heterozygous is marked het.',
      give:'A dataset, a region or gene, and a set of accessions.',
      get:'A genotype table and a VCF. From here, "Send selection to…" hands the same result to any other tool; SNPImpact, SNPCompare and SNPTree also offer it on their empty page (Load data from SNPVersity).' },
    { id:'snpgwas', name:'GWAS Explorer', icon:'gwas', color:'#cf8a12', status:LIVE,
      tag:'Explore curated GWAS results to identify trait-associated SNPs',
      what:'Interactive Manhattan plots of curated published GWAS results from accessions included in SNPVersity. Pick a trait or trait type (Plant Architecture, Yield Component, Flowering Time), publication, or population to explore. Pan and zoom the canvas plot across all 10 chromosomes, jump to a chromosome or search by SNP ID or coordinate range, and toggle between the study’s published significance threshold and your own custom threshold. Drag-select a peak to open a sortable table of its SNPs, filterable to significant-only, and hand the region and/or the accessions off to SNPVersity. Download all or selected results with P<0.001.',
      give:'A trait, population, and publication result to view, or arrive with none selected and browse available results.',
      get:'An interactive Manhattan plot, results downloads, a sortable per-region SNP table, and a hand-off of the selected region and/or NAM accessions to SNPVersity.' },
    { id:'snptrait', name:'SNPTrait', icon:'leaf', color:'#1f8a4c', status:LIVE,
      tag:'Select lines by passport and trait metadata',
      what:'A line selector for datasets with germplasm metadata — on maize, the 933 lines of the GRIN-linked 2026 release (v1.4; all 26 NAM parents) (NAM, Ames 282, WiDiv and other GRIN accessions). Lines are grouped by primary panel and can be filtered with facets built from GRIN passport and evaluation data: panel membership, Grzybowski et al. (2023) subpopulation, country and state of origin, improvement status and era, kernel type and colour, cob colour, Stewart’s wilt score, sample heterozygosity (Inbred, Elevated heterozygosity, Heterozygous sample), and tertile bins of 1000-kernel weight, plant height, GDU to silk and ear length. Numeric GRIN traits (per-accession means) can also be filtered by range and shown as table columns. Datasets without metadata show their project grouping only. Tick several values in one filter section to keep lines with any of them (OR); different sections combine (AND), and an unticked value shows +n, the lines it would add.',
      give:'A dataset, then facet choices, a text search (name, GRIN accession, pedigree, state, subpopulation) and optional trait ranges.',
      get:'A selection you can export (CSV/JSON, including the trait columns and each line’s sample heterozygosity) or send to SNPVersity, choosing to replace its current accession list or add to it.' },
    { id:'snpimpact', name:'SNPImpact', icon:'star', color:'#7c3aed', status:LIVE,
      tag:'Rank candidate variants',
      what:'Prioritizes the variants in a region regardless of which accessions you picked. It orders them by an AI-based score (PlantCAD DNA-model and ESM protein-model predictions) combined with predicted consequence and Pfam domain annotation, so likely causal changes rise to the top. Filter by consequence, priority, score, domain, or Site QC, and flag a shortlist. Variants at flagged sites (het only, het excess) and sites with no carrier are hidden by default; a note above the table counts them, and the Site QC filter shows them. A curated allele is never hidden and carries its SNPCurate mark.',
      give:'A region sent from SNPVersity, or loaded on the empty page with Load data from SNPVersity (accessions are ignored here).',
      get:'A ranked, filterable variant table and a shortlist of candidate alleles.' },
    { id:'snpfunction', name:'SNPFunction', icon:'func', color:'#2563eb', status:LIVE,
      tag:'Gene function & allele mining',
      what:'A gene-scoped dossier, independent of any one region. It opens with a functional-annotation record — gene symbol and names, a sourced description, protein domain architecture drawn to scale, and external cross-references — and MaizeGDB\u2019s own Gene Ontology and Pathways views, read live from the MaizeGDB API: GO terms placed in the ontology with their evidence, and the gene\u2019s metabolic pathways step by step beside its KEGG maps. Below that it computes the gene\u2019s variant burden across the whole panel and lists which accessions carry each damaging or knockout allele. The burden, the knockout count and the allele list use usable sites; alleles at flagged sites and alleles with no carrier in this release sit in two collapsed groups under the list, and a banner warns when most of a gene’s protein-changing sites are flagged. The gene’s curated alleles head the allele card, and each allele sends its homozygous carriers, heterozygous carriers, all carriers or reference lines to SNPVersity.',
      give:'A gene model ID and a dataset.',
      get:'A functional-annotation dossier (identity, domains, cross-references), MaizeGDB\u2019s GO and pathway views, a variant-burden breakdown, and a damaging-allele catalog with carrier and reference lines.' },
    { id:'snpcompare', name:'SNPCompare', icon:'compare', color:'#0e7490', status:LIVE,
      tag:'Similarity to a focal accession',
      what:'Ranks every accession by identity-by-state similarity to one focal accession. Similarity is the fraction of co-called sites with the same genotype. Works genome-wide (precomputed), for the current region (computed in-browser from a SNPVersity result), or both side-by-side with a delta that surfaces region-specific relatedness such as introgression. The region scope uses usable sites by default (Region sites: Usable only or All sites); the genome-wide values are precomputed over all sites.',
      give:'A focal accession, and optionally a region handed off from SNPVersity or loaded with Load data from SNPVersity.',
      get:'A ranked similarity table with project, SRA ID, accession name, and PI number (the line’s GRIN accession).' },
    { id:'snptree', name:'SNPTree', icon:'tree', color:'#15803d', status:LIVE,
      tag:'Local phylogeny',
      what:'Builds a local phylogenetic tree from the genotype matrix already in memory, using identity-by-state distances (UPGMA / neighbour-joining). Useful for reading haplotype structure, introgression, and how lines cluster in a region. Distances use usable sites by default (Sites: Usable only or All).',
      give:'A SNPVersity result set, sent from SNPVersity or loaded with Load data from SNPVersity (warns first if the region × accessions would be slow to compute).',
      get:'An interactive tree, downloadable as Newick, MEGA, or PHYLIP.' },
    { id:'snpmatrix', name:'SNPMatrix', icon:'grid', color:'#b45309', status:LIVE,
      tag:'Pairwise distance matrix',
      what:'Computes the pairwise identity-by-state distance among your selected accessions and draws it as a heatmap. Reorder by input order or by clustering, switch between IBS distance and % identity, and color rows by bioproject. Like SNPTree, it uses usable sites by default (Sites: Usable only or All sites).',
      give:'A SNPVersity result set (warns first above ~450 accessions, or if the region × accessions would be slow to compute).',
      get:'A heatmap plus downloads: CSV distance matrix, PHYLIP, PNG, and SVG.' },
    { id:'snpimpute', name:'SNPImpute', icon:'impute', color:'#0891b2', status:SOON,
      tag:'Impute sequence & function',
      what:'Pan-genome\u2013guided imputation that fills missing genotypes and carries functional predictions through, across light and deep sequencing depths.',
      give:'A genotype set with missing calls.',
      get:'Imputed genotypes with function predictions.' },
    { id:'snpcurate', name:'SNPCurate', icon:'check', color:'#a16207', status:LIVE,
      tag:'Published alleles and what this release shows of them',
      what:'A registry of published maize alleles, from causal changes validated by molecular work to published markers and sites that tag a haplotype. Each allele the release genotypes shows its carriers (by subpopulation and country of origin), Site QC class, language-model scores and priority, its rank among the gene\u2019s missense variants, and the GRIN trait means of the lines that carry it. An allele the release cannot show (not a site, a structural variant, unreliable calls) is listed with the reason. Marks: gold for a validated causal change genotyped here, outline for any other allele genotyped here, grey for one the release cannot show; the same marks appear on these alleles in SNPVersity, SNPImpact, SNPFunction, SNPFold and SNPGeo. The table opens sorted by mark (gold, outline, grey), then by gene name.',
      give:'Nothing: browse the registry, or arrive from a mark in another tool.',
      get:'A record per allele, with hand-offs of its carriers to SNPVersity, its gene to SNPFunction and the site to SNPGeo.' },
    { id:'snpfold', name:'SNPFold', icon:'fold', color:'#be185d', status:LIVE,
      tag:'Variants on protein structure',
      what:'Maps coding variants onto predicted protein structure. A linear protein browser aligns variants with Pfam domains, secondary structure, and per-residue pLDDT confidence; an on-demand 3D viewer shows the fold colored by confidence, domain, or impact; and a per-variant readout interprets each change (domain, local confidence, secondary structure, predicted ΔΔG). Usable alleles only (on by default) hides coding variants at flagged or no-carrier sites from the track, the table and the 3D view. A curated allele is ringed on the track and marked in the table, and each variant sends its carriers or reference lines to SNPVersity.',
      give:'A gene with an available structure model.',
      get:'A structure-aware, per-variant interpretation of coding changes.' },
    { id:'snpgeo', name:'SNPGeo', icon:'map', color:'#059669', status:LIVE,
      tag:'Map where each allele is found',
      what:'Shows the geographic distribution of every variant in a gene (search a B73 v5 gene model) or in a region sent from SNPVersity, using each line\u2019s GRIN country and state of origin. Countries — or, in the North America view, U.S. states, Canadian provinces and Mexican states — are coloured by one of three statistics: the carrier fraction over all lines of the dataset known from that place, carriers among the lines genotyped there (the default: a reference \u2194 alternative scale, blue to orange), or the alternative-allele frequency over called alleles (heterozygotes count one allele). Click a country for its per-state table and the list of carrier lines; step through variants with the arrow keys; export the map as PNG or SVG. Each variant shows its Site QC class (SNPGeo hides nothing) and a curated allele its SNPCurate mark.',
      give:'A gene model ID and a dataset, or a SNPVersity result (SNPGeo re-queries the region for all lines of the dataset when the result covered only part of it).',
      get:'A choropleth per variant, per-country and per-state counts with allele frequencies, the carrier lines, and a variant table with the DNA (PlantCAD1, PlantCAD2, Evo2) and protein (ESM1, ESM2, ESM3, ESM-C) language-model scores.' },
    { id:'paneffect', name:'PanEffect', icon:'effect', color:'#b45309', status:LIVE,
      tag:'Missense effects across the pan-genome',
      what:'Visualizes the predicted effect of every possible amino-acid substitution across a protein using ESM protein-language-model scores. A B73 reference view and a pan-genome view each show a full-length substitution heatmap with a zoomable window, aligned to Pfam domains and predicted secondary structure; the pan-genome view adds the natural variation seen across the maize assemblies, colored by heterotic group. Choose a gene model and an ESM model, or arrive from SNPVersity, SNPFold, SNPFunction or SNPImpact on a specific missense call and that substitution is highlighted.',
      give:'A gene model and an ESM model (ESM1 / ESM2 / ESM3); optionally a missense variant handed off from another tool.',
      get:'B73 and pan-genome substitution heatmaps with domain and secondary-structure context, plus a downloadable per-variant effects file.' },
    { id:'snpgermplasm', name:'SNPGermplasm', icon:'germ', color:'#16a34a', status:SOON,
      tag:'Collection management',
      what:'Applies genotype-driven analytics to germplasm management, identifying redundancy, uniqueness, and priority materials for curation and deployment.',
      give:'A collection of accessions.',
      get:'Redundancy, uniqueness, and priority flags.' },
  ];

  /* ---- glossary ---- */
  const GLOSSARY = [
    ['B73 RefGen v5', 'Every coordinate, gene model, and annotation in the suite is anchored to the B73 version 5 maize reference genome.'],
    ['Genotype dosage', 'Each call is read as 0 (0/0, reference), 1 (heterozygous), or 2 (1/1, alternate homozygous). Missing calls (./.) are left out of a comparison rather than counted as a match.'],
    ['IBS distance / % identity', 'Identity-by-state compares two accessions site by site over the calls they share. Distance is the mean allele difference; % identity is 100 − distance. Used by SNPTree, SNPMatrix, and SNPCompare.'],
    ['Predicted effect & impact', 'Each variant carries a predicted consequence (missense, LOF, splice, indel, synonymous, …) rolled up into an impact tier: HIGH, MODERATE, LOW, or MODIFIER, most severe wins when a site lists several.'],
    ['PlantCAD score', 'A DNA language-model prediction of how disruptive a change is. The GRIN-linked 2026 dataset carries PlantCAD1 and the second-generation PlantCAD2 for every SNP, and Evo2 for SNPs within 1 kb of a gene.'],
    ['ESM score', 'A protein language-model prediction of the effect of an amino-acid change. The GRIN-linked 2026 dataset carries ESM1, ESM2, ESM3 and ESM-C for missense sites.'],
    ['MAF', 'Minor allele frequency, the frequency of the less common allele. SNPVersity can filter a region by a minimum MAF.'],
    ['Sample heterozygosity', 'The share of clean sites where a line is heterozygous. A clean site passes Site QC and has at least 3 homozygous carriers (16,463,157 sites). Inbred lines are near 0.8% (the median is 0.78%). A line above 5% is a Heterozygous sample (20 lines, marked het in SNPVersity and SNPFunction); above 2% it has Elevated heterozygosity (26 lines); the other 887 are Inbred. SNPTrait filters by it.'],
    ['Site QC', 'Each site of the GRIN-linked 2026 release is classed from its heterozygous (het) and homozygous (hom) carriers among all 933 lines, first match wins: No carrier (no line carries the allele), Het only (no line is homozygous), Het excess (more het than hom carriers), Het elevated (het carriers are a quarter or more of all carriers), Pass (the rest). The lines are inbreds, so het-dominated sites usually mean reads from another copy of the sequence map there. Flagged means het only or het excess; usable means pass or het elevated. Of the release\u2019s 46,054,265 sites, 18,177,811 pass, 3,756,376 are het elevated, 14,354,590 het excess, 4,752,396 het only, and 5,013,092 have no carrier. Defaults: SNPVersity and SNPGeo show every site with its class; SNPImpact, SNPFunction and SNPFold hide flagged and no-carrier sites (each has a control to show them); SNPTree, SNPMatrix and SNPCompare (region scope) use usable sites only, with a switch to all sites. Downloaded VCFs keep every site, with NHET, NHOM and SITEQC in INFO.'],
    ['Curated allele marks', 'SNPCurate records published maize alleles. A gold star (★) marks a validated causal change genotyped in this release; an outline star (☆), a published marker, associated change or tagging site genotyped here; a grey circle (○), a known allele this release cannot show. Click a mark to open the allele’s record in SNPCurate.'],
    ['Reference lines', 'For an allele, the lines homozygous for the reference allele: the panel minus the carriers and the lines with a missing call. SNPFunction and SNPFold send them to SNPVersity beside the carriers.'],
    ['Pfam domain', 'When a variant falls inside a Pfam domain of the gene\u2019s canonical protein, that domain is shown and linked to InterPro; \u2014 means the site is outside every domain, which is true of most sites (Pfam domains cover about 1% of each chromosome).'],
    ['VCF', 'The Variant Call Format file SNPVersity generates for your query. It is the exact matrix the other tools reuse when you send a selection. For the GRIN-linked 2026 release its INFO field also carries NHET, NHOM and SITEQC: the Site QC counts and class of each site.'],
    ['Gene Ontology (GO)', 'Standardized terms describing a gene product\u2019s biological process, molecular function, and cellular component. SNPFunction shows MaizeGDB\u2019s GO for the gene, each term with its evidence code: most maize gene-model terms are computational (PANNZER, the NAM annotation, UniProt imports); experimental ones come from curation of the gene\u2019s locus.'],
    ['GO evidence', 'How a GO term was assigned, from its evidence code: experimental (IDA, IMP, IGI, IPI, IEP), by similarity (ISS and kin), author statement (TAS, NAS, IC), or computational (IEA, COMP). Terms only an InterPro domain suggests (InterPro2GO) are shown apart as domain-implied, not as annotations.'],
    ['Heterotic group', 'A maize breeding classification (stiff-stalk, non-stiff-stalk, Iodent, Lancaster, tropical, teosinte, and others). PanEffect colors the pan-genome rows by the assembly\u2019s heterotic group.'],
  ];


  /* ---- complete definitions by tool ----
   * Keep table-header wording identical to the live tools so users can
   * search this page for the label they see on screen.
   */
  const DEFINITIONS = [
    { tool:'Shared terms', color:'#64748b', items:[
      ['Accession', 'A named maize line, sample, or sequencing run represented by one genotype column.'],
      ['Allele', 'One observed DNA state at a genomic position. REF is the B73 v5 reference allele; ALT is an alternate allele.'],
      ['Variant', 'A genomic position where at least one accession differs from the B73 v5 reference.'],
      ['SNP', 'Single-nucleotide polymorphism: a one-base substitution.'],
      ['INDEL', 'Insertion or deletion relative to the reference sequence.'],
      ['Gene model', 'The B73 v5 identifier and annotated exon, intron, CDS, and transcript structure assigned to a gene.'],
      ['Consequence / Effect', 'The predicted molecular result of a variant, such as synonymous, missense, splice-site, frameshift, stop gained, or intronic.'],
      ['Impact', 'A broad severity class assigned from the consequence: HIGH, MODERATE, LOW, or MODIFIER.'],
      ['Priority', 'The SNPMaize ranking tier (TOP, HIGH, MODERATE, or LOW) produced by combining consequence, model scores, and domain context.'],
      ['Domain', 'A Pfam-annotated protein domain overlapping the affected residue; — means the residue is outside every domain.'],
      ['Het', 'Heterozygous: the accession carries one reference and one alternate allele, usually displayed as 0/1 or 1/0.'],
      ['Hom', 'Alternate homozygous: the accession carries two alternate alleles, displayed as 1/1.'],
      ['Missing / ./.', 'No usable genotype call. Missing calls are excluded from pairwise similarity calculations.'],
      ['AF', 'Alternate-allele frequency: the frequency of the ALT allele among called chromosomes in the analyzed panel.'],
      ['MAF', 'Minor-allele frequency: the frequency of the less common allele, constrained to 0–0.5.'],
      ['Carrier', 'An accession with at least one copy of the alternate allele.'],
      ['Co-called sites', 'Sites where both accessions in a pair have non-missing genotype calls.'],
      ['B73 RefGen v5', 'The maize reference assembly used for coordinates, REF alleles, gene models, and annotations throughout SNPMaize.'],
      ['Site QC', 'The class of a site from its heterozygous and homozygous carriers among all lines of the release: Pass, Het elevated, Het excess, Het only, or No carrier (see Scores & annotations).'],
      ['Usable / flagged site', 'Usable: Site QC Pass or Het elevated. Flagged: Het only or Het excess. A site with no carrier is neither.'],
      ['het (line mark)', 'Beside a line: its sample is heterozygous at more than 5% of clean sites (a Heterozygous sample); hover it for the share. Inbred lines are near 0.8%.'],
      ['Curated allele', 'A published allele recorded in SNPCurate, marked gold, outline or grey wherever it appears.'],
      ['Reference lines', 'Lines homozygous for the reference allele at a site: the panel minus the carriers and the lines with a missing call.'],
      ['Load data from SNPVersity', 'A button on the empty page of SNPImpact, SNPCompare and SNPTree when SNPVersity holds a table result. It names the region, lines and variants, and runs the same hand-off as Send selection to….'],
    ]},
    { tool:'SNPVersity', color:'#2563eb', items:[
      ['CHR', 'Reference chromosome containing the variant.'],
      ['POS', 'One-based genomic coordinate on B73 v5.'],
      ['REF', 'Reference allele in B73 v5.'],
      ['ALT', 'Alternate allele represented by the row.'],
      ['Gene model', 'B73 v5 gene model overlapping or associated with the variant.'],
      ['+N (several consequences)', 'A site can affect more than one gene, for example intronic in one and downstream of its neighbor. The Gene model and Effect columns show the most severe consequence; +N lists the others on hover. SNPFunction, SNPFold and SNPGeo (gene search) read each gene\u2019s own consequence, and the ESM scores go with the substitution they were computed for.'],
      ['Effect', 'Predicted variant consequence from the annotation source.'],
      ['Impact', 'Predicted severity category: HIGH, MODERATE, LOW, or MODIFIER.'],
      ['Domain', 'Pfam protein domain overlapping the affected coding residue, when available.'],
      ['MQ', 'Mapping quality: a phred-scaled measure of confidence that reads were aligned to the correct genomic location; higher is better.'],
      ['COMP', 'Completeness: the proportion of accessions with a non-missing genotype call at that site.'],
      ['maxR²', 'Maximum linkage-disequilibrium r² with any variant 400-5,000 bp away (PLINK 1.9); values closer to 1 indicate stronger correlation. For the GRIN-linked 2026 set it is computed from the 933 release genotypes without filtering (blank = no partner variant in range, or monomorphic). The table rounds it to two decimals.'],
      ['MAF', 'Minor-allele frequency among the selected or source accessions, depending on the returned record.'],
      ['Site QC', 'The class of the site from its heterozygous and homozygous carriers among all lines of the release; hover a cell for its counts. The filter above the table shows All sites (the default), Hide flagged and no-carrier, or Passing only.'],
      ['PlantCAD1 / PlantCAD2', 'DNA language-model variant scores, on every SNP of the GRIN-linked 2026 set (indels have none). More extreme disruptive scores are prioritized according to the score convention used by the data pipeline. In the GRIN-linked 2026 set PlantCAD1/PlantCAD2 (and Evo2) are rounded to 0.1, like the ESM scores. Their headers break over two lines (Plant / CAD1) so the columns are as narrow as the other scores.'],
      ['Evo2', 'Evo2 7B DNA language-model log-likelihood ratio for the allele (256-bp left context), scored for SNPs within 1 kb of a gene; more negative is more disruptive. Filled on every chromosome of the GRIN-linked 2026 set.'],
      ['ESM1 / ESM2 / ESM3', 'Protein language-model scores for amino-acid substitutions, filled for missense sites.'],
      ['ESM-C', 'ESM C 600M protein language-model log-likelihood ratio for the amino-acid substitution, missense sites only; more negative is more disruptive. Shown where the dataset carries it.'],
      ['Accession genotype columns', 'Each accession column shows its genotype at the site: 0/0 reference homozygous, 0/1 heterozygous, 1/1 alternate homozygous, or ./. missing.'],
      ['het', 'Beside an accession name: the line is a Heterozygous sample (more than 5% of clean sites heterozygous).'],
      ['Curated row', 'A row whose site is a curated allele: its SNPCurate mark beside the position, and its site and annotation cells shaded light gold. The curated alleles of the interval, grey ones included, are named above the table.'],
      ['Gene models in this region', 'The B73 v5 gene models the queried interval overlaps, each linked to its MaizeGDB gene page.'],
      ['Dataset', 'A defined variant collection with its own accession panel, filters, included variant types, and score columns.'],
      ['Sites', 'Number of variant positions in the complete dataset, not necessarily the number returned by the current query.'],
      ['Imputed', 'Whether missing genotypes were statistically inferred in that dataset.'],
    ]},
    { tool:'GWAS Explorer', color:'#cf8a12', items:[
      ['SNP', 'The variant’s identifier.'],
      ['Chr', 'Chromosome carrying the SNP.'],
      ['BP', 'Genomic coordinate on B73 v5.'],
      ['A1/A2', 'The two alleles tested at the marker; A1 is the effect allele the beta is estimated against.'],
      ['Freq', 'Allele frequency of A1 in the GWAS panel.'],
      ['Beta', 'Estimated allele-substitution effect size.'],
      ['SE', 'Standard error of the beta estimate.'],
      ['P', 'Raw p-value for the marker–trait association.'],
      ['−log₁₀P', 'The p-value transformed as −log₁₀(P) so smaller p-values plot higher; used for the y-axis and the significance threshold.'],
      ['Trait', 'The phenotype being explored, grouped into Plant Architecture, Yield Component, or Flowering Time.'],
      ['Measure: Trait value (intercept)', 'GWAS run on the trait value in the mean environment.'],
      ['Measure: Linear plasticity (slope)', 'GWAS run on the plasticity of a given trait, as modeled by accession-specific linear slope across environments.'],
      ['GWAS method', 'The association-testing model used to generate the results. These include: MLM (GCTA --mlma), see https://yanglab.westlake.edu.cn/software/gcta/#MLMA.'],
      ['Total markers', 'The number of markers tested in the complete GWAS, independent of how many are plotted after memory-saving filtering.'],
      ['Significant markers', 'The number of tested markers at or above the active significance threshold.'],
      ['Published significance threshold', 'The p-value cutoff reported by the source publication, shown by default on the plot.'],
      ['Custom threshold', 'A user-chosen alternative to the published threshold, computed live from a chosen method and α.'],
      ['SimpleM', 'A multiple-testing correction using α divided by the effective number of independent tests, recommended for GWAS.'],
      ['Bonferroni', 'A multiple-testing correction using α divided by the total marker count, more conservative than SimpleM.'],
      ['α (alpha)', 'The genome-wide false-positive rate used to compute a significance threshold.'],
    ]},
    { tool:'SNPTrait', color:'#1f8a4c', items:[
      ['Sample heterozygosity', 'A facet from the line’s share of clean sites that are heterozygous: Inbred, Elevated heterozygosity (above 2%), or Heterozygous sample (above 5%).'],
      ['sampleQC / hetShare', 'Export columns: the sample heterozygosity class, and the heterozygous share of clean sites (0 to 1).'],
      ['+n', 'Beside an unticked value: the lines ticking it would add to the selection.'],
    ]},
    { tool:'SNPImpact', color:'#7c3aed', items:[
      ['Gene', 'Gene model associated with the candidate variant.'],
      ['Variant', 'Genomic change, generally shown as position and REF→ALT alleles.'],
      ['Consequence', 'Specific predicted molecular consequence of the change.'],
      ['Domain', 'Pfam domain containing the affected amino acid, when present.'],
      ['PlantCAD1 / PlantCAD2', 'DNA language-model scores used to estimate regulatory or sequence disruption.'],
      ['Evo2 / ESM-C', 'A further DNA score (Evo2, SNPs within 1 kb of a gene) and protein score (ESM-C, missense); empty where the store has no score yet.'],
      ['ESM / ESM2 / ESM3', 'Protein language-model scores used to estimate the effect of an amino-acid substitution.'],
      ['Priority', 'Integrated candidate tier. TOP is the strongest prioritization, followed by HIGH, MODERATE, and LOW.'],
      ['Site QC', 'The class of the variant’s site. The filter shows All sites, Hide flagged and no-carrier (the default), or Passing only; a note above the table counts the variants it hides. A curated allele is never hidden.'],
      ['Shortlist / flag', 'A user-selected marker for retaining a candidate variant for later review or export.'],
      ['Gene-model diagram', 'A compact display of exons, introns, coding sequence, strand, and the variant position.'],
      ['Exon', 'A transcript segment retained in the mature RNA; coding portions contribute to the protein sequence.'],
      ['Intron', 'A transcribed segment removed during RNA splicing.'],
      ['CDS', 'Coding sequence: the portion of exons translated into protein.'],
      ['max LD r²', 'In a variant’s details: the highest linkage-disequilibrium r² with a variant 400-5,000 bp away (maxR² in SNPVersity), to two decimals.'],
      ['MQ', 'Mapping quality for the variant site.'],
      ['MAF', 'Minor-allele frequency for the candidate variant.'],
    ]},
    { tool:'SNPFunction', color:'#2563eb', items:[
      ['Non-syn : syn', 'The number or ratio of nonsynonymous coding variants to synonymous coding variants. A higher value indicates more amino-acid-changing variation relative to silent variation; it is descriptive and is not by itself a formal dN/dS estimate.'],
      ['Exon : intron', 'The number or ratio of variants in annotated exons to variants in introns of the gene model. Infinity (∞) means exon variants were observed but no intron variants were counted.'],
      ['Domain-disrupting', 'Coding variants that alter an amino acid located inside an annotated Pfam domain.'],
      ['Knockout alleles', 'Alleles predicted to strongly disrupt gene function, such as frameshift, stop-gained, essential splice, or other loss-of-function changes.'],
      ['Mean PlantCAD1 / Mean PlantCAD2', 'Average DNA language-model score across the gene variants included in the burden summary.'],
      ['Mean ESM / Mean ESM2 / Mean ESM3', 'Average protein language-model score across scored amino-acid-changing variants in the gene.'],
      ['Allele', 'The specific genomic REF→ALT change represented by a damaging-allele row.'],
      ['Consequence', 'Predicted molecular effect of that allele.'],
      ['Domain', 'Pfam domain overlapping the affected residue.'],
      ['PlantCAD1 / PlantCAD2', 'DNA language-model score for the allele.'],
      ['Evo2', 'Evo2 DNA language-model score for the allele (SNPs within 1 kb of a gene), with ESM-C beside the ESM scores; empty where the store has no score yet.'],
      ['ESM / ESM2 / ESM3', 'Protein language-model score for the resulting amino-acid change.'],
      ['Priority', 'Integrated SNPMaize evidence tier for the allele.'],
      ['Site QC', 'The class of the allele’s site from its heterozygous and homozygous carriers across the panel.'],
      ['Het', 'Number of accessions carrying the allele heterozygously.'],
      ['Hom', 'Number of accessions carrying the allele as alternate homozygous.'],
      ['AF', 'Alternate-allele frequency across the whole analyzed panel.'],
      ['Variant burden', 'The count and composition of variants assigned to the gene across the full dataset panel, over usable sites when the release has Site QC.'],
      ['Damaging allele', 'An allele selected because its consequence and/or prediction scores indicate a potentially important functional effect.'],
      ['Flagged calls', 'A collapsed group under the allele list: damaging alleles whose site is flagged (het only or het excess), with the same columns.'],
      ['No carrier in this release', 'A collapsed group of damaging alleles no line of the release carries.'],
      ['Curated alleles', 'The gene’s SNPCurate alleles at the top of the allele card, each with its site, carriers, Site QC and priority, or why the release cannot show it.'],
      ['Open in SNPVersity', 'Per allele: Homozygous carriers, Heterozygous carriers, All carriers, or Reference lines (homozygous for the reference allele) become SNPVersity’s accession list.'],
      ['Functional annotation', 'A per-gene dossier assembled for the 39,756 canonical B73 v5 gene models: identity, description, protein domains, and cross-references. Its GO and KEGG lists are the fallback shown when MaizeGDB cannot be reached.'],
      ['Gene symbol / aliases', 'The primary gene symbol and any additional names or synonyms recorded for the model.'],
      ['Description source', 'Provenance of the functional description, shown as a badge: MaizeGDB (a curated gene product or locus name), UniProt, InterPro (predicted from domains), or no informative source.'],
      ['Evidence chips', 'An at-a-glance row of what is annotated for the gene: GO, Pfam, KEGG KO, Pathway, UniProt, Symbol.'],
      ['Protein domain architecture', 'A to-scale diagram of the protein with Pfam domains as positioned blocks, plus a list giving Pfam and InterPro IDs, residue span, percent of protein covered, and InterProScan E-value.'],
      ['Gene Ontology view', 'MaizeGDB\u2019s GO view, one row per aspect (Biological process, Molecular function, Cellular component): the ancestry of the gene\u2019s terms up to the root, a fingerprint of the plant GO-slim categories they fall under (one square per category, always in the same order, so genes can be compared at a glance), and the terms with their evidence beneath. Hover a square, name, term or node to highlight what belongs together; click to pin.'],
      ['Plant GO slim', 'A fixed subset of GO categories chosen for plants (94 shown: 45 biological process, 25 molecular function, 24 cellular component). A lit square means at least one of the gene\u2019s terms falls under that category; a dashed one, only a domain-implied term.'],
      ['Obsolete term', 'A GO term the current GO release has retired; it is struck through and placed through its replacement.'],
      ['Pathways view', 'MaizeGDB\u2019s pathway view. Metabolic pathways come from the pan-genome pathway explorer (E2P2 on the 26 NAM founders): each pathway is a strip of reaction steps, this gene\u2019s step in gold, steps with another B73 gene outlined, steps with no B73 gene dashed, and a bar under each step for how many founders have a gene there. KEGG and EggNOG come from EnTAP: the EggNOG description and COG category, KEGG orthologs, and the KEGG maps KEGG lists for maize.'],
      ['Fallback annotation', 'When MaizeGDB cannot be reached, SNPFunction lists the GO terms and KEGG pathways from its own annotation file instead, without the ontology placement. Its \u201cMaizeGDB\u201d GO is MaizeGDB\u2019s computational gene-model GO, and its KEGG pathways cover only genes with UniProt/Entrez cross-references.'],
      ['Cross-references', 'External identifiers for the gene: UniProt, NCBI Gene, and B73 v4 / v3 gene models.'],
      ['Annotation build', 'The build date, assembly, and annotation version the functional record was generated from.'],
    ]},
    { tool:'SNPCurate', color:'#a16207', items:[
      ['Mark', 'Gold: validated causal change, genotyped here. Outline: published marker, associated change or tag, genotyped here. Grey: known allele this release cannot show.'],
      ['Gene', 'Gene symbol and B73 v5 model.'],
      ['Change', 'The published change.'],
      ['Trait', 'The trait the allele affects.'],
      ['Status in this release', 'Whether the release genotypes the allele, and why not when it does not (not a site, a structural variant, unreliable calls).'],
      ['Carriers (het / hom)', 'Heterozygous / homozygous carriers among the release lines.'],
      ['Site QC', 'The class of the site from its heterozygous and homozygous carriers.'],
      ['Priority', 'The SNPMaize priority of the change, as SNPImpact computes it.'],
      ['Reference', 'The first publication; the record lists all.'],
      ['Default order', 'By mark (gold, outline, grey), then by gene name; click Gene or Carriers to sort.'],
    ]},
    { tool:'SNPCompare', color:'#0e7490', items:[
      ['#', 'Current rank after sorting and filtering.'],
      ['Project', 'BioProject or dataset project associated with the accession.'],
      ['SRA ID', 'Sequence Read Archive run identifier associated with the accession.'],
      ['Accession name', 'Human-readable line or germplasm name.'],
      ['PI number', 'GRIN accession number of the line (PI, Ames or NSL); empty for a set without GRIN links.'],
      ['Global sim', 'Genome-wide identity-by-state similarity between the focal accession and the comparison accession: matching genotypes divided by co-called sites.'],
      ['Local sim', 'Identity-by-state similarity calculated only from variants in the current SNPVersity region.'],
      ['Δ (local−global)', 'Local similarity minus global similarity. Positive values indicate the pair is more similar in the selected region than genome-wide; negative values indicate less similarity in the region.'],
      ['Global miss%', 'Percentage of genome-wide comparison sites where at least one accession lacks a genotype call.'],
      ['Local miss%', 'Percentage of sites in the selected region where at least one accession lacks a genotype call.'],
      ['Similarity', 'Fraction of co-called sites at which the two accessions have the same genotype.'],
      ['Missing%', 'Percentage of evaluated sites where either member of the pair is missing a genotype.'],
      ['Co-called sites', 'Number of sites in the local region with genotype calls for both accessions.'],
      ['Focal accession', 'The accession against which every other accession is ranked.'],
      ['Global', 'Precomputed genome-wide comparison scope.'],
      ['This region / Local', 'Comparison calculated from the current SNPVersity genotype matrix.'],
      ['Both', 'Side-by-side display of global and local values plus their difference.'],
      ['Region sites', 'This region scope: Usable only (the default) leaves out flagged sites and sites with no carrier; All sites uses every site. The genome-wide values are precomputed over all sites.'],
    ]},
    { tool:'SNPTree', color:'#15803d', items:[
      ['IBS allele distance', 'Mean pairwise allele-dosage difference across co-called sites. Identical genotypes contribute 0, opposite homozygotes 1, and a homozygote-versus-heterozygote comparison 0.5.'],
      ['Informative sites', 'Variant sites that contain more than one observed genotype state among the selected accessions.'],
      ['Shared sites', 'Sites with non-missing calls for both accessions in a pair.'],
      ['Sites', 'Usable only (the default) leaves out flagged sites (het only, het excess) and sites with no carrier; All uses every site.'],
      ['UPGMA', 'Unweighted Pair Group Method with Arithmetic Mean, an agglomerative clustering method that assumes an ultrametric tree.'],
      ['Neighbour-Joining (NJ)', 'A distance-based tree-building method that does not require equal evolutionary rates among branches.'],
      ['Branch length', 'Distance assigned to a tree edge from the pairwise IBS distance calculation.'],
      ['Newick', 'Compact parenthetical text format for tree topology and branch lengths.'],
      ['MEGA pairwise format', 'Distance-matrix text format that can be imported by MEGA software.'],
      ['PHYLIP', 'Fixed-layout distance-matrix format used by many phylogenetic programs.'],
      ['Local phylogeny', 'A distance tree describing similarity in the selected genomic interval; it should not automatically be interpreted as a whole-genome species tree.'],
    ]},
    { tool:'SNPMatrix', color:'#b45309', items:[
      ['IBS distance', 'Pairwise mean allele-dosage difference across co-called sites; 0 means identical across compared calls and larger values indicate more difference.'],
      ['% identity', 'Similarity view derived from the distance matrix and displayed as a percentage.'],
      ['Shared-site count', 'Number of non-missing sites used for a particular pairwise matrix cell.'],
      ['Sites', 'Usable only (the default) or All sites, as in SNPTree.'],
      ['Input order', 'Rows and columns remain in the same order as the accession selection.'],
      ['Clustered order', 'Rows and columns are reordered by UPGMA clustering to place similar accessions near one another.'],
      ['Bioproject bars', 'Color strips indicating the project or BioProject associated with each accession.'],
      ['Heatmap cell', 'The pairwise distance or identity value for the row accession versus the column accession.'],
      ['CSV distance matrix', 'Comma-separated pairwise IBS distance table.'],
      ['PNG / SVG', 'Raster and scalable-vector image exports of the displayed heatmap.'],
      ['PHYLIP', 'Distance-matrix export for compatible phylogenetic software.'],
    ]},
    { tool:'SNPFold', color:'#be185d', items:[
      ['Variant', 'Coding DNA change mapped to a protein residue.'],
      ['Consequence', 'Predicted molecular effect of the coding change.'],
      ['Residue', 'One-based amino-acid position in the displayed protein sequence.'],
      ['Domain', 'Pfam domain overlapping the residue.'],
      ['Local pLDDT', 'Confidence of the structure model shown (AlphaFold2, Boltz2 or ESMFold) at the affected residue, from 0 to 100. Higher values indicate greater confidence in the local predicted structure.'],
      ['Structure', 'DSSP-style secondary-structure assignment at the residue: α-helix, β-strand, or loop/coil.'],
      ['PlantCAD / PlantCAD2', 'DNA language-model scores for the underlying nucleotide variant.'],
      ['ESM1 / ESM2 / ESM3', 'Protein language-model scores for the amino-acid substitution.'],
      ['Evo2', 'Evo2 DNA language-model score for the variant (SNPs within 1 kb of a gene); empty where the store has no score yet.'],
      ['ESM-C', 'ESM C 600M protein language-model log-likelihood ratio for the amino-acid substitution, missense sites only; more negative is more disruptive. Shown where the dataset carries it.'],
      ['IUPred2', 'Predicted intrinsic disorder at the residue, from 0 to 1 (IUPred2). Higher values mean the residue is more likely to be intrinsically disordered than to adopt a fixed fold.'],
      ['Anchor2', 'Predicted probability (0 to 1) that a disordered residue lies in a protein-binding region (ANCHOR2). Higher values flag disordered segments likely to become ordered upon binding a partner.'],
      ['Activity', 'Annotated functional site overlapping the residue, drawn from InterProScan member databases (for example an active site, binding site, or conserved functional feature), when present.'],
      ['Priority', 'Integrated evidence tier for the structure-mapped variant.'],
      ['Carriers', 'Number of accessions carrying the alternate allele; expanded details separate heterozygous and homozygous carriers.'],
      ['Site QC', 'The class of the variant’s site from its heterozygous and homozygous carriers across the panel.'],
      ['Usable alleles only', 'On by default: hides coding variants at flagged or no-carrier sites from the track, the table and the 3D view.'],
      ['Curated ring', 'A ring around a variant on the protein track, and a mark in the table, for a curated allele (gold for a validated causal change).'],
      ['Reference lines', 'In a variant’s details: the lines homozygous for the reference allele, which Open in SNPVersity can send beside the carriers.'],
      ['pLDDT', 'Predicted Local Distance Difference Test: the per-residue confidence each structure predictor reports (AlphaFold2, Boltz2 and ESMFold all do); the legend names the model shown. Common interpretation: ≥90 very high confidence, 70–89 confident, 50–69 low, and <50 very low.'],
      ['Secondary structure', 'Local protein conformation classified as helix, strand, or coil/loop, generated from the structure model.'],
      ['ΔΔG', 'Predicted change in protein folding free energy after mutation. Positive and negative interpretations depend on the scoring convention used by the source model; magnitude reflects predicted structural effect.'],
      ['Protein browser', 'Linear alignment of protein residues, domains, secondary structure, confidence, and variant markers.'],
      ['3D viewer color: confidence', 'Colors residues by pLDDT confidence band.'],
      ['3D viewer color: domain', 'Colors residues by Pfam-domain membership.'],
      ['3D viewer color: impact', 'Colors or highlights variants according to predicted functional severity.'],
      ['3D', 'A per-row control in the variant table that force-shows that residue in the 3D structure viewer, independent of the blanket "Variant residues" toggle.'],
      ['InterProScan sites (by program)', 'Residue-level functional-site and feature annotations from InterProScan, drawn as lanes in the protein browser and grouped by the member database (program) that produced them — CDD, PIRSR, and SFLD. All program lanes are shown even when the current protein has no hits, so the annotation is discoverable.'],
      ['CDD (Conserved Domains Database)', "NCBI's Conserved Domain Database. In the InterProScan site lanes, CDD contributes conserved-domain footprints and the functional-site residues within them."],
      ['PIRSR (PIR Site Rules)', 'PIR Site Rules from the Protein Information Resource. Curated site rules propagate functional-site residues (such as active or binding sites) to proteins that match a site-rule signature.'],
      ['SFLD (Structure–Function Linkage Database)', 'The Structure–Function Linkage Database links protein structural features to specific chemical functions and contributes conserved functional-site residues, chiefly for enzyme superfamilies.'],
    ]},
    { tool:'PanEffect', color:'#b45309', items:[
      ['B73 Reference View', 'The substitution heatmap computed on the B73 v5 protein: every residue (columns) against all 20 possible amino acids (rows).'],
      ['Pan-genome View', 'The same substitution heatmap projected onto a protein multiple-sequence alignment across the maize assemblies, so natural variation and its predicted effect can be read together. Available for canonical transcripts.'],
      ['ESM score', 'A protein-language-model estimate of how tolerated an amino-acid substitution is; more negative is more disruptive.'],
      ['ESM model (ESM1 / ESM2 / ESM3)', 'The protein language model used to score substitutions, selectable in the panel.'],
      ['Substitution heatmap', 'A per-residue grid colored by predicted effect. Hover a cell to read the position, wild-type → substitution, and score.'],
      ['Zoomed region', 'A 50-residue window of the full heatmap, positioned with the slider, showing per-cell substitution letters and the wild-type residue track.'],
      ['Show all variant effects', 'Colors every possible substitution in the B73 heatmap. Hand-offs from the other tools open in this view with their substitution highlighted. (The MaizeGDB 2026 High-Coverage view is turned off in this release: PanEffect\u2019s files do not yet flag the GRIN-linked 2026 substitutions.)'],
      ['PFAM Domains track', 'Pfam domains drawn to scale along the protein and linked to InterPro, shared with the SNPFold and SNPFunction views.'],
      ['Secondary structure', 'Predicted per-residue secondary structure (helix / strand / coil) drawn above the heatmap for context.'],
      ['Species / Gene model', 'In the pan-genome zoomed view, switches the row labels between assembly (species) names and their gene-model IDs.'],
      ['Heterotic group color', 'Pan-genome rows are colored by the assembly\u2019s heterotic group (stiff-stalk, non-stiff-stalk, Iodent, Lancaster, tropical, teosinte, and others).'],
      ['Variant effects file', 'The per-substitution ESM score table for the gene, downloadable from the summary.'],
    ]},
    { tool:'Dataset table', color:'#475569', items:[
      ['Dataset', 'Named variant collection and filtering configuration used for a query.'],
      ['Reference', 'Genome assembly to which reads and variants were aligned.'],
      ['Accessions', 'Number of samples or accession columns available in the dataset.'],
      ['Sites', 'Approximate or reported number of variant sites in the complete dataset.'],
      ['Filters', 'Quality and inclusion rules used to construct the dataset.'],
      ['Het', 'Whether heterozygous genotype calls are retained.'],
      ['INDELs', 'Whether insertion and deletion variants are included.'],
      ['Imputed', 'Whether missing genotype calls have been statistically inferred.'],
    ]},
  ];


  /* ---- FAQ ---- */
  const FAQ = [
    ['Where does the data come from?',
     'A query sends your region and accession list to the server, which reads the real HDF5 variant store, writes a VCF for exactly that slice, and returns it. The tools parse that VCF into the tables and matrices you see, so everything downstream is one consistent result.'],
    ['How do I move a selection between tools?',
     'Run a query in SNPVersity or highlight a region in GWAS Explorer, then use "Send selection to…" in the top bar (or the buttons on a result). The same genotype matrix is handed to SNPImpact, SNPCompare, SNPTree, SNPMatrix, and SNPGeo without re-querying. If you open SNPImpact, SNPCompare or SNPTree first, its empty page offers the current SNPVersity result: press Load data from SNPVersity. SNPMatrix and SNPTree can also pass their set on to each other. The selection chip in the top bar (a list icon and a count) shows the accessions SNPVersity will query, from any tool: click it to browse, filter and edit them in the line selector; outside SNPVersity the change applies to SNPVersity\u2019s next query, with Undo.'],
    ['Why did a query give me a download instead of a table?',
     'When the result is large — roughly when variants × selected accessions passes ~40 million, or the region alone would exceed ~400,000 variant sites — parsing and rendering it in the browser would freeze the page, so SNPVersity returns a downloadable VCF instead. Choose a smaller region, fewer accessions, or a lower-density SNP set to get the interactive table back. The run bar predicts which you’ll get before you run. The server builds at most about 2 billion genotype calls (variants × accessions) in one request — roughly 100 Mb for all 933 GRIN-linked lines — and says so instead of building anything larger. A built VCF stays downloadable for 24 hours.'],
    ['Why are some cells \u2014 or N/A?',
     'A \u2014 in the Domain column means the site lies outside every Pfam domain of the canonical proteins: introns, UTRs, intergenic sites and coding stretches between domains, so most sites. Language-model columns read N/A where a score does not apply: ESM and ESM-C for non-missense sites, PlantCAD for indels, Evo2 beyond 1 kb of a gene. MQ and coverage are not recorded for the GRIN-linked 2026 call set, so SNPVersity leaves those two columns out for it.'],
    ['Is there a limit on how many accessions I can compare?',
     'SNPTree, SNPMatrix and SNPCompare warn before a long computation — the cost grows with variants × accessions², so a large region with many accessions is what gets slow. You can build anyway. Sending a result to those tools is blocked outright only when it would be certain to crash the tab. SNPImpact renders up to 1,500 variants at a time.'],
    ['What can I download?',
     'A VCF from SNPVersity; a CSV distance matrix, PHYLIP, PNG, and SVG from SNPMatrix; Newick, MEGA, and PHYLIP trees from SNPTree. A genome-wide significant-SNP CSV, or a CSV of a selected region, from GWAS Explorer. Comparison and impact tables can be exported from their own pages, and SNPTrait exports its selection as CSV or JSON. The VCF of the GRIN-linked 2026 release keeps every site, with NHET, NHOM and SITEQC in INFO, whatever the tools hide.'],
    ['Where do the SNPTrait and SNPGeo metadata come from?',
     'For the MaizeGDB GRIN-linked 2026 dataset, from the USDA-ARS GRIN-Global passport and evaluation records of each line\u2019s accession, joined to the Grzybowski et al. (2023) sample names. Origin is the GRIN country and state of origin (developed, donated or collected); no coordinates are inferred. Trait values are per-accession means (numeric) or modes (coded) over the GRIN observations. Base maps: Natural Earth (public domain).'],
    ['Why does SNPGeo warn that carrier percentages understate frequencies?',
     'The carrier percentage for a country is over every line of the dataset known from it. If only some of those lines were genotyped (a partial SNPVersity selection), the rest count as non-carriers. Use the gene search, which queries every line, or the allele-composition and allele-frequency modes, which use only the called lines.'],
    ['Why do SNPImpact, SNPFunction and SNPFold show fewer variants than SNPVersity?',
     'Site QC. The lines of the GRIN-linked 2026 release are inbreds, so a site where heterozygous calls dominate usually means reads from another copy of the sequence map there. SNPImpact, SNPFunction and SNPFold hide flagged sites (het only, het excess) and sites with no carrier by default; SNPVersity and SNPGeo show every site with its class. To see them, use the Site QC filter in SNPImpact, open the Flagged calls and No carrier groups in SNPFunction, or untick Usable alleles only in SNPFold. Curated alleles stay in view. SNPTree, SNPMatrix and SNPCompare (region scope) likewise compute distances from usable sites, with a switch to all sites.'],
    ['What do the gold, outline and grey marks mean?',
     'They are SNPCurate marks on published alleles: gold (★) a validated causal change genotyped in this release, outline (☆) a published marker, associated change or tagging site genotyped here, grey (○) a known allele the release cannot show. In SNPVersity a curated row is also shaded light gold. Click a mark to open the allele’s record in SNPCurate.'],
    ['Why is a line marked het?',
     'Its sample is heterozygous at more than 5% of clean sites, where inbred lines are near 0.8%, which can mean residual heterozygosity, an outcross or a mixed sample. Read its calls with care. 20 lines of the release are marked; SNPTrait’s Sample heterozygosity facet selects or excludes them.'],
    ['Which tools are ready to use now?',
     'SNPVersity, GWAS Explorer, SNPTrait, SNPImpact, SNPFunction, SNPCurate, SNPFold, SNPGeo, SNPCompare, SNPTree, SNPMatrix, and PanEffect are live. SNPTrait and SNPGeo currently have metadata for the MaizeGDB GRIN-linked 2026 dataset only. SNPImpute, SNPDensity, and SNPGermplasm are on the roadmap.'],
  ];

  /* ---- what's new (the changes of October 2026) ---- */
  const NEWS = [
    ['SNPMaize', 'SNPTools is now SNPMaize: a SNP toolkit to explore variant diversity across maize germplasm.'],
    ['Site QC', 'Every site of the GRIN-linked 2026 release is classed from its heterozygous and homozygous carriers (Pass, Het elevated, Het excess, Het only, No carrier). SNPVersity and SNPGeo show the class on every site; SNPImpact, SNPFunction and SNPFold hide flagged and no-carrier sites by default; SNPTree, SNPMatrix and SNPCompare compute distances from usable sites, with a switch to all sites. Downloaded VCFs carry NHET, NHOM and SITEQC.'],
    ['Sample heterozygosity', 'Each line is classed Inbred, Elevated heterozygosity or Heterozygous sample from its share of heterozygous clean sites. SNPTrait filters by it and exports it; the 20 heterozygous samples are marked het in SNPVersity and SNPFunction.'],
    ['SNPCurate', 'A new registry of published maize alleles and what this release shows of them. Its gold, outline and grey marks appear in SNPVersity (where a curated row is shaded light gold), SNPImpact, SNPFunction, SNPFold and SNPGeo.'],
    ['Reference lines', 'SNPFunction and SNPFold can send the lines homozygous for the reference allele to SNPVersity, beside the homozygous, heterozygous and all carriers.'],
    ['Load data from SNPVersity', 'SNPImpact, SNPCompare and SNPTree offer the current SNPVersity result on their empty page, in one click.'],
    ['SNPCompare', 'A PI number column gives each line’s GRIN accession.'],
    ['SNPVersity table', 'maxR² is shown to two decimals, the PlantCAD columns are as narrow as the other scores, and the gene models of the region link to their MaizeGDB pages.'],
    ['Phones', 'SNPTrait and SNPFunction fit a 375 px screen.'],
  ];

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  /* Turns a bare URL inside already-escaped text into a clickable link,
     trimming trailing sentence punctuation (periods, commas, closing
     parens, ...) out of the link itself. */
  function linkify(escaped){
    return escaped.replace(/https?:\/\/[^\s<]+/g, function(url){
      let trail = '';
      while (/[.,;:)]$/.test(url)) { trail = url.slice(-1) + trail; url = url.slice(0, -1); }
      return '<a href="'+url+'" target="_blank" rel="noopener">'+url+'</a>'+trail;
    });
  }
  function ico(k){ return (typeof ICONS!=='undefined' && ICONS[k]) || ''; }

  /* ---- dataset table (live from data.js when available) ---- */
  function datasets(){
    try { if (typeof Data!=='undefined' && Data.datasets) return Data.datasets(); } catch(e){}
    return [];
  }
  function flag(on){ return on ? '<span class="hp-yes">Yes</span>' : '<span class="hp-no">\u2013</span>'; }
  function datasetTable(){
    const ds = datasets();
    if (!ds.length) return '';
    const rows = ds.map(d=>`<tr>
        <td><b>${esc(d.name)}</b><div class="hp-sub">${esc(d.sub||'')}</div></td>
        <td class="hp-mono">${esc(d.ref)}</td>
        <td class="hp-num">${esc(d.acc)}</td>
        <td class="hp-num">${esc(d.sites)}</td>
        <td>${(d.filters||[]).map(f=>`<span class="hp-chip">${esc(f)}</span>`).join(' ')}</td>
        <td class="hp-c">${flag(d.het)}</td>
        <td class="hp-c">${flag(d.indel)}</td>
        <td class="hp-c">${flag(d.impute)}</td>
      </tr>`).join('');
    return `<div class="hp-tablewrap"><table class="hp-table">
      <thead><tr>
        <th>Dataset</th><th>Reference</th><th>Accessions</th><th>Sites</th>
        <th>Filters</th><th>Het</th><th>INDELs</th><th>Imputed</th>
      </tr></thead><tbody>${rows}</tbody></table></div>
      <p class="hp-fine">This release offers the <b>MaizeGDB GRIN-linked 2026</b> dataset (release v1.4): 926 lines of Grzybowski et al. (2023), Beagle-imputed, plus 7 lines called separately at the same sites, including NAM founder CML103, each joined to its GRIN accession. Every SNPMaize tool reads the same set. Each site also carries a Site QC class from its heterozygous and homozygous carriers, and each line a sample heterozygosity class; both are computed from the release’s genotypes and leave them unchanged.</p>`;
  }

  /* ---- in-depth MaizeGDB 2026 dataset description ---- */
  const MGDB2026_PROJECTS = [
    [539,'PRJCA009749','WGS resequencing of 1,604 maize inbred lines','10.1038/s41477-022-01190-2'],
    [521,'PRJNA531553','Deep DNA resequencing of the association mapping panel','10.1038/s41588-019-0427-6'],
    [453,'PRJNA609577','Zea mays genome sequencing','10.1038/s41588-020-0671-9'],
    [340,'PRJNA783885','Maize landrace whole-genome resequencing','10.1038/s41467-022-32180-9'],
    [232,'PRJEB56320','Maize Wisconsin Diversity Panel resequencing','10.1111/tpj.16123'],
    [183,'PRJEB56320','Zea mays sequences (teosinte)','10.1038/s41588-022-01184-y'],
    [77,'PRJNA641489','Maize Nested Association Mapping (NAM)','10.1126/science.abg5289'],
    [76,'PRJEB31061','Maize Haplotype Map version 3','10.1093/gigascience/gix134'],
    [67,'PRJNA399729','Maize landraces from six highland and lowland populations','10.1093/gigascience/gix134'],
    [67,'PRJNA300309','European maize diversity','10.1371/journal.pgen.1006666'],
    [57,'PRJNA783885','Genetic diversity of Zea (teosinte)','10.1038/ng.2313'],
    [49,'PRJNA389800','Whole-genome sequencing of the maize 282 panel','10.1093/gigascience/gix134'],
    [33,'PRJNA260788','European maize genomes','10.1038/s41588-020-0671-9'],
    [7,'PRJNA479960','South American maize genome sequencing','10.1126/science.aav0207'],
    [4,'PRJEB32225','Zm-B73-REFERENCE-NAM-5.0','10.1126/science.abg5289'],
    [3,'PRJEB56265','Resequencing of three Polish maize inbred lines','10.1111/tpj.16123'],
    [1,'PRJEB61159','Coastal preceramic maize from Paredones, Peru','10.7554/eLife.83149'],
    [1,'PRJNA352392','A 5,310-year-old maize cob from the Tehuacan Valley, Mexico','10.1016/j.cub.2016.09.036']
  ];

  const MGDB2026_EFFECTS = [
    ['Intergenic',274890726,90284097],
    ["5' UTR",836712,398940],
    ['Synonymous',1409639,670709],
    ['Missense',1767459,747451],
    ['Stop-related',84653,26405],
    ['Frameshift',222401,71174],
    ['Intron',12043697,5916121],
    ['Non-coding',6777,3133],
    ["3' UTR",1231422,640830],
    ['Other',120889,46629],
    ['Total',290043644,97572347]
  ];

  const MGDB2026_CHR = [
    ['Chr1',308452471,51382911,41928591,14262099],
    ['Chr2',243675191,40457547,32643386,11008991],
    ['Chr3',238017767,39598148,32693292,10861052],
    ['Chr4',250330460,42142813,35830675,12496251],
    ['Chr5',226353449,36679188,30747268,10298954],
    ['Chr6',181357234,28914586,23742034,8001535],
    ['Chr7',185808916,30391975,24862777,8231456],
    ['Chr8',182411202,29883638,24964435,8054636],
    ['Chr9',163004744,26922259,21747152,7278487],
    ['Chr10',152435371,25319421,20883764,7078616],
    ['Total',2131846805,351692486,290043374,97572077]
  ];

  function fmtInt(n){ return Number(n).toLocaleString('en-US'); }

  function maize2026Details(){
    const projects = MGDB2026_PROJECTS.map(r=>`<tr>
      <td class="hp-num">${fmtInt(r[0])}</td>
      <td class="hp-mono">${esc(r[1])}</td>
      <td>${esc(r[2])}</td>
      <td class="hp-mono">${esc(r[3])}</td>
    </tr>`).join('');

    const effects = MGDB2026_EFFECTS.map((r,i)=>`<tr${i===MGDB2026_EFFECTS.length-1?' class="hp-total"':''}>
      <td>${esc(r[0])}</td><td class="hp-num">${fmtInt(r[1])}</td><td class="hp-num">${fmtInt(r[2])}</td>
    </tr>`).join('');

    const chr = MGDB2026_CHR.map((r,i)=>`<tr${i===MGDB2026_CHR.length-1?' class="hp-total"':''}>
      <td class="hp-mono">${esc(r[0])}</td><td class="hp-num">${fmtInt(r[1])}</td>
      <td class="hp-num">${fmtInt(r[2])}</td><td class="hp-num">${fmtInt(r[3])}</td><td class="hp-num">${fmtInt(r[4])}</td>
    </tr>`).join('');

    return `
        <br>
      <details class="hp-data-detail" open>
        <summary><b>MaizeGDB 2026 dataset: composition, processing, and scale</b><span class="hp-chev">${ico('caret')}</span></summary>
        <div class="hp-data-body">
          <p>The MaizeGDB 2026 resource combines public whole-genome resequencing from <b>2,710 maize accessions</b>, including diverse inbred lines, landraces, teosintes, association panels, NAM founders and related materials, and several historically important samples. All reads were processed through a standardized variant-calling workflow against <b>B73 RefGen_v5</b>, so coordinates, reference alleles, gene models, and downstream annotations use one common reference system.</p>

          <h3>How the two MaizeGDB 2026 datasets differ</h3>
          <div class="hp-comparegrid">
            <div class="hp-dcard"><b>High Coverage</b><p>Contains approximately 290 million loci that passed mapping-quality and genotypic-coverage requirements. This version retains broader variation and is useful when sensitivity and variant discovery are the main goals.</p></div>
            <div class="hp-dcard"><b>High Quality</b><p>Contains approximately 98 million loci. It applies the same mapping-quality and coverage filters plus an additional high-confidence linkage-disequilibrium criterion. This more conservative set is useful when specificity and confidence are priorities.</p></div>
          </div>
          <p class="hp-fine">“High Coverage” describes the broader filtered set; “High Quality” is the stricter subset. A locus count refers to a genomic variant position in the complete dataset, not the number of rows returned for a particular region or accession selection.</p>

          <h3>What is stored and displayed</h3>
          <p>The underlying resource is maintained in VCF and HDF5 forms. SNPVersity sends the selected dataset, genomic interval, and accession list to the server, which extracts the requested slice and returns a VCF. SNPMaize then reuses that same genotype matrix in SNPImpact, SNPCompare, SNPTree, SNPMatrix, and related pages. Each site may include REF and ALT alleles, accession genotypes, predicted molecular consequence, gene association, mapping quality, genotype completeness, linkage-disequilibrium support, allele frequency, Pfam domain context, and DNA- or protein-language-model scores when available.</p>

          <h3>Source projects</h3>
          <p>The panel was assembled from multiple public projects rather than one experiment. This increases biological and geographic diversity, but it also means sequencing depth, library preparation, and project design can differ among accessions. The standardized alignment, calling, and filtering workflow reduces—though does not completely remove—these study-to-study differences.</p>
          <div class="hp-tablewrap"><table class="hp-table hp-compact">
            <thead><tr><th>Accessions</th><th>BioProject</th><th>Project or population</th><th>DOI</th></tr></thead>
            <tbody>${projects}</tbody>
          </table></div>

          <h3>Variant-effect composition</h3>
          <p>Most loci are intergenic because much of the maize genome lies outside annotated coding regions. Coding and gene-associated categories are much smaller but are especially important for SNPImpact, SNPFunction, and SNPFold. “Stop-related” summarizes variants annotated in the source statistics as stop effects; the exact transcript-level consequence shown in a result may be more specific.</p>
          <div class="hp-tablewrap"><table class="hp-table hp-compact">
            <thead><tr><th>Effect category</th><th>High Coverage</th><th>High Quality</th></tr></thead>
            <tbody>${effects}</tbody>
          </table></div>

          <h3>Distribution by chromosome</h3>
          <p>Variant counts broadly track chromosome length, although local diversity, repetitive sequence, mappability, selection, and the composition of the accession panel also affect density. “Raw” is the pre-filter total reported in the source summary; the two filtered columns show the successive retained sets.</p>
          <div class="hp-tablewrap"><table class="hp-table hp-compact">
            <thead><tr><th>Chromosome</th><th>Length (bp)</th><th>Raw variants</th><th>High Coverage</th><th>High Quality</th></tr></thead>
            <tbody>${chr}</tbody>
          </table></div>

          <h3>Important interpretation notes</h3>
          <ul class="hp-notes">
            <li><b>Reference-relative calls:</b> REF and ALT are defined relative to B73 v5; “alternate” does not mean rare, harmful, or derived.</li>
            <li><b>Accessions are not all independent:</b> some projects include replicates, related lines, founders, or multiple runs. Project and accession metadata should be considered when interpreting similarity.</li>
            <li><b>Missingness varies:</b> genotype completeness can differ by site and accession. Pairwise tools exclude sites where either member lacks a usable call.</li>
            <li><b>Consequence is transcript dependent:</b> one genomic variant may receive several annotations when it overlaps multiple transcripts; interfaces generally display the most severe or most relevant consequence.</li>
            <li><b>Model scores are predictions:</b> PlantCAD and ESM scores prioritize candidates but do not establish biological causality. Use them together with frequency, consequence, domain, structure, phenotype, and experimental evidence.</li>
            <li><b>Counts can differ slightly among summaries:</b> totals generated at different pipeline stages or from differently normalized records may vary by a small number of sites. The live dataset metadata and downloadable files are authoritative for an analysis.</li>
          </ul>
        </div>
      </details>`;
  }

  function statusPill(s){
    const live = s===LIVE;
    return `<span class="hp-status ${live?'on':'off'}">${live?'Live':'In development'}</span>`;
  }
  function pageCard(p){
    const live = p.status===LIVE;
    const openBtn = live
      ? `<button class="hp-open" onclick="go('${p.id}')">Open ${esc(p.name)} ${ico('caret')}</button>`
      : `<span class="hp-open muted" title="On the roadmap">Coming soon</span>`;
    return `<details class="hp-tool">
      <summary>
        <span class="hp-ti" style="background:${p.color}">${ico(p.icon)}</span>
        <span class="hp-tt"><b>${esc(p.name)}</b><span class="hp-tag">${esc(p.tag)}</span></span>
        ${statusPill(p.status)}
        <span class="hp-chev">${ico('caret')}</span>
      </summary>
      <div class="hp-tbody">
        <p>${esc(p.what)}</p>
        <div class="hp-io">
          <div><span class="hp-lbl">You give it</span>${esc(p.give)}</div>
          <div><span class="hp-lbl">You get back</span>${esc(p.get)}</div>
        </div>
        <div class="hp-actions">${openBtn}</div>
      </div>
    </details>`;
  }


  function definitionGroup(g){
    const items = g.items.map(x=>`<div class="hp-def"><dt>${esc(x[0])}</dt><dd>${linkify(esc(x[1]))}</dd></div>`).join('');
    return `<details class="hp-defgroup">
      <summary><span class="hp-defdot" style="background:${g.color}"></span><b>${esc(g.tool)}</b><span class="hp-count">${g.items.length} definitions</span><span class="hp-chev">${ico('caret')}</span></summary>
      <dl class="hp-deflist">${items}</dl>
    </details>`;
  }

  function render(page){
    injectCSS();
    page = page || document.getElementById('page');
    page.className = 'page fade';
    const crumb=document.getElementById('crumbTool'); if(crumb) crumb.innerHTML='<b>Help &amp; FAQ</b>';

    const toolCards = PAGES.map(pageCard).join('');
    const gloss = GLOSSARY.map(g=>`<div class="hp-gl"><dt>${esc(g[0])}</dt><dd>${esc(g[1])}</dd></div>`).join('');
    const definitions = DEFINITIONS.map(definitionGroup).join('');
    const news = NEWS.map(n=>`<li><b>${esc(n[0])}</b> ${esc(n[1])}</li>`).join('');
    const faq = FAQ.map(f=>`<details class="hp-faq"><summary>${esc(f[0])}<span class="hp-chev">${ico('caret')}</span></summary><div>${esc(f[1])}</div></details>`).join('');

    page.innerHTML = `
      <section class="hp-hero">
        <div class="hp-eyebrow">SNPMaize · SNPVersity 2.1 · B73 RefGen v5</div>
        <h1>Help &amp; FAQ</h1>
        <p>SNPMaize is a SNP toolkit to explore variant diversity across maize germplasm. Everything starts from a genomic
           query and flows between tools without re-running it — this page explains each tool, the datasets behind
           them, and the vocabulary you\u2019ll meet along the way.</p>
        <div class="hp-jump">
          <a href="#hp-new">What’s new</a>
          <a href="#hp-flow">How it works</a>
          <a href="#hp-tools">The tools</a>
          <a href="#hp-data">Datasets</a>
          <a href="#hp-definitions">Definitions</a>
          <a href="#hp-gloss">Scores &amp; annotations</a>
          <a href="#hp-faq">FAQ</a>
        </div>
      </section>

      <section id="hp-new" class="hp-sec">
        <div class="hp-h"><span class="hp-n">01</span><h2>What’s new</h2></div>
        <p class="hp-lead">Changes of October 2026.</p>
        <ul class="hp-news">${news}</ul>
      </section>

      <section id="hp-flow" class="hp-sec">
        <div class="hp-h"><span class="hp-n">02</span><h2>How the suite fits together</h2></div>
        <div class="hp-flow">
          <div class="hp-step"><span class="hp-si" style="background:#2563eb">${ico('search')}</span>
            <b>Query</b><p>In SNPVersity, provide a region or gene and the accessions you care about, or choose them from published results using GWAS Explorer.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#1f8a4c">${ico('table')}</span>
            <b>Result</b><p>You get a genotype table and a VCF for exactly that slice of the genome.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#b45309">${ico('compare')}</span>
            <b>Send onward</b><p>Hand the same matrix to any other tool with "Send selection to…", or press Load data from SNPVersity in the tool — no re-query.</p></div>
        </div>
        <p class="hp-fine">A query is run against the real variant store and returned as a VCF; every other tool reuses that
          one result, so a set you build once stays consistent as you rank it, compare it, cluster it, or draw it.</p>
      </section>

      <section id="hp-tools" class="hp-sec">
        <div class="hp-h"><span class="hp-n">03</span><h2>The tools</h2></div>
        <p class="hp-lead">Fourteen tools. Twelve are live today; SNPImpute and SNPGermplasm are on the roadmap and will
          share the same data and coordinates. Expand any tool for what it does, what it takes, and what it returns.</p>
        <div class="hp-tools">${toolCards}</div>
      </section>

      <section id="hp-data" class="hp-sec">
        <div class="hp-h"><span class="hp-n">04</span><h2>Datasets</h2></div>
        <p class="hp-lead">Each query runs against one dataset. All are called against the B73 v5 reference; they differ
          in how they were filtered, how many accessions and sites they hold, and which score columns they carry.</p>
        ${datasetTable()}
        ${''/* maize2026Details(): the MaizeGDB 2026 composition section, kept below for when that set is offered again */}
      </section>

      <section id="hp-definitions" class="hp-sec">
        <div class="hp-h"><span class="hp-n">05</span><h2>Definitions &amp; table columns</h2></div>
        <p class="hp-lead">Definitions are grouped by tool and use the same labels shown in the interfaces. Expand a group or use your browser's find command to locate a column heading.</p>
        <div class="hp-defgroups">${definitions}</div>
      </section>

      <section id="hp-gloss" class="hp-sec">
        <div class="hp-h"><span class="hp-n">06</span><h2>Scores &amp; annotations</h2></div>
        <dl class="hp-gloss">${gloss}</dl>
      </section>

      <section id="hp-faq" class="hp-sec">
        <div class="hp-h"><span class="hp-n">07</span><h2>Frequently asked</h2></div>
        <div class="hp-faqs">${faq}</div>
      </section>

      <section class="hp-foot">
        <div>
          <b>Still stuck?</b>
          <p>SNPMaize is part of MaizeGDB, the Maize Genetics and Genomics Database.</p>
        </div>
        <div class="hp-foot-btns">
          <button class="hp-open" onclick="go('snpversity')">Start in SNPVersity ${ico('caret')}</button>
          <a class="hp-open ghost" href="https://maizegdb.org" target="_blank" rel="noopener">Visit MaizeGDB ${ico('caret')}</a>
        </div>
      </section>`;

    // jump links: smooth scroll within the .scroll container
    page.querySelectorAll('.hp-jump a').forEach(a=>{
      a.addEventListener('click', e=>{
        const t=document.querySelector(a.getAttribute('href')); if(!t) return;
        e.preventDefault(); t.scrollIntoView({behavior:'smooth', block:'start'});
      });
    });
  }

  /* Entry point used by the MaizeGDB masthead link. Mirrors go() but doesn't
     assume 'help' is present in the TOOLS map, so core.js needs no changes. */
  function open(){
    if (typeof S!=='undefined') S.tool='help';
    if (typeof closeMenu==='function') closeMenu();
    if (typeof toggleRail==='function') toggleRail(false);
    const sb=document.getElementById('sendBtn'); if(sb) sb.style.display='none';
    if (typeof renderNav==='function') renderNav();
    document.querySelectorAll('.mg-help').forEach(el=>el.classList.add('on'));
    window.scrollTo(0,0);
    const sc=document.querySelector('.scroll'); if(sc) sc.scrollTop=0;
    render(document.getElementById('page'));
  }

  function injectCSS(){
    if (document.getElementById('snphelp-css')) return;
    const s=document.createElement('style'); s.id='snphelp-css';
    s.textContent = `
      .hp-hero{padding:6px 0 18px;border-bottom:1px solid var(--line,#e6e9ef);margin-bottom:26px}
      .hp-eyebrow{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--blue-600,#2563eb);margin-bottom:10px}
      .hp-hero h1{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:40px;line-height:1.05;letter-spacing:-.02em;margin:0 0 12px;color:var(--ink,#141922)}
      .hp-hero p{max-width:70ch;color:var(--muted,#5b6b83);font-size:15px;line-height:1.6;margin:0 0 18px}
      .hp-jump{display:flex;gap:8px;flex-wrap:wrap}
      .hp-jump a{font-size:12.5px;font-weight:600;color:var(--ink,#141922);text-decoration:none;padding:7px 13px;border:1px solid var(--line,#e6e9ef);border-radius:999px;background:#fff;transition:.15s}
      .hp-jump a:hover{border-color:var(--blue-600,#2563eb);color:var(--blue-600,#2563eb)}

      .hp-sec{margin:0 0 42px;scroll-margin-top:16px}
      .hp-h{display:flex;align-items:baseline;gap:12px;margin-bottom:16px}
      .hp-h h2{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:23px;letter-spacing:-.01em;margin:0;color:var(--ink,#141922)}
      .hp-n{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12px;font-weight:600;color:#c8a24a;border:1px solid #e7d6a8;background:#fdf8ec;border-radius:6px;padding:3px 7px}
      .hp-lead{max-width:74ch;color:var(--muted,#5b6b83);font-size:14px;line-height:1.6;margin:0 0 16px}
      .hp-fine{max-width:78ch;color:var(--muted,#5b6b83);font-size:12.5px;line-height:1.6;margin:12px 0 0}
      .hp-fine b,.hp-lead b{color:var(--ink,#141922)}
      .hp-news{max-width:80ch;margin:0;padding:0 0 0 18px;font-size:13.5px;line-height:1.6;color:var(--muted,#5b6b83)}
      .hp-news li{margin:0 0 7px}
      .hp-news b{color:var(--ink,#141922)}

      /* flow */
      .hp-flow{display:flex;align-items:stretch;gap:10px;flex-wrap:wrap}
      .hp-step{flex:1 1 210px;border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff;padding:16px}
      .hp-step b{display:block;font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15px;color:var(--ink,#141922);margin:10px 0 5px}
      .hp-step p{margin:0;font-size:13px;line-height:1.55;color:var(--muted,#5b6b83)}
      .hp-si{display:inline-flex;width:34px;height:34px;border-radius:9px;align-items:center;justify-content:center;color:#fff}
      .hp-si svg{width:19px;height:19px}
      .hp-arrow{display:flex;align-items:center;color:var(--faint,#aab4c4)}
      .hp-arrow svg{width:22px;height:22px}

      /* tool accordions */
      .hp-tools{display:flex;flex-direction:column;gap:10px}
      .hp-tool{border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff;overflow:hidden}
      .hp-tool[open]{border-color:#cdd6e6;box-shadow:0 1px 0 rgba(20,25,34,.03)}
      .hp-tool summary{display:flex;align-items:center;gap:13px;padding:14px 16px;cursor:pointer;list-style:none}
      .hp-tool summary::-webkit-details-marker{display:none}
      .hp-ti{display:inline-flex;width:34px;height:34px;border-radius:9px;align-items:center;justify-content:center;color:#fff;flex:0 0 auto}
      .hp-ti svg{width:19px;height:19px}
      .hp-tt{display:flex;flex-direction:column;gap:2px;min-width:0}
      .hp-tt b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15.5px;color:var(--ink,#141922)}
      .hp-tag{font-size:12.5px;color:var(--muted,#5b6b83)}
      .hp-status{margin-left:auto;font-size:10.5px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;padding:3px 9px;border-radius:999px;border:1px solid}
      .hp-status.on{color:#1f8a4c;border-color:#bfe3cc;background:#f0faf3}
      .hp-status.off{color:#b06a12;border-color:#f0dcbb;background:#fdf6ec}
      .hp-chev{color:var(--faint,#aab4c4);transition:transform .18s;flex:0 0 auto}
      .hp-chev svg{width:18px;height:18px;transform:rotate(90deg)}
      .hp-tool[open] .hp-chev svg{transform:rotate(-90deg)}
      .hp-tbody{padding:2px 16px 18px 63px}
      .hp-tbody>p{margin:0 0 14px;font-size:13.5px;line-height:1.62;color:#3a465a;max-width:76ch}
      .hp-io{display:flex;gap:24px;flex-wrap:wrap;padding:13px 15px;background:var(--blue-50,#f2f6ff);border:1px solid #e4ecfb;border-radius:10px}
      .hp-io>div{flex:1 1 240px;font-size:13px;line-height:1.5;color:#3a465a}
      .hp-lbl{display:block;font-size:10px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--blue-600,#2563eb);margin-bottom:4px}
      .hp-actions{margin-top:14px}
      .hp-open{display:inline-flex;align-items:center;gap:6px;font-family:inherit;font-size:13px;font-weight:600;color:#fff;background:var(--ink,#141922);border:1px solid var(--ink,#141922);border-radius:9px;padding:9px 15px;cursor:pointer;text-decoration:none;transition:.15s}
      .hp-open:hover{background:#000}
      .hp-open svg{width:15px;height:15px}
      .hp-open.ghost{background:#fff;color:var(--ink,#141922)}
      .hp-open.ghost:hover{border-color:var(--blue-600,#2563eb);color:var(--blue-600,#2563eb)}
      .hp-open.muted{background:#f4f6fa;border-color:var(--line,#e6e9ef);color:var(--muted,#5b6b83);cursor:default}

      /* dataset table */
      .hp-tablewrap{overflow-x:auto;border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff}
      .hp-table{border-collapse:collapse;width:100%;font-size:13px;min-width:720px}
      .hp-table th{text-align:left;font-size:10.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted,#5b6b83);font-weight:700;padding:11px 14px;border-bottom:1px solid var(--line,#e6e9ef);background:#fafbfd;white-space:nowrap}
      .hp-table td{padding:12px 14px;border-bottom:1px solid #eef1f6;vertical-align:top;color:#3a465a}
      .hp-table tr:last-child td{border-bottom:none}
      .hp-sub{font-size:11.5px;color:var(--muted,#5b6b83);margin-top:2px;font-weight:400}
      .hp-mono{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12px;white-space:nowrap}
      .hp-num{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12.5px;white-space:nowrap;color:var(--ink,#141922)}
      .hp-c{text-align:center}
      .hp-chip{display:inline-block;font-size:11px;color:#3a465a;background:#f1f4f9;border:1px solid #e4e9f1;border-radius:6px;padding:2px 7px;margin:2px 2px 0 0;white-space:nowrap}
      .hp-yes{color:#1f8a4c;font-weight:600}
      .hp-no{color:var(--faint,#aab4c4)}


      /* complete definitions */
      .hp-defgroups{display:flex;flex-direction:column;gap:9px}
      .hp-defgroup{border:1px solid var(--line,#e6e9ef);border-radius:11px;background:#fff;overflow:hidden}
      .hp-defgroup[open]{border-color:#cdd6e6}
      .hp-defgroup summary{display:flex;align-items:center;gap:10px;padding:13px 15px;cursor:pointer;list-style:none}
      .hp-defgroup summary::-webkit-details-marker{display:none}
      .hp-defgroup summary b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15px;color:var(--ink,#141922)}
      .hp-defdot{width:10px;height:10px;border-radius:50%;flex:0 0 auto}
      .hp-count{font-size:11.5px;color:var(--muted,#5b6b83);margin-left:2px}
      .hp-defgroup summary .hp-chev{margin-left:auto}
      .hp-deflist{margin:0;padding:0 15px 8px 35px;display:grid;grid-template-columns:repeat(auto-fit,minmax(min(330px,100%),1fr));gap:0 28px}
      .hp-def{padding:11px 0;border-top:1px solid #eef1f6}
      .hp-def dt{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:13.5px;font-weight:600;color:var(--ink,#141922);margin-bottom:3px}
      .hp-def dd{margin:0;font-size:12.8px;line-height:1.55;color:var(--muted,#5b6b83);overflow-wrap:anywhere}


      /* glossary */
      .hp-gloss{margin:0;display:grid;grid-template-columns:repeat(auto-fit,minmax(min(310px,100%),1fr));gap:2px 26px}
      .hp-gl{padding:13px 0;border-bottom:1px solid #eef1f6}
      .hp-gl dt{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:14px;font-weight:600;color:var(--ink,#141922);margin-bottom:4px}
      .hp-gl dd{margin:0;font-size:13px;line-height:1.55;color:var(--muted,#5b6b83)}

      /* faq */
      .hp-faqs{display:flex;flex-direction:column;gap:8px}
      .hp-faq{border:1px solid var(--line,#e6e9ef);border-radius:10px;background:#fff}
      .hp-faq summary{display:flex;align-items:center;gap:12px;padding:14px 16px;cursor:pointer;font-family:var(--disp,'Space Grotesk',sans-serif);font-size:14.5px;font-weight:600;color:var(--ink,#141922);list-style:none}
      .hp-faq summary::-webkit-details-marker{display:none}
      .hp-faq .hp-chev{margin-left:auto}
      .hp-faq>div{padding:0 16px 16px;font-size:13.5px;line-height:1.62;color:#3a465a;max-width:80ch}

      /* footer */
      .hp-foot{display:flex;align-items:center;gap:18px;flex-wrap:wrap;justify-content:space-between;border:1px solid var(--line,#e6e9ef);border-radius:14px;background:linear-gradient(180deg,#fff,#fbfcfe);padding:20px 22px;margin-bottom:24px}
      .hp-foot b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:16px;color:var(--ink,#141922)}
      .hp-foot p{margin:4px 0 0;font-size:13px;color:var(--muted,#5b6b83)}
      .hp-foot-btns{display:flex;gap:10px;flex-wrap:wrap}

      @media (max-width:620px){
        .hp-hero h1{font-size:31px}
        .hp-tbody{padding-left:16px}
        .hp-arrow{transform:rotate(90deg);align-self:center}
      }`;
    document.head.appendChild(s);
  }

  /* Public lookup so other tools can pull their own definitions subset
     into an in-page popover instead of duplicating the text — e.g. GWAS
     Explorer's significance-threshold and region-table "?" buttons. */
  function definitionsFor(toolName) {
    const g = DEFINITIONS.find(d => d.tool === toolName);
    return g ? g.items : [];
  }

  if (typeof SNPTools!=='undefined') SNPTools.register('help', { render });
  return { render, open, definitionsFor };
})();
if (typeof window!=='undefined'){ window.SNPHelp = SNPHelp; window.openHelp = SNPHelp.open; }
