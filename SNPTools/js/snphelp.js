/* =====================================================================
 *  snphelp.js — Help & FAQ page for the SNPTools suite.
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
      what:'The front door of the suite and the starting point for most work. Choose a dataset, type a genomic interval (or an FGSG / FVEG / FVERT4 gene model ID), and pick the accessions you want. SNPVersity queries the variant store and returns a color-coded genotype table plus a downloadable VCF — allele states, predicted effects, and DNA/protein language-model scores included.',
      give:'A dataset, a region or gene, and a set of accessions.',
      get:'A genotype table and a VCF. From here, "Send selection to…" hands the same result to any other tool.' },
    { id:'snptrait', name:'SNPTrait', icon:'leaf', color:'#1f8a4c', status:LIVE,
      tag:'Strain Selector — pick isolates by metadata',
      what:'An interactive catalogue of the isolates in a dataset, grouped by population structure. Filter by population, species, host, chemotype, or country (rich metadata is available for F. graminearum), batch-select, then send the selected isolates straight to SNPVersity to build a VCF.',
      give:'A dataset and metadata filters over its isolates.',
      get:'A selected set of isolates handed to SNPVersity, also exportable as CSV / JSON.' },
    { id:'snpimpact', name:'SNPImpact', icon:'star', color:'#7c3aed', status:LIVE,
      tag:'Rank candidate variants',
      what:'Prioritizes the variants in a region regardless of which accessions you picked. It orders them by AI variant-effect scores (DNA models FunDLM and EVO2, protein models ESM1/ESM2/ESM3/ESM-C) combined with predicted consequence and Pfam domain annotation, so likely causal changes rise to the top. Filter by consequence, priority, score, or domain, and flag a shortlist.',
      give:'A region — queried directly in SNPImpact (reference + interval or a gene model), or sent over from SNPVersity. Accessions are ignored here.',
      get:'A ranked, filterable variant table and a shortlist of candidate alleles.' },
    { id:'snpfunction', name:'SNPFunction', icon:'func', color:'#2563eb', status:LIVE,
      tag:'Gene function & allele mining',
      what:'A gene-scoped dossier, independent of any one region. It opens with a functional-annotation record — gene symbol and names, a sourced description, protein domain architecture drawn to scale, Gene Ontology grouped by aspect, KEGG pathways and orthology, and external cross-references, each tagged by provenance (curated vs predicted). Below that it computes the gene\u2019s variant burden across the whole panel and lists which accessions carry each damaging or knockout allele.',
      give:'A gene model ID and a dataset.',
      get:'A functional-annotation dossier (identity, domains, GO, pathways, cross-references), a variant-burden breakdown, and a damaging-allele catalog with carrier lines.' },
    { id:'snpcompare', name:'SNPCompare', icon:'compare', color:'#0e7490', status:LIVE,
      tag:'Similarity to a focal accession',
      what:'Ranks every accession by identity-by-state similarity to one focal accession. Similarity is the fraction of co-called sites with the same genotype. Works genome-wide (precomputed), for the current region (computed in-browser from a SNPVersity result), or both side-by-side with a delta that surfaces region-specific relatedness such as introgression.',
      give:'A focal accession, and optionally a region handed off from SNPVersity.',
      get:'A ranked similarity table with project, SRA ID, and accession name.' },
    { id:'snptree', name:'SNPTree', icon:'tree', color:'#15803d', status:LIVE,
      tag:'Local phylogeny',
      what:'Builds a local phylogenetic tree from the genotype matrix already in memory, using identity-by-state distances (UPGMA / neighbour-joining). Useful for reading haplotype structure, introgression, and how lines cluster in a region.',
      give:'A SNPVersity result set (up to 250 accessions before a warning).',
      get:'An interactive tree, downloadable as Newick, MEGA, or PHYLIP.' },
    { id:'snpfold', name:'SNPFold', icon:'fold', color:'#be185d', status:LIVE,
      tag:'Variants on protein structure',
      what:'Maps coding variants onto predicted protein structure. A linear protein browser aligns variants with Pfam domains, secondary structure, and per-residue pLDDT confidence; an on-demand 3D viewer shows the fold colored by confidence, domain, or impact; and a per-variant readout interprets each change (domain, local confidence, secondary structure, predicted \u0394\u0394G).',
      give:'A gene with an available structure model.',
      get:'A structure-aware, per-variant interpretation of coding changes.' },
    { id:'paneffect', name:'PanEffect', icon:'effect', color:'#b45309', status:LIVE,
      tag:'Missense effects across a protein',
      what:'Visualizes the predicted effect of every possible amino-acid substitution across a protein using ESM protein-language-model scores, aligned to Pfam domains and predicted secondary structure. The PanEffect engine has not yet been ported to the Fusarium references and is on the roadmap.',
      give:'A gene model and an ESM model.',
      get:'A full-length substitution-effect heatmap with domain and secondary-structure context.' },
  ];

  /* ---- glossary ---- */
  const GLOSSARY = [
    ['Reference genomes', 'Coordinates, gene models, and annotations are anchored to one of three Fusarium references: F. graminearum PH-1 (FGSG, chr1-chr4, 36.4 Mb, 512 isolates), F. verticillioides 7600 (FVEG, chr1-chr11, 41.1 Mb, 113 isolates), and F. verticillioides MRC826 (FVERT4, chr1-chr12, 42.9 Mb, 113 isolates).'],
    ['Genotype (haploid)', 'Fusarium isolates are haploid, so each call is a single allele: 0 (reference), a non-zero index (an alternate allele), or missing (shown as N). There is no heterozygous state. Missing calls are left out of a comparison rather than counted as a match.'],
    ['IBS distance / % identity', 'Identity-by-state compares two accessions site by site over the calls they share. Distance is the mean allele difference; % identity is 100 − distance. Used by SNPTree, SNPMatrix, and SNPCompare.'],
    ['Predicted effect & impact', 'Each variant carries a predicted consequence (missense, LOF, splice, indel, synonymous, …) rolled up into an impact tier: HIGH, MODERATE, LOW, or MODIFIER, most severe wins when a site lists several.'],
    ['DNA scores (FunDLM, EVO2)', 'Two DNA language-model predictions of how disruptive a nucleotide change is, from the VCF FUNDLM_SCORE and EVO2_SCORE fields. Both abstain on INDELs (shown as N/A). Present for all three references.'],
    ['Protein scores (ESM1, ESM2, ESM3, ESM-C)', 'Four protein language-model predictions of the effect of an amino-acid change, from the VCF ESM1_SCORE / ESM2_SCORE / ESM3_SCORE / ESMC_SCORE fields. Populated for coding (missense) variants across all three references.'],
    ['MAF', 'Minor allele frequency, the frequency of the less common allele. SNPVersity can filter a region by a minimum MAF.'],
    ['Pfam domain', 'When a variant falls inside a known protein domain, that domain is shown and linked to InterPro. Domain annotation is still being loaded for some regions, where it reads as \u2014.'],
    ['VCF', 'The Variant Call Format file SNPVersity generates for your query. It is the exact matrix the other tools reuse when you send a selection.'],
    ['Gene Ontology (GO)', 'Standardized terms describing a gene product\u2019s biological process, molecular function, and cellular component. SNPFunction groups them by aspect and tags each by source: UniProt is curated; InterPro2GO is predicted from domains.'],
    ['Curated vs predicted', 'Annotation provenance. Curated evidence (UniProt and curated sources) is human-reviewed; predicted evidence (InterPro / InterPro2GO) is inferred from protein domains. SNPFunction labels each description and GO term accordingly and can hide predicted-only terms.'],
    ['Population structure', 'F. graminearum isolates are grouped into population-structure clusters: NA1, NA2, NA3, Admixture, Outgroups, and Unknown. SNPTrait and the accession picker group and colour isolates by these.'],
  ];


  /* ---- complete definitions by tool ----
   * Keep table-header wording identical to the live tools so users can
   * search this page for the label they see on screen.
   */
  const DEFINITIONS = [
    { tool:'Shared terms', color:'#64748b', items:[
      ['Accession', 'A Fusarium isolate (or sequencing run) represented by one genotype column.'],
      ['Allele', 'One observed DNA state at a genomic position. REF is the reference-genome allele; ALT is an alternate allele.'],
      ['Variant', 'A genomic position where at least one isolate differs from the reference genome.'],
      ['SNP', 'Single-nucleotide polymorphism: a one-base substitution.'],
      ['INDEL', 'Insertion or deletion relative to the reference sequence.'],
      ['Gene model', 'The reference gene-model identifier (FGSG / FVEG / FVERT4) and its annotated exon, intron, CDS, and transcript structure.'],
      ['Consequence / Effect', 'The predicted molecular result of a variant, such as synonymous, missense, splice-site, frameshift, stop gained, or intronic.'],
      ['Impact', 'A broad severity class assigned from the consequence: HIGH, MODERATE, LOW, or MODIFIER.'],
      ['Priority', 'The SNPTools ranking tier (TOP, HIGH, MODERATE, or LOW) produced by combining consequence, model scores, and domain context.'],
      ['Domain', 'A Pfam-annotated protein domain overlapping the affected residue; — means no loaded domain hit.'],
      ['Het', 'Heterozygous state — not applicable to Fusarium. Isolates are haploid, so each genotype is a single allele (0 or an alternate index). The SNPFunction Het column reads 0 for these datasets.'],
      ['Hom', 'For haploid Fusarium isolates this is simply the count of isolates carrying an alternate allele (there is no true homozygous state).'],
      ['Missing / ./.', 'No usable genotype call (shown as N). Missing calls are excluded from pairwise similarity calculations.'],
      ['AF', 'Alternate-allele frequency: the frequency of the ALT allele among called chromosomes in the analyzed panel.'],
      ['MAF', 'Minor-allele frequency: the frequency of the less common allele, constrained to 0–0.5.'],
      ['Carrier', 'An accession with at least one copy of the alternate allele.'],
      ['Co-called sites', 'Sites where both accessions in a pair have non-missing genotype calls.'],
      ['Reference genome', 'The Fusarium reference assembly (PH-1, 7600, or MRC826) used for coordinates, REF alleles, gene models, and annotations throughout SNPTools.'],
    ]},
    { tool:'SNPVersity', color:'#2563eb', items:[
      ['CHR', 'Reference chromosome containing the variant.'],
      ['POS', 'One-based genomic coordinate on the reference genome.'],
      ['REF', 'Reference allele in the reference genome.'],
      ['ALT', 'Alternate allele represented by the row.'],
      ['Gene model', 'Reference gene model overlapping or associated with the variant.'],
      ['Effect', 'Predicted variant consequence from the annotation source.'],
      ['Impact', 'Predicted severity category: HIGH, MODERATE, LOW, or MODIFIER.'],
      ['Domain', 'Pfam protein domain overlapping the affected coding residue, when available.'],
      ['MQ', 'Mapping quality: a phred-scaled measure of confidence that reads were aligned to the correct genomic location; higher is better.'],
      ['COMP', 'Completeness: the proportion of accessions with a non-missing genotype call at that site.'],
      ['maxR²', 'Maximum linkage-disequilibrium r² used by the dataset filter or imputation-quality workflow; values closer to 1 indicate stronger correlation.'],
      ['MAF', 'Minor-allele frequency among the selected or source accessions, depending on the returned record.'],
      ['FunDLM / EVO2', 'The two DNA language-model variant scores (VCF FUNDLM_SCORE / EVO2_SCORE). More extreme scores are more disruptive; both abstain on INDELs.'],
      ['ESM1 / ESM2 / ESM3 / ESM-C', 'The four protein language-model scores for amino-acid substitutions (VCF ESM1_SCORE / ESM2_SCORE / ESM3_SCORE / ESMC_SCORE), populated for coding variants.'],
      ['Isolate genotype columns', 'Each isolate column shows its haploid genotype at the site: 0 reference (green), a non-zero allele index (orange), or N missing (grey).'],
      ['Dataset', 'A defined variant collection with its own accession panel, filters, included variant types, and score columns.'],
      ['Sites', 'Number of variant positions in the complete dataset, not necessarily the number returned by the current query.'],
      ['Imputed', 'Whether missing genotypes were statistically inferred in that dataset.'],
    ]},
    { tool:'SNPImpact', color:'#7c3aed', items:[
      ['Gene', 'Gene model associated with the candidate variant.'],
      ['Variant', 'Genomic change, generally shown as position and REF→ALT alleles.'],
      ['Consequence', 'Specific predicted molecular consequence of the change.'],
      ['Domain', 'Pfam domain containing the affected amino acid, when present.'],
      ['FunDLM / EVO2', 'The two DNA language-model scores used to estimate sequence disruption (abstain on INDELs).'],
      ['ESM1 / ESM2 / ESM3 / ESM-C', 'The four protein language-model scores used to estimate the effect of an amino-acid substitution.'],
      ['Priority', 'Integrated candidate tier. TOP is the strongest prioritization, followed by HIGH, MODERATE, and LOW.'],
      ['Shortlist / flag', 'A user-selected marker for retaining a candidate variant for later review or export.'],
      ['Gene-model diagram', 'A compact display of exons, introns, coding sequence, strand, and the variant position.'],
      ['Exon', 'A transcript segment retained in the mature RNA; coding portions contribute to the protein sequence.'],
      ['Intron', 'A transcribed segment removed during RNA splicing.'],
      ['CDS', 'Coding sequence: the portion of exons translated into protein.'],
      ['Imputation r²', 'Estimated squared correlation between imputed and true genotypes; values nearer 1 indicate greater confidence.'],
      ['MQ', 'Mapping quality for the variant site.'],
      ['MAF', 'Minor-allele frequency for the candidate variant.'],
    ]},
    { tool:'SNPFunction', color:'#2563eb', items:[
      ['Non-syn : syn', 'The number or ratio of nonsynonymous coding variants to synonymous coding variants. A higher value indicates more amino-acid-changing variation relative to silent variation; it is descriptive and is not by itself a formal dN/dS estimate.'],
      ['Exon : intron', 'The number or ratio of variants in annotated exons to variants in introns of the gene model. Infinity (∞) means exon variants were observed but no intron variants were counted.'],
      ['Domain-disrupting', 'Coding variants that alter an amino acid located inside an annotated Pfam domain.'],
      ['Knockout alleles', 'Alleles predicted to strongly disrupt gene function, such as frameshift, stop-gained, essential splice, or other loss-of-function changes.'],
      ['Mean FunDLM / Mean EVO2', 'Average DNA language-model score across the gene variants included in the burden summary.'],
      ['Mean ESM1 / ESM2 / ESM3 / ESM-C', 'Average protein language-model score across scored amino-acid-changing variants in the gene.'],
      ['Allele', 'The specific genomic REF→ALT change represented by a damaging-allele row.'],
      ['Consequence', 'Predicted molecular effect of that allele.'],
      ['Domain', 'Pfam domain overlapping the affected residue.'],
      ['FunDLM / EVO2', 'DNA language-model scores for the allele.'],
      ['ESM1 / ESM2 / ESM3 / ESM-C', 'Protein language-model scores for the resulting amino-acid change.'],
      ['Priority', 'Integrated SNPTools evidence tier for the allele.'],
      ['Het', 'Number of accessions carrying the allele heterozygously.'],
      ['Hom', 'Number of accessions carrying the allele as alternate homozygous.'],
      ['AF', 'Alternate-allele frequency across the whole analyzed panel.'],
      ['Variant burden', 'The count and composition of variants assigned to the gene across the full dataset panel.'],
      ['Damaging allele', 'An allele selected because its consequence and/or prediction scores indicate a potentially important functional effect.'],
      ['Functional annotation', 'A per-gene dossier: identity, description, protein domains, Gene Ontology, pathways, and cross-references, where the reference annotation provides them.'],
      ['Gene symbol / aliases', 'The primary gene symbol and any additional names or synonyms recorded for the model.'],
      ['Description source', 'Provenance of the functional description, shown as a badge: curated (UniProt and curated sources), InterPro (predicted from domains), or no informative source.'],
      ['Evidence chips', 'An at-a-glance row of what is annotated for the gene: GO, Pfam, KEGG KO, Pathway, UniProt, Symbol.'],
      ['Protein domain architecture', 'A to-scale diagram of the protein with Pfam domains as positioned blocks, plus a list giving Pfam and InterPro IDs, residue span, percent of protein covered, and InterProScan E-value.'],
      ['Gene Ontology (BP / MF / CC)', 'GO terms grouped by aspect — Biological process, Molecular function, Cellular component — each linking out to AmiGO.'],
      ['GO source (UniProt / InterPro2GO)', 'Provenance of each GO term. UniProt and curated sources are human-reviewed; InterPro2GO is predicted from protein domains.'],
      ['Curated only', 'A toggle that hides GO terms supported only by prediction (InterPro2GO), leaving curated evidence.'],
      ['Obsolete term', 'A GO term whose status is no longer current; it is retained but flagged.'],
      ['KEGG orthology (KO) / Pathways', 'KEGG orthology assignments, pathway memberships, and KEGG gene IDs, where the gene cross-references to UniProt or Entrez.'],
      ['Cross-references', 'External identifiers for the gene: UniProt, NCBI Gene, and FungiDB.'],
      ['Annotation build', 'The build date, assembly, and annotation version the functional record was generated from.'],
    ]},
    { tool:'SNPCompare', color:'#0e7490', items:[
      ['#', 'Current rank after sorting and filtering.'],
      ['Project', 'BioProject or dataset project associated with the accession.'],
      ['SRA ID', 'Sequence Read Archive run identifier associated with the accession.'],
      ['Accession name', 'Human-readable line or germplasm name.'],
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
    ]},
    { tool:'SNPTree', color:'#15803d', items:[
      ['IBS allele distance', 'Mean pairwise allele difference across co-called sites. For haploid isolates, two isolates with the same allele contribute 0 and different alleles contribute 1 (there is no half-step heterozygous comparison).'],
      ['Informative sites', 'Variant sites that contain more than one observed genotype state among the selected accessions.'],
      ['Shared sites', 'Sites with non-missing calls for both accessions in a pair.'],
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
      ['Local pLDDT', 'Predicted-structure confidence near the affected residue, from 0 to 100, from the loaded model (Boltz2 or ESMFold). Higher values indicate greater confidence in the local fold.'],
      ['Structure', 'DSSP-style secondary-structure assignment at the residue: α-helix, β-strand, or loop/coil.'],
      ['FunDLM / EVO2', 'DNA language-model scores for the underlying nucleotide variant.'],
      ['ESM1 / ESM2 / ESM3 / ESM-C', 'Protein language-model scores for the amino-acid substitution.'],
      ['IUPred2', 'Predicted intrinsic disorder at the residue, from 0 to 1 (IUPred2). Higher values mean the residue is more likely to be intrinsically disordered than to adopt a fixed fold.'],
      ['Anchor2', 'Predicted probability (0 to 1) that a disordered residue lies in a protein-binding region (ANCHOR2). Higher values flag disordered segments likely to become ordered upon binding a partner.'],
      ['Activity', 'Annotated functional site overlapping the residue, drawn from InterProScan member databases (for example an active site, binding site, or conserved functional feature), when present.'],
      ['Priority', 'Integrated evidence tier for the structure-mapped variant.'],
      ['Carriers', 'Number of accessions carrying the alternate allele; expanded details separate heterozygous and homozygous carriers.'],
      ['pLDDT', 'Predicted Local Distance Difference Test score from AlphaFold. Common interpretation: ≥90 very high confidence, 70–89 confident, 50–69 low, and <50 very low.'],
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
      ['Reference View', 'The substitution heatmap computed on the reference protein: every residue (columns) against all 20 possible amino acids (rows).'],
      ['Pan-genome View', 'The same substitution heatmap projected onto a protein multiple-sequence alignment across related assemblies, so natural variation and its predicted effect can be read together. Available for canonical transcripts.'],
      ['ESM score', 'A protein-language-model estimate of how tolerated an amino-acid substitution is; more negative is more disruptive.'],
      ['ESM model (ESM1 / ESM2 / ESM3)', 'The protein language model used to score substitutions, selectable in the panel.'],
      ['Substitution heatmap', 'A per-residue grid colored by predicted effect. Hover a cell to read the position, wild-type → substitution, and score.'],
      ['Zoomed region', 'A 50-residue window of the full heatmap, positioned with the slider, showing per-cell substitution letters and the wild-type residue track.'],
      ['Observed-variants view', 'Colors the reference heatmap using substitutions observed in the dataset. Hand-offs from SNPVersity or SNPFold open in this view; "Show all variant effects" colors every possible substitution.'],
      ['PFAM Domains track', 'Pfam domains drawn to scale along the protein and linked to InterPro, shared with the SNPFold and SNPFunction views.'],
      ['Secondary structure', 'Predicted per-residue secondary structure (helix / strand / coil) drawn above the heatmap for context.'],
      ['Species / Gene model', 'In a pan-genome view, switches the row labels between assembly (species) names and their gene-model IDs.'],
      ['Row grouping color', 'In a pan-genome view, rows can be colored by grouping metadata for the assemblies.'],
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
     'Run a query in SNPVersity, then use "Send selection to…" in the top bar (or the buttons on a result). The same genotype matrix is handed to SNPImpact, SNPCompare, SNPTree, SNPMatrix, and SNPTrait without re-querying. SNPMatrix and SNPTree can also pass their set on to each other.'],
    ['Why did a wide region give me a download instead of a table?',
     'Regions larger than one million bases skip the in-browser table and return a downloadable VCF instead, so the page stays responsive. Narrow the interval to get the interactive table back.'],
    ['Why are some cells \u2014 or N/A?',
     'A \u2014 in a domain column means Pfam annotation for that position hasn\u2019t been loaded yet. The DNA scores (FunDLM, EVO2) read N/A on INDELs, and the protein scores (ESM1/ESM2/ESM3/ESM-C) are only defined for coding (missense) changes, so they are blank on non-coding or intergenic sites. MQ and coverage read N/A when the source didn\u2019t record them.'],
    ['Is there a limit on how many accessions I can compare?',
     'The distance tools warn before doing very large computations in the browser: SNPTree above 250 accessions and SNPMatrix above 400. You can build anyway, it just may be slow. SNPImpact renders up to 1,500 variants at a time.'],
    ['What can I download?',
     'A VCF from SNPVersity; a CSV distance matrix, PHYLIP, PNG, and SVG from SNPMatrix; Newick, MEGA, and PHYLIP trees from SNPTree. Comparison and impact tables can be exported from their own pages.'],
    ['Which tools are ready to use now?',
     'SNPVersity, SNPImpact, SNPFunction, SNPCompare, SNPTree, SNPMatrix, SNPFold, and SNPTrait are live. PanEffect, SNPImpute, SNPDensity, and SNPGermplasm are on the roadmap and marked in development in the sidebar.'],
  ];

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
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
      <p class="hp-fine">Isolates are grouped by reference genome — F. graminearum (FGSG), F. verticillioides 7600 (FVEG), and F. verticillioides MRC826 (FVERT4) — and, for F. graminearum, by population structure (NA1, NA2, NA3, Admixture, Outgroups, Unknown). Every dataset carries six AI variant-effect scores: two DNA language models (FunDLM, EVO2) and four protein language models (ESM1, ESM2, ESM3, ESM-C). Site counts are the complete per-genome totals across all chromosomes.</p>`;
  }


  function fmtInt(n){ return Number(n).toLocaleString('en-US'); }

  function fusariumDatasetDetails(){
    return `
        <br>
      <details class="hp-data-detail" open>
        <summary><b>Fusarium datasets: references, populations, and processing</b><span class="hp-chev">${ico('caret')}</span></summary>
        <div class="hp-data-body">
          <p>SNPTools for Fusarium serves variant data for three reference genomes. Pick the reference that matches the
             isolates you want to study; gene-model IDs and coordinates are specific to each reference.</p>

          <div class="hp-comparegrid">
            <div class="hp-dcard"><b>F. graminearum PH-1 (2026)</b><p><span class="hp-mono">FGSG_</span> gene models across chr1&ndash;chr4 (36.4 Mb). <b>512 isolates</b>, grouped by population structure (NA1, NA2, NA3, Admixture, Outgroups, Unknown). Variants: <b>2,485,599</b> High Coverage (2,339,005 SNPs + 146,594 INDELs), <b>2,126,508</b> High Quality (111,754 SNPs + 8,636 INDELs).</p></div>
            <div class="hp-dcard"><b>F. verticillioides 7600</b><p><span class="hp-mono">FVEG_</span> gene models across chr1&ndash;chr11 (41.1 Mb). <b>113 isolates</b>. Variants: <b>1,227,771</b> High Coverage (1,139,829 SNPs + 87,942 INDELs), <b>937,497</b> High Quality (872,140 SNPs + 65,357 INDELs).</p></div>
            <div class="hp-dcard"><b>F. verticillioides MRC826</b><p><span class="hp-mono">FVERT4_</span> gene models across chr1&ndash;chr12 (42.9 Mb). <b>113 isolates</b>. Variants: <b>1,402,131</b> High Coverage (1,304,907 SNPs + 97,224 INDELs), <b>1,064,463</b> High Quality (990,183 SNPs + 74,280 INDELs).</p></div>
          </div>

          <h3>High Quality vs High Coverage</h3>
          <div class="hp-comparegrid">
            <div class="hp-dcard"><b>High Coverage (HC)</b><p>The broader filtered set (mapping-quality and genotype-coverage requirements). Retains more variation &mdash; useful for discovery and sensitivity.</p></div>
            <div class="hp-dcard"><b>High Quality (HQ)</b><p>The stricter subset &mdash; the same filters plus a high-confidence linkage-disequilibrium criterion (LD max R&sup2; &gt; 0.5). Useful when specificity and confidence matter most.</p></div>
          </div>
          <p class="hp-fine"><b>F. graminearum</b> retains 2,126,508 HQ of 2,485,599 HC sites (~85%), and the two <b>F. verticillioides</b> genomes retain slightly fewer variation (7600: 937,497 of 1,227,771, ~76%; MRC826: 1,064,463 of 1,402,131, ~76%).</p>

          <h3>Haploid genotypes</h3>
          <p>Fusarium isolates are <b>haploid</b>, so each genotype is a single allele: <b>0</b> is the reference allele
             (green), any non-zero index is an alternate allele (orange), and a missing or uncalled genotype reads as
             <b>N</b> (grey). There is no heterozygous state, so the identity-by-state tools (SNPTree, SNPMatrix,
             SNPCompare) score a site as a full difference when two isolates carry different alleles and as a match
             otherwise.</p>

          <h3>What is stored and displayed</h3>
          <p>The data is maintained as per-chromosome HDF5 stores. SNPVersity sends the selected dataset, genomic
             interval, and isolate list to the server, which extracts the requested slice and returns a VCF. SNPTools
             then reuses that same genotype matrix in SNPImpact, SNPCompare, SNPTree, SNPMatrix, and related pages. Each
             site may include REF and ALT alleles, isolate genotypes, predicted molecular consequence, gene association,
             mapping quality, genotype completeness, linkage-disequilibrium support (maxR&sup2;), minor-allele frequency,
             Pfam domain context when available, and six AI variant-effect scores (two DNA and four protein language models).</p>

          <h3>Language-model scores</h3>
          <p>Every dataset carries <b>six</b> AI variant-effect scores. Two are <b>DNA</b> language models &mdash;
             <b>FunDLM</b> (VCF <span class="hp-mono">FUNDLM_SCORE</span>) and <b>EVO2</b> (<span class="hp-mono">EVO2_SCORE</span>) &mdash;
             which estimate how disruptive a nucleotide change is and both abstain on INDELs (shown as N/A). Four are
             <b>protein</b> language models &mdash; <b>ESM1</b>, <b>ESM2</b>, <b>ESM3</b>, and <b>ESM-C</b>
             (<span class="hp-mono">ESM1_SCORE</span> … <span class="hp-mono">ESMC_SCORE</span>) &mdash; scored for coding
             (missense) substitutions. SNPImpact/SNPFold/SNPFunction combine the protein models into a composite used for
             ranking.</p>

          <h3>Important interpretation notes</h3>
          <ul class="hp-notes">
            <li><b>Reference-relative calls:</b> REF and ALT are defined relative to each reference genome; &ldquo;alternate&rdquo; does not mean rare, harmful, or derived.</li>
            <li><b>Isolates are not all independent:</b> collections can include closely related isolates. Population, host, chemotype, and country metadata (available for F. graminearum) should be considered when interpreting similarity.</li>
            <li><b>Missingness varies:</b> genotype completeness can differ by site and isolate. Pairwise tools exclude sites where either isolate lacks a usable call.</li>
            <li><b>Consequence is model dependent:</b> one genomic variant may receive several annotations when it overlaps multiple gene models; interfaces generally display the most severe or most relevant consequence.</li>
            <li><b>Model scores are predictions:</b> the DNA (FunDLM, EVO2) and protein (ESM1/ESM2/ESM3/ESM-C) scores prioritize candidates but do not establish biological causality. Use them together with frequency, consequence, domain, structure, and experimental evidence.</li>
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
    const items = g.items.map(x=>`<div class="hp-def"><dt>${esc(x[0])}</dt><dd>${esc(x[1])}</dd></div>`).join('');
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
    const faq = FAQ.map(f=>`<details class="hp-faq"><summary>${esc(f[0])}<span class="hp-chev">${ico('caret')}</span></summary><div>${esc(f[1])}</div></details>`).join('');

    page.innerHTML = `
      <section class="hp-hero">
        <div class="hp-eyebrow">SNPTools · Fusarium · SNPVersity</div>
        <h1>Help &amp; FAQ</h1>
        <p>SNPTools is an integrated suite for exploring Fusarium sequence variation. Everything starts from a genomic
           query and flows between tools without re-running it — this page explains each tool, the datasets behind
           them, and the vocabulary you\u2019ll meet along the way.</p>
        <div class="hp-jump">
          <a href="#hp-flow">How it works</a>
          <a href="#hp-tools">The tools</a>
          <a href="#hp-data">Datasets</a>
          <a href="#hp-definitions">Definitions</a>
          <a href="#hp-gloss">Scores &amp; annotations</a>
          <a href="#hp-faq">FAQ</a>
        </div>
      </section>

      <section id="hp-flow" class="hp-sec">
        <div class="hp-h"><span class="hp-n">01</span><h2>How the suite fits together</h2></div>
        <div class="hp-flow">
          <div class="hp-step"><span class="hp-si" style="background:#2563eb">${ico('search')}</span>
            <b>Query</b><p>In SNPVersity, choose a dataset, a region or gene, and the accessions you care about.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#1f8a4c">${ico('table')}</span>
            <b>Result</b><p>You get a genotype table and a VCF for exactly that slice of the genome.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#b45309">${ico('compare')}</span>
            <b>Send onward</b><p>Hand the same matrix to any other tool with "Send selection to…" — no re-query.</p></div>
        </div>
        <p class="hp-fine">A query is run against the real variant store and returned as a VCF; every other tool reuses that
          one result, so a set you build once stays consistent as you rank it, compare it, cluster it, or draw it.</p>
      </section>

      <section id="hp-tools" class="hp-sec">
        <div class="hp-h"><span class="hp-n">02</span><h2>The tools</h2></div>
        <p class="hp-lead">Twelve tools across five stages. Eight are live today; the rest are on the roadmap and share
          the same data and coordinates. Expand any tool for what it does, what it takes, and what it returns.</p>
        <div class="hp-tools">${toolCards}</div>
      </section>

      <section id="hp-data" class="hp-sec">
        <div class="hp-h"><span class="hp-n">03</span><h2>Datasets</h2></div>
        <p class="hp-lead">Each query runs against one dataset. Each is called against its Fusarium reference genome; they differ
          in how they were filtered, how many accessions and sites they hold, and which score columns they carry.</p>
        ${datasetTable()}
        ${fusariumDatasetDetails()}
      </section>

      <section id="hp-definitions" class="hp-sec">
        <div class="hp-h"><span class="hp-n">04</span><h2>Definitions &amp; table columns</h2></div>
        <p class="hp-lead">Definitions are grouped by tool and use the same labels shown in the interfaces. Expand a group or use your browser's find command to locate a column heading.</p>
        <div class="hp-defgroups">${definitions}</div>
      </section>

      <section id="hp-gloss" class="hp-sec">
        <div class="hp-h"><span class="hp-n">05</span><h2>Scores &amp; annotations</h2></div>
        <dl class="hp-gloss">${gloss}</dl>
      </section>

      <section id="hp-faq" class="hp-sec">
        <div class="hp-h"><span class="hp-n">06</span><h2>Frequently asked</h2></div>
        <div class="hp-faqs">${faq}</div>
      </section>

      <section class="hp-foot">
        <div>
          <b>Still stuck?</b>
          <p>SNPTools is developed by MaizeGDB. This build serves Fusarium variant data.</p>
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
      .hp-deflist{margin:0;padding:0 15px 8px 35px;display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:0 28px}
      .hp-def{padding:11px 0;border-top:1px solid #eef1f6}
      .hp-def dt{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:13.5px;font-weight:600;color:var(--ink,#141922);margin-bottom:3px}
      .hp-def dd{margin:0;font-size:12.8px;line-height:1.55;color:var(--muted,#5b6b83)}


      /* glossary */
      .hp-gloss{margin:0;display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:2px 26px}
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

  if (typeof SNPTools!=='undefined') SNPTools.register('help', { render });
  return { render, open };
})();
if (typeof window!=='undefined'){ window.SNPHelp = SNPHelp; window.openHelp = SNPHelp.open; }
