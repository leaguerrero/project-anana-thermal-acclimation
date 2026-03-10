# project-anana-thermal-acclimation

Some cells in the gene_methylation_matrix contain NA values. These indicate cases where gene-level DNA methylation could not be calculated for a given gene in a given sample. Gene-level methylation estimates were derived from CpG methylation calls obtained from whole-genome bisulfite sequencing (WGBS). CpG sites were included only if they met coverage filtering criteria (minimum read depth ≥10 and additional quality filters) as described in the Methods.

For some genes, no CpG sites passed these filtering thresholds in a given sample. In these cases, gene-level methylation could not be estimated and the value is recorded as NA.

NA therefore represents missing measurements due to insufficient CpG coverage, not zero methylation.

Missing values were retained rather than imputed to preserve the original structure of the merged dataset used in downstream analyses.
