# project-anana-thermal-acclimation

Missing values (NA)
Some columns in the gene_methylation_matrix contain NA values. These values represent genes for which data were unavailable for one of the integrated data types used in downstream analyses. RNA-seq gene expression and whole-genome bisulfite sequencing (WGBS) methylation data were generated and processed independently. 

Gene-level methylation estimates were calculated only for genes that met coverage thresholds for CpG sites (minimum read depth ≥10 and coverage requirements across samples). As a result, some genes that had sufficient RNA-seq read counts did not have sufficient CpG coverage for reliable methylation estimation and therefore contain NA values in methylation columns. Conversely, some genes with methylation estimates may not have passed expression filtering thresholds.

NA values therefore indicate missing measurements rather than zero methylation or zero expression, and typically arise from differences in sequencing coverage or filtering criteria between the RNA-seq and WGBS datasets.

Missing values were retained rather than imputed to preserve the original structure of the merged dataset used in downstream analyses.
