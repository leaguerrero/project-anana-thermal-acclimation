# Overlap between gene sets

library(tidyverse)
setwd("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/")


## Total set HRGs and Transcriptionally Modulated Genes

total.gene.set <- read.csv("total.set.gene.names.csv", header = TRUE)

total.set.HRGs <- read.csv("nonaccl.HRG.gene.list.csv", header = TRUE)

amp.genes <- read.csv("amp.genes.csv", header = TRUE)
amp.gene.names <- data.frame(amp.genes[,1])
colnames(amp.gene.names) <- colnames(amp.genes)[1]

damp.genes <- read.csv("damp.genes.csv", header = TRUE)
damp.gene.names <- data.frame(damp.genes[,1])
colnames(damp.gene.names) <- colnames(damp.genes)[1]

tmgs <- rbind(amp.gene.names, damp.gene.names)

# Count the overlap
overlap.count <- sum(total.set.HRGs$gene %in% tmgs$gene)
print(overlap.count) #299

# Set up contingency table
total.genes <- nrow(total.gene.set)
total.set.HRG.count <- nrow(total.set.HRGs)
tmg.count <- nrow(tmgs)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, total.set.HRG.count - overlap.count,
                              tmg.count - overlap.count,
                              total.genes - total.set.HRG.count - tmg.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) #7.367248e-51

# Remove objects
rm(total.gene.set, total.set.HRGs, amp.genes, amp.gene.names, damp.genes, damp.gene.names, tmgs, overlap.count, total.genes, total.set.HRG.count, tmg.count, contingency.table, result, p.value)

## Core set HRGs and Transcriptionally Modulated Genes

total.gene.set <- read.csv("total.set.gene.names.csv", header = TRUE)

core.set.HRGs <- read.csv("core.HRG.list.csv", header = TRUE)

amp.genes <- read.csv("amp.genes.csv", header = TRUE)
amp.gene.names <- data.frame(amp.genes[,1])
colnames(amp.gene.names) <- colnames(amp.genes)[1]

damp.genes <- read.csv("damp.genes.csv", header = TRUE)
damp.gene.names <- data.frame(damp.genes[,1])
colnames(damp.gene.names) <- colnames(damp.genes)[1]

tmgs <- rbind(amp.gene.names, damp.gene.names)

# Count the overlap
overlap.count <- sum(core.set.HRGs$gene %in% tmgs$gene)
print(overlap.count) # 18

# Set up contingency table
total.genes <- nrow(total.gene.set)
core.set.HRG.count <- nrow(core.set.HRGs)
tmg.count <- nrow(tmgs)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, core.set.HRG.count - overlap.count,
                              tmg.count - overlap.count,
                              total.genes - core.set.HRG.count - tmg.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 0.8059585

# Remove objects
rm(total.gene.set, core.set.HRGs, amp.genes, amp.gene.names, damp.genes, damp.gene.names, tmgs, overlap.count, total.genes, core.set.HRG.count, tmg.count, contingency.table, result, p.value)

## Genes with diff meth due to acclimation and genes with diff meth due to heat stress
total.gene.acc.set <- read.csv("total.meth.set.acc.genes.csv", header = TRUE)
total.gene.hs.set <- read.csv("total.meth.set.hs.genes.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.acc.set, total.gene.hs.set)

dm.acc.gene.set <- read.csv("diff.meth.acc.df.csv", header = TRUE)

dm.hs.gene.set <- read.csv("diff.meth.hs.df.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.acc.gene.set$gene %in% dm.hs.gene.set$gene)
print(overlap.count) # 40

# Print the genes that overlap
dm.acc.dm.hs.gene.list <- intersect(dm.acc.gene.set$gene, dm.hs.gene.set$gene)
dm.acc.dm.hs.genes <- as.data.frame(dm.acc.dm.hs.gene.list)
colnames(dm.acc.dm.hs.genes) <- "gene"

sub.diff.meth.hs.df <- right_join(diff.meth.hs.df, dm.acc.dm.hs.genes, by = 'gene')
names(sub.diff.meth.hs.df)[5] <- "hs.meth.diff"

sub.diff.meth.acc.df <- right_join(diff.meth.acc.df, dm.acc.dm.hs.genes, by = 'gene')
names(sub.diff.meth.acc.df)[5] <- "acc.meth.diff"

sub.diff.meth <- full_join(sub.diff.meth.hs.df, sub.diff.meth.acc.df, by = 'gene')
sub.diff.meth <- sub.diff.meth[,c(1,5,9)]

sub.diff.meth.plot <- ggplot(sub.diff.meth, aes(x = acc.meth.diff, y = hs.meth.diff)) +
  geom_point() +
  xlab("hs.meth.diff") +
  ylab("acc.meth.diff") +
  ggtitle("Scatterplot of hs.meth.diff and acc.meth.diff")
sub.diff.meth.plot

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.acc.gene.count <- nrow(dm.acc.gene.set)
dm.hs.gene.count <- nrow(dm.hs.gene.set)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.acc.gene.count - overlap.count,
                              dm.hs.gene.count - overlap.count,
                              total.genes - dm.acc.gene.count - dm.hs.gene.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 4.024402e-12

# Remove objects
rm(total.gene.acc.set, total.gene.hs.set, total.gene.set, dm.acc.gene.set, dm.hs.gene.set, overlap.count, total.genes, dm.acc.gene.count, dm.hs.gene.count, contingency.table, result, p.value)


## Genes with diff meth due to acclimation and total set HRGs 
total.gene.acc.set <- read.csv("total.meth.set.acc.genes.csv", header = TRUE)
total.gene.exp.set <- read.csv("total.set.gene.names.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.acc.set, total.gene.exp.set)

dm.acc.gene.set <- read.csv("diff.meth.acc.df.csv", header = TRUE)

total.hrg.set <- read.csv("nonaccl.HRG.gene.list.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.acc.gene.set$gene %in% total.hrg.set$gene)
print(overlap.count) # 88

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.acc.gene.count <- nrow(dm.acc.gene.set)
total.hrg.count <- nrow(total.hrg.set)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.acc.gene.count - overlap.count,
                              total.hrg.count - overlap.count,
                              total.genes - dm.acc.gene.count - total.hrg.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 0.4961182

# Remove objects
rm(total.gene.acc.set, total.gene.exp.set, total.gene.set, dm.acc.gene.set, total.hrg.set, overlap.count, total.genes, dm.acc.gene.count, total.hrg.count, contingency.table, result, p.value)

## Genes with diff meth due to acclimation and total set HRGs 
total.gene.hs.set <- read.csv("total.meth.set.hs.genes.csv", header = TRUE)
total.gene.exp.set <- read.csv("total.set.gene.names.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.hs.set, total.gene.exp.set)

dm.hs.gene.set <- read.csv("diff.meth.hs.df.csv", header = TRUE)

total.hrg.set <- read.csv("nonaccl.HRG.gene.list.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.hs.gene.set$gene %in% total.hrg.set$gene)
print(overlap.count) # 101

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.hs.gene.count <- nrow(dm.hs.gene.set)
total.hrg.count <- nrow(total.hrg.set)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.hs.gene.count - overlap.count,
                              total.hrg.count - overlap.count,
                              total.genes - dm.hs.gene.count - total.hrg.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 0.171473

# Remove objects
rm(total.gene.hs.set, total.gene.exp.set, total.gene.set, dm.hs.gene.set, total.hrg.set, overlap.count, total.genes, dm.hs.gene.count, total.hrg.count, contingency.table, result, p.value)

## Genes with diff meth due to acclimation and Dampened genes 
total.gene.acc.set <- read.csv("total.meth.set.acc.genes.csv", header = TRUE)
total.gene.exp.set <- read.csv("total.set.gene.names.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.acc.set, total.gene.exp.set)

dm.acc.gene.set <- read.csv("diff.meth.acc.df.csv", header = TRUE)

damp.gene.set <- read.csv("damp.genes.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.acc.gene.set$gene %in% damp.gene.set$gene)
print(overlap.count) # 14 

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.acc.gene.count <- nrow(dm.acc.gene.set)
damp.gene.count <- nrow(damp.gene.set)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.acc.gene.count - overlap.count,
                              damp.gene.count - overlap.count,
                              total.genes - dm.acc.gene.count - damp.gene.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 0.008822156

# Remove objects
rm(total.gene.acc.set, total.gene.exp.set, total.gene.set, dm.acc.gene.set, damp.gene.set, overlap.count, total.genes, dm.acc.gene.count, damp.gene.count, contingency.table, result, p.value)

## Genes with diff meth due to acclimation and amplified genes 
total.gene.acc.set <- read.csv("total.meth.set.acc.genes.csv", header = TRUE)
total.gene.exp.set <- read.csv("total.set.gene.names.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.acc.set, total.gene.exp.set)

dm.acc.gene.set <- read.csv("diff.meth.acc.df.csv", header = TRUE)

amp.gene.set <- read.csv("amp.genes.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.acc.gene.set$gene %in% amp.gene.set$gene)
print(overlap.count) # 3 

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.acc.gene.count <- nrow(dm.acc.gene.set)
amp.gene.count <- nrow(amp.gene.set)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.acc.gene.count - overlap.count,
                              amp.gene.count - overlap.count,
                              total.genes - dm.acc.gene.count - amp.gene.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 1

# Remove objects
rm(total.gene.acc.set, total.gene.exp.set, total.gene.set, dm.acc.gene.set, amp.gene.set, overlap.count, total.genes, dm.acc.gene.count, amp.gene.count, contingency.table, result, p.value)

## Genes with diff meth due to hs and amplified genes 
total.gene.hs.set <- read.csv("total.meth.set.hs.genes.csv", header = TRUE)
total.gene.exp.set <- read.csv("total.set.gene.names.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.hs.set, total.gene.exp.set)

dm.hs.gene.set <- read.csv("diff.meth.hs.df.csv", header = TRUE)

amp.gene.set <- read.csv("amp.genes.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.hs.gene.set$gene %in% amp.gene.set$gene)
print(overlap.count) # 5 

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.hs.gene.count <- nrow(dm.hs.gene.set)
amp.gene.count <- nrow(amp.gene.set)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.hs.gene.count - overlap.count,
                              amp.gene.count - overlap.count,
                              total.genes - dm.hs.gene.count - amp.gene.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 0.3929648

# Remove objects
rm(total.gene.hs.set, total.gene.exp.set, total.gene.set, dm.hs.gene.set, amp.gene.set, overlap.count, total.genes, dm.hs.gene.count, amp.gene.count, contingency.table, result, p.value)

## Genes with diff meth due to hs and dampened genes 
total.gene.hs.set <- read.csv("total.meth.set.hs.genes.csv", header = TRUE)
total.gene.exp.set <- read.csv("total.set.gene.names.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.hs.set, total.gene.exp.set)

dm.hs.gene.set <- read.csv("diff.meth.hs.df.csv", header = TRUE)

damp.gene.set <- read.csv("damp.genes.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.hs.gene.set$gene %in% damp.gene.set$gene)
print(overlap.count) # 9 

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.hs.gene.count <- nrow(dm.hs.gene.set)
damp.gene.count <- nrow(damp.gene.set)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.hs.gene.count - overlap.count,
                              damp.gene.count - overlap.count,
                              total.genes - dm.hs.gene.count - damp.gene.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 0.4469999

# Remove objects
rm(total.gene.hs.set, total.gene.exp.set, total.gene.set, dm.hs.gene.set, damp.gene.set, overlap.count, total.genes, dm.hs.gene.count, damp.gene.count, contingency.table, result, p.value)

## Core HRGs and Diff Methylation following heat stress:
total.gene.exp.set <- read.csv("total.set.gene.names.csv", header = TRUE)
total.gene.hs.set <- read.csv("total.meth.set.hs.genes.csv", header = TRUE)
total.gene.set <- dplyr::union(total.gene.exp.set, total.gene.hs.set)

dm.hs.gene.set <- read.csv("diff.meth.hs.df.csv", header = TRUE)
core.set.HRGs <- read.csv("core.HRG.list.csv", header = TRUE)

# Count the overlap
overlap.count <- sum(dm.hs.gene.set$gene %in% core.set.HRGs$gene)
print(overlap.count) # 24

# Set up contingency table
total.genes <- nrow(total.gene.set)
dm.hs.gene.count <- nrow(dm.hs.gene.set)
core.gene.count <- nrow(core.set.HRGs)
overlap.count <- overlap.count

# Create a 2x2 contingency table
contingency.table <- matrix(c(overlap.count, dm.hs.gene.count - overlap.count,
                              core.gene.count - overlap.count,
                              total.genes - dm.hs.gene.count - core.gene.count + overlap.count),
                            nrow = 2, ncol = 2, byrow = TRUE)

# Perform Fisher's exact test
result <- fisher.test(contingency.table)

# Extract p-value
p.value <- result$p.value

# Print the p-value
print(p.value) # 0.001445349

# Remove objects
rm(total.gene.exp.set, total.gene.hs.set, total.gene.set, dm.hs.gene.set, core.set.HRGs, overlap.count, total.genes, dm.hs.gene.count, core.gene.count, contingency.table, result, p.value)


