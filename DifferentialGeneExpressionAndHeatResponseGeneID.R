# Differential gene expression and heat response gene identification
# Libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)

# File Set-up
my.dir <- "/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/HTSeq_out"
my.files <- grep(".txt", list.files(my.dir), value=TRUE)
my.metadata <- read.csv("/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/HTSeq_out/a_nana_metadata.csv")

# Create sample table
sampleNames <- c("D0_T10_C_1", "D0_T10_C_2", "D0_T10_C_3", "D0_T10_H_1", "D0_T10_H_2", "D0_T10_H_3",
                 "D0_T3_C_1", "D0_T3_C_2", "D0_T3_C_3", "D0_T3_H_1", "D0_T3_H_2", "D0_T3_H_3",
                 "D0_T4_C_1", "D0_T4_C_2", "D0_T4_C_3", "D0_T4_H_1", "D0_T4_H_2", "D0_T4_H_3",
                 "D0_T9_C_1", "D0_T9_C_2", "D0_T9_C_3", "D0_T9_H_1", "D0_T9_H_2", "D0_T9_H_3",
                 "D11_T10_C_1", "D11_T10_C_2", "D11_T10_C_3", "D11_T10_H_1", "D11_T10_H_2", "D11_T10_H_3",
                 "D11_T3_C_1", "D11_T3_C_2", "D11_T3_C_3", "D11_T3_H_1", "D11_T3_H_2", "D11_T3_H_3",
                 "D11_T4_C_1", "D11_T4_C_2", "D11_T4_C_3", "D11_T4_H_1", "D11_T4_H_2", "D11_T4_H_3",
                 "D11_T9_C_1", "D11_T9_C_2", "D11_T9_C_3", "D11_T9_H_1", "D11_T9_H_2", "D11_T9_H_3")
my.sampleTable <- data.frame(sampleName = sampleNames, 
                             fileName = my.files, 
                             condition = my.metadata)

# Convert variables to factors
factorVars <- c("condition.accl.temp", "condition.accl.time", "condition.stress.treatment", 
                "condition.tank", "condition.colony")

my.sampleTable[, factorVars] <- lapply(my.sampleTable[, factorVars], factor)
colnames(my.sampleTable)

# Create DESeqDataSet Object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable,
  directory = my.dir,
  design = ~ condition.accl.temp +
    condition.stress.treatment +
    condition.stress.treatment:condition.accl.temp)
# Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Differential Expression Analysis
a.nana.dds <- DESeq(ddsHTSeq) 
resultsNames(a.nana.dds)
#saveRDS(a.nana.dds, file = "/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/DataOutput/a.nana.dds.RData")

# Report: Number of transcripts
num_transcripts <- nrow(a.nana.dds)
num_transcripts

# Save for topGO analysis
total.set.gene.names <- rownames(assay(a.nana.dds))
total.set.gene.names.df <- data.frame(gene = total.set.gene.names)
# write.csv(total.set.gene.names.df, file = '~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.set.gene.names.csv', row.names = FALSE)

# Variance stabilizing transformation
vtd <- vst(a.nana.dds)

# PCA
deseq_PCA <- plotPCA(vtd,intgroup=c( "condition.stress.treatment"), returnData = TRUE) 
deseq_PCA$Acclimation<- my.metadata$Acclimation[match((deseq_PCA$name), my.metadata$ID)]

# Figure 2A:
# Plot PCA
HeatStress_ge_PCA_plot <- ggplot(deseq_PCA,aes(x=PC1,y=PC2, shape = condition.stress.treatment,  color = Acclimation)) +
  theme_classic(base_size = 20) +
  theme(legend.position="none") + 
  scale_color_manual(values = c("#d77943", "#7943d7")) +
  labs(title="Gene Expression", shape ="Stress Treatment") +
  geom_point(size = 5, alpha = 0.7)

HeatStress_ge_PCA_plot
#jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/Figures/ge_pca.jpeg', units="in", width=5, height=5, res=300)
#print(HeatStress_ge_PCA_plot)     
#dev.off() 

non.acclimated.files <- grep(".txt", list.files(my.dir), value = TRUE)[1:36]

# Filter metadata and create sample table
fil.metadata <- my.metadata[1:36,]
non.acclimated.sampleTable <- data.frame(
  sampleName = c("D0_T10_C_1", "D0_T10_C_2", "D0_T10_C_3", "D0_T10_H_1", "D0_T10_H_2", "D0_T10_H_3", "D0_T3_C_1", "D0_T3_C_2",
                 "D0_T3_C_3", "D0_T3_H_1", "D0_T3_H_2", "D0_T3_H_3", "D0_T4_C_1", "D0_T4_C_2", "D0_T4_C_3", "D0_T4_H_1",
                 "D0_T4_H_2", "D0_T4_H_3", "D0_T9_C_1", "D0_T9_C_2", "D0_T9_C_3", "D0_T9_H_1", "D0_T9_H_2", "D0_T9_H_3",
                 "D11_T10_C_1", "D11_T10_C_2", "D11_T10_C_3", "D11_T10_H_1", "D11_T10_H_2", "D11_T10_H_3", "D11_T3_C_1",
                 "D11_T3_C_2", "D11_T3_C_3", "D11_T3_H_1", "D11_T3_H_2", "D11_T3_H_3"),
  fileName = non.acclimated.files,
  condition = fil.metadata
)

# Convert factors in the sample table
columns_to_factorize <- c("condition.accl.temp", "condition.accl.time", "condition.stress.treatment", "condition.tank", "condition.colony")

for (column in columns_to_factorize) {
  non.acclimated.sampleTable[[column]] <- factor(non.acclimated.sampleTable[[column]])
}


# Create DESeqDataSet
non_acclimated_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = non.acclimated.sampleTable,
  directory = my.dir,
  design = ~ condition.accl.temp +
    condition.stress.treatment +
    condition.stress.treatment:condition.accl.temp
)

# Filter low count rows
nonacc.keep <- rowSums(counts(non_acclimated_ddsHTSeq)) >= 10
non_acclimated_ddsHTSeq <- non_acclimated_ddsHTSeq[nonacc.keep,]

# Differential Expression Analysis
non.acc.a.nana.dds <- DESeq(non_acclimated_ddsHTSeq)
resultsNames(non.acc.a.nana.dds)

# Subsetting dds object and performing DE analysis for different conditions
non.acc.a.nana.dds$group <- factor(paste0(
  non.acc.a.nana.dds$condition.stress.treatment,
  non.acc.a.nana.dds$condition.accl.time,
  non.acc.a.nana.dds$condition.accl.temp
))

design(non.acc.a.nana.dds) <- ~ group
group.dds <- DESeq(non.acc.a.nana.dds)
resultsNames(group.dds)

# Day 0 29
D0.29.response.results <- lfcShrink(group.dds,
                                    contrast = c("group", "stress029", "control029"),
                                    type = "ashr")
D0.29.response.results <- D0.29.response.results[complete.cases(D0.29.response.results),]

D0.29C.HRGs <- subset(D0.29.response.results,
                      padj < 0.01 &
                        abs(log2FoldChange) >= 2)

D0.29C.HRGs.tb <- D0.29C.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

# Day 0 31
D0.31.response.results <- lfcShrink(group.dds,
                                    contrast = c("group", "stress031", "control031"),
                                    type = "ashr")
D0.31.response.results <- D0.31.response.results[complete.cases(D0.31.response.results),]

D0.31C.HRGs <- subset(D0.31.response.results,
                      padj < 0.01 &
                        abs(log2FoldChange) >= 2)

D0.31C.HRGs.tb <- D0.31C.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

# Day 11 29
D11.29.response.results <- lfcShrink(group.dds,
                                     contrast = c("group", "stress1129", "control1129"),
                                     type = "ashr")
D11.29.response.results <- D11.29.response.results[complete.cases(D11.29.response.results),]

D11.29C.HRGs <- subset(D11.29.response.results,
                       padj < 0.01 &
                         abs(log2FoldChange) >= 2)

D11.29C.HRGs.tb <- D11.29C.HRGs %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
dim(D11.29C.HRGs.tb)

nonaccl.HRG <- union(union(D0.31C.HRGs.tb$gene, 
                           D11.29C.HRGs.tb$gene), 
                     D0.29C.HRGs.tb$gene)
nonaccl.HRG.gene.list <- as.data.frame(nonaccl.HRG)
colnames(nonaccl.HRG.gene.list) <- 'gene'
# write.csv(nonaccl.HRG.gene.list, file = '~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/nonaccl.HRG.gene.list.csv', row.names = FALSE)


core.HRGs <- inner_join(
  inner_join(D0.31C.HRGs.tb, 
             D11.29C.HRGs.tb, 
             by = "gene"), 
  D0.29C.HRGs.tb, by = "gene")


core.HRGs.gene.list <- as.data.frame(core.HRGs$gene)
colnames(core.HRGs.gene.list) <- 'gene'
# write.csv(core.HRGs.gene.list, file = '~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/core.HRG.list.csv', row.names = FALSE)

# Gene Ontology (GO) terms for core heat response genes

# Set-up
# Read in annotation table and reformat GO column
annos <- read.csv("~/LabNotebook/Chapter1/Amil_v2.01/Amillepora_trinotate_annotation_report.csv",na.strings=".")
annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
GOmap <- annos[,c("transcript_id","GOlist")]
# write.table(GOmap,"~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
geneID2GO <- readMappings(file="~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",sep="\t",IDsep=",")

# Read in set of all genes and make 'background' gene set
all <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.set.gene.names.csv")
allgenes <- all$gene

# Modify core heat response gene set for topGO analysis
## Note: Used VLOOKUP in Excel to modify core gene data set

core.genes <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/core.HRG.list.csv")
core.genesIG <- factor(as.numeric(allgenes%in%core.genes$gene))
names(core.genesIG) <- allgenes

##BP
GOcore <- new("topGOdata",ontology="BP",allGenes=core.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
coreFisher <- getSigGroups(GOcore,test.stat)
core.res <- GenTable(GOcore,classic=coreFisher,topNodes=length(coreFisher@score),numChar=100)
core.filt <- core.res[core.res$classic<0.01 & core.res$Significant>=10,]
#write.csv(core.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/core.genes.GO_BP.csv",sep=""))

##MF
GOcore <- new("topGOdata",ontology="MF",allGenes=core.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
coreFisher <- getSigGroups(GOcore,test.stat)
core.res <- GenTable(GOcore,classic=coreFisher,topNodes=length(coreFisher@score),numChar=100)
core.filt <- core.res[core.res$classic<0.01 & core.res$Significant>=10,]
#write.csv(core.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/core.genes.GO_MF.csv",sep=""))

##CC
GOcore <- new("topGOdata",ontology="CC",allGenes=core.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
coreFisher <- getSigGroups(GOcore,test.stat)
core.res <- GenTable(GOcore,classic=coreFisher,topNodes=length(coreFisher@score),numChar=100)
core.filt <- core.res[core.res$classic<0.01 & core.res$Significant>=10,]
#write.csv(core.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/core.genes.GO_CC.csv",sep=""))

# Gene Ontology (GO) terms for total set heat response genes

# Set-up
# Read in annotation table and reformat GO column
annos <- read.csv("~/LabNotebook/Chapter1/Amil_v2.01/Amillepora_trinotate_annotation_report.csv",na.strings=".")
annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
GOmap <- annos[,c("transcript_id","GOlist")]
# write.table(GOmap,"~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
geneID2GO <- readMappings(file="~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",sep="\t",IDsep=",")

# Read in set of all genes and make 'background' gene set
all <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.set.gene.names.csv")
allgenes <- all$gene

# Modify core heat response gene set for topGO analysis
## Note: Used VLOOKUP in Excel to modify core gene data set

total.hr.genes <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/nonaccl.HRG.gene.list.csv")
total.hr.genesIG <- factor(as.numeric(allgenes%in%total.hr.genes$gene))
names(total.hr.genesIG) <- allgenes

##BP
GOtotal.hr <- new("topGOdata",ontology="BP",allGenes=total.hr.genesIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
total.hrFisher <- getSigGroups(GOtotal.hr,test.stat)
total.hr.res <- GenTable(GOtotal.hr,classic=total.hrFisher,topNodes=length(total.hrFisher@score),numChar=100)
total.hr.filt <- total.hr.res[total.hr.res$classic<0.01 & total.hr.res$Significant>=10,]
write.csv(total.hr.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.hr.genes.GO_BP.csv",sep=""))

##MF
GOtotal.hr <- new("topGOdata",ontology="MF",allGenes=total.hr.genesIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
total.hrFisher <- getSigGroups(GOtotal.hr,test.stat)
total.hr.res <- GenTable(GOtotal.hr,classic=total.hrFisher,topNodes=length(total.hrFisher@score),numChar=100)
total.hr.filt <- total.hr.res[total.hr.res$classic<0.01 & total.hr.res$Significant>=10,]
write.csv(total.hr.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.hr.genes.GO_MF.csv",sep=""))

##CC
GOtotal.hr <- new("topGOdata",ontology="CC",allGenes=total.hr.genesIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
total.hrFisher <- getSigGroups(GOtotal.hr,test.stat)
total.hr.res <- GenTable(GOtotal.hr,classic=total.hrFisher,topNodes=length(total.hrFisher@score),numChar=100)
total.hr.filt <- total.hr.res[total.hr.res$classic<0.01 & total.hr.res$Significant>=10,]
write.csv(total.hr.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.hr.genes.GO_CC.csv",sep=""))


