# Quantifying transcriptional response modulation

# Libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(topGO)
library(stringr)


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

# transcriptional response modulation
a.nana.dds.d11.sub <- a.nana.dds[, a.nana.dds@colData@listData$condition.accl.time == "11"]

coeff.df <- as.data.frame(coef(a.nana.dds.d11.sub), colnames = TRUE)
wald.stat.df <- as.data.frame(cbind(row.names(coeff.df), as.numeric(coeff.df$condition.accl.temp31.condition.stress.treatmentstress)), colnames = TRUE)
colnames(wald.stat.df) <- c("gene", "coefficient")

significant.results <- results(a.nana.dds.d11.sub, name = "condition.accl.temp31.condition.stress.treatmentstress")
significant.results <- significant.results[complete.cases(significant.results),]
significant.results <- significant.results[significant.results$padj < 0.1,]
significant.results.df <- significant.results %>% 
                          data.frame() %>% 
                          rownames_to_column(var = "gene") 

significant.interaction.df <- dplyr::left_join(significant.results.df, 
                                               wald.stat.df, 
                                               by = "gene") %>% 
                              mutate(coefficient = as.numeric(coefficient))



# Dampened Genes
neg.wald.stat <- nrow(significant.interaction.df[significant.interaction.df$stat < 0,])
dampened.genes <- significant.interaction.df[significant.interaction.df$stat < 0,]
# write.csv(dampened.genes, file = "~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/damp.genes.csv", row.names = FALSE)

# Amplified Genes
pos.wald.stat <- nrow(significant.interaction.df[significant.interaction.df$stat > 0,])
amp.genes <- significant.interaction.df[significant.interaction.df$stat > 0,]
# write.csv(amp.genes, file = "~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/amp.genes.csv", row.names = FALSE)

# Difference in the number of Dampened and Amplified Genes
total.rows <- nrow(significant.interaction.df)
wald.stat.mat <- matrix(c(neg.wald.stat, pos.wald.stat), ncol = 2, byrow = TRUE)
colnames(wald.stat.mat) <- c("gene count with neg wald stat", "gene count with pos wald stat")
wald.stat.table <- as.table(wald.stat.mat)
wald.stat.table <- round(prop.table(wald.stat.table), 2)

test <- prop.test(x = neg.wald.stat, n = total.rows, p = 0.5)
test # p = 0.6840491 

# Figure 3C: 
significant.interaction.df <- significant.interaction.df %>%
                              mutate(modulation.cat = ifelse(stat > 0,
                                                             "amplified", 
                                                             "dampened"))
modulation.counts <- table(significant.interaction.df$modulation.cat)
modulation.counts

significant.interaction.df$logbm <- log10(significant.interaction.df$baseMean)

MyColor2 <- c("#d8b365", "#5ab4ac")
names(MyColor2) <- c("amp","damp")
compare_means(logbm ~ modulation.cat, data = significant.interaction.df, method = "t.test")

amp.damp.basemean <- ggplot(significant.interaction.df, aes(x = modulation.cat, y = logbm, fill = modulation.cat)) +
  theme_classic(base_size = 20) +
  theme(legend.position="none") +
  xlab(NULL) + 
  ylab("log10 mean expression") +
  # labs(title="Base mean of amplified and dampened genes") + 
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = c(0.5, 1), alpha = 0.8) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA, alpha = 0.8) +
  scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  stat_compare_means(label = "p.signif", size = 6, method = 't.test') 

amp.damp.basemean

# jpeg('~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/basemeans_amp_damp.jpeg', units="in", width=6, height=5, res=300)
# print(amp.damp.basemean)     
# dev.off() 


# Gene Ontology (GO) terms for amplified and dampened genes

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

# Modify amplified and dampened gene sets for topGO analysis
## Note: Used VLOOKUP in Excel to modify amp.genes and dampened.genes data sets

# Amplified genes
amp.genes <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/amp.genes.csv")
amp.genesIG <- factor(as.numeric(allgenes%in%amp.genes$gene))
names(amp.genesIG) <- allgenes

##BP
GOamp <- new("topGOdata",ontology="BP",allGenes=amp.genesIG, nodeSize=10,
            annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ampFisher <- getSigGroups(GOamp,test.stat)
amp.res <- GenTable(GOamp,classic=ampFisher,topNodes=length(ampFisher@score),numChar=100)
amp.filt <- amp.res[amp.res$classic<0.01 & amp.res$Significant>=10,]
# write.csv(amp.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/amp.genes.GO_BP.csv",sep=""))

##MF
GOamp <- new("topGOdata",ontology="MF",allGenes=amp.genesIG, nodeSize=10,
            annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ampFisher <- getSigGroups(GOamp,test.stat)
amp.res <- GenTable(GOamp,classic=ampFisher,topNodes=length(ampFisher@score),numChar=100)
amp.filt <- amp.res[amp.res$classic<0.01 & amp.res$Significant>=10,]
#write.csv(amp.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/amp.genes.GO_MF.csv",sep=""))

##CC
GOamp <- new("topGOdata",ontology="CC",allGenes=amp.genesIG, nodeSize=10,
            annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
ampFisher <- getSigGroups(GOamp,test.stat)
amp.res <- GenTable(GOamp,classic=ampFisher,topNodes=length(ampFisher@score),numChar=100)
amp.filt <- amp.res[amp.res$classic<0.01 & amp.res$Significant>=10,]
#write.csv(amp.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/amp.genes.GO_CC.csv",sep=""))

# Dampened Genes
damp.genes <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/damp.genes.csv")
damp.genesIG <- factor(as.numeric(allgenes%in%damp.genes$gene))
names(damp.genesIG) <- allgenes

##BP
GOdamp <- new("topGOdata",ontology="BP",allGenes=damp.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dampFisher <- getSigGroups(GOdamp,test.stat)
damp.res <- GenTable(GOdamp,classic=dampFisher,topNodes=length(dampFisher@score),numChar=100)
damp.filt <- damp.res[damp.res$classic<0.01 & damp.res$Significant>=10,]
#write.csv(damp.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/damp.genes.GO_BP.csv",sep=""))

##MF
GOdamp <- new("topGOdata",ontology="MF",allGenes=damp.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dampFisher <- getSigGroups(GOdamp,test.stat)
damp.res <- GenTable(GOdamp,classic=dampFisher,topNodes=length(dampFisher@score),numChar=100)
damp.filt <- damp.res[damp.res$classic<0.01 & damp.res$Significant>=10,]
#write.csv(damp.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/damp.genes.GO_MF.csv",sep=""))

##CC
GOdamp <- new("topGOdata",ontology="CC",allGenes=damp.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dampFisher <- getSigGroups(GOdamp,test.stat)
damp.res <- GenTable(GOdamp,classic=dampFisher,topNodes=length(dampFisher@score),numChar=100)
damp.filt <- damp.res[damp.res$classic<0.01 & damp.res$Significant>=10,]
#write.csv(damp.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/damp.genes.GO_CC.csv",sep=""))



















