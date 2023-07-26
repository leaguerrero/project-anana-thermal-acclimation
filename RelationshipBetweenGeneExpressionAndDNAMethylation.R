# Relationship between gene expression and DNA methylation
library(DESeq2)
library(methylKit)
library(rtracklayer)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

# Gene Expression: File Set-up
my.dir <- "/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/HTSeq_out"
my.files <- grep(".txt", list.files(my.dir), value=TRUE)
my.metadata <- read.csv("/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/HTSeq_out/a_nana_metadata.csv")

# Gene Expression: Create sample table
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

# Gene Expression: Convert variables to factors
factorVars <- c("condition.accl.temp", "condition.accl.time", "condition.stress.treatment", 
                "condition.tank", "condition.colony")

my.sampleTable[, factorVars] <- lapply(my.sampleTable[, factorVars], factor)
colnames(my.sampleTable)

# Gene Expression: Create DESeqDataSet Object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = my.sampleTable,
  directory = my.dir,
  design = ~ condition.accl.temp +
    condition.stress.treatment +
    condition.stress.treatment:condition.accl.temp)
# Gene Expression: Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Gene Expression: Differential Expression Analysis
a.nana.dds <- DESeq(ddsHTSeq) 

# Gene Expression: Making gene expression data frame, contrast is heat stress
heat.stress.DE.results <- results(a.nana.dds, 
                                  contrast = c("condition.stress.treatment", 
                                               "stress", 
                                               "control"))
heat.stress.DE.results <- as.data.frame(heat.stress.DE.results[complete.cases(heat.stress.DE.results), ])
heat.stress.DE.results$gene <- row.names(heat.stress.DE.results)
heat.stress.DE.results$log10_baseMean <- log10(heat.stress.DE.results$baseMean)
row.names(heat.stress.DE.results) <- NULL
heat.stress.DE.results <- heat.stress.DE.results[, c("gene", setdiff(names(heat.stress.DE.results), "gene"))]

# Methylation: File Set-up and meta data set-up
my.dir <- "/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/"
meta.file.name <- "meth_meta.csv"
meta.path <- file.path(my.dir, meta.file.name)
meta <- read.csv(meta.path, header = TRUE)
meta$Acclimation <- factor(meta$Acclimation,levels=c("Control","Accl"))
meta$Colony <- as.factor(meta$Colony)

# Methylation: Files for methylKit
setwd("~/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/MethylDackel_updated/")
file.list=list("D0T10C1C_CKDL200156542_rmdup_CpG.methylKit", 
               "D0T10C1H_CKDL200156543_rmdup_CpG.methylKit", 
               "D0T10C2C_CKDL200156544_rmdup_CpG.methylKit",
               "D0T10C2H_CKDL200156545_rmdup_CpG.methylKit", 
               "D0T10C3C_CKDL200156546_rmdup_CpG.methylKit",
               "D0T10C3H_CKDL200156547_rmdup_CpG.methylKit", 
               "D0T3C1C_CKDL200156536_rmdup_CpG.methylKit",
               "D0T3C1H_CKDL200156537_rmdup_CpG.methylKit",
               "D0T3C2C_CKDL200156538_rmdup_CpG.methylKit",
               "D0T3C2H_CKDL200156539_rmdup_CpG.methylKit",
               "D0T3C3C_CKDL200156540_rmdup_CpG.methylKit",
               "D0T3C3H_CKDL200156541_rmdup_CpG.methylKit",
               "D0T4C1C_CKDL200156560_rmdup_CpG.methylKit",
               "D0T4C1H_CKDL200156561_rmdup_CpG.methylKit",
               "D0T4C2C_CKDL200156562_rmdup_CpG.methylKit",
               "D0T4C2H_CKDL200156563_rmdup_CpG.methylKit",
               "D0T4C3C_CKDL200156564_rmdup_CpG.methylKit",
               "D0T4C3H_CKDL200156565_rmdup_CpG.methylKit",
               "D0T9C1C_CKDL200156566_rmdup_CpG.methylKit",
               "D0T9C1H_CKDL200156567_rmdup_CpG.methylKit",
               "D0T9C2C_CKDL200156568_rmdup_CpG.methylKit",
               "D0T9C2H_CKDL200156569_rmdup_CpG.methylKit",
               "D0T9C3C_CKDL200156570_rmdup_CpG.methylKit",
               "D0T9C3H_CKDL200156571_rmdup_CpG.methylKit",
               "D11T0C1C_CKDL200156554_rmdup_CpG.methylKit",
               "D11T0C1H_CKDL200156555_rmdup_CpG.methylKit",
               "D11T0C2C_CKDL200156556_rmdup_CpG.methylKit",
               "D11T0C2H_CKDL200156557_rmdup_CpG.methylKit",
               "D11T0C3C_CKDL200156558_rmdup_CpG.methylKit",
               "D11T0C3H_CKDL200156559_rmdup_CpG.methylKit",
               "D11T3C1C_CKDL200156548_rmdup_CpG.methylKit",
               "D11T3C1H_CKDL200156549_rmdup_CpG.methylKit",
               "D11T3C2C_CKDL200156550_rmdup_CpG.methylKit",
               "D11T3C2H_CKDL200156551_rmdup_CpG.methylKit",
               "D11T3C3C_CKDL200156552_rmdup_CpG.methylKit",
               "D11T3C3H_CKDL200156553_rmdup_CpG.methylKit",
               "D11T4C1C_CKDL200156572_rmdup_CpG.methylKit",
               "D11T4C1H_CKDL200156573_rmdup_CpG.methylKit",
               "D11T4C2C_CKDL200156574_rmdup_CpG.methylKit",
               "D11T4C2H_CKDL200156575_rmdup_CpG.methylKit",
               "D11T4C3C_CKDL200156576_rmdup_CpG.methylKit",
               "D11T4C3H_CKDL200156577_rmdup_CpG.methylKit",
               "D11T9C1C_CKDL200156578_rmdup_CpG.methylKit",
               "D11T9C1H_CKDL200156579_rmdup_CpG.methylKit",
               "D11T9C2C_CKDL200156580_rmdup_CpG.methylKit",
               "D11T9C2H_CKDL200156581_rmdup_CpG.methylKit",
               "D11T9C3C_CKDL200156582_rmdup_CpG.methylKit",
               "D11T9C3H_CKDL200156583_rmdup_CpG.methylKit")


sample.id=list("D0_T10_C_1",  
               "D0_T10_H_1",  
               "D0_T10_C_2",  
               "D0_T10_H_2",  
               "D0_T10_C_3",  
               "D0_T10_H_3",  
               "D0_T3_C_1",   
               "D0_T3_H_1",   
               "D0_T3_C_2",   
               "D0_T3_H_2",   
               "D0_T3_C_3",   
               "D0_T3_H_3",   
               "D0_T4_C_1",   
               "D0_T4_H_1",   
               "D0_T4_C_2",   
               "D0_T4_H_2",   
               "D0_T4_C_3",   
               "D0_T4_H_3",   
               "D0_T9_C_1",   
               "D0_T9_H_1",   
               "D0_T9_C_2",   
               "D0_T9_H_2",   
               "D0_T9_C_3",   
               "D0_T9_H_3",   
               "D11_T10_C_1", 
               "D11_T10_H_1", 
               "D11_T10_C_2", 
               "D11_T10_H_2", 
               "D11_T10_C_3", 
               "D11_T10_H_3", 
               "D11_T3_C_1",  
               "D11_T3_H_1",  
               "D11_T3_C_2",  
               "D11_T3_H_2",
               "D11_T3_C_3",  
               "D11_T3_H_3",  
               "D11_T4_C_1",
               "D11_T4_H_1",  
               "D11_T4_C_2",  
               "D11_T4_H_2",  
               "D11_T4_C_3",  
               "D11_T4_H_3",  
               "D11_T9_C_1",  
               "D11_T9_H_1",  
               "D11_T9_C_2",  
               "D11_T9_H_2",  
               "D11_T9_C_3",  
               "D11_T9_H_3")

acclimation.and.time.treatment <- rep(c(0, 1, 2, 3), each = 12)

# Methylation: methylKit object creation
myobj=methRead(file.list,                   
               sample.id=sample.id,       
               assembly="Amil.v2.01",
               treatment=acclimation.and.time.treatment,
               context="CpG",
               dbtype = "tabix",
               dbdir = "methylDB")

filtered.myobj=filterByCoverage(myobj,
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count=NULL,
                                hi.perc=99.9)

# Methylation: Combine gene data with methylation object
gr.obj <- import("/Users/LeslieDO/LabNotebook/Chapter1/Amil_v2.01/Amil.genes.gff") # This is the Fuller file.
gr.df <- data.frame(gr.obj)
names(gr.df)
# Methylation: Make a column with chr.start.end
gr.df$concat <- paste(gr.df$seqnames,gr.df$start,gr.df$end,sep=".")

# Methylation: Region counts
genes <- regionCounts(filtered.myobj,gr.obj,save.db=FALSE) # Gene counts for individuals

# Methylation: Make methylBase object
genes.unite <- methylKit::unite(genes,
                                destrand = FALSE, #Combine
                                min.per.group = 3L, # in at least 3 replicates per treatment
                                save.db = F)

# Methylation: Create methylation matrix
gene.pm <- percMethylation(genes.unite, rowids=TRUE)

# Methylation: Making data frame
gene.names <- gr.df$ID[match(rownames(gene.pm),gr.df$concat)]
rownames(gene.pm) <- gene.names
gene.pm.means.df <- data.frame(gene = gene.names, 
                               pm.means = rowMeans(gene.pm, na.rm = TRUE), 
                               row.names = NULL)

# Gene Expression and Methylation: Merge data frames
dim(heat.stress.DE.results) # 26965 
dim(gene.pm.means.df) # 12631

exp.and.perc.meth <- inner_join(heat.stress.DE.results,gene.pm.means.df, by = 'gene')
dim(exp.and.perc.meth) 

# Gene expression and methylation: linear regression of the relationship
mod.ge<-lm(log10_baseMean ~ pm.means, data = exp.and.perc.meth)
summary(mod.ge) # significant relationship between the percent methylation mean and the log10 base means


# Gene expression and methylation: Plot

plot_l10bMs_methylation <- ggplot(exp.and.perc.meth, aes(x=pm.means, y=log10_baseMean)) + 
  geom_point(alpha = 0.3) + geom_smooth(method='lm', formula= y~x) +
  theme_classic(base_size = 20) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle(expression(paste("Relationship between log10 transformed expression means\nand percent methylation"))) +
  xlab("Average % methylation") + 
  ylab("log10 transformed expression means") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 8)  
plot_l10bMs_methylation 


 # jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/rel.bet.log10bm.methylation.jpeg', units="in", width=13, height=10, res=300)
 # print(plot_l10bMs_methylation)     
 # dev.off()   

# Heat stress gene expression plasticity (log2FC contrast: heat stress vs control) and
# percent methylation
plot_hs_l2fcSD_methylation <- ggplot(exp.and.perc.meth, aes(x=pm.means, y=lfcSE)) + 
  geom_point(alpha = 0.3) + geom_smooth(method='lm', formula= y~x) +
  theme_classic(base_size = 20) +
  ggtitle("Relationship between heat stress Log2FC SD and percent methylation") +
  xlab("Average % Methylation") + 
  ylab("heat stress Log2FC SD")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3, label.y = 5, size = 8 )
plot_hs_l2fcSD_methylation
# jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/rel.bet.hslog2FC.methylation.jpeg', units="in", width=13, height=10, res=300)
# print(plot_hs_l2fcSD_methylation)
# dev.off()









