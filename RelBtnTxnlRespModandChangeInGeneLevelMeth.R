# Relationship between transcriptional response modulation and change
# change in gene level methylation

# Set Up:
library(DESeq2)
library(methylKit)
library(rtracklayer)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggdist)

# Gene Expression: File set up
my.dir <- "/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/HTSeq_out"
my.files <- grep(".txt", list.files(my.dir), value=TRUE)
my.metadata <- read.csv("/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/HTSeq_out/a_nana_metadata.csv")

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
# Pre-filter
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Gene Expression: Differential Expression Analysis
a.nana.dds <- DESeq(ddsHTSeq) 

# Gene Expression: Transcriptional response modulation data frame
a.nana.dds.d11.sub <- a.nana.dds[, a.nana.dds@colData@listData$condition.accl.time == "11"]

coeff.df <- as.data.frame(coef(a.nana.dds.d11.sub), colnames = TRUE)
wald.stat.df <- as.data.frame(cbind(row.names(coeff.df), 
                                    as.numeric(coeff.df$condition.accl.temp31.condition.stress.treatmentstress)),
                              colnames = TRUE)
wald.stat.df <- wald.stat.df[complete.cases(wald.stat.df),]
colnames(wald.stat.df) <- c("gene", "coefficient")

GE.results <- results(a.nana.dds.d11.sub, name = "condition.accl.temp31.condition.stress.treatmentstress")
GE.results <- GE.results[complete.cases(GE.results),]
GE.results.df <- GE.results %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") 

dim(GE.results.df)

GE.results.interaction.df <- dplyr::left_join(GE.results.df,
                                              wald.stat.df,
                                              by = "gene") %>% 
  mutate(coefficient = as.numeric(coefficient))

GE.results.interaction.df$coefficient <- as.numeric(GE.results.interaction.df$coefficient)
GE.results.interaction.df$category <- ifelse(GE.results.interaction.df$padj < 0.1 & GE.results.interaction.df$coefficient > 0, "Amplified",
                        ifelse(GE.results.interaction.df$padj < 0.1 & GE.results.interaction.df$coefficient < 0, "Dampened",
                               ifelse(GE.results.interaction.df$padj > 0.1, "neither", NA)))


amp.damp.results <- GE.results.interaction.df[GE.results.interaction.df$padj < 0.1,]
amp.damp.results <- amp.damp.results[,c(1,9)]

# Methylation: File set up
rm(my.dir)
my.dir <- "/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/"
meta.file.name <- "meth_meta.csv"
meta.path <- file.path(my.dir, meta.file.name)
meta <- read.csv(meta.path, header = TRUE)
meta$Acclimation <- factor(meta$Acclimation,levels=c("Control","Accl"))
meta$Colony <- as.factor(meta$Colony)

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

acc.file.list <- file.list[c(25:48)]
acc.sample.id <- sample.id[c(25:48)]
acclimation.treatment <- c(rep(0, 12), rep(1, 12))

# Methyltion: Creating methylation object
myobj.acc=methRead(acc.file.list,                   
                   sample.id=acc.sample.id,       
                   assembly="Amil.v2.01",
                   treatment=acclimation.treatment,
                   context="CpG",
                   dbtype = "tabix",
                   dbdir = "methylDB")

filtered.myobj.acc=filterByCoverage(myobj.acc,
                                    lo.count=10,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)

# Methylation: Integrate gene information
gr.obj <- import("/Users/LeslieDO/LabNotebook/Chapter1/Amil_v2.01/Amil.genes.gff") # This is the Fuller file.
gr.df <- data.frame(gr.obj)
names(gr.df)
# Make a column with chr.start.end
gr.df$concat <- paste(gr.df$seqnames,gr.df$start,gr.df$end,sep=".")

# Methylation: Region counts
genes.acc <- regionCounts(filtered.myobj.acc,
                          gr.obj,
                          save.db=F) # Gene counts for individuals


genes.unite.acclimation <- methylKit::unite(genes.acc,
                                            destrand = FALSE, #Combine
                                            min.per.group = 3L, # in at least 3 of the samples per treatment
                                            save.db = F)
pooled.meth.genes.acclimation <- pool(genes.unite.acclimation, 
                                      sample.ids=c("Control","Acclimation"))

dm.pooled.genes.acc=calculateDiffMeth(pooled.meth.genes.acclimation, 
                                      adjust="fdr")

# Methylation: Attach gene names to the data
diff.meth.acc.df  <- getData(dm.pooled.genes.acc)
diff.meth.acc.df$concat <-paste(diff.meth.acc.df$chr,diff.meth.acc.df$start,diff.meth.acc.df$end,sep=".")
diff.meth.acc.df <- diff.meth.acc.df %>%
  dplyr::select(concat, everything())

gene.names.dm <- gr.df$ID[match((diff.meth.acc.df$concat),gr.df$concat)]

rownames(diff.meth.acc.df) <- gene.names.dm

diff.meth.acc.df <- diff.meth.acc.df %>%
  dplyr::select(-concat, - chr, -start, -end)

diff.meth.acc.df <- rownames_to_column(diff.meth.acc.df, var = "rowname")
colnames(diff.meth.acc.df)[1] <- "gene"

# Gene Expression and Methylation: Merge the data frames
txn.mod.response.and.dname.change <- dplyr::inner_join(GE.results.interaction.df,
                                                       diff.meth.acc.df,
                                                       by = 'gene')


txn.mod.response.and.dname.change <- txn.mod.response.and.dname.change[,c('gene',
                                                                          'padj',
                                                                          'coefficient',
                                                                          'category',
                                                                          'meth.diff')]
# sort the rows
sorted.txn.mod.response.and.dname.change <- sorted_df <- txn.mod.response.and.dname.change %>%
  arrange(factor(category, levels = c("neither", "amplified", "dampened")))

# Plot the reordered data frame

MyColor <- c( "#dddddd", "#d8b365",  "#5ab4ac")
names(MyColor) <- c("neither", "Amplified", "Dampened")
diff.meth.vs.plasticity.interaction <- ggplot(sorted.txn.mod.response.and.dname.change, aes(x = meth.diff, y = coefficient, color = category)) +
  geom_hline(yintercept = 0, lty = 2, color = '#585858', lwd = 1) +
  geom_vline(xintercept = 0, lty = 2, color = '#585858', lwd = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = MyColor) +
  xlab("Change in DNA methylation following acclimation") + 
  ylab("Change in gene expression\n response to heat stress") +
  labs(color = NULL) +
  geom_point(size = 4, alpha = 0.5) + 
  stat_ellipse(level = 0.99, linetype = 2, color = "black") 


diff.meth.vs.plasticity.interaction
#jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/diff.meth.vs.plasticity.jpeg', units="in", width=8, height=5, res=300)
#print(diff.meth.vs.plasticity.interaction)     
#dev.off() 

# Testing this relationship with anova
anova(lm(txn.mod.response.and.dname.change$coefficient~txn.mod.response.and.dname.change$meth.diff))
# p val:  0.1812

# Plotting %methylation mean of amplified and dampened genes
# Creating data frame with percent methylation mean of genes

gene.pm.acc <- percMethylation(genes.unite.acclimation, rowids = TRUE)
gene.names <- gr.df$ID[match(rownames(gene.pm.acc),gr.df$concat)]
rownames(gene.pm.acc) <- gene.names # Replace the rownames with the real gene names
pm.means <- rowMeans(gene.pm.acc, na.rm =TRUE)

pm.means <-data.frame(gene = names(pm.means), 
                      mean.perc.meth = as.numeric(pm.means),
                      row.names = NULL)


# Merging data frames with transcriptionally modified genes
amp.damp.pm.means <- dplyr::inner_join(amp.damp.results, 
                                       pm.means,
                                       by = 'gene')

compare_means(mean.perc.meth ~ category, data = amp.damp.pm.means)
amp_damp_meth <- ggplot(amp.damp.pm.means, aes(x = category, y = mean.perc.meth, fill = category)) +
  theme_classic(base_size = 20) +
  theme(legend.position="none") +
  xlab(NULL) + 
  ylab("% methylation mean") +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = c(0.5, 1), alpha = 0.8) + 
  ggdist::stat_dots(side = "left", dotsize = .4, justification = 1.05, binwidth = .1) +
  #stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA, alpha = 0.8) +
  scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  stat_compare_means(label = "p.signif", size = 6)


amp_damp_meth
# jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/meth_amp_damp.jpeg', units="in", width=6, height=5, res=300)
# print(amp_damp_meth)
# dev.off()

# Percent methylation of damp and amp genes 
# of acclimated and control samples.
gene.pm.acc

# Control treatment percent methylation
D11.29C.treat.pm <- gene.pm.acc[,1:12]
D11.29C.treat.pm.means <- rowMeans(D11.29C.treat.pm, na.rm =TRUE)
D11.29C.treat.pm.means <- as.data.frame(D11.29C.treat.pm.means)
D11.29C.treat.pm.means <- rownames_to_column(D11.29C.treat.pm.means, var = "rowname")
colnames(D11.29C.treat.pm.means) <- c("gene", "cont.pm.mean")

# Acclimation treatment percent methylation
D11.31C.treat.pm <- gene.pm.acc[,13:24]
D11.31C.treat.pm.means <- rowMeans(D11.31C.treat.pm, na.rm =TRUE)
D11.31C.treat.pm.means <- as.data.frame(D11.31C.treat.pm.means)
D11.31C.treat.pm.means <- rownames_to_column(D11.31C.treat.pm.means, var = "rowname")
colnames(D11.31C.treat.pm.means) <- c("gene", "acc.pm.mean")

# Merge the data frames 
D11.methylation.means <- dplyr::inner_join(D11.29C.treat.pm.means,
                                           D11.31C.treat.pm.means, 
                                           by = 'gene')

# Merge methylation and transcriptional mod genes
acc.gene.perc.meth.df <- inner_join(amp.damp.results,
                                       D11.methylation.means,
                                       by = 'gene')

long.trans.mod.perc.met<- acc.gene.perc.meth.df %>% pivot_longer(
  cols = cont.pm.mean:acc.pm.mean,
  names_to = c("Acc.treatment"),
  values_to = "Perc.Meth"
)

# Plot
trans.mod.boxplot <- long.trans.mod.perc.met %>% 
  ggplot(aes(x = category, y = Perc.Meth, fill = Acc.treatment)) + 
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_blank())+ 
  scale_fill_brewer(name = "Acclimation Treatment", labels = c("Acclimation", "Control")) +
  geom_boxplot(color="black", alpha = 0.7) + 
  labs(y = "% Gene Body Methylation") + 
  guides(color="none") + 
  stat_compare_means(label = "p.signif", size = 6)

trans.mod.boxplot 
jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/SuppFigure4.jpeg', units="in", width=5, height=5, res=300)
print(trans.mod.boxplot)     
dev.off() 






