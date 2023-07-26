# Effects of Treatment on DNA methylation
library(methylKit)
library(vegan)
library(ggplot2)
library(plotly)
library(tidyr)
library(rtracklayer)
library(topGO)

# Load project meta data
my.dir <- "/Users/LeslieDO/LabNotebook/Chapter1/Mech_thermal_tol_plasticity/"
meta.file.name <- "meth_meta.csv"
meta.path <- file.path(my.dir, meta.file.name)
meta <- read.csv(meta.path, header = TRUE)
meta$Acclimation <- factor(meta$Acclimation,levels=c("Control","Accl"))
meta$Colony <- as.factor(meta$Colony)


# Files for methylKit
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

# Methylation object creation
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

meth=methylKit::unite(filtered.myobj, destrand=TRUE)
pm <- percMethylation(meth)


# PCoA
pcoa.data <- t(pm)
pcoa.data <- data.frame(pcoa.data)
names(pcoa.data) <- 1:ncol(pcoa.data)
pcoa.data <- cbind(pcoa.data, meta[, 2:5])

pcoa.bray <- vegdist(pcoa.data[,1:18716], method = "bray")
opar <- par(mfrow = c(1, 2), mar = c(3.1, 3.1, 3.1, 5.1), mgp = c(2, 0.5, 0), oma = c(0, 0, 0, 4))
pcoa.result <- aPCoA(pcoa.bray ~ Acclimation, pcoa.data, Acclimation)
pcoa.result.df <- cbind(as.data.frame(pcoa.result), meta[, 2:5])
pcoa.result.df$row_names <- rownames(pcoa.result.df)

# Plot PCoA Results
Acclimation.me.PCoA.plot <- ggplot(pcoa.result.df,aes(x=pcoa.result.df[,1],y=pcoa.result.df[,2], color = Acclimation, shape = Heat, text = row_names)) +
  theme_classic(base_size = 20) +
  theme(legend.position="none") + 
  scale_color_manual(values = c("#ef8a62", "#67a9cf"), labels = c("Control", "Acclimation")) +
  xlab("PCoA1") + 
  ylab("PCoA2") +
  labs(title="PCoA of DNA methylation", color = "Acclimation treatment", shape = "Assay Treatment") +
  geom_point(size = 5, alpha = 0.7)

Acclimation.me.PCoA.plot
# jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/gmeth_pcoa.jpeg', units="in", width=10, height=10, res=300)
# print(Acclimation.me.PCoA.plot)     
# dev.off() 

# Interactive plot to ID outliers
Acclimation.me.PCoA.plot.interactive <- ggplotly(Acclimation.me.PCoA.plot,  tooltip = "text")
Acclimation.me.PCoA.plot.interactive

# Removing outliers
outliers <- c("D11_T10_H_1", "D11_T3_C_1", "D11_T9_H_3", "D11_T4_H_1")
pcoa.data.no.outliers <- pcoa.data[!(row.names(pcoa.data) %in% outliers),]
meta.no.outlier <- meta[!(meta$ID %in% outliers),]

pcoa.bray.no.outliers <- vegdist(pcoa.data.no.outliers[,1:18716], method = "bray")

opar<-par(mfrow=c(1,2),
          mar=c(3.1, 3.1, 3.1, 5.1),
          mgp=c(2, 0.5, 0),
          oma=c(0, 0, 0, 4))
pcoa.result.no.outliers<-aPCoA(pcoa.bray.no.outliers~Acclimation,pcoa.data.no.outliers, Acclimation)
pcoa.result.no.outliers.df <- cbind(as.data.frame(pcoa.result.no.outliers), meta.no.outlier[,2:5])

Acclimation.me.PCoA.plot.removed.outliers <- ggplot(pcoa.result.no.outliers.df,aes(x=pcoa.result.no.outliers.df[,1],y=pcoa.result.no.outliers.df[,2], color = Acclimation, shape = Heat)) +
  theme_classic(base_size = 20) +
  theme(legend.position="none") + 
  scale_color_manual(values = c("#7943d7", "#d77943")) +
  xlab("PCoA1") + 
  ylab("PCoA2") +
  labs(title="DNA methylation",
       color = "Acclimation Treatment") +
  geom_point(size = 5, alpha = 0.7)

Acclimation.me.PCoA.plot.removed.outliers
# jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/gmeth_pcoa_no_outliers.jpeg', units="in", width=5, height=5, res=300)
# print(Acclimation.me.PCoA.plot.removed.outliers)     
# dev.off() 

# Percent methylation across CpG sites
meth.min.3.samps=methylKit::unite(filtered.myobj, 
                      destrand=TRUE,
                      min.per.group = 3L, # in at least 3 of replicates per treatment
                      save.db = F)

meth.min.3.samps.pm <- percMethylation(meth.min.3.samps)
min.3.pm.df <- as.data.frame(meth.min.3.samps.pm)
min.3.pm.df.long <- min.3.pm.df %>%
                   gather(key = "sample", value = "percent", na.rm = TRUE)


hist.cpg.meth.perc <- ggplot(min.3.pm.df.long, aes(x = percent, y = stat(count))) +
                      geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
                      theme_classic(base_size = 20) +
                      labs(title = "Histogram of percent methylation of CpG sites", 
                           x = "Percent methylation",
                           y = "Frequency")
hist.cpg.meth.perc
#jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/hist.cpg.meth.perc.jpeg', units="in", width=12, height=10, res=300)
#print(hist.cpg.meth.perc)     
#dev.off()

# Percent methylation across genes
gr.obj <- import("/Users/LeslieDO/LabNotebook/Chapter1/Amil_v2.01/Amil.genes.gff") # This is the Fuller file.
gr.df <- data.frame(gr.obj)
names(gr.df)
# Make a column with chr.start.end
gr.df$concat <- paste(gr.df$seqnames,gr.df$start,gr.df$end,sep=".")

# Region counts
genes <- regionCounts(filtered.myobj,gr.obj,save.db=FALSE) # Gene counts for individuals

genes.unite <- methylKit::unite(genes,
                                destrand = FALSE, #Combine
                                min.per.group = 3L, # in at least 3 replicates per treatment
                                save.db = F)

gene.pm <- percMethylation(genes.unite)
gene.pm.df <- as.data.frame(gene.pm)
gene.pm.df.long <- gene.pm.df %>%
  gather(key = "sample", value = "percent", na.rm = TRUE)


hist.gene.meth.perc <- ggplot(gene.pm.df.long, aes(x = percent, y = stat(count))) +
  geom_histogram(binwidth = 5, fill = "#555599" , color = "white") +
  theme_classic(base_size = 20) +
  labs(title = "Histogram of percent methylation of genes", 
       x = "Percent methylation",
       y = "Frequency")
hist.gene.meth.perc
# jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/hist.gene.meth.perc.jpeg', units="in", width=12, height=10, res=300)
# print(hist.gene.meth.perc)     
# dev.off()

# Effect of treatment on DNA methylation at the global level
## Acclimation temperature and time treatment
ind.meth.acc.temp.and.time <- colMeans(meth.min.3.samps.pm,na.rm=TRUE)
anova(lm(ind.meth.acc.temp.and.time~as.factor(meth.min.3.samps@treatment))) # p val: 0.1652

## Heat stress
heat.stress.treatment <- rep(c(0, 1), length.out = 48)

myobj.hs=methRead(file.list,                   
               sample.id=sample.id,       
               assembly="Amil.v2.01",
               treatment=heat.stress.treatment,
               context="CpG",
               dbtype = "tabix",
               dbdir = "methylDB")

filtered.myobj.hs=filterByCoverage(myobj.hs,
                                lo.count=10,
                                lo.perc=NULL,
                                hi.count=NULL,
                                hi.perc=99.9)

meth.sites.hs=methylKit::unite(filtered.myobj.hs, 
                               destrand=TRUE, 
                               min.per.group = 3L,
                               save.db = FALSE)

pm.hs <- percMethylation(meth.sites.hs)
ind.meth.hs <- colMeans(pm.hs,na.rm=T)
anova(lm(ind.meth.hs~as.factor(meth.sites.hs@treatment))) # p-value 0.378

## Acclimation (Day 11 29 samples vs Day 11 31)

acc.file.list <- file.list[c(25:48)]
acc.sample.id <- sample.id[c(25:48)]
acclimation.treatment <- c(rep(0, 12), rep(1, 12))

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

meth.sites.acc=methylKit::unite(filtered.myobj.acc, 
                               destrand=TRUE, 
                               min.per.group = 3L,
                               save.db = FALSE)

pm.acc <- percMethylation(meth.sites.acc)
ind.meth.acc <- colMeans(pm.acc,na.rm=TRUE)
anova(lm(ind.meth.acc~as.factor(meth.sites.acc@treatment))) # p val:0.4339

# Testing whether DNA methylation percent varies between genomic features 
# and if there are feature-specific changes in methylation following acclimation.
# To avoid confounding feature of day, using heat-stress treatment and control 
# samples from the acclimation control (D11 29C), and the acclimation treatment
# (Day 11 31C)

# Exon and Intron data from reference genome source: Fuller et al 2020 
Amil.coding.bed <- readTranscriptFeatures("~/LabNotebook/Chapter1/Amil_v2.01/Amil_coding.bed")
Amil.coding.bed$exons 
Amil.coding.bed$introns

# BED file source: https://github.com/Groves-Dixon-Matz-laboratory/benchmarking_coral_methylation/tree/master/windowStats
# Load in BED files 
line.repeats <- read.table("~/LabNotebook/Chapter1/Amil_v2.01/groves_2020_beds/LINE_repeats_rmdups.bed")
colnames(line.repeats) <- c("chr", "start", "end", "name")

sine.repeats <- read.table("~/LabNotebook/Chapter1/Amil_v2.01/groves_2020_beds/SINE_repeats_rmdups.bed")
colnames(sine.repeats) <- c("chr", "start", "end", "name")

tss <- read.table("~/LabNotebook/Chapter1/Amil_v2.01/groves_2020_beds/tssBoundaries.bed")
colnames(tss) <- c("chr", "start", "end", "name")

rc.repeats <- read.table("~/LabNotebook/Chapter1/Amil_v2.01/groves_2020_beds/RC_repeats_rmdups.bed")
colnames(rc.repeats) <- c("chr", "start", "end", "name")

promoter <- read.table("~/LabNotebook/Chapter1/Amil_v2.01/groves_2020_beds/promoterBoundaries.bed")
colnames(promoter) <- c("chr", "start", "end", "name")


# Assign coordinates based on window inputs
line.seqnames = line.repeats$chr
line.start=line.repeats$start
line.end=line.repeats$end
line.names=line.repeats$name
line.ranges <- IRanges(start = line.repeats$start, end = line.repeats$end)
line.repeats$concat <- paste(line.seqnames,line.start,line.end,sep=".")

sine.seqnames = sine.repeats$chr
sine.start=sine.repeats$start
sine.end=sine.repeats$end
sine.names=sine.repeats$name
sine.ranges <- IRanges(start = sine.repeats$start, end = sine.repeats$end)
sine.repeats$concat <- paste(sine.seqnames,sine.start,sine.end,sep=".")

tss.seqnames = tss$chr
tss.start=tss$start
tss.end=tss$end
tss.names=tss$name
tss.ranges <- IRanges(start = tss$start, end = tss$end)
tss$concat <- paste(tss.seqnames,tss.start,tss.end,sep=".")

rc.seqnames = rc.repeats$chr
rc.start=rc.repeats$start
rc.end=rc.repeats$end
rc.names=rc.repeats$name
rc.ranges <- IRanges(start = rc.repeats$start, end = rc.repeats$end)
rc.repeats$concat <- paste(rc.seqnames,rc.start,rc.end,sep=".")

promoter.seqnames = promoter$chr
promoter.start=promoter$start
promoter.end=promoter$end
promoter.names=promoter$name
promoter.ranges <- IRanges(start = promoter$start, end = promoter$end)
promoter$concat <- paste(promoter.seqnames,promoter.start,promoter.end,sep=".")

#create the Granges object specifying the regions
line.win <- GRanges(seqnames = line.seqnames, ranges = line.ranges)
sine.win <- GRanges(seqnames = sine.seqnames, ranges = sine.ranges)
tss.win <- GRanges(seqnames = tss.seqnames, ranges = tss.ranges)
rc.win <- GRanges(seqnames = rc.seqnames, ranges = rc.ranges)
promoter.win <- GRanges(seqnames = promoter.seqnames, ranges = promoter.ranges)

# Create counts object for regions
exon.regions = regionCounts(object = filtered.myobj,
                            regions = Amil.coding.bed$exons,
                            save.db = FALSE)

intron.regions = regionCounts(object = filtered.myobj,
                              regions = Amil.coding.bed$introns,
                              save.db = FALSE)

line.regions = regionCounts(object = filtered.myobj,
                            regions = line.win,
                            save.db = FALSE)

sine.regions = regionCounts(object = filtered.myobj,
                            regions = sine.win,
                            save.db = FALSE)

tss.regions = regionCounts(object = filtered.myobj,
                           regions = tss.win,
                           save.db = FALSE)

rc.regions = regionCounts(object = filtered.myobj,
                          regions = rc.win,
                          save.db = FALSE)

promoter.regions = regionCounts(object = filtered.myobj,
                                regions = promoter.win,
                                save.db = FALSE)


# Create MethylKit objects, minimum 3 samples per treatment
exon.unite.min = methylKit::unite(exon.regions,
                                  min.per.group = 3L,
                                  save.db = FALSE)

intron.unite.min = methylKit::unite(intron.regions,
                                    min.per.group = 3L,
                                    save.db = FALSE)

line.unite.min = methylKit::unite(line.regions, 
                                  min.per.group = 3L,
                                  save.db = FALSE)

sine.unite.min = methylKit::unite(sine.regions,
                                  min.per.group = 3L,
                                  save.db = FALSE)

tss.unite.min = methylKit::unite(tss.regions,
                                 min.per.group = 3L,
                                 save.db = FALSE)

rc.unite.min = methylKit::unite(rc.regions,
                                min.per.group = 3L,
                                save.db = FALSE)

promoter.unite.min = methylKit::unite(promoter.regions,
                                      min.per.group = 3L,
                                      save.db = FALSE)

# Make a percent methylation matrix for each genomic feature
exon.percent.meth <- percMethylation(exon.unite.min, rowids=TRUE)
dim(exon.percent.meth) #4728

intron.percent.meth <- percMethylation(intron.unite.min)
dim(intron.percent.meth) #13931

line.percent.meth <- percMethylation(line.unite.min, rowids=TRUE)
dim(line.percent.meth) #4760

sine.percent.meth <- percMethylation(sine.unite.min, rowids=TRUE)
dim(sine.percent.meth) #28

tss.percent.meth <- percMethylation(tss.unite.min,rowids=TRUE)
dim(tss.percent.meth) #1825

rc.percent.meth <- percMethylation(rc.unite.min, rowids=TRUE)
dim(rc.percent.meth) #142

promoter.percent.meth <- percMethylation(promoter.unite.min, rowids=TRUE)
dim(promoter.percent.meth) #3638

# Genomic features data frames with the average %methylation calculated across 
# samples

# Exon
exon.df <- as.data.frame(exon.percent.meth)
exon.df$mean <- rowMeans(exon.df, na.rm = TRUE)
exon.df$feature <- c(rep("Exon", times = nrow(exon.df)))
dim(exon.df) #4728

# Intron
intron.df <- as.data.frame(intron.percent.meth)
intron.df$mean <- rowMeans(intron.df, na.rm = TRUE)
intron.df$feature <- c(rep("Intron", times = nrow(intron.df)))
dim(intron.df) #13931

# Lines df
lines.df <- as.data.frame(line.percent.meth)
lines.df$mean <- rowMeans(lines.df, na.rm = TRUE)
lines.df$feature <- c(rep("LINE repeat", times = nrow(lines.df)))
dim(lines.df) #4760

# Sines df
sines.df <- as.data.frame(sine.percent.meth)
sines.df$mean <- rowMeans(sines.df, na.rm = TRUE)
sines.df$feature <- c(rep("SINE repeat", times = nrow(sines.df)))
dim(sines.df) #28

# TSS df
tss.df <- as.data.frame(tss.percent.meth)
tss.df$mean <- rowMeans(tss.df, na.rm = TRUE)
tss.df$feature <- c(rep("TSS", times = nrow(tss.df)))
dim(tss.df) #1825

# RC df
rc.df <- as.data.frame(rc.percent.meth, )
rc.df$mean <- rowMeans(rc.df, na.rm=TRUE)
rc.df$feature <- c(rep("RC repeat", times = nrow(rc.df)))
dim(rc.df) #142

# Promoter df
promoter.df <- as.data.frame(promoter.percent.meth)
promoter.df$mean <- rowMeans(promoter.df, na.rm = TRUE)
promoter.df$feature <- c(rep("Promoter", times = nrow(promoter.df)))
dim(promoter.df) #3638 

# Make a single data frame from the percent methylation matrices
features.methylation.df <- rbind(promoter.df,
                                 tss.df,
                                 exon.df, 
                                 intron.df, 
                                 lines.df, 
                                 sines.df, 
                                 rc.df)

# Plot the mean methylation of the features
features.colors <- c("#8476dd","#5170c2", "#2d66a0", "#fa6ef5", "#2a4858", "#23587b", "#bd75ef")
methylation.features.violin.plot <- ggplot(features.methylation.df, aes(x=feature, y=mean, color=feature)) + 
                                    geom_violin(linewidth = 1.5) +
                                    theme_classic(base_size = 20) +
                                    geom_boxplot(width=0.1) +
                                    labs(title="Percent methylation per feature",
                                         x="Genomic Feature", 
                                         y = "% Methylation") +
                                    theme(legend.position = "none") +
                                    scale_color_manual(values = features.colors) +
                                    scale_x_discrete(limits=c("Promoter", 
                                                              "TSS",
                                                              "Exon",
                                                              "Intron",
                                                              "LINE repeat",
                                                              "SINE repeat",
                                                              "RC repeat")) +
                                    theme( axis.text.x = element_text(size=14)) +
                                    geom_hline(yintercept = median(features.methylation.df$mean),
                                               linetype = "dashed", color = "lightgray")

methylation.features.violin.plot 

#jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/methylation.features.violin.plot.jpeg', units="in", width=15, height=10, res=300)
#print(methylation.features.violin.plot)     
#dev.off()

# Test if DNA methylation varies between genomic features
# Reshape data
long.features.df<- features.methylation.df %>% 
  pivot_longer(
  cols = D0_T10_C_1:D11_T9_H_3,
  names_to = c("Day","Tank", "Treatment", "Replicate"), 
  names_pattern = "D(.*)_T(.*)_(.*)_(.*)",
  values_to = "Perc.Meth"
) %>%
  mutate(Acclimation = ifelse(Tank %in% c("4", "9") &
                                Day %in% "11",
                              "Acclimation", "Control"
  ))

# perform ANOVA on the dataset
feature.model <- aov(mean ~ feature, data = long.features.df)

# print the ANOVA table
summary(feature.model) # <2e-16 ***

# We see there are differences between gene features
# and % methylation. Perform post-hoc test to tell 
# us which specific groups are significantly different
# from each other
# run Tukey's HSD post-hoc test
posthoc <- TukeyHSD(feature.model)

# print the post-hoc test results
posthoc

# Summary: 
# There are significant differences in the mean values between most pairs of features.
# Intron has significantly higher mean values compared to Exon, LINE repeat, Promoter, RC repeat, SINE repeat, and TSS.
# Exon has significantly higher mean values compared to TSS, but not significantly different from LINE repeat or Promoter.

# Test if acclimation affects methylation in the genomic features
acclimation.affect.feature.methylation <- ggplot(long.features.df, aes(fill = Acclimation,x=feature, y=Perc.Meth)) + 
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(size=10)) +
  geom_boxplot(width=0.4) +
  labs(title="Affect of acclimation on percent methylation of features",x="Genomic Feature", y = "% Methylation") +
  scale_x_discrete(limits=c("Promoter", 
                            "TSS",
                            "Exon",
                            "Intron",
                            "LINE repeat",
                            "SINE repeat",
                            "RC repeat")) +
  scale_fill_brewer(palette = "BuGn") +
  stat_compare_means(label = "p.signif", method = 'wilcox.test') +
  guides(fill = guide_legend(title = NULL))

acclimation.affect.feature.methylation
# jpeg('/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/Figures/acclimation.affect.feature.methylation.jpeg', units="in", width=15, height=10, res=300)
# print(acclimation.affect.feature.methylation)     
# dev.off()

# Differential methylation (dm) of genes 
## Heat Stress
genes.hs <- regionCounts(filtered.myobj.hs,
                         gr.obj,
                         save.db=F) # Gene counts for individuals


genes.unite.heat.stress <- methylKit::unite(genes.hs,
                                            destrand = FALSE, #Combine
                                            min.per.group = 3L, # in at least 3 of the samples per treatment
                                            save.db = F)

# Number of genes where methylation can be summarized in at least 3 replicates per treatment
dim(genes.unite.heat.stress) # 16531

pooled.meth.genes.heat.stress <- pool(genes.unite.heat.stress, 
                                      sample.ids=c("Control","Heat Stress"))

dm.pooled.genes.hs=calculateDiffMeth(pooled.meth.genes.heat.stress, 
                                     adjust="fdr")

# get differentially methylated regions with specific cutoffs
all.diff.gene.hs = getMethylDiff(dm.pooled.genes.hs, 
                                 difference=25,
                                 qvalue=0.05,
                                 type="all") #451 genes

hyper.diff.gene.hs = getMethylDiff(dm.pooled.genes.hs, 
                                   difference=25,
                                   qvalue=0.05,
                                   type="hyper") #205 genes

hypo.diff.gene.hs = getMethylDiff(dm.pooled.genes.hs, 
                                  difference=25,
                                  qvalue=0.05,
                                  type="hypo") #246

## Attach gene names to output
diff.meth.hs.df  <- getData(all.diff.gene.hs)
diff.meth.hs.df$concat <-paste(diff.meth.hs.df$chr,diff.meth.hs.df$start,diff.meth.hs.df$end,sep=".")
diff.meth.hs.df <- diff.meth.hs.df %>%
  dplyr::select(concat, everything())

gene.names.dm <- gr.df$ID[match((diff.meth.hs.df$concat),gr.df$concat)]

rownames(diff.meth.hs.df) <- gene.names.dm

diff.meth.hs.df <- diff.meth.hs.df %>%
  dplyr::select(-concat, - chr, -start, -end)

diff.meth.hs.df <- rownames_to_column(diff.meth.hs.df, var = "rowname")
colnames(diff.meth.hs.df)[1] <- "gene"

#write.csv(diff.meth.hs.df ,"/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/diff.meth.hs.df.csv", row.names = FALSE)

# topGO of differentially methylated genes after heat stress
# Set-up
# Read in annotation table and reformat GO column
#annos <- read.csv("~/LabNotebook/Chapter1/Amil_v2.01/Amillepora_trinotate_annotation_report.csv",na.strings=".")
#annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
#annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
#GOmap <- annos[,c("transcript_id","GOlist")]
# write.table(GOmap,"~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
geneID2GO <- readMappings(file="~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",sep="\t",IDsep=",")

# Make 'background' gene set from `genes.unite.heat.stress`
background.hs.genes <- getData(genes.unite.heat.stress)
background.hs.genes$concat <-paste(background.hs.genes$chr,background.hs.genes$start,background.hs.genes$end,sep=".")
background.hs.genes <- background.hs.genes %>%
  dplyr::select(concat, everything())

rm(gene.names.dm)
gene.names.dm <- gr.df$ID[match((background.hs.genes$concat),gr.df$concat)]

rownames(background.hs.genes) <- gene.names.dm

background.hs.genes <- background.hs.genes %>%
  dplyr::select(-concat, - chr, -start, -end)

background.hs.genes<- rownames_to_column(background.hs.genes, var = "rowname")
colnames(background.hs.genes)[1] <- "gene"

background.hs.genes <- background.hs.genes[,1]
#write.csv(background.hs.genes ,"/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.meth.set.hs.genes.csv", row.names = FALSE)

all <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.meth.set.hs.genes.csv")
allgenes <- all$gene

# Modify dm heat stress gene set for topGO analysis
## Note: Used VLOOKUP in Excel to modify dm heat stress data set

dm.hs.genes <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/diff.meth.hs.df.csv")
dm.hs.genesIG <- factor(as.numeric(allgenes%in%dm.hs.genes$gene))
names(dm.hs.genesIG) <- allgenes

##BP
GOdm.hs <- new("topGOdata",ontology="BP",allGenes=dm.hs.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dm.hs.Fisher <- getSigGroups(GOdm.hs,test.stat)
dm.hs.res <- GenTable(GOdm.hs,classic=dm.hs.Fisher,topNodes=length(dm.hs.Fisher@score),numChar=100)
dm.hs.filt <- dm.hs.res[dm.hs.res$classic<0.01 & dm.hs.res$Significant>=10,]
#write.csv(dm.hs.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/dm.hs.genes.GO_BP.csv",sep=""))

##MF
GOdm.hs <- new("topGOdata",ontology="MF",allGenes=dm.hs.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dm.hs.Fisher <- getSigGroups(GOdm.hs,test.stat)
dm.hs.res <- GenTable(GOdm.hs,classic=dm.hs.Fisher,topNodes=length(dm.hs.Fisher@score),numChar=100)
dm.hs.filt <- dm.hs.res[dm.hs.res$classic<0.01 & dm.hs.res$Significant>=10,]
#write.csv(dm.hs.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/dm.hs.genes.GO_MF.csv",sep=""))

##CC
GOdm.hs <- new("topGOdata",ontology="CC",allGenes=dm.hs.genesIG, nodeSize=10,
             annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dm.hs.Fisher <- getSigGroups(GOdm.hs,test.stat)
dm.hs.res <- GenTable(GOdm.hs,classic=dm.hs.Fisher,topNodes=length(dm.hs.Fisher@score),numChar=100)
dm.hs.filt <- dm.hs.res[dm.hs.res$classic<0.01 & dm.hs.res$Significant>=10,]
#write.csv(dm.hs.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/dm.hs.genes.GO_CC.csv",sep=""))

# Differential methylation (dm) of genes 
## Acclimation

genes.acc <- regionCounts(filtered.myobj.acc,
                         gr.obj,
                         save.db=F) # Gene counts for individuals


genes.unite.acclimation <- methylKit::unite(genes.acc,
                                            destrand = FALSE, #Combine
                                            min.per.group = 3L, # in at least 3 of the samples per treatment
                                            save.db = F)

# Number of genes where methylation can be summarized in at least 3 replicates per treatment
dim(genes.unite.acclimation) #13761

pooled.meth.genes.acclimation <- pool(genes.unite.acclimation, 
                                      sample.ids=c("Control","Acclimation"))

dm.pooled.genes.acc=calculateDiffMeth(pooled.meth.genes.acclimation, 
                                     adjust="fdr")

# get differentially methylated regions with specific cutoffs
all.diff.gene.acc = getMethylDiff(dm.pooled.genes.acc, 
                                 difference=25,
                                 qvalue=0.05,
                                 type="all") # 416 genes

hyper.diff.gene.acc = getMethylDiff(dm.pooled.genes.acc, 
                                   difference=25,
                                   qvalue=0.05,
                                   type="hyper") # 205 genes

hypo.diff.gene.acc = getMethylDiff(dm.pooled.genes.acc, 
                                  difference=25,
                                  qvalue=0.05,
                                  type="hypo") # 211 genes

## Attach gene names to output
diff.meth.acc.df  <- getData(all.diff.gene.acc)
diff.meth.acc.df$concat <-paste(diff.meth.acc.df$chr,diff.meth.acc.df$start,diff.meth.acc.df$end,sep=".")
diff.meth.acc.df <- diff.meth.acc.df %>%
  dplyr::select(concat, everything())

rm(gene.names.dm)
gene.names.dm <- gr.df$ID[match((diff.meth.acc.df$concat),gr.df$concat)]

rownames(diff.meth.acc.df) <- gene.names.dm

diff.meth.acc.df <- diff.meth.acc.df %>%
  dplyr::select(-concat, - chr, -start, -end)

diff.meth.acc.df <- rownames_to_column(diff.meth.acc.df, var = "rowname")
colnames(diff.meth.acc.df)[1] <- "gene"

#write.csv(diff.meth.acc.df ,"/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/diff.meth.acc.df.csv", row.names = FALSE)

# topGO of differentially methylated genes after acclimation
# Set-up
# Read in annotation table and reformat GO column
#annos <- read.csv("~/LabNotebook/Chapter1/Amil_v2.01/Amillepora_trinotate_annotation_report.csv",na.strings=".")
#annos$GOlist <- str_replace_all(annos$gene_ontology_blast,"(GO:[0-9]*)\\^.*?\\`","\\1,")
#annos$GOlist <- str_replace_all(annos$GOlist,"\\^.*?$","")
#GOmap <- annos[,c("transcript_id","GOlist")]
# write.table(GOmap,"~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")
geneID2GO <- readMappings(file="~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/Amil_GOmap.txt",sep="\t",IDsep=",")

# Make 'background' gene set from `genes.unite.acclimation`
background.acc.genes <- getData(genes.unite.acclimation)
background.acc.genes$concat <-paste(background.acc.genes$chr,background.acc.genes$start,background.acc.genes$end,sep=".")
background.acc.genes <- background.acc.genes %>%
  dplyr::select(concat, everything())

rm(gene.names.dm)
gene.names.dm <- gr.df$ID[match((background.acc.genes$concat),gr.df$concat)]

rownames(background.acc.genes) <- gene.names.dm

background.acc.genes <- background.acc.genes %>%
  dplyr::select(-concat, - chr, -start, -end)

background.acc.genes<- rownames_to_column(background.acc.genes, var = "rowname")
colnames(background.acc.genes)[1] <- "gene"

background.acc.genes <- background.acc.genes[,1]
#write.csv(background.acc.genes ,"/Users/LeslieDO/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.meth.set.acc.genes.csv", row.names = FALSE)

rm(all,allgenes)
all <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/total.meth.set.acc.genes.csv")
allgenes <- all$gene

# Modify dm acclimation gene set for topGO analysis
## Note: Used VLOOKUP in Excel to modify dm acclimation data set

dm.acc.genes <- read.csv("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/diff.meth.acc.df.csv")
dm.acc.genesIG <- factor(as.numeric(allgenes%in%dm.acc.genes$gene))
names(dm.acc.genesIG) <- allgenes

##BP
GOdm.acc <- new("topGOdata",ontology="BP",allGenes=dm.acc.genesIG, nodeSize=10,
               annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dm.acc.Fisher <- getSigGroups(GOdm.acc,test.stat)
dm.acc.res <- GenTable(GOdm.acc,classic=dm.acc.Fisher,topNodes=length(dm.acc.Fisher@score),numChar=100)
dm.acc.filt <- dm.acc.res[dm.acc.res$classic<0.01 & dm.acc.res$Significant>=10,]
#write.csv(dm.acc.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/dm.acc.genes.GO_BP.csv",sep=""))

##MF
GOdm.acc <- new("topGOdata",ontology="MF",allGenes=dm.acc.genesIG, nodeSize=10,
               annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dm.acc.Fisher <- getSigGroups(GOdm.acc,test.stat)
dm.acc.res <- GenTable(GOdm.acc,classic=dm.acc.Fisher,topNodes=length(dm.acc.Fisher@score),numChar=100)
dm.acc.filt <- dm.acc.res[dm.acc.res$classic<0.01 & dm.acc.res$Significant>=10,]
#write.csv(dm.acc.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/dm.acc.genes.GO_MF.csv",sep=""))

##CC
GOdm.acc <- new("topGOdata",ontology="CC",allGenes=dm.acc.genesIG, nodeSize=10,
               annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
dm.acc.Fisher <- getSigGroups(GOdm.acc,test.stat)
dm.acc.res <- GenTable(GOdm.acc,classic=dm.acc.Fisher,topNodes=length(dm.acc.Fisher@score),numChar=100)
dm.acc.filt <- dm.acc.res[dm.acc.res$classic<0.01 & dm.acc.res$Significant>=10,]
#write.csv(dm.acc.filt,paste("~/LabNotebook/Chapter1/Patterns_of_methylation_and_transcriptional_plasticity/DataOutput/dm.acc.genes.GO_CC.csv",sep=""))




