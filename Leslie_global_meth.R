library(methylKit)
library(beeswarm)
library(genomation)
library(GRanges)

setwd("/Users/LeslieDO/OneDrive - UC Davis/A_nana_DE_DM/")

meta <- read.csv("meth_meta.csv", header = TRUE)
meta$Acclimation <- factor(meta$Acclimation,levels=c("Control","Accl"))

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

myobj=methRead(file.list,                    # Corresponding File prefix:
                    sample.id=list("D0_T10_C_1",  # D0T10C1C
                                   "D0_T10_H_1",  # D0T10C1H
                                   "D0_T10_C_2",  # D0T10C2C
                                   "D0_T10_H_2",  # D0T10C2H
                                   "D0_T10_C_3",  # D0T10C3C
                                   "D0_T10_H_3",  # D0T10C3H
                                   "D0_T3_C_1",   # D0T3C1C
                                   "D0_T3_H_1",   # D0T3C1H
                                   "D0_T3_C_2",   # D0T3C2C 
                                   "D0_T3_H_2",   # D0T3C2H
                                   "D0_T3_C_3",   # D0T3C3C
                                   "D0_T3_H_3",   # D0T3C3H
                                   "D0_T4_C_1",   # D0T4C1C
                                   "D0_T4_H_1",   # D0T4C1H
                                   "D0_T4_C_2",   # D0T4C2C
                                   "D0_T4_H_2",   # D0T4C2H 
                                   "D0_T4_C_3",   # D0T4C3C
                                   "D0_T4_H_3",   # D0T4C3H
                                   "D0_T9_C_1",   # D0T9C1C
                                   "D0_T9_H_1",   # D0T9C1H
                                   "D0_T9_C_2",   # D0T9C2C
                                   "D0_T9_H_2",   # D0T9C2H
                                   "D0_T9_C_3",   # D0T9C3C
                                   "D0_T9_H_3",   # D0T9C3H
                                   "D11_T10_C_1", # D11T0C1C
                                   "D11_T10_C_2", # D11T0C2C
                                   "D11_T10_H_2", # D11T0C2H
                                   "D11_T10_C_3", # D11T0C3C
                                   "D11_T10_H_3", # D11T0C3H
                                   "D11_T3_C_1",  # D11T3C1C
                                   "D11_T3_H_1",  # D11T3C1H
                                   "D11_T3_C_2",  # D11T3C2C
                                   "D11_T3_H_2",  # D11T3C2H
                                   "D11_T3_C_3",  # D11T3C3C
                                   "D11_T3_H_3",  # D11T3C3H
                                   "D11_T4_C_1",  # D11T4C1C
                                   "D11_T4_H_1",  # D11T4C1H
                                   "D11_T4_C_2",  # D11T4C2C
                                   "D11_T4_H_2",  # D11T4C2H
                                   "D11_T4_C_3",  # D11T4C3C
                                   "D11_T4_H_3",  # D11T4C3H
                                   "D11_T9_C_1",  # D11T9C1C
                                   "D11_T9_H_1",  # D11T9C1H
                                   "D11_T9_C_2",  # D11T9C2C
                                   "D11_T9_H_2",  # D11T9C2H
                                   "D11_T9_C_3",  # D11T9C3C
                                   "D11_T9_H_3"), # D11T9C3H         
                    assembly="Amil.v2.01",
                    treatment=c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3),
                    context="CpG",
                    dbtype = "tabix",
                    dbdir = "methylDB")

filtered.myobj=filterByCoverage(myobj,
                                   lo.count=10,
                                   lo.perc=NULL,
                                   hi.count=NULL,
                                   hi.perc=99.9)


meth.sites=methylKit::unite(filtered.myobj, 
                            destrand=FALSE, 
                            min.per.group = 4L,
                            save.db = TRUE,
                            suffix = "meth.sites",
                            dbtype = "tabix",
                            dbdir = "methylBaseDB")                                   


# Number of sites: 377,672

perc <- percMethylation(meth.sites)
ind.meth <- colMeans(perc,na.rm=T)  # na.rm = logical parameter that tells the 
                                    # function whether or not to remove NA 
                                    # values from the calculation.

hist(ind.meth,breaks=20)  # distribution of average total percent methylation 
                          # per individual.

beeswarm(ind.meth~meth.sites@treatment)  # The bee swarm plot is a one-dimensional
                                         # scatter plot like "stripchart", but with 
                                         # closely-packed, non-overlapping points.
                                         # Treatment 0: Day 0 29C
                                         # Treatment 1: Day 0 31C                                
                                         # Treatment 2: Day 11 29C
                                         # Treatment 3: Day 11 31C 

anova(lm(ind.meth~as.factor(meth.sites@treatment))) # No affect of treatment on global methylation

###Some exploratory analyses
clusterSamples(meth.sites,dist="correlation",method="ward",plot=T)
PCASamples(meth.sites,screeplot=T)
PCASamples(meth.sites)

###Calculate differential methylation
#First, see if there are differences in heat stress
heat.meth <- reorganize(meth.sites,
                        sample.ids=meta$ID,
                        treatment=as.numeric(as.factor(meta$Heat)))
heat.diff <- calculateDiffMeth(meth.sites,adjust="fdr") # lenient q
heat.bases <- getMethylDiff(heat.diff,qvalue=0.05) # just 2 positions - I'll lump H/C for now

## Identical to see if two vectors are identical.
## Reorganize to subset
## pull the results dataframe out of the dm to do other things

#Now subset Day 11
accl.meth <- reorganize(meth.sites,
                        sample.ids=meta$ID[meta$Day==11],
                        treatment=as.numeric(meta$Acclimation[meta$Day==11]))

accl.diff <- calculateDiffMeth(accl.meth,adjust="fdr")
accl.bases <- getMethylDiff(accl.diff,qvalue=0.01,meth.cutoff=25)
diffMethPerChr(accl.diff,plot=T,qvalue.cutoff=0.01)

#Now subset Day 0 (this is a kind of control)
cont.meth <- reorganize(meth.sites,
                        sample.ids=meta$ID[meta$Day==0], #select all samples from day 0
                        treatment=as.numeric(meta$Acclimation[meta$Day==0]))
cont.diff <- calculateDiffMeth(cont.meth,adjust="fdr")
cont.bases <- getMethylDiff(cont.diff,qvalue=0.01,meth.cutoff=25)
diffMethPerChr(cont.diff,plot=T,qvalue.cutoff=0.01) # Some significant positions - individual variation? False positives?

#Acclimation only Day 0 - Day 11
time.meth <- reorganize(meth.sites,
                        sample.ids=meta$ID[meta$Acclimation=="Accl"],
                        treatment=as.numeric(as.factor(meta$Day[meta$Acclimation=="Accl"])))
time.diff <- calculateDiffMeth(time.meth,adjust="fdr")
time.bases <- getMethylDiff(time.diff,qvalue=0.01,meth.cutoff=25)
diffMethPerChr(time.diff,plot=T,qvalue.cutoff=0.01)

# Get Annotations of Differentially methylated bases
coding.obj=readTranscriptFeatures("Amil.coding.bed")
diff.ann <- annotateWithGeneParts(as(time.bases,"GRanges"),
                                  coding.obj)
plotTargetAnnotation(diff.ann,precedence=TRUE,
                     main="differential methylation annotation")
genes.with.dm <- getAssociationWithTSS(diff.ann)
class(genes.with.dm)
genes.with.dm$feature.name

#exons <- regionCounts(filtered.myobj,coding.obj$exons)

### What regions are differentially methylated between 11 days at 29 
### and 11 days at 32??

# Subset Acclimation treatment 29
cont2.meth <- reorganize(meth.sites,
                        sample.ids=meta$ID[meta$Acclimation=="Control"], # Select all samples held at 29
                        treatment=as.numeric(as.factor(meta$Day[meta$Acclimation=="Control"]))) #treatment is day 0 vs day 11 at 29. Expect to see no or few sites.
cont2.diff <- calculateDiffMeth(cont2.meth,adjust="fdr")
cont2.bases <- getMethylDiff(cont2.diff,qvalue=0.01,meth.cutoff=25)
diffMethPerChr(cont2.diff,plot=T,qvalue.cutoff=0.01) # Some significant positions - individual variation? False positives?

# Day 11 only Control - Accl
accl2.meth <- reorganize(meth.sites,
                        sample.ids=meta$ID[meta$Day==11],
                        treatment=as.numeric(as.factor(meta$Acclimation[meta$Day=="11"])))
accl2.diff <- calculateDiffMeth(accl2.meth,adjust="fdr")
accl2.bases <- getMethylDiff(accl2.diff,qvalue=0.1,meth.cutoff=25)
diffMethPerChr(accl2.diff,plot=T,qvalue.cutoff=0.01)
