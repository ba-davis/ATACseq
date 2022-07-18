#!/usr/bin/Rscript

# use DiffBind on the macs2 peaks
# called peaks on each individual file
# pvalue 0.05 cutoff

library(DiffBind)

# inpath to bam files
bamPath <- "/home/groups/hoolock/u1/bd/Projects/ECP38/bowtie2_hg19_max50/tmp/5_noMT_mapq30_proper_dedup_bams"
# inpath to peak files
peakpath <- "/home/groups/hoolock/u1/bd/Projects/ECP38/macs2/hg19/individual_callPeaks/p05/output/narrowPeak/canonical"

# make the "sample sheet" data frame
# note: diffBind dba function has set colnames for a sample sheet
#  in this case:
#     Factor = Time
#     Replicate = Pair
#     Tissue = group
samples <- data.frame(SampleID=c("C_I_0","C_II_0","C_III_0","P_I_0","P_II_0","P_III_0",
				 "C_I_3","C_II_3","C_III_3","P_I_3","P_II_3","P_III_3",
				 "C_I_12","C_II_12","C_III_12","P_I_12","P_II_12","P_III_12"),
		      Condition=c(rep("control", 3), rep("PCOS", 3), rep("control", 3), rep("PCOS", 3), rep("control", 3), rep("PCOS", 3)),
		      Factor=c(rep("0", 6), rep("3", 6), rep("12", 6)),
		      Replicate=c(rep(c("p1","p2","p3"), 6)),
		      Tissue=c(rep("C_0", 3), rep("P_0", 3), rep("C_3", 3), rep("P_3", 3), rep("C_12", 3), rep("P_12", 3)),
		      bamReads=c(paste0(bamPath, "/C_I_0.x2000.noMT.mapq30.proper.dedup.bam"),
		      		 paste0(bamPath, "/C_II_0.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/C_III_0.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_I_0.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_II_0.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_III_0.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/C_I_3.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/C_II_3.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/C_III_3.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_I_3.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_II_3.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_III_3.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/C_I_12.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/C_II_12.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/C_III_12.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_I_12.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_II_12.x2000.noMT.mapq30.proper.dedup.bam"),
				 paste0(bamPath, "/P_III_12.x2000.noMT.mapq30.proper.dedup.bam")),
		      Peaks=c(paste0(peakpath, "/C_I_0.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
		      	      paste0(peakpath, "/C_II_0.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/C_III_0.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_I_0.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_II_0.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_III_0.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/C_I_3.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/C_II_3.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/C_III_3.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_I_3.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_II_3.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_III_3.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/C_I_12.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/C_II_12.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/C_III_12.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_I_12.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_II_12.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak"),
			      paste0(peakpath, "/P_III_12.x2000.noMT.mapq30.proper.dedup.bam.p05_peaks.canonical.narrowPeak")),
		      PeakCaller=c(rep("narrow", 18))
)


###########################
# 1) Read in the Peaksets #
###########################

mydat <- dba(sampleSheet=samples, config=data.frame(RunParallel=TRUE, reportInit="DBA", DataType=DBA_DATA_GRANGES, AnalysisMethod=DBA_DESEQ2, minQCth=30, th=0.05, bUsePval="FALSE", bRemoveRandom=TRUE, bCorPlot=FALSE))
mydat # 18 Samples, 1,081,604 sites in matrix (1,827,853 total)

# using peak caller data, produce correlation heatmap using cross-correlations of each row in binding matrix
png("corplot.png")
plot(mydat)
dev.off()


###########################################
# 2) Counting Reads in Consensus Peak Set #
###########################################

# The next step is to calculate a binding matrix with scores based on read counts for every
#  sample (affinity scores), rather than confidence scores for only those peaks called in a specific
#  sample (occupancy scores).
# also, choose to re-center the peaks at the summit and extend 150bp on either side so all peaks will be the same width
# by default, minOverlap parameter is set to 2, meaning only include peaks in at least 2 peaksets
mydat2 <- dba.count(mydat, summits=150)
mydat2

# we can make a new correlation heatmap based on the affinity scores
png("corplot2.png")
plot(mydat2)
dev.off()

# make PCA of the consensus peak counts
png("pca.mydat2.png")
dba.plotPCA(mydat2, attributes=DBA_TISSUE, label=DBA_ID)
dev.off()


###########################
# 3) Establish a Contrast #
###########################

# PCOS_0 vs Control_0 (BLOCK for Pair Match aka "Replicate")
mydat3.comp0 <- dba.contrast(mydat2, group1=mydat2$masks$P_0, group2=mydat2$masks$C_0, name1="PCOS_0", name2="Control_0", block=DBA_REPLICATE)

# PCOS_3 vs Control_3 (BLOCK for Pair Match aka "Replicate")
mydat3.comp3 <- dba.contrast(mydat2, group1=mydat2$masks$P_3, group2=mydat2$masks$C_3, name1="PCOS_3", name2="Control_3", block=DBA_REPLICATE)

# PCOS_12 vs Control_12 (BLOCK for Pair Match aka "Replicate")
mydat3.comp12 <- dba.contrast(mydat2, group1=mydat2$masks$P_12, group2=mydat2$masks$C_12, name1="PCOS_12", name2="Control_12", block=DBA_REPLICATE)


####################################
# 4) Perform Differential Analysis #
####################################

mydat4.0 <- dba.analyze(mydat3.comp0)
mydat4.3 <- dba.analyze(mydat3.comp3)
mydat4.12 <- dba.analyze(mydat3.comp12)

mydat4.0
mydat4.3
mydat4.12

# Can make new correlation heatmap based on the diff bound sites
png("corrplot.diff.PCOS.0_vs_Control.0.png")
plot(mydat4.0, contrast=1, method=DBA_DESEQ2_BLOCK)
dev.off()

png("corrplot.diff.PCOS.3_vs_Control.3.png")
plot(mydat4.3, contrast=1, method=DBA_DESEQ2_BLOCK)
dev.off()

png("corrplot.diff.PCOS.12_vs_Control.12.png")
plot(mydat4.12, contrast=1, method=DBA_DESEQ2_BLOCK)
dev.off()

# Can make PCA plots based on the diff bound sites
png("PCOS.0_vs_Control.0.diff.pca.png")
dba.plotPCA(mydat4.0, method=DBA_DESEQ2_BLOCK, contrast=1, attributes=DBA_TISSUE, label=DBA_ID)
dev.off()

png("PCOS.3_vs_Control.3.diff.pca.png")
dba.plotPCA(mydat4.3, method=DBA_DESEQ2_BLOCK, contrast=1, attributes=DBA_TISSUE, label=DBA_ID)
dev.off()

png("PCOS.12_vs_Control.12.diff.pca.png")
dba.plotPCA(mydat4.12, method=DBA_DESEQ2_BLOCK, contrast=1, attributes=DBA_TISSUE, label=DBA_ID)
dev.off()

#########################################################
# 5) PValue Distributions of Differentially Bound Sites #
#########################################################

# obtain stats for all binding sites in order to produce pvalue distribution histograms
mydat.DB.0 <- dba.report(mydat4.0, th=1, method=DBA_DESEQ2_BLOCK, bCounts=TRUE)
mydat.DB.3 <- dba.report(mydat4.3, th=1, method=DBA_DESEQ2_BLOCK, bCounts=TRUE)
mydat.DB.12 <- dba.report(mydat4.12, th=1, method=DBA_DESEQ2_BLOCK, bCounts=TRUE)

mydat.DB.0
mydat.DB.3
mydat.DB.12

# make pval distributions
png("PCOS.0_vs_Control.0.pval.hist.png")
hist(as.data.frame(mydat.DB.0)$p.value, breaks=40, col="grey")
dev.off()

png("PCOS.3_vs_Control.3.pval.hist.png")
hist(as.data.frame(mydat.DB.3)$p.value, breaks=40, col="grey")
dev.off()

png("PCOS.12_vs_Control.12.pval.hist.png")
hist(as.data.frame(mydat.DB.12)$p.value, breaks=40, col="grey")
dev.off()


#########################
# 6) Export ALL results #
#########################

res1 <- as.data.frame(mydat.DB.0)
res1$PeakID <- paste("Peak", rownames(res1), sep="_")
#foo <- foo[ ,c(12,1,2,3,4,5,6,7,8,9,10,11)]

write.table(res1, "PCOS.0_vs_Control.0.ALL.results.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

res2 <- as.data.frame(mydat.DB.3)
res2$PeakID <- paste("Peak", rownames(res2), sep="_")
#foo <- foo[ ,c(12,1,2,3,4,5,6,7,8,9,10,11)]

write.table(res2, "PCOS.3_vs_Control.3.ALL.results.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

res3 <- as.data.frame(mydat.DB.12)
res3$PeakID <- paste("Peak", rownames(res3), sep="_")
#foo <- foo[ ,c(12,1,2,3,4,5,6,7,8,9,10,11)]

write.table(res3, "PCOS.12_vs_Control.12.ALL.results.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


####################################
# 7) Export SIG Results FDR < 0.05 #
####################################

# Obtain and export sig regions (FDR < 0.05, both gains and losses)
mydat.DB.0.sig <- res1[res1$FDR < 0.05, ]
write.table(mydat.DB.0.sig, "PCOS.0_vs_Control.0.SIG.results.both.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mydat.DB.3.sig <- res2[res2$FDR < 0.05, ]
write.table(mydat.DB.3.sig, "PCOS.3_vs_Control.3.SIG.results.both.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mydat.DB.12.sig <- res3[res3$FDR < 0.05, ]
write.table(mydat.DB.12.sig, "PCOS.12_vs_Control.12.SIG.results.both.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


# Obtain and export sig regions (FDR < 0.05) gains only (Fold > 0)
mydat.DB.0.sig.gains <- mydat.DB.0.sig[mydat.DB.0.sig$Fold > 0, ]
write.table(mydat.DB.0.sig.gains, "PCOS.0_vs_Control.0.SIG.results.gains.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mydat.DB.3.sig.gains <- mydat.DB.3.sig[mydat.DB.3.sig$Fold > 0, ]
write.table(mydat.DB.3.sig.gains, "PCOS.3_vs_Control.3.SIG.results.gains.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mydat.DB.12.sig.gains <- mydat.DB.12.sig[mydat.DB.12.sig$Fold > 0, ]
write.table(mydat.DB.12.sig.gains, "PCOS.12_vs_Control.12.SIG.results.gains.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


# Obtain and export sig regions (FDR < 0.05) losses only (Fold < 0)
mydat.DB.0.sig.losses <- mydat.DB.0.sig[mydat.DB.0.sig$Fold < 0, ]
write.table(mydat.DB.0.sig.losses, "PCOS.0_vs_Control.0.SIG.results.losses.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mydat.DB.3.sig.losses <- mydat.DB.3.sig[mydat.DB.3.sig$Fold < 0, ]
write.table(mydat.DB.3.sig.losses, "PCOS.3_vs_Control.3.SIG.results.losses.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mydat.DB.12.sig.losses <- mydat.DB.12.sig[mydat.DB.12.sig$Fold < 0, ]
write.table(mydat.DB.12.sig.losses, "PCOS.12_vs_Control.12.SIG.results.losses.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


###########
# Boxplot #
###########

png("PCOS.0_vs_Control.0.boxplot.png")
dba.plotBox(mydat4.0, method=DBA_DESEQ2_BLOCK)
dev.off()

png("PCOS.3_vs_Control.3.boxplot.png")
dba.plotBox(mydat4.3, method=DBA_DESEQ2_BLOCK)
dev.off()

png("PCOS.12_vs_Control.12.boxplot.png")
dba.plotBox(mydat4.12, method=DBA_DESEQ2_BLOCK)
dev.off()


###########
# MA Plot #
###########

png("MA.PCOS.0_vs_Control.0.png")
dba.plotMA(mydat4.0, th=0.05, bUsePval=FALSE, method=DBA_DESEQ2_BLOCK, contrast=1)
dev.off()

png("MA.PCOS.3_vs_Control.3.png")
dba.plotMA(mydat4.3, th=0.05, bUsePval=FALSE, method=DBA_DESEQ2_BLOCK, contrast=1)
dev.off()

png("MA.PCOS.12_vs_Control.12.png")
dba.plotMA(mydat4.12, th=0.05, bUsePval=FALSE, method=DBA_DESEQ2_BLOCK, contrast=1)
dev.off()

sessionInfo()
