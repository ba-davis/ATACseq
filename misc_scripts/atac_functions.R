










# FUNCTION to remove noncanonical contig entries from a narrowPeak file from macs2
# these new canonical narrowPeak files will be read in to diffBind for differential analysis
cleanNarrowPeak <- function(narrowPeak, name, cans=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')) {
  mytab <- read.delim(narrowPeak, header=F)
  
  mytab2 <- mytab[mytab$V1 %in% cans, ]

  write.table(mytab2, paste(name, "canonical.narrowPeak", sep="."), sep="\t", col.names=F, row.names=F, quote=F)
}


# FUNCTION to make a sig BED file track of differential peaks from diffBind sig results df (both gains and losses)
# use "chr" format, compatible with UCSC or IGV
# use canonical contigs only
# input args:
#   res: sig results df
#   name: desired prefix name of outfile
#   cans: list of canonical contigs to keep, default is human
#   intermediate: default FALSE, if TRUE will only export the original df but with canonical contigs only (no other files produced)
#   track: default FALSE, if TRUE will include a track header and 8 col bed file and ucsc contig format
#   ensembl: default FALSE, if TRUE will keep ensembl contig names (1,2,etc)
# so default of this function will produce a 4col bed file with canonical contigs and ucsc format (used for GREAT)
makeDiffPeakBed <- function(res, name, cans=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'), intermediate=FALSE, track=FALSE, ensembl=FALSE) {
  res <- read.delim(res, header=T)
  res2 <- res[res$seqnames %in% cans, ]

  if (intermediate == TRUE) {
    write.table(res2, paste(name, "canonical.txt", sep="."), sep="\t", col.names=T, row.names=F, quote=F)
    stop("Only producing a new canonical contig results df")
  }

  # keep only the chr, start, end, and PeakID columns
  res3 <- res2[ ,c(1,2,3,ncol(res2))]
  res3$start <- res3$start - 1

  if (ensembl == TRUE) {
    write.table(res3, paste(name, "ensembl.bed", sep="."), sep="\t", col.names=F, row.names=F, quote=F)
  }

  # change seqnames to ucsc format
  res3$seqnames <- paste0("chr", res3$seqnames)

  write.table(res3, paste(name, "bed", sep="."), sep="\t", col.names=F, row.names=F, quote=F)

  # add track line to exported file
  if (track == TRUE) {
    mytab <- res2[ ,c(1,2,3,ncol(res2),9)]
    mytab$seqnames <- paste0("chr", mytab$seqnames)
    mytab$start <- mytab$start - 1

    mytab[ ,6] <- "."
    mytab[ ,7] <- mytab$start
    mytab[ ,8] <- mytab$end
    mytab[ ,9] <- ifelse(mytab[ ,5] > 0, '255,0,0', ifelse(mytab[ ,5] < 0, '0,0,255', '0,0,0'))

    # add track line to exported file
    cat(paste0("track type=bed name=", name), file=paste(name, "track.bed", sep="."))
    cat("\n", file=paste(name, "track.bed", sep="."), append=TRUE)

    # export bed info
    write.table(mytab,
                paste(name, "track.bed", sep="."),
                sep="\t",
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE,
                append=TRUE
    )
  }
}


# FUNCTION to use ChIPseeker to annotate differential peaks
#library(GenomicFeatures)
#library(ChIPseeker)
# first, make a TxDb object from a gtf file
annoMyRegions <- function(peaks, biomart, name) {
  # read in differential peaks file
  mydat <- makeGRangesFromDataFrame(read.delim(peaks, header=T), keep.extra.columns=TRUE)

  # annotate with ChIPseeker function
  res <- annotatePeak(mydat,
         tssRegion = c(-3000, 3000),
         TxDb = txdb,
         level = "transcript",
         assignGenomicAnnotation = TRUE,
         overlap="TSS"
  )

  res
  # convert to data frame
  foo <- as.data.frame(res)

  # read in biomart info
  biomartinfo <- read.delim(biomart, header=T)
  # remove repeat rows based on geneID
  biomartinfo <- biomartinfo[!(duplicated(biomartinfo$Gene.stable.ID)), ]
  colnames(biomartinfo)[1] <- "geneId"

  # merge
  mydat <- merge(foo, biomartinfo, by="geneId", all.x=T)
  mydat2 <- mydat[ ,c(2:ncol(foo),1,(ncol(foo)+1):ncol(mydat))]

  write.table(mydat2, paste(name, "annot.txt", sep='.'), sep="\t", col.names=T, row.names=F, quote=F)

  return(res)
}
