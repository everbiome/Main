#Required R packages
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
#Load fastq files
read_path <- "/xxx/" # CHANGE to the directory containing the fastq files after unzipping.
list.files(read_path)
#Filter and Trim
# Sort ensures forward/reverse reads are in same order, change the fastq patterns if needed
fnFs <- sort(list.files(read_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(read_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1]
#Quality check
pdf("/Volumes/SDB1/Dropbox/Dropbox/Bee2/Quality_Forward.pdf",5,5)
plotQualityProfile(fnFs[1])
dev.off()
pdf("/Volumes/SDB1/Dropbox/Dropbox/Bee2/Quality_Reverse.pdf",5,5)
plotQualityProfile(fnRs[1])
dev.off()
#Filter and trim based on read length and quality
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
