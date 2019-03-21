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

miseq_path <- "/Volumes/SDB1/Dropbox/Dropbox/Bee2/" # CHANGE to the directory containing the fastq files after unzipping.
list.files(miseq_path)
library("openxlsx")
merge=read.xlsx("/Volumes/SDB1/Dropbox/Dropbox/Bee2/100ngdna-送样准备-细菌测序(1).xlsx", sheet=4)
Barcode1=substring(merge[50:54,3],1,9)
fileR=cbind(paste(Barcode2,".fq",sep=""),paste("F",1:5,".fq",sep=""))
Barcode2=substring(merge[62:73,3],1,9)
write.table(Barcode2,"/Volumes/SDB1/Dropbox/Dropbox/Bee2/Rcode.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(fileF,"/Volumes/SDB1/Dropbox/Dropbox/Bee2/Fname.txt",row.names=F,col.names=F,sep="\t",quote=F)
#Filter and Trim
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1]
pdf("/Volumes/SDB1/Dropbox/Dropbox/Bee2/Quality_Forward.pdf",5,5)
plotQualityProfile(fnFs[1])
dev.off()
pdf("/Volumes/SDB1/Dropbox/Dropbox/Bee2/Quality_Reverse.pdf",5,5)
plotQualityProfile(fnRs[1])
dev.off()

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
