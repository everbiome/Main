#!bin/bash

#Extract barcodes for pair-end read1 and read2
umi_tools extract --stdin=Bee2/merged_filtered.fq.gz --bc-pattern=NNNNNNNNN --log=processedF.log --stdout processedF.fastq.gz
#Trim the primers from read1
./bbmap/bbduk.sh in=Bee2/processedF.fastq.gz out=Bee2/TrimmedbarcodedF.fastq.gz ftl=18
./bbmap/reformat.sh in=Bee2/TrimmedbarcodedF.fastq.gz out=Bee2/renamed_trimedF.fastq.gz trd
#Merge Trimed_barcoded_read1 with untrimed read2
./bbmap/reformat.sh in=Bee2/renamed_trimedF.fastq.gz out=Bee2/renamed_trimedF_rev.fastq.gz rcomp
umi_tools extract --stdin=Bee2/renamed_trimedF_rev.fastq.gz --bc-pattern=NNNNNNNNN --log=processedR.log --stdout Bee2/processedFR.fastq.gz
./bbmap/bbduk.sh in=Bee2/processedFR.fastq.gz out=Bee2/TrimmedbarcodedFR.fastq.gz ftl=18

#Demultiplex reads
./bbmap/demuxbyname.sh in=Bee2/TrimmedbarcodedFR.fastq.gz out=Bee2/out_%_#.fq prefixmode=f names=AGCGTGTCC_AGGAACGCG outu=filename

Conda activate qiime2

qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path fq33.csv --output-path qiime2/single-end-demux.qza --input-format SingleEndFastqManifestPhred33
mv qiime2/single-end-demux.qza qiime2/demux.qza;cd qiime2
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
qiime tools view demux.qzv

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 390\
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza


qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file ../Bee2_metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1500 \
  --m-metadata-file ../Bee2_metadata.tsv \
  --output-dir core-metrics-results


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ../Bee2_metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file ../Bee2_metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier ../../gg-13-8-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file ../Bee2_metadata.tsv  \
  --o-visualization taxa-bar-plots.qzv
