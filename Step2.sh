#!bin/bash

#!bin/bash
#Folders: filtered merged split qiime2
#Extract barcodes for pair-end read1 and read2
echo 'Merge filtered R1 and R2 pair-end reads'
./bbmap/bbmerge.sh in1=filtered/R1.fq.gz in2=filtered/R2.fq.gz out=merged/merged.fq.gz
echo 'Extract barcodes from merged pair-end reads'
umi_tools extract --stdin=merged/merged.fq.gz --bc-pattern=NNNNNNNNN  --stdout=split/split1.fq.gz
./bbmap/bbduk.sh in=split/split1.fq.gz out=split/split12.fq.gz ftl=18
./bbmap/reformat.sh n=split/split12.fq.gz out=split/split123.fq.gz trd
./bbmap/reformat.sh in=split/split123.fq.gz out=split/split1234.fq.gz rcomp
umi_tools extract --stdin=split/split1234.fq.gz --bc-pattern=NNNNNNNNN --stdout split/split12345.fq.gz
./bbmap/bbduk.sh in=split/split12345.fq.gz out=split/R1R2.fq.gz ftl=18
rm split/s*.*
./bbmap/reformat.sh in=merged/merged.fq.gz out=merged/mergedr.fq.gz rcomp
umi_tools extract --stdin=merged/mergedr.fq.gz --bc-pattern=NNNNNNNNN  --stdout=split/split1.fq.gz
./bbmap/bbduk.sh in=split/split1.fq.gz out=split/split12.fq.gz ftl=18
./bbmap/reformat.sh n=split/split12.fq.gz out=split/split123.fq.gz trd
./bbmap/reformat.sh in=split/split123.fq.gz out=split/split1234.fq.gz rcomp
umi_tools extract --stdin=split/split1234.fq.gz --bc-pattern=NNNNNNNNN --stdout split/split12345.fq.gz
./bbmap/bbduk.sh in=split/split12345.fq.gz out=split/R2R1.fq.gz ftl=18
rm split/s*.*
cat split/R2R1.fq.gz split/R1R2.fq.gz > split.fq.gz
rm split/R*.*
#Demultiplex reads
./bbmap/demuxbyname.sh in=split.fq.gz out=split/%.fq prefixmode=f names=barcode.txt
rm split.fq.gz
#load qiime2
echo 'Start QIIME2'
Conda activate qiime2
echo 'Import demux filtered merged reads'
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path fq33.csv --output-path qiime2/demux.qza --input-format SingleEndFastqManifestPhred33
cd qiime2
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
echo 'Dada2 denoise reads'
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 390\
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv
echo 'Summarize feature-table'
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file ../metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
echo 'phylogeny tree'
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
echo 'Core diversity at depth at 1500'
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1500 \
  --m-metadata-file ../metadata.tsv \
  --output-dir core-metrics-results
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ../metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file ../metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
echo 'Taxonomy clssification by Greengene'
qiime feature-classifier classify-sklearn \
  --i-classifier ../../gg-13-8-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file ../metadata.tsv  \
  --o-visualization taxa-bar-plots.qzv
