#!/bin/bash
#SBATCH -J qiime_import_individuals -t 04:00:00 -n 1 --ntasks=1 -o qiime1_individuals.out --mem=10G
module load QIIME2/2023.2

## importing reads into qiime2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ../metadata/manifest_individuals.tsv \
--input-format PairedEndFastqManifestPhred33V2 --output-path demux-paired-end-individuals.qza

## creating visualization to pick trim and trunc lengths
## visualize demuxed reads
qiime demux summarize --i-data demux-paired-end-individuals.qza --o-visualization demux-paired-end-individuals.qzv 

