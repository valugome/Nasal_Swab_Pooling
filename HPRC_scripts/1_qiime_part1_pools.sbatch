#!/bin/bash
#SBATCH -J qiime_import_pools -t 02:00:00 -n 1 --ntasks=1 -o qiime1_pools.out --mem=10G
module load QIIME2/2023.2

## importing reads into qiime2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ../metadata/manifest_pools.tsv \
--input-format PairedEndFastqManifestPhred33V2 --output-path demux-paired-end-pools.qza

## creating visualization to pick trim and trunc lengths
## visualize demuxed reads
qiime demux summarize --i-data demux-paired-end-pools.qza --o-visualization demux-paired-end-pools.qzv 

