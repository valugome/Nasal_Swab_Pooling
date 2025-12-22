#!/bin/bash
#SBATCH -J dada2_pools -t 48:00:00 -n 1 --ntasks=12 --mem=48G -o dada2-pools.out
module load QIIME2/2023.2

## DADA2
qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end-pools.qza --p-n-threads 0 --o-table table-pools.qza \
--p-trim-left-f 18 --p-trunc-len-f 248 \
--p-trim-left-r 22 --p-trunc-len-r 241 \
--o-representative-sequences rep-seqs-pools.qza --output-dir DADA2/Pools

## export DADA2 stats
qiime tools export --input-path DADA2/Pools/denoising_stats.qza --output-path DADA2/Pools/exported_stats_pools




 

