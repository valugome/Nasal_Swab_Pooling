#!/bin/bash
#SBATCH -J dada2_individuals -t 10:00:00 -n 1 --ntasks=12 --mem=48G -o dada2-individuals.out
module load QIIME2/2023.2

## DADA2
qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end-individuals.qza --p-n-threads 0 --o-table table-individuals.qza \
--p-trim-left-f 18 --p-trunc-len-f 248 \
--p-trim-left-r 22 --p-trunc-len-r 241 \
--o-representative-sequences rep-seqs-individuals.qza --output-dir DADA2/Individuals

## export DADA2 stats
qiime tools export --input-path DADA2/Individuals/denoising_stats.qza --output-path DADA2/Individuals/exported_stats




 

