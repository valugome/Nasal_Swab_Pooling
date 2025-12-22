#!/bin/bash
#SBATCH -J merging -t 04:00:00 -n 1 --cpus-per-task=4 --mem=48G -o merge_poolandindividuals.out --account=132729268120
module load QIIME2/2023.2

#merging the ASV count tables
qiime feature-table merge \
--i-tables table-pools.qza table-individuals.qza \
--o-merged-table merged-table.qza

#merging the representative sequences for the ASVs
qiime feature-table merge-seqs \
--i-data rep-seqs-pools.qza rep-seqs-individuals.qza \
--o-merged-data merged-repseqs.qza
