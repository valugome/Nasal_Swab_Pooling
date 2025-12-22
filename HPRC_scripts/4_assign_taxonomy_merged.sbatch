#!/bin/bash
#SBATCH -J classify_merged -t 48:00:00 -n 1 --ntasks=12 --mem=200G -o classify.out
module load QIIME2/2023.2

## classifying
qiime feature-classifier classify-sklearn --i-classifier /scratch/group/vero_research/valerialugo/PooledSwabs/QIIME2/QIIME_Classifiers_Lee/silva-138.1-ssu-nr99-341f-785r-classifier.qza \
--i-reads merged-repseqs.qza --o-classification taxonomy.qza --p-n-jobs -1

## copy the taxonomy file into created no-mito-no-chloro directory
cp taxonomy.qza no-mito-no-chloro/

## filtering ASVs
qiime taxa filter-table --i-table merged-table.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast --o-filtered-table no-mito-no-chloro/table.qza

## filtering rep-seqs
qiime taxa filter-seqs --i-sequences merged-repseqs.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast --o-filtered-sequences no-mito-no-chloro/rep-seqs.qza

