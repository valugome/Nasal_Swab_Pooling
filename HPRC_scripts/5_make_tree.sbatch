#!/bin/bash
#SBATCH -J tree -t 48:00:00 -n 1 --cpus-per-task=12 --mem=200G -o tree.out --account=132729268120
module load QIIME2/2023.2

# making a mid-point rooted tree (4 steps)
##1.Multiple sequence allignment on our repseqs using MAFFT 
qiime alignment mafft --i-sequences rep-seqs.qza  --p-n-threads 'auto' \
--o-alignment aligned-rep-seqs.qza
##2. Mask regions of the allignment that arent useful (gaps, highly variable regions, etc.)
qiime alignment mask --i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza
## 3. Make the tree with FastTree (it will use multile threads so we use 12 cpus per task)
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --p-n-threads 'auto' \
--o-tree unrooted-tree.qza
## 4. Midpoint rooting the tree (takes the lenght of the longest branch in two terminal nodes (two ASVs). So whatever the longest branch lenght is, calculates the midpoint by finding the equal length between that point and both the terminal nodes)
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza
