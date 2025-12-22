module load QIIME2/2023.2

# make exported dir
mkdir exported

## export to phyloseq dir
qiime tools export --input-path rep-seqs.qza --output-path exported/
qiime tools export --input-path taxonomy.qza --output-path exported/
qiime tools export --input-path rooted-tree.qza --output-path exported/
qiime tools export --input-path table.qza --output-path exported/
