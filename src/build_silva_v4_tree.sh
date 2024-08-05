#!/bin/bash

echo "alignment mafft"
# build reference rooted phylogenetic tree for usage with alpha diversity -
# requires having run src/get_silva_reads_n_classifier.sh before
qiime alignment mafft \
  --i-sequences "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza" \
  --p-n-threads 8 \
  --p-parttree \
  --verbose \
  --o-alignment "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned.qza"

echo "Ben's alignment masking"
python ../src/bens_masking.py --seqs "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned.qza"

echo "phylogeny fasttree"
qiime phylogeny fasttree \
  --i-alignment "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned-ben-masked.qza" \
  --p-n-threads $2 \
  --o-tree "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-unrooted-tree-ben.qza" \
  --verbose

echo "phylogeny midpoint-root"
qiime phylogeny midpoint-root \
  --i-tree "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-unrooted-tree-ben.qza" \
  --o-rooted-tree "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-rooted-tree-ben.qza" \
  --verbose

echo "Remove taxonomic information from leaves"
python ../src/process_rooted_tree.py --tree "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-rooted-tree-ben.qza"
