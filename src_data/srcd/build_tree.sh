#!/bin/bash

# build reference rooted phylogenetic tree for usage with alpha diversity -
# requires having run get_silva_tree.sh before
qiime alignment mafft \
  --i-sequences "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza" \
  --p-n-threads $2 \
  --o-alignment "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned.qza"
qiime alignment mask \
  --i-alignment "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned.qza" \
  --o-masked-alignment "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned-masked.qza"
qiime phylogeny fasttree \
  --i-alignment "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned-masked.qza" \
  --o-tree "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-unrooted-tree.qza"
qiime phylogeny midpoint-root \
  --i-tree "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-unrooted-tree.qza" \
  --o-rooted-tree "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-rooted-tree.qza"
