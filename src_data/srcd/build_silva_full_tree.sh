#!/bin/bash
# get aligned full length RNA sequences
curl -o "$1/SILVA_138_1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz" \
    https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz

gunzip "$1/SILVA_138_1_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"

qiime tools import \
  --input-path "$1/SILVA_138_1_SSURef_NR99_tax_silva_full_align_trunc.fasta" \
  --output-path "$1/silva-138-1-ssu-nr99-aln-rna.qza" \
  --type 'FeatureData[AlignedRNASequence]'

# reverse transcribe to DNA
qiime rescript reverse-transcribe \
  --i-rna-sequences "$1/silva-138-1-ssu-nr99-aln-rna.qza" \
  --o-dna-sequences "$1/silva-138-1-ssu-nr99-seqs-aligned.qza"

# mask out gaps using Ben Kaehler's approach from:
# https://gist.github.com/BenKaehler/d9291d59bce5cd3d2a90c73b822b3a21 linked
# here:
# https://forum.qiime2.org/t/estimating-time-and-memory-for-masking-large-dataset/21166/9
# only keeping columns with less than ~99.56% gaps
echo "Ben's alignment masking"
python srcd/bens_masking.py --seqs "$1/silva-138-1-ssu-nr99-seqs-aligned.qza"

# build a tree
echo "phylogeny fasttree"
qiime phylogeny fasttree \
  --i-alignment "$1/silva-138-1-ssu-nr99-seqs-aligned-ben-masked.qza" \
  --p-n-threads $2 \
  --o-tree "$1/silva-138-1-ssu-nr99-seqs-unrooted-tree-ben.qza" \
  --verbose

echo "phylogeny midpoint-root"
qiime phylogeny midpoint-root \
  --i-tree "$1/silva-138-1-ssu-nr99-seqs-unrooted-tree-ben.qza" \
  --o-rooted-tree "$1/silva-138-1-ssu-nr99-seqs-rooted-tree-ben.qza" \
  --verbose

echo "Remove taxonomic information from leaves"
python srcd/process_rooted_tree.py --tree "$1/silva-138-1-ssu-nr99-seqs-rooted-tree-ben.qza"