# aligned RNA sequence
#!/bin/bash
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

# mask out gaps
qiime alignment mask \
  --i-alignment "$1/silva-138-1-ssu-nr99-seqs-aligned.qza" \
  --o-masked-alignment "$1/silva-138-1-ssu-nr99-seqs-aligned-masked.qza"

# build a tree
qiime phylogeny fasttree \
  --i-alignment "$1/silva-138-1-ssu-nr99-seqs-aligned-masked.qza" \
  --p-n-threads $2 \ 
  --o-tree "$1/silva-138-1-ssu-nr99-seqs-unrooted-tree.qza"

qiime phylogeny midpoint-root \
  --i-tree "$1/silva-138-1-ssu-nr99-seqs-unrooted-tree.qza" \
  --o-rooted-tree "$1/silva-138-1-ssu-nr99-seqs-rooted-tree.qza"