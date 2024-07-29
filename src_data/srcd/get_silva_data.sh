#!/bin/bash


# following this tutorial we are fetching SILVA reference reads + tax. classifiers
# https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences "$1/silva-138.1-ssu-nr99-rna-seqs.qza" \
    --o-silva-taxonomy "$1/silva-138.1-ssu-nr99-tax.qza"

qiime rescript reverse-transcribe \
    --i-rna-sequences "$1/silva-138.1-ssu-nr99-rna-seqs.qza" \
    --o-dna-sequences "$1/silva-138.1-ssu-nr99-seqs.qza"

qiime rescript cull-seqs \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs.qza" \
    --o-clean-sequences "$1/silva-138.1-ssu-nr99-seqs-cleaned.qza"

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-cleaned.qza" \
    --i-taxonomy "$1/silva-138.1-ssu-nr99-tax.qza" \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs "$1/silva-138.1-ssu-nr99-seqs-filt.qza" \
    --o-discarded-seqs "$1/silva-138.1-ssu-nr99-seqs-discard.qza"

qiime rescript dereplicate \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-filt.qza" \
    --i-taxa "$1/silva-138.1-ssu-nr99-tax.qza" \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "$1/silva-138.1-ssu-nr99-seqs-derep-uniq.qza" \
    --o-dereplicated-taxa "$1/silva-138.1-ssu-nr99-tax-derep-uniq.qza"

qiime feature-classifier extract-reads \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-derep-uniq.qza" \
    --p-f-primer GTGYCAGCMGCCGCGGTAA \
    --p-r-primer GGACTACNVGGGTWTCTAAT \
    --p-n-jobs $2 \
    --p-read-orientation 'forward' \
    --o-reads "$1/silva-138.1-ssu-nr99-seqs-515f-806r.qza"

qiime rescript dereplicate \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-515f-806r.qza" \
    --i-taxa "$1/silva-138.1-ssu-nr99-tax-derep-uniq.qza" \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza" \
    --o-dereplicated-taxa  "$1/silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza"
# ! are the `o-dereplicated-sequences` the SILVA reference v4 reads that can be used for closed reference clustering?

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza" \
    --i-reference-taxonomy "$1/silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza" \
    --o-classifier "$1/silva-138.1-ssu-nr99-515f-806r-classifier.qza"

# build reference rooted phylogenetic tree for usage with alpha diversity
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

# remove unneeded files
rm -r "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-unrooted-tree.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned-masked.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq-aligned.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-515f-806r.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-derep-uniq.qza" 
rm -r "$1/silva-138.1-ssu-nr99-tax-derep-uniq.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-filt.qza" 
rm -r "$1/silva-138.1-ssu-nr99-seqs-discard.qza"
rm -r "$1/silva-138.1-ssu-nr99-rna-seqs.qza"
rm -r "$1/silva-138.1-ssu-nr99-tax.qza"
rm -r "$1/silva-138.1-ssu-nr99-rna-seqs.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-cleaned.qza"