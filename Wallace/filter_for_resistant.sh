#!/bin/bash

VAR=$(awk ' BEGIN {ORS=""} /^>/ {print (NR==1 ? "" : "\n") $0 " "; next } 
    {print $0}' wallace_genes.fasta)

awk 'BEGIN {FS="[_\t ]"; ORS="\n"} NR==FNR{a[$1]; next} $2 in a' \
    wallace_exclusive_resistant_centroids.txt <(echo "$VAR") | \
    awk '{sub(" ", "\n"); print}' > wallace_genes_exclusive_resistant.fasta