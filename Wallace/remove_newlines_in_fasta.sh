#!/bin/bash

awk 'BEGIN {ORS=""} /^>/ { sub("\r", "\n", $0); print (NR==1 ? "" : "\n") $0; next } {sub("\r", "", $0); print}' wallace_genes.fasta > wallace_genes_2.fasta
mv wallace_genes_2.fasta wallace_genes.fasta