#!/bin/bash

# Merge all genemania genetic interaction files 
awk '
    FNR==1 && NR!=1 { while (/^Gene_A/) getline; }
    1 {print}
' data_db/genemania/genetic_int/*.txt > results/genemania/genemania_genint_all.txt

# Remove third column from genemania_genint_all.txt
awk '{print $1,$2}' results/genemania/genemania_genint_all.txt > results/genemania/genemania_genint_12.txt

# Remove tmp file
rm -f results/genemania/genemania_genint_all.txt

## Give permission to execute script
chmod +x  scripts/sortlines_linux_amd64

# Sort alphabetucally genemania genetic interactions
./scripts/sortlines_linux_amd64 -i results/genemania/genemania_genint_12.txt > results/genemania/genemania_genint_sorted.txt
rm -f results/genemania/genemania_genint_12.txt

# Remove duplicates in genemania genetic interactions
awk '!a[$1$2]++' results/genemania/genemania_genint_sorted.txt > results/genemania/genemania_genint_sorted_noduplicated.txt
rm -f results/genemania/genemania_genint_sorted.txt

# Count entries in genemania preprocessed
wc -l results/genemania/genemania_genint_sorted_noduplicated.txt > results/comparisons/gm_genint_count.txt

