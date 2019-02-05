#!/bin/bash

# Merge all genemania physical interaction files 
awk '
    FNR==1 && NR!=1 { while (/^Gene_A/) getline; }
    1 {print}
' data_db/genemania/ppi/*.txt > results/genemania/genemania_ppi_all.txt

# Remove third column from genemania_ppi_all.txt
awk '{print $1,$2}' results/genemania/genemania_ppi_all.txt > results/genemania/genemania_ppi_12.txt

# Remove tmp file
rm -f results/genemania/genemania_ppi_all.txt

## Give permission to execute script
chmod +x  scripts/sortlines_linux_amd64

# Sort alphabetucally genemania ppi
./scripts/sortlines_linux_amd64 -i results/genemania/genemania_ppi_12.txt > results/genemania/genemania_ppi_sorted.txt
rm -f results/genemania/genemania_ppi_12.txt

# Remove duplicates in genemania ppi
awk '!a[$1$2]++' results/genemania/genemania_ppi_sorted.txt > results/genemania/genemania_ppi_sorted_noduplicated.txt
rm -f results/genemania/genemania_ppi_sorted.txt

# Count entries in genemania preprocessed
 wc -l results/genemania/genemania_ppi_sorted_noduplicated.txt > results/comparisons/gm_ppi_count.txt
