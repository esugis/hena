#!/bin/bash

## Give permission to execute script
chmod +x  scripts/sortlines_linux_amd64

# PPI
echo Sorting interacting PPI pairs from integrated dataset
mv results/integrated_int/integrated_int_ppi.txt  results/integrated_int/integrated_int_ppi_12.txt

# Sort alphabetucally ppi from integrated ds
./scripts/sortlines_linux_amd64 -i results/integrated_int/integrated_int_ppi_12.txt  > results/integrated_int/integrated_int_ppi_sorted.txt
rm -f results/integrated_int/integrated_int_ppi_12.txt

# Remove duplicates in ppi from integrated ds
awk '!a[$1$2]++' results/integrated_int/integrated_int_ppi_sorted.txt  > results/integrated_int/integrated_int_ppi_noduplicated.txt

rm -f results/integrated_int/integrated_int_ppi_sorted.txt

# Count entries in co-expression file
 wc -l results/integrated_int/integrated_int_ppi_noduplicated.txt > results/comparisons/integrated_int_ppi_noduplicated_count.txt

# Co-expression
echo Sorting interacting co-expressed pairs from integrated dataset

mv results/integrated_int/integrated_int_coexp.txt  results/integrated_int/integrated_int_coexp_12.txt

# Sort alphabetucally co-expression from integrated ds
./scripts/sortlines_linux_amd64 -i results/integrated_int/integrated_int_coexp_12.txt  > results/integrated_int/integrated_int_coexp_sorted.txt

rm -f results/integrated_int/integrated_int_coexp_12.txt

# Remove duplicates in co-expression from integrated ds
awk '!a[$1$2]++' results/integrated_int/integrated_int_coexp_sorted.txt  > results/integrated_int/integrated_int_coexp_noduplicated.txt

rm -f results/integrated_int/integrated_int_coexp_sorted.txt

# Count entries in co-expression file
wc -l results/integrated_int/integrated_int_coexp_noduplicated.txt > results/comparisons/integrated_int_coexp_noduplicated_count.txt


# Epistasis
echo Sorting interacting epistatic pairs from integrated dataset

mv results/integrated_int/integrated_int_epi.txt  results/integrated_int/integrated_int_epi_12.txt

# Sort alphabetucally epistatic pairs from integrated ds
./scripts/sortlines_linux_amd64 -i results/integrated_int/integrated_int_epi_12.txt  > results/integrated_int/integrated_int_epi_sorted.txt

rm -f results/integrated_int/integrated_int_epi_12.txt

# Remove duplicates in epistatic pairs from integrated ds
awk '!a[$1$2]++' results/integrated_int/integrated_int_epi_sorted.txt  > results/integrated_int/integrated_int_epi_noduplicated.txt

rm -f results/integrated_int/integrated_int_epi_sorted.txt

# Count entries in epistasis file
wc -l results/integrated_int/integrated_int_epi_noduplicated.txt > results/comparisons/integrated_int_epi_noduplicated_count.txt





