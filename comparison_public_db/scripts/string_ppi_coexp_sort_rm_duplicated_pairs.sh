#!/bin/bash

## Give permission to execute script
chmod +x  scripts/sortlines_linux_amd64

echo Sorting interacting PPI pairs from STRING

# Sort alphabetucally string ppi
scripts/sortlines_linux_amd64 -i results/string_db/string_ppi_filt.txt > results/string_db/string_ppi_filt_sorted.txt

# Remove duplicates in string ppi
awk '!a[$1$2]++' results/string_db/string_ppi_filt_sorted.txt > results/string_db/string_ppi_filt_noduplicated.txt
rm -f results/string_db/string_ppi_filt_sorted.txt

# Count entries in string preprocessed
wc -l results/string_db/string_ppi_filt_noduplicated.txt > results/comparisons/string_ppi_count.txt


echo Sorting co-expressed pairs from STRING

# Sort alphabetucally string coexpression
scripts/sortlines_linux_amd64 -i results/string_db/string_coexp_filt.txt > results/string_db/string_coexp_filt_sorted.txt

# Remove duplicates in string coexpression
awk '!a[$1$2]++' results/string_db/string_coexp_filt_sorted.txt > results/string_db/string_coexp_filt_noduplicated.txt

rm -f results/string_db/string_coexp_filt_sorted.txt

# Count entries in string preprocessed coexpression
wc -l results/string_db/string_coexp_filt_noduplicated.txt > results/comparisons/string_coexp_count.txt

