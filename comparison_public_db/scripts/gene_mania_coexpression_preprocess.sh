#!/bin/bash

awk '
FNR==1 && NR!=1 { while (/^Gene_A/) getline; }
1 {print}
' data_db/genemania/coexp/*.txt > results/genemania/genemania_coexp_all.txt

# remove third column
# sed -i -r 's/(\s+)?\S+//3' genemania_coexp_all.txt > genemania_coexp_12.txt
awk '$3="";1' results/genemania/genemania_coexp_all.txt > results/genemania/genemania_coexp_12.txt

rm -f results/genemania/genemania_coexp_all.txt
# give permission to execute script
#cd scripts/
chmod +x  scripts/sortlines_linux_amd64

# sort alphabetucally
./scripts/sortlines_linux_amd64 -i results/genemania/genemania_coexp_12.txt > results/genemania/genemania_coexp_sorted.txt

rm -f results/genemania/genemania_coexp_12.txt

# remove duplicates
awk '!a[$1$2]++' results/genemania/genemania_coexp_sorted.txt > results/genemania/genemania_coexp_sorted_noduplicated.txt

rm -f results/genemania/genemania_coexp_sorted.txt

# Count entries in genemania preprocessed
wc -l results/genemania/genemania_coexp_sorted_noduplicated.txt > results/comparisons/gm_coexp_count.txt
