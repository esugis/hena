#!/bin/bash

####### Find overlaps in PPIs
echo Compare PPIs from STRING and integrated dataset
tail -n +2 results/string_db/string_ppi_filt_noduplicated.txt  > results/string_db/string_ppi_filt_noduplicated_tmp.txt
#rm -f results/string_db/string_ppi_filt_noduplicated.txt
mv results/string_db/string_ppi_filt_noduplicated_tmp.txt results/string_db/string_ppi_filt_noduplicated.txt

# Sort both files before comparisons
sort  results/integrated_int/integrated_int_ppi_noduplicated.txt > results/integrated_int/integrated_int_ppi_noduplicated_sorted.txt
sort results/string_db/string_ppi_filt_noduplicated.txt > results/string_db/string_ppi_filt_noduplicated_sorted.txt


comm -12 results/integrated_int/integrated_int_ppi_noduplicated_sorted.txt results/string_db/string_ppi_filt_noduplicated_sorted.txt > results/comparisons/string_intds_ppi_overlap_sort.txt


#sort  results/integrated_int/integrated_int_ppi_noduplicated.txt results/string_db/string_ppi_filt_noduplicated.txt| uniq -d > results/comparisons/string_intds_ppi_overlap_sort.txt

# Count overlaps of the gene pairs in ppi from integrated dataset and string physical interactions.
wc -l results/comparisons/string_intds_ppi_overlap_sort.txt > results/comparisons/string_intds_ppi_overlap_sort_count.txt


echo Compare PPIs from GeneMania and integrated dataset
tail -n +2 results/genemania/genemania_ppi_sorted_noduplicated.txt > results/genemania/genemania_ppi_sorted_noduplicated_tmp.txt
#rm -f results/genemania/genemania_ppi_sorted_noduplicated.txt
mv results/genemania/genemania_ppi_sorted_noduplicated_tmp.txt results/genemania/genemania_ppi_sorted_noduplicated.txt

# Count overlaps of the gene pairs in ppi from integrated dataset and genemania physical interactions.
sort  results/integrated_int/integrated_int_ppi_noduplicated.txt > results/integrated_int/integrated_int_ppi_noduplicated_sorted.txt

sort results/genemania/genemania_ppi_sorted_noduplicated.txt > results/genemania/genemania_ppi_sorted_noduplicated_sorted.txt

comm -12 results/integrated_int/integrated_int_ppi_noduplicated_sorted.txt results/genemania/genemania_ppi_sorted_noduplicated_sorted.txt > results/comparisons/gm_intds_ppi_overlap_sort.txt
#sort  results/integrated_int/integrated_int_ppi_noduplicated.txt results/genemania/genemania_ppi_sorted_noduplicated.txt | uniq -d > results/comparisons/gm_intds_ppi_overlap_sort.txt

wc -l results/comparisons/gm_intds_ppi_overlap_sort.txt > results/comparisons/gm_intds_ppi_overlap_sort_count.txt


####### Find overlaps in co-expression
echo Compare co-expression from STRING and integrated dataset
tail -n +2 results/string_db/string_coexp_filt_noduplicated.txt  > results/string_db/string_coexp_filt_noduplicated_tmp.txt
#rm -f results/string_db/string_coexp_filt_noduplicated.txt
mv results/string_db/string_coexp_filt_noduplicated_tmp.txt results/string_db/string_coexp_filt_noduplicated.txt

# Count overlaps of the gene pairs in co-expression from integrated dataset and string coexpression interactions.
sort  results/integrated_int/integrated_int_coexp_noduplicated.txt results/string_db/string_coexp_filt_noduplicated.txt| uniq -d > results/comparisons/string_intds_coexp_overlap_sort.txt

wc -l results/comparisons/string_intds_coexp_overlap_sort.txt > results/comparisons/string_intds_coexp_overlap_sort_count.txt


echo Compare co-expression from GeneMania and integrated dataset
tail -n +2 results/genemania/genemania_coexp_sorted_noduplicated.txt > results/genemania/genemania_coexp_sorted_noduplicated_tmp.txt
#rm -f results/genemania/genemania_coexp_sorted_noduplicated.txt
mv results/genemania/genemania_coexp_sorted_noduplicated_tmp.txt results/genemania/genemania_coexp_sorted_noduplicated.txt

# Count overlaps of the gene pairs in co-expression from integrated dataset and genemania co-expression interactions.
sort  results/integrated_int/integrated_int_coexp_noduplicated.txt results/genemania/genemania_coexp_sorted_noduplicated.txt | uniq -d > results/comparisons/gm_intds_coexp_overlap_sort.txt

wc -l results/comparisons/gm_intds_coexp_overlap_sort.txt > results/comparisons/gm_intds_coexp_overlap_sort_count.txt

####### Find overlaps in genetic interactions

echo Compare PPIs from GeneMania and integrated dataset
tail -n +2 results/genemania/genemania_genint_sorted_noduplicated.txt > results/genemania/genemania_genint_sorted_noduplicated_tmp.txt
#rm -f results/genemania/genemania_genint_sorted_noduplicated.txt
mv results/genemania/genemania_genint_sorted_noduplicated_tmp.txt results/genemania/genemania_genint_sorted_noduplicated.txt

# Count overlaps of the gene pairs in ppi from integrated dataset and genemania genetic interactions.
sort  results/integrated_int/integrated_int_epi_noduplicated.txt results/genemania/genemania_genint_sorted_noduplicated.txt | uniq -d > results/comparisons/gm_intds_genint_overlap_sort.txt

wc -l results/comparisons/gm_intds_genint_overlap_sort.txt > results/comparisons/gm_intds_genint_overlap_sort_count.txt



####### Find overlaps in PPI STRING GeneMania
# Count overlaps of the gene pairs in co-expression from STRING and GeneMania physical interactions.
sort  results/string_db/string_coexp_filt_noduplicated.txt results/genemania/genemania_coexp_sorted_noduplicated.txt | uniq -d > results/comparisons/gm_string_coexp_overlap_sort.txt
wc -l results/comparisons/gm_string_coexp_overlap_sort.txt > results/comparisons/gm_string_coexp_overlap_sort_count.txt


####### Find overlaps in co-expression STRING GeneMania
# Count overlaps of the gene pairs in co-expression from STRING and GeneMania co-expression interactions.
sort  results/string_db/string_ppi_filt_noduplicated.txt results/genemania/genemania_ppi_sorted_noduplicated.txt | uniq -d > results/comparisons/gm_string_ppi_overlap_sort.txt
wc -l results/comparisons/gm_string_ppi_overlap_sort.txt > results/comparisons/gm_string_ppi_overlap_sort_count.txt
