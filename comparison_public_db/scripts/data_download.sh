#!/bin/bash

# We are in the comparison catalogue
# Create data folder
mkdir data_db

# Data folder for genemania dataset
mkdir data_db/genemania
mkdir data_db/genemania/coexp
mkdir data_db/genemania/ppi
mkdir data_db/genemania/genetic_int
wget -A "Co-expression*" -R http://genemania.org/data/current/Homo_sapiens/ --no-check-certificate -O data_db/genemania/coexp

cd data_db/genemania/coexp/genemania.org/data/current/Homo_sapiens
rm -f "Predicted*.txt"

# Move files to the dedicated folders
find . -name '*Co-expression*' -exec mv -i {} ~/hena/comparison/data_db/genemania/coexp/ \;
find . -name '*Physical_Interactions*' -exec mv -i {} ~/hena/comparison/data_db/genemania/ppi/ \;
find . -name '*Genetic_Interactions*' -exec mv -i {} ~/hena/comparison/data_db/genemania/genetic_int/ \;


# Data folder folder for string_db
mkdir data_db/string_db
cd data_db/string_db/

#wget https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz
#wget https://stringdb-static.org/download/protein.links.detailed.v10.5/9606.protein.links.detailed.v10.5.txt.gz
#wget https://stringdb-static.org/download/protein.links.full.v10.5/9606.protein.links.full.v10.5.txt.gz
#wget https://stringdb-static.org/download/protein.actions.v10.5/9606.protein.actions.v10.5.txt.gz

wget https://ndownloader.figshare.com/files/14629310 -O 9606.protein.actions.v10.5.txt
wget https://ndownloader.figshare.com/files/14629331 -O 9606.protein.links.detailed.v10.5.txt
wget https://ndownloader.figshare.com/files/14629337 -O 9606.protein.links.v10.5.txt
wget https://ndownloader.figshare.com/files/14629475 -O 9606.protein.links.full.v10.5.txt

# Unzip
#gunzip 9606.protein.actions.v10.5.txt.gz
#gunzip 9606.protein.links.detailed.v10.5.txt.gz
#gunzip 9606.protein.links.full.v10.5.txt.gz
#gunzip 9606.protein.links.v10.5.txt.gz

# Create forlder for results
mkdir  results
mkdir results/genemania
mkdir results/string_db
mkdir results/integrated_int
mkdir results/comparisons
