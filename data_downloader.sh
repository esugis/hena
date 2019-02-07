#!/bin/bash
echo "Downloading raw data sets for integration"
mkdir data

echo "Downloading epistasis data."
mkdir data/epistasis
wget https://ndownloader.figshare.com/files/9918205?private_link=e6f47fd4b00e39b4fe2b -O data/epistasis/ADNI_CT_epistasis.txt
wget https://ndownloader.figshare.com/files/9918220?private_link=b2c816eaa22c0e7273e8 -O data/epistasis/ADNI_VER_epistasis.tsv
wget https://ndownloader.figshare.com/files/9918214?private_link=ce0b39cccbe26214de12 -O data/epistasis/HBTRC_epistasis.tsv
wget https://ndownloader.figshare.com/files/9918217?private_link=a2c42f9aaeb04758bdff -O data/epistasis/TGEN_epistasis.tsv
echo "Done!"

echo "Downloading positive selection data."
mkdir data/ps
wget https://ndownloader.figshare.com/files/8943232?private_link=e53c14d0b033c5e09f53 -O data/ps/positive_darwinian_selection.csv
echo "Done!"

echo "Downloading GWAS data."
mkdir data/gwas
wget https://ndownloader.figshare.com/files/9919579?private_link=9d1e67084f4b6f98dc96 -O data/gwas/IGAP_stage_1_2_combined.txt
echo "Done!"

echo "Downloading PPI related to brain ageing data."
mkdir data/pba
wget https://ndownloader.figshare.com/files/13993100 -O data/pba/PBA_PPI_HS.txt
echo "Done!"

echo "Downloading gene expression in disease and healthy samples."
mkdir data/adn
wget https://ndownloader.figshare.com/files/8753989?private_link=187aad2ba5bdf5775cb4 -O data/adn/E-GEOD-18309.nc
wget https://ndownloader.figshare.com/files/8753992?private_link=187aad2ba5bdf5775cb4 -O data/adn/E-GEOD-28146.nc
wget https://ndownloader.figshare.com/files/8753983?private_link=187aad2ba5bdf5775cb4 -O data/adn/E-GEOD-4757.nc
wget https://ndownloader.figshare.com/files/8753977?private_link=187aad2ba5bdf5775cb4 -O data/adn/E-MEXP-2280.nc
wget https://ndownloader.figshare.com/files/8753986?private_link=187aad2ba5bdf5775cb4 -O data/adn/E-GEOD-29652.nc
wget https://ndownloader.figshare.com/files/8753980?private_link=187aad2ba5bdf5775cb4 -O data/adn/E-GEOD-5281.nc
echo "Done!"

echo "Downloading PPI data sets from IntAct."
mkdir data/intact
wget https://ndownloader.figshare.com/files/9924304?private_link=45d221ae3860a96d35f9 -O data/intact/intact_hs_v_4_2_6.txt
wget https://ndownloader.figshare.com/files/9924376?private_link=ffe149df62963ab7e1bb -O data/intact/alzheimers_intact_v_4_2_6.txt
wget https://ndownloader.figshare.com/files/9924400?private_link=d23c343174f7bda8e345 -O data/intact/synapse_intact_v_4_2_6.txt
echo "Done!"

echo "Downloading gene expression data sets from Human Allen Brain Atlas."
mkdir data/allenbrain
mkdir data/allenbrain/178236545_ds
mkdir data/allenbrain/178238266_ds
mkdir data/allenbrain/178238316_ds
mkdir data/allenbrain/178238359_ds
mkdir data/allenbrain/178238387_ds
mkdir data/allenbrain/178238373_ds

wget http://human.brain-map.org/api/v2/well_known_file_download/178236545 -O data/allenbrain/normalized_microarray_donor15697.zip
unzip data/allenbrain/normalized_microarray_donor15697.zip -d data/allenbrain/178236545_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238359 -O data/allenbrain/normalized_microarray_donor12876.zip
unzip data/allenbrain/normalized_microarray_donor12876.zip -d data/allenbrain/178238359_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238316 -O data/allenbrain/normalized_microarray_donor14380.zip
unzip data/allenbrain/normalized_microarray_donor14380.zip -d data/allenbrain/178238316_ds 
wget http://human.brain-map.org/api/v2/well_known_file_download/178238266 -O data/allenbrain/normalized_microarray_donor15496.zip
unzip data/allenbrain/normalized_microarray_donor15496.zip -d data/allenbrain/178238266_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238387 -O data/allenbrain/normalized_microarray_donor9861.zip
unzip data/allenbrain/normalized_microarray_donor9861.zip -d data/allenbrain/178238387_ds
wget http://human.brain-map.org/api/v2/well_known_file_download/178238373 -O data/allenbrain/normalized_microarray_donor10021.zip
unzip data/allenbrain/normalized_microarray_donor10021.zip -d data/allenbrain/178238373_ds
echo "Done!"


# Download precomputed data
echo "Starting to download precomputed data sets. Data will be downloaded in .txt and .Rdata formats."
mkdir results
echo "Downloading precomputed data set of co-expression in disease and normal samples."
mkdir results/adn
mkdir results/adn/integration
wget https://ndownloader.figshare.com/files/14184956?private_link=42ff273d27a5bb40eaeb -O results/adn/integration/alzcoexp_int.RData
wget https://ndownloader.figshare.com/files/13993085?private_link=42ff273d27a5bb40eaeb -O results/adn/integration/alzcoexp_int.txt
echo "Done!"

mkdir results/allenbrain
echo "Downloading precomputed data set of aggregated gene expression in brain regions."
wget https://ndownloader.figshare.com/files/13178600?private_link=6bd5651259759a77b67a -O results/allenbrain/brain_tissue_zscores_aggreg.RData
echo "Done!"

echo "Downloading precomputed data sets of co-expression in disease related brain regions."
wget https://ndownloader.figshare.com/files/14184959?private_link=06f0d1a6d7d74b1bf4e0 -O results/allenbrain/CA1_coexp_int.RData
wget https://ndownloader.figshare.com/files/14000897?private_link=06f0d1a6d7d74b1bf4e0 -O results/allenbrain/CA1_coexp_int.txt

wget https://ndownloader.figshare.com/files/14184962?private_link=d3e83e8f1559b87e7b28 -O results/allenbrain/CA2_coexp_int.RData
wget https://ndownloader.figshare.com/files/13998587?private_link=d3e83e8f1559b87e7b28 -O results/allenbrain/CA2_coexp_int.txt

wget https://ndownloader.figshare.com/files/14184965?private_link=2ea2fb119b31526219c3 -O results/allenbrain/CA3_coexp_int.RData
wget https://ndownloader.figshare.com/files/14000012?private_link=2ea2fb119b31526219c3 -O results/allenbrain/CA3_coexp_int.txt

wget https://ndownloader.figshare.com/files/14184974?private_link=90ad68063ece2c0e854d -O results/allenbrain/CA4_coexp_int.RData
wget https://ndownloader.figshare.com/files/14000702?private_link=90ad68063ece2c0e854d -O results/allenbrain/CA4_coexp_int.txt

wget https://ndownloader.figshare.com/files/14184977?private_link=62db540f6f03e763cd42 -O results/allenbrain/DG_coexp_int.RData
wget https://ndownloader.figshare.com/files/14000852?private_link=62db540f6f03e763cd42 -O results/allenbrain/DG_coexp_int.txt

wget https://ndownloader.figshare.com/files/14184980?private_link=1a9f9e64c63134aa2b2c -O results/allenbrain/SptN_coexp_int.RData
wget https://ndownloader.figshare.com/files/14000885?private_link=1a9f9e64c63134aa2b2c -O results/allenbrain/SptN_coexp_int.txt

wget https://ndownloader.figshare.com/files/14184983?private_link=05d0c33f8a7265f3911a -O results/allenbrain/subiculum_coexp_int.RData
wget https://ndownloader.figshare.com/files/14000891?private_link=05d0c33f8a7265f3911a -O results/allenbrain/subiculum_coexp_int.txt
echo "Done!"


echo "Downloading precomputed epistasis data sets."
mkdir results/epistasis
wget https://ndownloader.figshare.com/files/14185007?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_tgen_int.txt
wget https://ndownloader.figshare.com/files/14184992?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_tgen_int.RData

wget https://ndownloader.figshare.com/files/14184995?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_hbtrc_int.txt
wget https://ndownloader.figshare.com/files/14184998?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_hbtrc_int.RData

wget https://ndownloader.figshare.com/files/14185001?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_adni_ver_int.txt
wget https://ndownloader.figshare.com/files/14185004?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_adni_ver_int.RData

wget https://ndownloader.figshare.com/files/14185010?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_adni_cog_int.txt
wget https://ndownloader.figshare.com/files/14185013?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epi_adni_cog_int.RData

wget https://ndownloader.figshare.com/files/14184986?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epistasis_all_int.txt
wget https://ndownloader.figshare.com/files/14184989?private_link=d6fdd4f6394e45f23f00 -O results/epistasis/epistasis_all_int.RData
echo "Done!"

echo "Downloading preprocessed PPI data sets from IntAct."
mkdir results/intact
wget https://ndownloader.figshare.com/files/14185028?private_link=76c1bf613581b62699f7 -O results/intact/syn_intact_int.txt
wget https://ndownloader.figshare.com/files/14185031?private_link=76c1bf613581b62699f7 -O results/intact/syn_intact_int.RData

wget https://ndownloader.figshare.com/files/14185034?private_link=76c1bf613581b62699f7 -O results/intact/intact_int.txt
wget https://ndownloader.figshare.com/files/14185037?private_link=76c1bf613581b62699f7 -O results/intact/intact_int.RData

wget https://ndownloader.figshare.com/files/14185040?private_link=76c1bf613581b62699f7 -O results/intact/alz_intact_int.txt
wget https://ndownloader.figshare.com/files/14185043?private_link=76c1bf613581b62699f7 -O results/intact/alz_intact_int.RData
echo "Done!"

echo "Downloading preprocessed PPI data set associated with brain ageing."
mkdir results/pba
wget https://ndownloader.figshare.com/files/14185052?private_link=45b14507dd2e3a0807fd -O results/pba/pba_int.txt
wget https://ndownloader.figshare.com/files/14185055?private_link=45b14507dd2e3a0807fd -O results/pba/pba_int.RData
echo "Done!"

echo "Downloading positive selection data."
mkdir results/ps
wget https://ndownloader.figshare.com/articles/7635881?private_link=153d29faff4c4a60d0c2 -O results/ps/ps.txt
wget https://ndownloader.figshare.com/articles/7635881?private_link=153d29faff4c4a60d0c2 -O results/ps/ps.RData
echo "Done!"

echo "Downloading preprocessed GWAS data."
mkdir results/gwas
wget https://ndownloader.figshare.com/files/14185082?private_link=e95ab740073249913d51 -O results/gwas/gwas_ensg.txt
wget https://ndownloader.figshare.com/files/14185085?private_link=e95ab740073249913d51 -O results/gwas/gwas_ensg.RData
echo "Done!"

echo "Downloading data set of integrated interactions and node attributes."
mkdir results/integration

wget https://ndownloader.figshare.com/files/14185094?private_link=cf5a98092dd0c10e9790 -O results/integration/integrated_int.RData
wget https://ndownloader.figshare.com/files/14185109?private_link=cf5a98092dd0c10e9790 -O results/integration/integrated_int.txt

wget https://ndownloader.figshare.com/files/14185124?private_link=6ef6836d311a90ccce49 -O results/integration/node_attributes.RData
wget https://ndownloader.figshare.com/files/14185127?private_link=6ef6836d311a90ccce49 -O results/integration/node_attributes.txt
echo "Done!"

echo "All precomputed data is now located in the corresponding individual folders in the catalogue hena/results/"

mkdir comparison_public_db
mkdir case_study
mkdir case_study/datasets
mkdir case_study/datasets/genes_data
