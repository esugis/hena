# HENA: Heterogeneous network-based data set related to Alzheimer's disease.
The scripts in this repository were used to create integrated dataset described in the article:
*HENA: Heterogeneous network-based data set related to Alzheimer's disease, Sügis et. al., 2019*.      
*The integrated dataset consists of two parts: edge list of interactions deposited in one file and the file with node(gene) attributes.*
*For detailed description of the methods please see the publication.*

## IMPORTANT. BEFORE WE START

**NB! This file contains information about the structure of the files and some useful notes on how to repeat the analysis.
Please be aware that the analysis of the collected heterogeneous data is highly memory and computational time consuming.
It is highly recommended to execute it on the server!**

To exicute only *case study* analysis with the application of Graph Convolutional Networks you can download nesessary files from figshare repository (see Section Case Study below).


## INSTALLATION REQUIREMENTS
 R (version > 3.4.2) that can be installed from [here](https://cran.r-project.org/). 


## GETTING STARTED
1. Clone repository to your local machine or to the server use the following command:  
```
git clone https://github.com/esugis/hena.git  
```
Or download and unpack project .zip archive.
By default the project will be cloned or unpacked to the folder `hena`.  

## STRUCTURE
Navigate to the the folder `hena` main default catalogue of the project.

In the terminal you can do it by typing:

``` 
cd hena
```
All the following commands will be executed assuming you are located in the folder  `hena`.

Project is divided into 3 major parts:  
- **data**  a place to store downloaded datasets  
- **scripts** a place to store analysis scripts  
- **results** a place to store the results  

Additionally, there 2 more parts describing the comaprison with the public databases and the case study using GCNs. The details are provided in the corresponding sections (see below).

All instructions provided below assume the following project structure:
   
```
hena
|── README.md
|── data_downloader.sh
|── libraries.R
|── analysis_adn.R
|── analysis_allen_brain.R
|── analysis_epistasis.R
|── analysis_gwas.R
|── analysis_positive_selection.R
|── analysis_ppi_intact.R
|── analysis_ppi_related_to_brain_ageing.R
|
|── data
│   |── adn
│   |── allenbrain
│   |── epistasis
│   |── gwas
│   |── intact
│   |── pba
│   └── ps
|
|── scripts
|   |── libraries.R
|   |── adn
|   │   |── E_GEOD_18309.R
|   │   |── E_GEOD_28146.R
|   │   |── E_GEOD_29652.R
|   │   |── E_GEOD_4757.R
|   │   |── E_GEOD_5281.R
|   │   |── E_MEXP_2280.R
|   │   |── RRA_probesets.R
|   │   |── affy2ensg.R
|   │   |── coexp2undirrected_selfloops_rm.R
|   │   |── coexp_int.R
|   │   └── final_adn.R
|   |── epistasis
|   │   |── combine.R
|   │   |── epi_adni.R
|   │   |── epi_adni_cog.R
|   │   |── epi_hbtrc.R
|   │   └── epi_tgen.R
|   |── intact
|   │   |── intact.R
|   │   |── alz_intact.R
|   │   └── synapse_intact.R
|   |── pba
|   │   └── pba_int.R
|   |── gwas
|   │   └── gwas2ensg.R
|   |── allenbrain
|   │   |── coexp_per_tissue_preprocessing.R
|   │   |── coexp_in_selected_brain_regions.R
|   │   |── expression_all_regions_zscores.R
|   │   └── p2ensg.R
|   |── ps
|   │   └── positive_selection.R
|   └── integration
|       |── integrate.R
|       |── integrate_node_attributes.R
|       └── plot_stats.R
|
|
└── results
|   |── adn
|   |   |── all_probes
|   |   └── integration
|   |── allenbrain 
|   |   |── all_probes
|   |   └── integration
|   |── epistasis 
|   |── gwas 
|   |── intact 
|   |── pba 
|   |── ps 
|   └── integration
|
|
└── comparison_public_db
|
└── case_study
```

## DATA

Data is stored in the folder  `hena/data`
For automatic data download run `data_downloader.sh` script from the main project directory in your terminal:  
```
./data_downloader.sh
```
After executing the script folder `data` will contain sub-folders for individual data types used in the analysis.  
All necessary datasets will be downloaded into the corresponding locations:

```
|── data
│   |── adn #contains microarray datasets in .nc format
│   |── allenbrain #contains microarray datasets downloaded from Allen Brain Atlas
│   |── epistasis #contains epistatic datasets computed based on ADNI, TGEN and HBTRC cohorts data
│   |── gwas #contains GWAS dataset
│   |── intact #contains PPI datasets downloaded from IntAct database
│   |── pba #contains PPI dataset related to brain ageing
│   |── ps #contains positive selection dataset
```

## SCRIPTS. USER MANUAL "REPEATING THE ANALYSIS".

### General information.
The whole analysis is divided into separate parts corresponding to theseparate analyses of the individual data sets.
Each datatype analisys is wrapped in to the corresponding R script.  
Each analysis (except for the integration part) cab be executed independently in any order.
To repeat these individual analyses please run the commands shown below from the root project catalogue:

First you would need to install all required R libraries by running the following command:

```
Rscript libraries.R
```
NB! In case you use Ubuntu 18.04 please install the following packages before running `script libraries.R`

`r-base build-essential netcdf-bin libudunits2-dev libnetcdf-dev gfortran liblapack-dev libopenblas-dev libxml2-dev libssl-dev`

Run the analysis of epistatic interactions:
```
Rscript analysis_epistasis.R
```
Run the analysis of GWAS data related to Alzheimer's disease:
```
Rscript analysis_gwas.R
```
Run the analysis of positive selection data related to Alzheimer's disease:
```
Rscript analysis_positive_selection.R
```
Run the analysis of protein-protein interactions related to Azheimer's disease, synaptic interactions and highly confident interactions in human:
```
Rscript analysis_ppi_intact.R
```
Run the analysis of protein-protein interactions related to brain ageing:
```
Rscript analysis_ppi_related_to_brain_ageing.R
```
Run the analysis of gene expression in six whole-brain microarrays from Allen Brain Atlas:
```
Rscript analysis_allen_brain.R
```
Run gene co-expression analysis in Alzheimer's disease and healthy samples. 
**We do not recommend to run this co-expression analysis on your laptop.** Precomputed file is located in folder  `results/adn/integration/`.
```
Rscript analysis_adn.R
```
To integrate all individual results  and plot interaction statistics run the following script:
```
Rscript integrate.R
```
As the result of the integration you will receive two final files describing all the interactions `integrated_int.txt` and  the corresponding nodes (genes) attributes  `node_attributes.txt`.

Additionally, after HENA has been assembled, script  `hena/scripts/integration/plot_stats.R` outputs a pdf file with the histogram describing the number of interactions of each individual interaction type coming from individual data sources and a table (as txt file) with the corresponding number of nodes and edges.
You execute it as following:
```
Rscript plot_stats.R
```

In case you want to understand more details about the individual steps of the project, please see the next section "SCRIPTS. INDIVIDUAL ANALYSIS STEPS EXPLAINED." NB! This section is mostly aimed for building a better intuition about the  analyses you have just performed.


## SCRIPTS. INDIVIDUAL ANALYSIS STEPS EXPLAINED.
This section provides the description of individual analysis steps and the corresponding scripts.
As a user you don't need to run these scripts. Instead please follow the instructions described in the section "SCRIPTS. USER MANUAL "REPEATING THE ANALYSIS".

#### Calculation of gene co-expression using microarray data. 
**Computations of all against co-expression between 50kx50k probes in every dataset with intermediate step of saving the results requires substantial operative memory and file storage resources. We do not recommend to run this co-expression analysis on your laptop.** Precomputed file is located in folder  `results/adn/integration/
`
1. Filtering out probesets with SD < 0.29. Extracting Alzheimer’s related and healthy samples from the data sets. Calculating the co-expression between  all probesets in each of the data sets using  Spearman correlation coefficient. For each individual probeset in each of the datasets script creates 2 separate files in .txt and .RData formats. Files are named after the probeset. Created files contain the names of the correlated probesets and the corresponding Spearman coefficient.
```
|── adn
│   |── E_GEOD_18309.R
│   |── E_GEOD_28146.R
│   |── E_GEOD_29652.R
│   |── E_GEOD_4757.R
│   |── E_GEOD_5281.R
│   |── E_MEXP_2280.R
```
Results will be written into files for individual probesets to the folder corresponding to the dataset name, i.e  `AgedBrainSYSBIO/results/adn/all_probes/rdata/E_GEOD_18309`

2. Ranking the co-expressed values for each probeset in each of the datasets. Aggregating the ranks in all the datasets.
```
    │   |── RRA_probesets.R  
```
3. Converting the values to ensg ids. The script also handles the cases when affymetrix probeset id was not recognised by the gconvert() function due to the presence of the unrecognised symbols in the probeset name.
```
    │   |── affy2ensg.R  
```
4. Assembling together calculated RRA scores for all the probes.
```
    │   |── coexp_int.R 
```
5. Removing “self loops”(co-expression of gene with itself) from the co-expression dataset.
```
    │   |── coexp2undirrected_selfloops_rm.R  
```
6. Adding columns interaction_type, data_source.
```
    │   └──  final_adn.R  
```

Script adds missing columns to match the integrated dataset format interaction_type="co-expression" and data_source="ADN”. The filnal co-expression part of the integrated dataset is stored in the following files:
```
`AgedBrainSYSBIO/results/adn/integration/adn_coexp_int.RData` 
`AgedBrainSYSBIO/results/adn/integration/adn_coexp_int.txt` 
```
#### Processing epistatic interactions.

1. epi_adni.R Script processes epistatic interactions associated with change in ventrical volume in ADNI cohort.
```
    |── epistasis             
    │   |── epi_adni.R  
```
2. Preprocessing epistatic interactions in HBTRC cohort dataset.
```
│   |── epi_hbtrc.R  
```
3. Preprocessing epistatic interactions in TGEN cohort dataset.
```
│   └── epi_tgen.R    
```
4. Preprocessing epistatic interactions associated with a set of cognitive traits in ADNI cohort dataset.
```
│   |── epi_adni_cog.R  
```
5. Coombining all epistatic data into one dataset.
```
│   |── combine.R 
```
#### Reading in and processing datasets dowloaded from IntAct database.
1. Preprocessing all human protein-protein interactions from IntAct database. 
```
|── intact  
│   |── intact.R  
```
2. Preprocessing protein-protein interactions from IntAct database where one interacting partner is associated with Alzheimre's disease.
```
│   |── alz_intact.R   
```
3. Preprocessing synaptic protein-protein interactions from IntAct database.
```
│   └── synapse_intact.R   
```
#### Reading in and processing of the protein-protein interaction dataset realted to brainageing generated by AgedBrainSYSBIO consortium.
```
|── pba  
│   └── pba_int.R     
``` 

#### Processing of GWAS data.   
```
|── gwas  
│   └── gwas2ensg.R      
```
#### Analysing gene expression datasets from Human Allen Brain Atlas. 

1. Convert probe names to Ensembl namespace.   
```
|   |── p2ensg.R  
```
2. Compute z-scores based on gene expression in brain regions over all six brains.
```
|   |── expression_all_regions_zscores.R
```
3.  Datasets preprocessing for computing co-expression in Alzheimer's disease associated regions.
```
|   |── coexp_per_tissue_preprocessing.R
```
4. Computing gene co-expression in CA1, CA2, CA3, CA4, DG, SptN and subiculum brain regions.
Brain region annotations are considered as provided by Allen Brain Atlas.
```
|   |── coexp_in_selected_brain_regions.R
```
As a results the `results/allenbrain` folder will contain the following files in RData and txt format:
```
|   |── CA1_coexp_int_aba.RData
|   |── CA2_coexp_int_aba.RData
|   |── CA3_coexp_int_aba.RData
|   |── CA4_coexp_int_aba.RData
|   |── SptN_coexp_int_aba.RData
|   |── subiculum_coexp_int_aba.RData

|   |── CA1_coexp_int_aba.txt
|   |── CA2_coexp_int_aba.txt
|   |── CA3_coexp_int_aba.txt
|   |── CA4_coexp_int_aba.txt
|   |── SptN_coexp_int_aba.txt
|   |── subiculum_coexp_int_aba.txt
|   |── DG_coexp_int_aba.txt
```

#### Reading in and processing of the dataset of positively selected genes generated by AgedBrainSYSBIO consortium.
```
└── ps  
    └── positive_selection.R     
```
### Integrate individual results for interactions and interactors.
```
└── integration  
```
1. Combine interactions.
```
    |── integrate_interactions.R  
```
2. Combine node (gene) attributes.
```
    └── integrate_node_attributes.R  
```
## RESULTS
1. The results of the individual data type analyses are save as .RData and .txt formats. They are placed into the newly created folder `hena/results`.      

```
└── results
    |── adn #contains resultig co-expression interactions
    |── allenbrain #contains co-expression intercation in brain regions related to Alzheimer's disease and aggregated expression in 231 brain regions
    |── epistasis #contains preprocessed epistatic interactions from ADNI, TGEN and HBTRC cohorts
    |── gwas #contains preprocessed data from GWAS studies
    |── intact #contains preprcessed PPI from IntAct database
    |── pba #contains preperocessed PPI related to brain ageing
    |── ps #contains preprocessed positive selection p-values
```
2. The resulted integrated dataset cosists of two parts: integrated interactions and node attributes. These parts are stored in the sub-folder `integration` correspondingly in the files `integrated_int.txt` and `integrated_int_attributes.txt`.  Additionally, ahistogram of interactions in the data set is save to the file `integrated_int_stats_stacked_hist.pdf`.
```
    └──integration #contains integrated interactions of individual data types and node attributes
```

## COMPARISON WITH THE PUBLIC DATABASES
Since there are no alternative Alzheimer's disease scpecific datasets available for comparison, we have benchmarked integrated dataset using resources that collect non-disease specific datasets. For that purpose GENEMANIA and STRING data were downloaded and arranged into groups of interactions of class co-expression, physical protein interactions (PPI) and genetic interactions (the latter was compared with epistatic intecations). 

The folder  `comparison_public_db` is located in the project root catalogue and has the following structure:

``` 
comparison_public_db
|
|── data_db
│    |── genemania
│    │   |── coexp
│    │   |── ppi
│    │   └── genetic_int
│    └── string_db
|
|── scripts
│        |── data_download.sh
│        |── sortlines_linux_amd64 
│        |── gene_mania_coexpression_preprocess.sh 
│        |── gene_mania_genint_preprocess.sh
│        |── gene_mania_int_ppi_preprocess.sh
│        |── extract_interactions_integrated_int.R
│        |── integrated_ppi_coexp_epi_sort_rm_duplicated_pairs.sh
│        |── string_coexp_ppi_preprocess.R
│        |── string_ppi_coexp_sort_rm_duplicated_pairs.sh
│        |── compare_datasets.sh
│        └── plot_comparisons.R
│
└── results
        |── genemania
        |── string_db
        |── integrated_int
        └── comparisons

```

To repeat the comparison, please navigate to the folder  `comparison_public_db`  from the project root directory in your terminal.
```
cd comparison_public_db
```
All the following commands will be executed assuming you are located in the folder  `comparison_public_db`.

### Repeating the comparison. User manual.

1. To download datasets from GeneMania and STRING databases; place co-expression, ppi, and genetic interactions into the corresponding folders please run the following command:
```
./scripts/data_download.sh
```
2. To prepare datasets of individual interaction types from GeneMania for the comparison run the following set of commands: 
```
./scripts/gene_mania_coexpression_preprocess.sh 
./scripts/gene_mania_genint_preprocess.sh
./scripts/gene_mania_int_ppi_preprocess.sh
```
3.  To extract interactions for the comaprison with the public resources from the integrated dataset run the following commands:
Selected interaction as a result will include positive co-expression, protein-protein interactions and epistatic interactions excluding intergenic region interactions. 
```
Rscript scripts/extract_interactions_integrated_int.R
./scripts/integrated_ppi_coexp_epi_sort_rm_duplicated_pairs.sh
```
4. To prepare STRING database data for the comparison please execute the following commands:
```
Rscript string_coexp_ppi_preprocess.R
./string_ppi_coexp_sort_rm_duplicated_pairs.sh
```
5. Run the comparison by executing the following:
```
./scripts/compare_datasets.sh
```
### Interpreting the results
The results folder contains sub-directories related to the individual comparison steps.
1. The preprocessed data for individual comaprisons are located in the folloeing directories:
```
|── genemania
|── string_db
|── integrated_int
```

2. We first fing the number of ppi, coexpression and genetic interactions in each individual data source, i.e. GeneMania, STRING, integrated dataset.
The counts are written in the corresponding files. Then the overlaps between the ineractions of each types in GeneMania, STRING and integated dataset are counted. These results are save to the ``comparisons`` direstory.

3. Visualizing the results.

The visualization of the overlapping interactions are created using UpsetterR tool wrapper for R that allows to plot vienn-diagram-like visualizations for the large number of datasets.

To create the plots run the following script:

```
Rscript scripts/plot_comparisons.R
```
As a result two files will be created:
co-expression.pdf and ppi.pdf that demostrate the number of interactions overlap in GeneMania, STRING and integrated dataset.


## CASE STUDY. APPLICATION OF GRAPH CONVOLUTIONAL NETWORKS TO HETEROGENEOUS GRAPHS.

## Alzheimer's gene classification use-case of graph analytics

## Overview
The goal of this case study is to identify genes that are associated with Alzheimer's disease using information about genes and interactions of different types between pairs of genes. The findings of this case study indicate that graph structure is a rich data source that helps to capture complex relationships and find the distinctive patterns in biological data that are not easily detectable otherwise. 

This repository contains the code to reproduces the results presented in the article by Sügis et al "Heterogeneous network-based dataset for Alzheimer's disease". 


## Data

### Data acquisition 
If you have repeated the assembley of HENA data set in the previous steps,  please copy two files `integrated_int.RData` and `node_attributes.RData` from  the project folder `hena/results/integration/` to `hena/case_study/data/`.

Alternatively you can download HENA interactions and node attributes from figshare repository and place them into the folder `hena/case_study/data/` .
HENA interactions file: `https://ndownloader.figshare.com/files/14185094`
HENA node attributes file: `https://ndownloader.figshare.com/files/14185124`

### Data set preparation for GCN application
For this case study we have excluded IGRI to keep interactions between nodes that were mapped directly to genes, and co-expression interactions in disease related brain regions with low co-expression values (< 0.5).

In order to perform classification task we have used gene status reagarding to its realtion to any disease in human based on the evolutionary studies review by Spataro et al. to define a negative class (see the article for details). This study identified ~1500 genes that have not shown  realtion to any of the diseases based on he evolutionary studies. In the node classification task these 1500 genes will constitute a negative class. The file with the the gene status is located in here `case_study/data/negative_class/gene_status.csv`.

Data set preparation for the application of GCN  is performed  by executing the following R scripts from the root project directory `hena`.
These scripts assume you have executed script `hena/data_downloader.sh`
```
Rscript case_study/scripts/highly_exp_genes_in_dis_brain_regions.R
Rscript case_study/scripts/prepare_data_GCN.R
```
This will result in two files `case_study/datasets/genes_data/interactions.csv` and  `case_study/datasets/genes_data/node_attributes.csv`. 

**NB! We recommend to follow the steps of the procedure as it might not run sucessfully otherwise.**

### Installation requirements
* Install Python 3.6.3. We recommend installing Anaconda that you can download from [here](https://conda.io/docs/user-guide/install/download.html) 
* You need to install [docker](https://docs.docker.com/)  
* Some of the steps require R version 3.4.2. It is highly recommended to install [RStudio](https://www.rstudio.com/products/rstudio/download/) as well.

### Running the experiment
1. It is recommended to create a conda environment using Anaconda:
```
conda create -n sgusecases python=3.6 anaconda
```
2. Activate the newly created environment `source activate sgusecases` and install the required python packages:
`pip install -r requirements.txt`

3. Open RStudio and copy the following lines to install R packages:
```
cran.packages <- c("igraph", "randomForest", "h2o", "ModelMetrics", "ggplot2", "tidyr", "dplyr", "data.table", "stringr", "gProfileR")
install.packages(cran.packages)
```
4. Navigate to `case_study` directory in your terminal and  clone GraphSAGE repository  `git clone https://github.com/williamleif/GraphSAGE.git`.
Install the requirements from `case_study/GraphSAGE` folder  `$ pip install -r requirements.txt`. Navigate back to `case_study`.
5. Clone the stellargraph repository `git clone https://github.com/stellargraph/stellargraph.git -b usecase/genes`
6. Navigate to `stellargraph` folder and install it `pip install -e .`
7. Once everything is installed navigate back to `case_study/genes/` and run `./genes_analysis.sh` that produces the files with all the results:

```
genes_data
├── all_features.csv
├── biological_features.csv
├── coexpression_graphsage_embs.csv
├── edgelist.csv
├── epistasis_graphsage_embs.csv
├── graph_features.csv
├── hinsage_all_features.csv
├── hinsage_biological_features.csv
├── hinsage_graph_features.csv
├── hops_coexpression_1.csv
├── hops_coexpression_2.csv
├── hops_epistasis_1.csv
├── hops_epistasis_2.csv
├── hops_ppi_1.csv
├── hops_ppi_2.csv
├── interactions.csv
├── node_attributes.csv
├── nodes_hops.csv
├── ppi_graphsage_embs.csv
├── reconstruction_mse.png
├── rf_predictions_graph_set.csv
├── rf_predictions_alzheimer.csv
└── hinsage_predictions_alzheimer.csv
```

Files  `case_study/datasets/genes_data/rf_predictions_alzheimer.csv` and  `case_study/datasets/genes_data/hinsage_predictions_alzheimer.csv`  contain the lists of genes shown by independent research to be associated with the disease, and the corresponding probabilities provided by random forest and HinSAGE. 


### Acknowledgements ###
The case study was developed together with Anna Leontjeva https://github.com/annitrolla.
We would like to acknowledge the funding support from the European Union’s Seventh Framework Programme for research, technological development and demonstration under grant agreement No 305299 (http://www.agedbrainsysbio.eu). We also give thanks to the Estonian Research Council grants (PSG59; IUT34-4); the European Regional Development Fund for CoE of Estonian ICT research EXCITE project. Additionally, we thank Yuriy Tyshetskiy and Kevin Jung who provided insight and expertise that greatly assisted the development of the case study, and Ivan Kuzmin for assisting in solving non-trivial technical issues.Additionally, we thank Yuriy Tyshetskiy and Kevin Jung who provided insight and expertise that greatly assisted the development of the case study, and Ivan Kuzmin for assisting in solving non-trivial technical issues.
