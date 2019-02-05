# This script installs required packages

#install.packages("r-base",repos = 'http://cran.us.r-project.org')
#install.packages("build-essential", repos = 'http://cran.us.r-project.org')
#install.packages("netcdf-bin", repos = 'http://cran.us.r-project.org')
#install.packages("libudunits2-dev", repos = 'http://cran.us.r-project.org')
#install.packages("libnetcdf-dev", repos = 'http://cran.us.r-project.org')
#install.packages("gfortran", repos = 'http://cran.us.r-project.org')
#install.packages("liblapack-dev", repos = 'http://cran.us.r-project.org')
#install.packages("libopenblas-dev", repos = 'http://cran.us.r-project.org')

install.packages("gProfileR",repos = 'http://cran.us.r-project.org')
library(gProfileR) 

install.packages("reshape",repos = 'http://cran.us.r-project.org')
library(reshape)

install.packages("reshape2",repos = 'http://cran.us.r-project.org')
library(reshape2)

install.packages("foreach",repos = 'http://cran.us.r-project.org')
library(foreach)

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

install.packages("stringr",repos = 'http://cran.us.r-project.org')
library(stringr)

install.packages("R.utils",repos = 'http://cran.us.r-project.org')
library(R.utils)

install.packages("ncdf4",repos = 'http://cran.us.r-project.org')
library(ncdf4)

install.packages("RobustRankAggreg",repos = 'http://cran.us.r-project.org')
library(RobustRankAggreg)

install.packages("doMC",repos = 'http://cran.us.r-project.org')
library(doMC)

install.packages("tidyr",repos = 'http://cran.us.r-project.org')
library(tidyr)

install.packages("dplyr",repos = 'http://cran.us.r-project.org')
library(dplyr)

source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
library(preprocessCore)

install.packages("psych",repos = 'http://cran.us.r-project.org')
library("psych")

install.packages("doParallel",repos = 'http://cran.us.r-project.org')
library(doParallel)

install.packages("plyr",repos = 'http://cran.us.r-project.org')
library(plyr)

install.packages("ggplot2",repos = 'http://cran.us.r-project.org')
library(ggplot2)

install.packages("UpSetR")
library(UpSetR)

install.packages("igraph")
library(igraph)
