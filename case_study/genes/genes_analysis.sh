#!/bin/bash

# assumes you have node_attributes.csv and interactions.csv in ../../datasets/genes_data/
# calculate graph features based on a gene neighborhood

check_datasets() {
  # check that files exist
  if [ ! -f "$1"/interactions.csv ]; then
    echo "Error: interaction dataset not found!"
    exit 1
  fi

  if [ ! -f "$1"/node_attributes.csv ]; then
    echo "Error: node attributes dataset not found!"
    exit 1
  fi
}

generate_neighborhood() {
  python two_hops_neighborhood.py
  Rscript hop_ratio_preproc.r
}

# Attempt to do a clean exit
trap cleanup SIGINT EXIT
cleanup() {
  exit 130
}

script_dir=$(dirname "$0")

# If no path is passed, assume dataset is two directories up
dataset_dir="$script_dir"/../datasets/genes_data/

# Sanity check if datasets exist
check_datasets "$dataset_dir"

# Generate neighborhood files
generate_neighborhood

base_dir="$(cd "$(dirname "$1")";cd ..;pwd -P)"
submodule_path="$base_dir/GraphSAGE"

# Make sure submoulde is cloned
(cd "$script_dir" && git submodule update --init --recursive)

# Build graphsage docker image if it doesn't exist
if [ ! "$(docker images graphsage -q)" ] ; then
  docker build -t graphsage -f "$submodule_path"/Dockerfile "$submodule_path"
fi

# Create folders for data
mkdir -p "$submodule_path"/genes/ppi
mkdir -p "$submodule_path"/genes/coexpression
mkdir -p "$submodule_path"/genes/epistasis

# Generate separate inputs for GraphSage in json
python genes_graphsage_input.py
graphsage_location="$submodule_path"

docker run -it -v "$graphsage_location":/notebooks graphsage bash -c "python graphsage/utils.py 'genes/ppi/genes-G.json' 'genes/ppi/genes-walks.txt' && python -m graphsage.unsupervised_train --train_prefix genes/ppi/genes --model graphsage_mean --max_total_steps 500 --validate_iter 10"

docker run -it -v "$graphsage_location":/notebooks graphsage bash -c "python graphsage/utils.py 'genes/coexpression/genes-G.json' 'genes/coexpression/genes-walks.txt' && python -m graphsage.unsupervised_train --train_prefix genes/coexpression/genes --model graphsage_mean --max_total_steps 500 --validate_iter 10"

docker run -it -v "$graphsage_location":/notebooks graphsage bash -c "python graphsage/utils.py 'genes/epistasis/genes-G.json' 'genes/epistasis/genes-walks.txt' && python -m graphsage.unsupervised_train --train_prefix genes/epistasis/genes --model graphsage_mean --max_total_steps 500 --validate_iter 10"

python genes_embs_preproc.py
Rscript modelling.r
python genes_hinsage.py > ../datasets/genes_data/results.txt 
Rscript interpret_results.r
