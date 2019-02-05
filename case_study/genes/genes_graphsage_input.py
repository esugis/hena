# -*- coding: utf-8 -*-
#
# Copyright 2018 Data61, CSIRO
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This script prepares input files for graphSage embedding calculation
Requires:
- gene ids and attributes in the csv format.
- gene interactions
Check README for more details.

"""
import pandas as pd
import networkx as nx
import numpy as np

import json
from networkx.readwrite import json_graph


def graphsage_input_files(nodes, edges, edge_type, path, output):
    interactions_type = edges[(edges.int_type == edge_type) | (edges.int_type == edge_type.upper())]
    edgelist = interactions_type.iloc[:, 0:2]
    print("Data loaded.")

    # ugly hack, saving json that is accepted by GraphSage only works when it is loaded from csv,
    # not from pandas
    edgelist.to_csv(path + '/edgelist.csv', header=None, index=None)
    G_edgelist = nx.read_edgelist(path + '/edgelist.csv', delimiter=',')
    print("Edgelist created.")

    # split the data into train/val/test
    np.random.seed(35456)
    nodes['rnd'] = np.random.uniform(0, 1, len(nodes))
    nodes['train'] = nodes['rnd'] <= .5
    nodes['val'] = nodes['rnd'].between(0.5, 1, inclusive=False)
    nodes['test'] = nodes['rnd'] >= 1

    # assign labels and split indicators to the graph as node attributes
    for n in G_edgelist.nodes():
        G_edgelist.node[n]['label'] = nodes.label[nodes.ensg == n].tolist()[0]
        G_edgelist.node[n]['val'] = nodes.val[nodes.ensg == n].bool()
        G_edgelist.node[n]['test'] = nodes.test[nodes.ensg == n].bool()
    print("Edgelist is modified.")

    # convert the graph into the networkx format required by GraphSage
    graph_dt = json_graph.node_link_data(G_edgelist)
    with open(output + edge_type + '/genes-G.json', 'w') as fp:
        json.dump(graph_dt, fp)

    id_map = nodes['ensg'].to_dict()
    inv_id_map = {v: k for k, v in id_map.items()}
    with open(output + edge_type + '/genes-id_map.json', 'w') as fp:
        json.dump(inv_id_map, fp)

    class_map = nodes.set_index('ensg')['label'].to_dict()
    with open(output + edge_type + '/genes-class_map.json', 'w') as fp:
        json.dump(class_map, fp)

    features = nodes.iloc[:, 2:20] # only graph-related features
    features.fillna(inplace=True, value=0)
    features.reset_index(inplace=True, drop=True)
    np.save(output + edge_type + '/genes-feats', features)

    print("All files saved.")

if __name__ == "__main__":
    path = "../datasets/genes_data"
    output = "../GraphSAGE/genes/"
    node_attributes = pd.read_csv(path + "/nodes_hops.csv", sep=",")
    interactions = pd.read_csv(path + "/interactions.csv", sep=" ")

    # calculate the files for each subgraph
    for t in ["ppi", "coexpression", "epistasis"]:
        graphsage_input_files(node_attributes, interactions, t, path, output)



