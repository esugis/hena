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
Calculation of two-hop neighborhood ratios and counts for the gene use case
Check README for more details.
Input: nodes and interactions
Output: csv files of alzheimer genes count in the radius of 1 and 2-hop neighborhood
for different edge types
"""
import numpy as np
import pandas as pd
import networkx as nx


def hop_calculation(nodes, edges, edge_type, depth, path):
    interactions_type = edges[(edges.int_type == edge_type) | (edges.int_type==edge_type.upper())]
    edgelist = interactions_type.iloc[:, 0:2]
    G_edgelist = nx.from_pandas_dataframe(edgelist, source="ensg1", target="ensg2")
    filename = '{}/hops_{}_{}.csv'.format(path, edge_type, depth)
    with open(filename, 'w', 10) as f:
        for i, ensg in enumerate(G_edgelist.nodes()):
            if i % 1000 == 0:
                print("Processed %d ENSGs" % i)
            hops = nx.single_source_shortest_path_length(G_edgelist, ensg, cutoff=depth)
            hop_k = {k: v for k, v in hops.items() if v == depth}
            hopk_labels = nodes[nodes.ensg.isin(hop_k)]
            f.write('"%s", %d, %d\n' % (ensg, np.sum(hopk_labels.label == 'disease'), len(hopk_labels)))


if __name__ == "__main__":
    path = "../datasets/genes_data"
    node_attributes = pd.read_csv(path + "/node_attributes.csv", sep=" ")
    interactions = pd.read_csv(path + "/interactions.csv", sep=" ")

    # run the hop calculation
    for t in ["ppi", "coexpression", "epistasis"]:
        for d in [1, 2]:
            hop_calculation(node_attributes, interactions, t, d, path)
