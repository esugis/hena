"""
Example usage:
python genes_embs_preproc.py

"""
import numpy as np
import pandas as pd

# load embeddings
def graphsage_preproc(edge_type, path):
    filename = "../GraphSage/unsup-{}/graphsage_mean_small_0.000010/".format(edge_type)
    embeds = np.load(filename + "val.npy")
    ids = pd.read_csv(filename + "val.txt", header=None)
    embeds_dt = pd.DataFrame(embeds)
    embeds_dt["ensg"] = ids
    embeds_dt.to_csv("{}/{}_graphsage_embs.csv".format(path, edge_type))

if __name__ == "__main__":
    path = "../datasets/genes_data"

    # calculate the files for each subgraph
    for t in ["ppi", "coexpression", "epistasis"]:
        graphsage_preproc(t, path)

    print("Embeddings are preprocessed.")

