import keras
import numpy as np
import pandas as pd
from typing import List
import pdb
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from stellargraph.core.graph import StellarGraph
from stellargraph.layer.hinsage import HinSAGE
import keras
from keras import backend as K
import numpy as np

from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    average_precision_score,
    roc_auc_score,
)
import os


class GeneGraph:
    def __init__(self, edge_data_file, gene_attr_file):
        # Read edge data
        edge_data = pd.read_csv(edge_data_file, delimiter=' ')

        edge_data["ensg1"] = edge_data["ensg1"].str.replace("ENSG", "").astype(int)
        edge_data["ensg2"] = edge_data["ensg2"].str.replace("ENSG", "").astype(int)

        # Read attribute data
        # gene_attr = pd.read_csv(gene_attr_file, sep="\t").drop_duplicates("ensg")
        gene_attr = pd.read_csv(gene_attr_file, sep=",").drop_duplicates("ensg")
        gene_attr = gene_attr[gene_attr.ensg.str.contains("LRG") == False]
        gene_attr["ensg"] = gene_attr["ensg"].str.replace("ENSG", "").astype(int)

        gene_attr.set_index("ensg", inplace=True)
        print(gene_attr["class"].value_counts())
        gene_attr_known = gene_attr[gene_attr["class"] != 2]
        print(gene_attr_known["class"].value_counts())

        # List of IDs
        rnd = np.random.uniform(0, 1, len(gene_attr_known))
        self.ids_train = gene_attr_known.loc[rnd <= 0.5].index.values
        self.ids_val = gene_attr_known.loc[(rnd > 0.5) & (rnd < 0.7)].index.values
        self.ids_test = gene_attr_known.loc[rnd >= 0.7].index.values
        self.ids_unknown = gene_attr[gene_attr["class"] == 2].index.values

        # Features
        # self.feats = gene_attr.drop(["node_type"], axis=1)
        self.feats = gene_attr.drop(["class"], axis=1)
        # create an all-zeros feature vector for the special node with ind -1, i.e., a non-existant node:
        self.feats.loc[-1] = [0] * self.feats.shape[1]

        # Labels
        # self.labels = gene_attr["node_type"].map(lambda x: x == "alz")
        self.labels = gene_attr["class"].map(lambda x: x == 1)

        # YT: create separate adjacency lists, one per edge type:
        edge_types = np.unique(edge_data["int_type"])
        self.adj = dict()
        for edge_type in edge_types:
            self.adj[edge_type] = (
                edge_data.loc[edge_data["int_type"] == edge_type]
                    .groupby(["ensg1"])["ensg2"]
                    .apply(list)
            )

            self.adj[edge_type][-1] = []

        for edge_type in edge_types:
            missing_genes = set(gene_attr.index).difference(
                set(self.adj[edge_type].index)
            )
            adj_ext = pd.Series(
                [list()] * len(missing_genes), index=missing_genes
            )  # form a series to append to self.adj[edge_type]
            self.adj[edge_type] = self.adj[edge_type].append(
                adj_ext, verify_integrity=True
            )  # append adj_ext

    def get_feats(self, indices: List[int]):
        """Get features of nodes whose list is given in indices"""
        return self.feats.loc[indices].fillna(0).as_matrix()

    def get_labels(self, indices: List[int]):
        return np.array(self.labels[indices], dtype=np.float64)

    def sample_neighs(self, indices: List[int], ns: int):
        """Neighbour sampling method"""

        def with_adj(adj_curr):
            if ns > 0:
                return [
                    [-1] * ns
                    if ((not isinstance(adj, list)) and pd.isnull(adj)) or len(adj) == 0
                    else [
                        adj[i] for i in np.random.randint(len(adj), size=ns)
                    ]  # YT: sample ns neighbours of each node in indices
                    for adj in adj_curr.loc[indices].values
                ]
            else:
                return [[-1] for adj in adj_curr.loc[indices].values]
        return tuple([with_adj(adj) for adj in self.adj.values()])

    def get_batch(
            self, indices: List[int], ns: List[int]
    ):
        nb = len(indices)
        flatten = lambda l: [item for sublist in l for item in sublist]

        neigh_1hop = self.sample_neighs(indices, ns[0])
        neigh_2hop = dict()
        for i, et in enumerate(self.adj.keys()):
            neigh_2hop[et] = self.sample_neighs(flatten(neigh_1hop[i]), ns[1])

        return (
            self.get_labels(indices),
            [
                self.get_feats(flatten(inds)).reshape([nb, -1, self.feats.shape[1]])
                for inds in [
                [indices],
                *neigh_1hop,
                *[n for n in neigh_2hop[et] for et in self.adj.keys()],
            ]
            ],
        )


class DataGenerator(keras.utils.Sequence):
    """Generates data for Keras"""

    def __init__(
            self,
            g: GeneGraph,
            nf: int,
            ns: List[int],
            batch_size: int = 1000,
            name: str = "train",
    ):
        """Initialization"""
        if isinstance(batch_size, int):
            self.batch_size = batch_size
        else:
            raise Exception(
                "DataGenerator: batch_size should be of type int, got {}".format(
                    type(batch_size)
                )
            )

        if isinstance(nf, int):
            self.nf = nf
        else:
            raise Exception(
                "DataGenerator: nf should be of type int, got {}".format(type(nf))
            )

        self.g = g
        self.ns = ns
        if name == "train":
            self.ids = g.ids_train
        elif name == "validate":
            self.ids = g.ids_val
        else:
            raise Exception(
                'DataGenerator: name is {}; should be either "train" or "validate"'.format(
                    name
                )
            )

        self.data_size = len(self.ids)
        self.idx = 0
        self.name = name
        self.on_epoch_end()  # shuffle the data entries

    def __len__(self):
        """Denotes the number of batches per epoch"""
        return int(np.ceil(self.data_size / self.batch_size))

    def __getitem__(self, index):
        """Generate one batch of data"""
        if self.idx > self.data_size:
            raise Exception(
                "DataGenerator: index {} exceeds data size {}. This shouldn't happen!".format(
                    self.idx, self.data_size
                )
            )
        elif self.idx == self.data_size:
            #print(
            #     "DataGenerator: index {} is equal to data size {}. Calling self.on_epoch_end()...".format(
            #         self.idx, self.data_size
            #     )
            # )
            self.on_epoch_end()
        else:
            pass

        end = min(self.idx + self.batch_size, self.data_size)
        indices = list(self.ids[range(self.idx, end)])
        tgt, inp = self.g.get_batch(indices, self.ns)
        self.idx = end

        return inp, tgt

    def on_epoch_end(self):
        """Updates indexes after each epoch"""
        self.idx = 0
        if self.name == "train":
            np.random.shuffle(self.ids)


class TestDataGenerator(keras.utils.Sequence):
    """Generates data for Keras"""

    def __init__(self, g: GeneGraph, nf: int, ns: List[int], batch_size: int = 1000):
        """Initialization"""
        self.batch_size = batch_size
        self.g = g
        self.ids = g.ids_test
        self.data_size = len(self.ids)
        self.nf = nf
        self.ns = ns
        self.idx = 0
        self.y_true = []

    def __len__(self):
        """Denotes the number of batches per epoch"""
        return int(np.ceil(self.data_size / self.batch_size))

    def __getitem__(self, index):
        """Generate one batch of data"""
        end = min(self.idx + self.batch_size, self.data_size)
        indices = list(self.ids[range(self.idx, end)])
        tgt, inp = self.g.get_batch(indices, self.ns)
        self.y_true += [tgt]
        self.idx = end
        return inp

    def on_epoch_end(self):
        """Updates indexes after each epoch"""
        self.idx = 0

class UnknownDataGenerator(keras.utils.Sequence):
    """Generates data for Keras"""

    def __init__(self, g: GeneGraph, nf: int, ns: List[int], batch_size: int = 1000):
        """Initialization"""
        self.batch_size = batch_size
        self.g = g
        self.ids = g.ids_unknown
        self.data_size = len(self.ids)
        self.nf = nf
        self.ns = ns
        self.idx = 0

    def __len__(self):
        """Denotes the number of batches per epoch"""
        return int(np.ceil(self.data_size / self.batch_size))

    def __getitem__(self, index):
        """Generate one batch of data"""
        end = min(self.idx + self.batch_size, self.data_size)
        indices = list(self.ids[range(self.idx, end)])
        tgt, inp = self.g.get_batch(indices, self.ns)
        #self.y_true += [tgt]
        self.idx = end

        # print('TestDataGenerator: index={}, batch size={}, self.idx={}, self.data_size={}'.format(index, len(indices),
        #                                                                                            self.idx,
        #                                                                                            self.data_size))

        return inp

    def on_epoch_end(self):
        """Updates indexes after each epoch"""
        self.idx = 0


class GeneHinSageClassifier(object):
    """Gene classification model with HinSAGE+logistic regression architecture"""

    def __init__(self, nf, n_samples, emb_dim=256):
        """

        :param nf: number of node features
        :param n_samples: list of numbers of node samples, per edge type, for each layer (hop)
        :param emb_dim: dimensionality of the hidden node representations (embeddings)
        """
        self.nf = nf
        self.n_samples = n_samples
        self.emb_dim = emb_dim
        n_samples_model = [ns if ns > 0 else 1 for ns in n_samples]
        self.build_model(n_samples_model)  # build the model

    def build_model(self, n_samples):
        def n_at(i):
            return np.product(n_samples[:i])

        def create_weighted_binary_crossentropy(zero_weight, one_weight):
            def weighted_binary_crossentropy(y_true, y_pred):
                b_ce = K.binary_crossentropy(y_true, y_pred)
                weight_vector = y_true * one_weight + (1. - y_true) * zero_weight
                weighted_b_ce = weight_vector * b_ce
                return K.mean(weighted_b_ce)

            return weighted_binary_crossentropy

        hs = HinSAGE(
            output_dims=[self.emb_dim] * len(n_samples),
            n_samples=n_samples,
            input_neigh_tree=[
                ("gene", [1, 2, 3]),
                ("gene", [4, 5, 6]),
                ("gene", [7, 8, 9]),
                ("gene", [10, 11, 12]),
                ("gene", []),
                ("gene", []),
                ("gene", []),
                ("gene", []),
                ("gene", []),
                ("gene", []),
                ("gene", []),
                ("gene", []),
                ("gene", []),
            ],
            input_dim={"gene": self.nf},
        )

        x_inp = [
            keras.Input(shape=(1, self.nf)),
            keras.Input(shape=(n_at(1), self.nf)),
            keras.Input(shape=(n_at(1), self.nf)),
            keras.Input(shape=(n_at(1), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
            keras.Input(shape=(n_at(2), self.nf)),
        ]

        x_out = keras.layers.Reshape((self.emb_dim,))(hs(x_inp))
        pred = keras.layers.Activation("sigmoid")(keras.layers.Dense(1)(x_out))

        # stack the layers together into a model, and compile the model with desired loss:
        self.model = keras.Model(inputs=x_inp, outputs=pred)

        self.model.compile(
            optimizer=keras.optimizers.Adam(lr=0.01),
            loss=create_weighted_binary_crossentropy(
                0.6, 4
            ),
            metrics=["accuracy"],
        )

    def train(self, g: GeneGraph, epochs=1):
        train_iter = DataGenerator(g, self.nf, self.n_samples, name="train")
        valid_iter = DataGenerator(g, self.nf, self.n_samples, name="validate")

        self.model.fit_generator(
            train_iter,
            epochs=epochs,
            validation_data=valid_iter,
            max_queue_size=10,
            shuffle=True,
            verbose=2,
        )

    def test(self, g: GeneGraph, threshold=0.5):
        test_iter = TestDataGenerator(g, self.nf, self.n_samples)
        y_preds_proba = self.model.predict_generator(test_iter)
        y_preds_proba = np.reshape(y_preds_proba, (-1,))
        y_preds = np.array(y_preds_proba >= threshold, dtype=np.float64)
        y_trues = np.concatenate(test_iter.y_true).ravel()[
                  : len(g.ids_test)
                  ]
        # Evaluate metrics (binary classification task):
        precision = precision_score(y_trues, y_preds)
        recall = recall_score(y_trues, y_preds)
        average_precision = average_precision_score(y_trues, y_preds_proba)
        f1score = f1_score(y_trues, y_preds)
        roc_auc = roc_auc_score(y_trues, y_preds_proba)
        idx = g.ids_test

        met = lambda f, k: sum(
            [
                f(y_trues[i], y_preds[i]) and y_preds[i] == k
                for i in range(len(g.ids_test))
            ]
        )
        conf_matrix = {
            "tn": met(lambda t, p: t == p, 0),
            "fp": met(lambda t, p: t != p, 1),
            "fn": met(lambda t, p: t != p, 0),
            "tp": met(lambda t, p: t == p, 1),
        }

        return precision, recall, f1score, average_precision, roc_auc, conf_matrix, y_preds_proba, idx

    def predict(self, g: GeneGraph, threshold=0.5):
        pred_iter = UnknownDataGenerator(g, self.nf, self.n_samples)
        y_preds_proba = self.model.predict_generator(pred_iter)
        y_preds_proba = np.reshape(y_preds_proba, (-1,))
        y_preds = np.array(y_preds_proba >= threshold, dtype=np.float64)
        idx = g.ids_unknown
        return y_preds_proba, y_preds, idx


if __name__ == "__main__":
    data_dir = "../datasets/genes_data/"
    edge_data_fname = os.path.join(data_dir, "interactions.csv")
    feature_sets = ["biological_features.csv", "graph_features.csv", "all_features.csv"]
    n_samples = [2, 5]

    for fs in feature_sets:
        print("Start processing " + fs)
        print("Reading graph...")
        gene_attr_fname = os.path.join(data_dir, fs)
        g = GeneGraph(edge_data_fname, gene_attr_fname)
        nf = g.feats.shape[1]
        # Create a model:
        gene_model = GeneHinSageClassifier(nf, n_samples, emb_dim=64)

        print(
             "Training the {}-layer model with n_samples={}...".format(
                 len(n_samples), n_samples
             )
        )
        gene_model.train(g, epochs=10)
        print("Evaluating the model on test set...")
        threshold = 0.5
        precision, recall, f1score, average_precision, roc_auc, conf_matrix, y_preds_proba_test, ids_test = gene_model.test(
            g, threshold
        )
        print("Precision score: {0:0.2f}".format(precision))
        print("Recall score: {0:0.2f}".format(recall))
        print("F1 score: {0:0.2f}".format(f1score))
        print("Average precision-recall score: {0:0.2f}".format(average_precision))
        print("ROC AUC score: {0:0.2f}".format(roc_auc))
        print("Confusion matrix for threshold={}:".format(threshold))
        print(conf_matrix)

        y_preds_proba, y_preds, ids = gene_model.predict(g, threshold)
        results = pd.DataFrame(ids)
        results.columns = ['id']
        results['pred'] = y_preds_proba
        results['ensg'] = results['id'].astype(str).apply(lambda x: x.zfill(11))
        results['ensg'] = 'ENSG' + results['ensg'].astype(str)

        #results['ensg'] = 'ENSG00000' + results['id'].astype(str)
        output_name = "hinsage_{}".format(fs)
        output_fname = os.path.join(data_dir, output_name)
        results.to_csv(output_fname)

        test_results = pd.DataFrame(ids_test)
        test_results.columns = ['id']
        test_results['pred'] = y_preds_proba_test
        test_results['ensg'] = test_results['id'].astype(str).apply(lambda x: x.zfill(11))
        test_results['ensg'] = 'ENSG' + test_results['ensg'].astype(str)

        output_name_test = "hinsage_test_{}".format(fs)
        output_fname_test = os.path.join(data_dir, output_name_test)
        test_results.to_csv(output_fname_test)




