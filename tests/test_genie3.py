import unittest
from app.genie3.genie3 import Genie3
import scanpy as sc
import os
import pandas as pd
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from app.model import GeneRegulatoryNetworkParams


class TestGenie3(unittest.TestCase):

    def test_genie3_small(self):
        adata = sc.read_h5ad(os.path.join(
            os.path.dirname(__file__), "data", "in", "h5ad", "GSE43580.h5ad"))
        with open(os.path.join(
                os.path.dirname(__file__), "data", "in", "tfs", "tf.txt"), "r") as f:
            tf_list = [line.strip() for line in f.read().splitlines()]
        genie3 = Genie3(
            adata, tf_list, GeneRegulatoryNetworkParams(n_trees=50))
        weights, links = genie3.run()
        self.assertTrue(weights.shape[0] == 1030)
        self.assertTrue(weights.shape[1] == adata.shape[1])
        self.assertTrue(links.shape[1] == 3)

    def test_genie3_large(self):
        adata = sc.read_h5ad(os.path.join(
            os.path.dirname(__file__), "data", "in", "h5ad", "SKCM.h5ad"))
        with open(os.path.join(
                os.path.dirname(__file__), "data", "in", "tfs", "tf.txt"), "r") as f:
            tf_list = [line.strip() for line in f.read().splitlines()]
        genie3 = Genie3(
            adata, tf_list, GeneRegulatoryNetworkParams(n_trees=50))
        weights, links = genie3.run()
        self.assertTrue(weights.shape[0] == 1030)
        self.assertTrue(weights.shape[1] == adata.shape[1])
        self.assertTrue(links.shape[1] == 3)
