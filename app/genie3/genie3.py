from anndata import AnnData
import numpy as np
import scipy.sparse as sp
import pandas as pd
import multiprocessing
from rpy2.robjects import pandas2ri, r
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from app.utils import get_R_script_base_path
from app.model import GeneRegulatoryNetworkParams


class Genie3:
    def __init__(self, adata: AnnData, tf_list: list[str], params: GeneRegulatoryNetworkParams):
        self.adata = adata
        self.tf_list = tf_list
        self.params = params

    def run(self):
        # Get gene and cell names
        genes = self.adata.var_names.tolist()
        cells = self.adata.obs_names.tolist()

        # Create expression dataframe (genes as rows, cells as columns)
        # For GENIE3, we need genes as rows and cells as columns
        # AnnData has cells as rows and genes as columns, so we need to transpose
        if isinstance(self.adata.X, np.ndarray):
            # For dense array, just transpose
            expr_matrix_t = self.adata.X.T
        else:
            # For sparse matrix, convert to CSR first if needed then transpose
            expr_matrix_t = sp.csr_matrix(self.adata.X).T.toarray()

        # Create DataFrame with appropriate index and columns
        expr_df = pd.DataFrame(expr_matrix_t, index=genes, columns=cells)

        # Filter for only TFs that exist in the dataset
        existing_tfs = [tf for tf in self.tf_list if tf in genes]
        print(
            f"Found {len(existing_tfs)} out of {len(self.tf_list)} TFs in the dataset")

        if len(existing_tfs) == 0:
            raise ValueError(
                "None of the provided TFs were found in the dataset")

        # Initialize R
        print("Initializing R environment...")
        ro.r.source(f"{get_R_script_base_path()}/genie3/genie3.R")

        # Convert expression data to R
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_expr_matrix = ro.conversion.py2rpy(expr_df)

            # Convert TF list to R
            r_tf_list = ro.StrVector(existing_tfs)

            # Run GENIE3
            n_trees = self.params["n_trees"]
            n_threads = self.params.get(
                "n_threads", multiprocessing.cpu_count())
            print(
                f"Running GENIE3 with {len(existing_tfs)} TFs, {n_trees} trees, and {n_threads} threads...")
            r_run_genie3 = r['run_genie3']
            r_weight_matrix = r_run_genie3(
                r_expr_matrix,
                r_tf_list,
                ntrees=n_trees,
                nthreads=n_threads
            )

            r_extract_network_links = r['extract_network_links']
            r_network_links = r_extract_network_links(r_weight_matrix)

        return r_weight_matrix, r_network_links[r_network_links["weight"] > 0].copy()
