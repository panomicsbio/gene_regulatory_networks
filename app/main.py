import os
import logging
import scanpy as sc
from anndata import AnnData
from app.utils import get_directory_files
from app.genie3.genie3 import Genie3
from app.model import GeneRegulatoryNetworkParams
import json

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_inputs(input_root_dir) -> tuple[AnnData, list[str], GeneRegulatoryNetworkParams]:
    adata_files = get_directory_files(os.path.join(input_root_dir, "adata"))
    if len(adata_files) != 1:
        raise ValueError(
            f"Expected exactly one AnnData file in {input_root_dir}/adata")
    adata = sc.read_h5ad(adata_files[0])
    tf_files = get_directory_files(os.path.join(input_root_dir, "tfs"))
    if len(tf_files) != 1:
        raise ValueError(
            f"Expected exactly one TF file in {input_root_dir}/tfs")
    with open(tf_files[0], 'r') as f:
        tf_list = [line.strip() for line in f.read().splitlines()]
    strings_json_path = os.path.join(input_root_dir, "strings.json")
    with open(strings_json_path, 'r') as f:
        raw_params = json.load(f)
    params = GeneRegulatoryNetworkParams(n_trees=int(
        raw_params["n_trees"]))
    return adata, tf_list, params


def main():
    logger.info("Starting GRN")

    input_root_dir = "/mount/data/in"
    output_root_dir = "/mount/data/out"

    os.makedirs(os.path.join(output_root_dir,
                             "genie_weights"), exist_ok=True)
    os.makedirs(os.path.join(output_root_dir,
                             "genie_links"), exist_ok=True)

    try:
        adata, tf_list, params = load_inputs(input_root_dir)
        genie3 = Genie3(adata, tf_list, params)
        weights, links = genie3.run()
        # Save the results
        weights.to_csv(os.path.join(output_root_dir,
                                    "genie_weights", "genie_weights.csv"))
        links.to_csv(os.path.join(output_root_dir,
                                  "genie_links", "genie_links.csv"))
    except Exception as e:
        logger.error(f"Failed to run GRN: {str(e)}")
        raise


if __name__ == "__main__":
    main()
