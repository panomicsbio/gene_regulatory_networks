"""
GENIE3 Gene Regulatory Network Analysis Implementation
This script performs Gene Regulatory Network analysis using GENIE3 on an AnnData object.
It bridges Python and R using rpy2.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from rpy2.robjects import pandas2ri, r
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
import scipy.sparse as sp

# Define the R script as a string
r_script = """
# Function to install packages if not already installed
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package == "GENIE3") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("GENIE3")
    } else {
      install.packages(package)
    }
  }
  library(package, character.only = TRUE)
}

# Install necessary packages
install_if_missing("GENIE3")
install_if_missing("igraph")
install_if_missing("RColorBrewer")

#' Run GENIE3 analysis on expression data
#'
#' @param expr_matrix Expression matrix with genes as rows and samples as columns
#' @param tf_list List of transcription factor gene names
#' @param targets_list Optional list of target genes (default: all genes)
#' @param regulators_per_target Maximum number of regulators per target (default: 10)
#' @param ntrees Number of trees in the random forest (default: 1000)
#' @param nthreads Number of parallel threads (default: 4)
#' @return A weight matrix of regulatory interactions
run_genie3 <- function(expr_matrix, 
                       tf_list,
                       targets_list = NULL,
                       regulators_per_target = 10,
                       ntrees = 1000,
                       nthreads = 4) {
  
  # Validate inputs
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("Expression data must be a matrix or data frame")
  }
  
  # Convert to matrix if data.frame (GENIE3 requires a matrix)
  if (is.data.frame(expr_matrix)) {
    gene_names <- rownames(expr_matrix)
    cell_names <- colnames(expr_matrix)
    expr_matrix <- as.matrix(expr_matrix)
    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- cell_names
    cat("Converted data frame to matrix for GENIE3\n")
  }
  
  if (is.null(tf_list) || length(tf_list) == 0) {
    stop("TF list cannot be empty")
  }
  
  # Ensure TFs are in the expression matrix
  valid_tfs <- intersect(tf_list, rownames(expr_matrix))
  if (length(valid_tfs) == 0) {
    stop("None of the provided TFs are in the expression matrix")
  }
  
  if (length(valid_tfs) < length(tf_list)) {
    warning(paste(length(tf_list) - length(valid_tfs), "TFs not found in expression matrix and will be ignored"))
  }
  
  # Run GENIE3
  cat("Running GENIE3 with", length(valid_tfs), "transcription factors...\n")
  
  weight_matrix <- GENIE3::GENIE3(
    exprMatrix = expr_matrix,
    regulators = valid_tfs,
    targets = targets_list,
    nTrees = ntrees,
    nCores = nthreads,
    verbose = TRUE
  )
  
  cat("GENIE3 analysis completed successfully.\n")
  return(weight_matrix)
}

#' Extract network links from GENIE3 weight matrix
#'
#' @param weight_matrix Weight matrix from GENIE3
#' @param regulators_per_target Number of top regulators to keep for each target
#' @param threshold Optional global threshold for weights (ignores regulators_per_target if provided)
#' @return Data frame with regulatory links (TF, target, weight)
extract_network_links <- function(weight_matrix, 
                                  regulators_per_target = 10,
                                  threshold = NULL) {
  
  # Get links
  if (!is.null(threshold)) {
    # Use global threshold
    cat("Extracting links with weight threshold:", threshold, "\n")
    link_list <- GENIE3::getLinkList(weight_matrix, threshold = threshold)
  } else {
    # Use top N regulators per target
    cat("Extracting top", regulators_per_target, "regulators per target\n")
    link_list <- GENIE3::getLinkList(weight_matrix, 
                                    reportMax = regulators_per_target)
  }
  
  return(link_list)
}

#' Create and save network visualization
#'
#' @param link_list Data frame with columns: regulator, target, weight 
#' @param weight_threshold Minimum weight to include in visualization
#' @param max_links Maximum number of links to include in visualization
#' @param output_file File path to save the visualization
#' @param layout Graph layout algorithm (default: "fr" for Fruchterman-Reingold)
#' @return igraph object of the network
visualize_network <- function(link_list,
                             weight_threshold = 0.05,
                             max_links = 200,
                             output_file = "grn_network.pdf",
                             layout = "fr") {
  
  # Filter links by weight
  filtered_links <- link_list[link_list$weight >= weight_threshold,]
  
  # Take top links if needed
  if (nrow(filtered_links) > max_links) {
    filtered_links <- filtered_links[order(filtered_links$weight, decreasing = TRUE),][1:max_links,]
    cat("Taking top", max_links, "links for visualization\n")
  }
  
  # Create graph
  cat("Creating network with", nrow(filtered_links), "edges and", 
      length(unique(c(filtered_links$regulatoryGene, filtered_links$targetGene))), "nodes\n")
  
  g <- igraph::graph_from_data_frame(
    d = filtered_links[, c("regulatoryGene", "targetGene", "weight")],
    directed = TRUE,
    vertices = NULL
  )
  
  # Get regulators (in-degree = 0) and targets
  regulators <- names(which(igraph::degree(g, mode = "in") == 0))
  
  # Set node attributes
  V(g)$type <- ifelse(V(g)$name %in% regulators, "regulator", "target")
  V(g)$color <- ifelse(V(g)$type == "regulator", "red", "blue")
  V(g)$shape <- ifelse(V(g)$type == "regulator", "square", "circle")
  V(g)$size <- ifelse(V(g)$type == "regulator", 10, 5)
  
  # Set edge attributes
  E(g)$width <- E(g)$weight * 10
  
  # Create layout
  layout_coords <- igraph::layout_(g, switch(layout,
                                           "fr" = igraph::with_fr(),
                                           "kk" = igraph::with_kk(),
                                           "circle" = igraph::in_circle(),
                                           igraph::with_fr()))
  
  # Plot and save the network
  pdf(output_file, width = 12, height = 10)
  igraph::plot.igraph(g,
                    layout = layout_coords,
                    vertex.label.cex = 0.7,
                    vertex.label.color = "black",
                    vertex.label.font = 2,
                    edge.arrow.size = 0.5,
                    edge.curved = 0.2,
                    main = "Gene Regulatory Network")
  dev.off()
  
  cat("Network visualization saved to", output_file, "\n")
  return(g)
}

#' Calculate network statistics
#'
#' @param graph igraph network object
#' @return List of network statistics
calculate_network_stats <- function(graph) {
  stats <- list(
    num_nodes = igraph::vcount(graph),
    num_edges = igraph::ecount(graph),
    density = igraph::edge_density(graph),
    avg_path_length = mean(igraph::shortest.paths(graph), na.rm = TRUE),
    diameter = igraph::diameter(graph, directed = TRUE, weights = NA),
    clustering_coefficient = igraph::transitivity(graph)
  )
  
  # Add degree distributions
  in_deg <- igraph::degree(graph, mode = "in")
  out_deg <- igraph::degree(graph, mode = "out")
  stats$max_in_degree <- max(in_deg)
  stats$max_out_degree <- max(out_deg)
  stats$avg_in_degree <- mean(in_deg)
  stats$avg_out_degree <- mean(out_deg)
  
  # Top hubs (nodes with highest out-degree)
  top_hubs <- sort(out_deg, decreasing = TRUE)
  stats$top_5_hubs <- names(top_hubs[1:min(5, length(top_hubs))])
  
  # Output stats
  cat("Network Statistics:\n")
  cat("Number of nodes:", stats$num_nodes, "\n")
  cat("Number of edges:", stats$num_edges, "\n")
  cat("Network density:", round(stats$density, 4), "\n")
  cat("Average path length:", round(stats$avg_path_length, 4), "\n")
  cat("Network diameter:", stats$diameter, "\n")
  cat("Clustering coefficient:", round(stats$clustering_coefficient, 4), "\n")
  cat("Top 5 hub nodes:", paste(stats$top_5_hubs, collapse = ", "), "\n")
  
  return(stats)
}

#' Save network as graphml file for Cytoscape
#'
#' @param graph igraph network object
#' @param file_path Output file path
#' @return None
export_to_graphml <- function(graph, file_path = "grn_network.graphml") {
  igraph::write_graph(graph, file_path, format = "graphml")
  cat("Network exported to", file_path, "for use in Cytoscape\n")
}
"""


def run_genie3_analysis(anndata_file, tf_list, output_dir=None, n_top_tfs=10,
                        weight_threshold=0.05, n_trees=50, n_threads=8):
    """
    Run GENIE3 Gene Regulatory Network analysis on an AnnData object

    Parameters:
    -----------
    anndata_file : str or AnnData
        Path to h5ad file or AnnData object containing normalized log-transformed counts
    tf_list : list
        List of transcription factor gene names
    output_dir : str, optional
        Directory to save output files (default: current directory)
    n_top_tfs : int, optional
        Number of top TFs to retain per target gene (default: 10)
    weight_threshold : float, optional
        Threshold for link weights (default: 0.05)
    n_trees : int, optional
        Number of trees in the random forest (default: 1000)
    n_threads : int, optional
        Number of threads to use (default: 4)

    Returns:
    --------
    dict
        Dictionary with GENIE3 results, including:
        - weight_matrix: Pandas DataFrame of regulatory weights
        - links: Pandas DataFrame of regulatory links
        - network: networkx Graph object
        - stats: Dictionary of network statistics
    """
    # Create output directory if it doesn't exist
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    # Load AnnData object if string is provided
    if isinstance(anndata_file, str):
        print(f"Loading AnnData from {anndata_file}")
        adata = sc.read_h5ad(anndata_file)
    else:
        adata = anndata_file

    print(f"AnnData object: {adata.shape[1]} genes, {adata.shape[0]} cells")

    # Get gene and cell names
    genes = adata.var_names.tolist()
    cells = adata.obs_names.tolist()

    # Create expression dataframe (genes as rows, cells as columns)
    # For GENIE3, we need genes as rows and cells as columns
    # AnnData has cells as rows and genes as columns, so we need to transpose
    if isinstance(adata.X, np.ndarray):
        # For dense array, just transpose
        expr_matrix_t = adata.X.T
    else:
        # For sparse matrix, convert to CSR first if needed then transpose
        expr_matrix_t = sp.csr_matrix(adata.X).T.toarray()

    # Create DataFrame with appropriate index and columns
    expr_df = pd.DataFrame(expr_matrix_t, index=genes, columns=cells)

    # Filter for only TFs that exist in the dataset
    existing_tfs = [tf for tf in tf_list if tf in genes]
    print(f"Found {len(existing_tfs)} out of {len(tf_list)} TFs in the dataset")

    if len(existing_tfs) == 0:
        raise ValueError("None of the provided TFs were found in the dataset")

    # Initialize R
    print("Initializing R environment...")
    r(r_script)

    # Convert expression data to R
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr_matrix = ro.conversion.py2rpy(expr_df)

    # Convert TF list to R
    r_tf_list = ro.StrVector(existing_tfs)

    # Run GENIE3
    print(
        f"Running GENIE3 with {len(existing_tfs)} TFs and {n_trees} trees...")
    r_run_genie3 = r['run_genie3']
    r_weight_matrix = r_run_genie3(
        r_expr_matrix,
        r_tf_list,
        regulators_per_target=n_top_tfs,
        ntrees=n_trees,
        nthreads=n_threads
    )

    # Extract the links
    print(f"Extracting top {n_top_tfs} regulators per target...")
    r_extract_network_links = r['extract_network_links']
    r_links = r_extract_network_links(
        r_weight_matrix, regulators_per_target=n_top_tfs)

    # Visualize the network
    output_file = os.path.join(output_dir, "grn_network.pdf")
    r_visualize_network = r['visualize_network']
    r_graph = r_visualize_network(
        r_links,
        weight_threshold=weight_threshold,
        output_file=output_file
    )

    # Calculate network statistics
    r_calculate_network_stats = r['calculate_network_stats']
    r_stats = r_calculate_network_stats(r_graph)

    # Export network for Cytoscape
    graphml_file = os.path.join(output_dir, "grn_network.graphml")
    r_export_to_graphml = r['export_to_graphml']
    r_export_to_graphml(r_graph, file_path=graphml_file)

    # Convert R results back to Python
    with localconverter(ro.default_converter + pandas2ri.converter):
        weight_matrix = ro.conversion.rpy2py(r_weight_matrix)
        links = ro.conversion.rpy2py(r_links)

    # Return results
    results = {
        "weight_matrix": weight_matrix,
        "links": links,
        "output_dir": output_dir,
        "graphml_file": graphml_file,
        "visualization_file": output_file
    }

    print("GRN analysis complete!")
    return results


def plot_top_regulators(results, n_top=50, output_file=None):
    """
    Plot top TF regulators based on their overall regulatory influence

    Parameters:
    -----------
    results : dict
        Results from run_genie3_analysis
    n_top : int, optional
        Number of top TFs to show (default: 20)
    output_file : str, optional
        Path to save the plot (default: None, displays the plot)
    """
    links = results["links"]

    # Aggregate weights by regulator
    reg_influence = links.groupby("regulatoryGene")[
        "weight"].sum().sort_values(ascending=False)

    # Plot top regulators
    plt.figure(figsize=(12, 8))
    top_regs = reg_influence.head(n_top)

    sns.barplot(x=top_regs.values, y=top_regs.index)
    plt.title(
        f"Top {n_top} Transcription Factor Regulators by Overall Influence")
    plt.xlabel("Sum of Regulatory Weights")
    plt.ylabel("Transcription Factor")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def plot_regulatory_network(results, n_edges=100, output_file=None):
    """
    Plot the regulatory network using networkx

    Parameters:
    -----------
    results : dict
        Results from run_genie3_analysis
    n_edges : int, optional
        Number of top edges to include (default: 100)
    output_file : str, optional
        Path to save the plot (default: None, displays the plot)
    """
    import networkx as nx

    links = results["links"]

    # Take top N edges
    top_links = links.sort_values("weight", ascending=False).head(n_edges)

    # Create network
    G = nx.DiGraph()

    # Add edges with weights
    for _, row in top_links.iterrows():
        G.add_edge(row["regulatoryGene"],
                   row["targetGene"], weight=row["weight"])

    # Identify TFs (nodes with out-degree > 0)
    tfs = [node for node, out_deg in G.out_degree() if out_deg > 0]

    # Set node attributes
    node_colors = ["red" if node in tfs else "blue" for node in G.nodes()]
    node_sizes = [300 if node in tfs else 100 for node in G.nodes()]

    # Set edge attributes
    edge_widths = [G[u][v]["weight"] * 5 for u, v in G.edges()]

    # Create layout
    pos = nx.spring_layout(G, seed=42)

    # Plot
    plt.figure(figsize=(12, 10))
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.8)
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.5, edge_color="gray",
                           connectionstyle="arc3,rad=0.1", arrowsize=10)
    nx.draw_networkx_labels(G, pos, font_size=8)

    plt.title(f"Top {n_edges} Regulatory Interactions")
    plt.axis("off")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Plot saved to {output_file}")
    else:
        plt.show()

    return G


# Example usage
if __name__ == "__main__":
    import scanpy as sc
    import numpy as np

    adata = sc.read_h5ad("tests/data/h5ad/GSE43580.h5ad")

    n_genes, n_cells = adata.shape[1], adata.shape[0]

    with open("tf.txt", "r") as f:
        tf_list = f.read().splitlines()

    # Run GENIE3 analysis
    results = run_genie3_analysis(adata, tf_list)

    # Plot results
    plot_top_regulators(results, output_file="top_regulators.png")
    plot_regulatory_network(results, output_file="regulatory_network.png")
