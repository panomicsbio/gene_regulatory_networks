#' Run GENIE3 analysis on expression data
#'
#' @param expr_matrix Expression matrix with genes as rows and samples as columns
#' @param tf_list List of transcription factor gene names
#' @param targets_list Optional list of target genes (default: all genes)
#' @param ntrees Number of trees in the random forest (default: 100)
#' @param nthreads Number of parallel threads (default: 4)
#' @return A weight matrix of regulatory interactions
run_genie3 <- function(expr_matrix,
                       tf_list,
                       targets_list = NULL,
                       ntrees = 100,
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

    gc()

    cat("GENIE3 analysis completed successfully.\n")
    return(as.data.frame(weight_matrix))
}

#' Extract network links from GENIE3 weight matrix
#'
#' @param weight_matrix Weight matrix from GENIE3
#' @param regulators_per_target Optional number of top regulators to keep for each target
#' @param threshold Optional global threshold for weights (ignores regulators_per_target if provided)
#' @return Data frame with regulatory links (TF, target, weight)
extract_network_links <- function(weight_matrix,
                                  regulators_per_target = NULL,
                                  threshold = NULL) {
    # Check if weight_matrix is a data frame and convert to matrix if needed
    if (is.data.frame(weight_matrix)) {
        message("Converting input data frame to matrix...")
        weight_matrix <- as.matrix(weight_matrix)
    }
    # Ensure it's a matrix now
    if (!is.matrix(weight_matrix)) {
        stop("Input must be a matrix or data frame")
    }
    if (!is.null(threshold) && !is.null(regulators_per_target)) {
        stop("Cannot provide both threshold and regulators_per_target")
    }
    if (!is.null(threshold)) {
        # Use global threshold
        cat("Extracting links with weight threshold:", threshold, "\n")
        link_list <- GENIE3::getLinkList(weight_matrix, threshold = threshold)
        return(link_list)
    }
    if (!is.null(regulators_per_target)) {
        # Use top N regulators per target
        cat("Extracting top", regulators_per_target, "regulators per target\n")
        link_list <- GENIE3::getLinkList(weight_matrix,
            reportMax = regulators_per_target
        )
        return(link_list)
    }
    link_list <- GENIE3::getLinkList(weight_matrix)
    return(link_list)
}
