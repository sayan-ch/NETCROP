#!/usr/bin/env Rscript

# Rebuild the cleaned Twitch network from the public raw edge and feature
# tables. This provenance script replaces the legacy inline Rcpp component
# search with netOP's public largest_connected_component() API.

suppressPackageStartupMessages(library(netOP))
if (packageVersion("netOP") != "0.1.1") stop("This script requires netOP v0.1.1.")

script_dir <- function() {
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(...) NULL)
  if (!is.null(source_file)) return(dirname(normalizePath(source_file)))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
  }
  normalizePath(getwd())
}

root <- script_dir()
edges <- utils::read.csv(file.path(root, "full_twitch_edges.csv"))
features <- utils::read.csv(file.path(root, "full_twitch_features.csv"))

node1 <- edges$numeric_id_1 + 1L
node2 <- edges$numeric_id_2 + 1L
keep_upper <- node1 < node2
node1 <- node1[keep_upper]
node2 <- node2[keep_upper]
features$node <- features$numeric_id + 1L

cutoff_views <- stats::quantile(features$views, 0.1)
cutoff_life_time <- stats::quantile(features$life_time, 0.1)
keep_features <- features$language != "EN" & features$dead_account == 0 &
  features$views >= cutoff_views & features$life_time >= cutoff_life_time
twitch_sub <- features[keep_features, , drop = FALSE]
twitch_nodes <- twitch_sub$node
keep_edges <- node1 %in% twitch_nodes & node2 %in% twitch_nodes
i <- match(node1[keep_edges], twitch_nodes)
j <- match(node2[keep_edges], twitch_nodes)
A <- Matrix::sparseMatrix(i = i, j = j, x = 1,
                          dims = c(length(twitch_nodes), length(twitch_nodes)))
A <- A + Matrix::t(A)

lcc <- largest_connected_component(A)
node_id <- twitch_nodes[lcc$nodes]
twitch_A <- lcc$submatrix
edge_summary <- as.data.frame(summary(twitch_A))
edge_summary <- edge_summary[edge_summary$i != edge_summary$j, c("i", "j")]
community <- twitch_sub[match(node_id, twitch_sub$node), "language"]
community_table <- data.frame(
  user = seq_along(community),
  lang = community,
  langnum = as.integer(factor(community))
)

utils::write.csv(edge_summary, file.path(root, "Twitch_edge_list.csv"),
                 row.names = FALSE)
utils::write.csv(community_table, file.path(root, "Twitch_community.csv"),
                 row.names = FALSE)
save(twitch_A, community_table, node_id, file = file.path(root, "Twitch.rda"))
cat(sprintf("Wrote the %d-node Twitch largest connected component.\n",
            nrow(twitch_A)))
