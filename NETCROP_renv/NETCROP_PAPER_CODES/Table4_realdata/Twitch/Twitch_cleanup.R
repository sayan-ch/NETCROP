setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table4_realdata", "Twitch")))

library(dplyr)
library(tidyr)
library(Matrix)

twitch_el <- read_csv(file.path("full_twitch_edges.csv"))
twitch_el <- twitch_el |>
  rename(node1 = numeric_id_1, node2 = numeric_id_2) |>
  mutate(node1 = node1 + 1, node2 = node2 + 1) |>
  filter(node1 < node2)

twitch_features <- read_csv(file.path("full_twitch_features.csv"))
twitch_features <- twitch_features |>
  rename(node = numeric_id) |>
  mutate(node = node + 1)

# clean to 32407 nodes
cutoff_views <- quantile(twitch_features$views, 0.1)
cutoff_life_time <- quantile(twitch_features$life_time, 0.1)

# remove EN speakers and dead accounts
# and nodes with views < cutoff views and life_time < cutoff life_time
twitch_sub <- twitch_features |>
  filter(language != 'EN') |>
  filter(
    dead_account == 0 &
      views >= cutoff_views &
      life_time >= cutoff_life_time
  )

nrow(twitch_sub)
table(twitch_sub$language)

twitch_nodes <- twitch_sub$node
twitch_el_sub <- twitch_el |>
  filter(node1 %in% twitch_nodes & node2 %in% twitch_nodes) |>
  mutate(node11 = match(node1, twitch_nodes), node22 = match(node2, twitch_nodes))

A <- sparseMatrix(i = twitch_el_sub$node11, j = twitch_el_sub$node22, x = 1,
                   dims = c(length(twitch_nodes), length(twitch_nodes)))

A <- A + t(A)

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// Treat NA / NaN as not TRUE for `adj != 0`, matching which(adj[i, ] != 0)
inline bool is_nonzero(double x) {
  return !ISNAN(x) && x != 0.0;
}

// Return largest component from CSR-style row adjacency
IntegerVector largest_component_from_csr(
    const std::vector<int>& row_ptr,
    const std::vector<int>& col_idx,
    int n
) {
  std::vector<unsigned char> visited(n, 0);
  std::vector<int> best;
  std::vector<int> comp;
  std::vector<int> queue;

  best.reserve(n);
  comp.reserve(n);
  queue.reserve(n);

  for (int v = 0; v < n; ++v) {
    if (visited[v]) continue;

    comp.clear();
    queue.clear();

    visited[v] = 1;
    queue.push_back(v);
    comp.push_back(v);

    std::size_t head = 0;

    while (head < queue.size()) {
      int u = queue[head++];

      for (int e = row_ptr[u]; e < row_ptr[u + 1]; ++e) {
        int w = col_idx[e];

        if (!visited[w]) {
          visited[w] = 1;
          queue.push_back(w);
          comp.push_back(w);
        }
      }
    }

    // which.max() in R returns the first maximum, so use strict >
    if (comp.size() > best.size()) {
      best = comp;
    }
  }

  IntegerVector out(best.size());
  for (std::size_t i = 0; i < best.size(); ++i) {
    out[i] = best[i] + 1; // convert back to R 1-based indexing
  }

  return out;
}

// [[Rcpp::export]]
IntegerVector cpp_largest_cc_dgC(SEXP adj_) {
  S4 adj(adj_);

  IntegerVector Dim = adj.slot("Dim");
  int n = Dim[0];
  int m = Dim[1];

  if (n != m) {
    stop("Input adjacency matrix must be square.");
  }

  if (n == 0) {
    return IntegerVector(0);
  }

  IntegerVector p = adj.slot("p");
  IntegerVector i = adj.slot("i");
  NumericVector x = adj.slot("x");

  // Build row adjacency from column-compressed dgCMatrix.
  // Original R code uses neighbors from rows:
  // neighbors[[row]] = which(adj[row, ] != 0)
  std::vector<int> row_counts(n, 0);

  for (int col = 0; col < m; ++col) {
    for (int idx = p[col]; idx < p[col + 1]; ++idx) {
      int row = i[idx];
      if (is_nonzero(x[idx])) {
        ++row_counts[row];
      }
    }
  }

  std::vector<int> row_ptr(n + 1, 0);
  for (int r = 0; r < n; ++r) {
    row_ptr[r + 1] = row_ptr[r] + row_counts[r];
  }

  std::vector<int> col_idx(row_ptr[n]);
  std::vector<int> cursor = row_ptr;

  // Fill in increasing column order, matching which(adj[row, ] != 0)
  for (int col = 0; col < m; ++col) {
    for (int idx = p[col]; idx < p[col + 1]; ++idx) {
      int row = i[idx];

      if (is_nonzero(x[idx])) {
        col_idx[cursor[row]++] = col;
      }
    }
  }

  return largest_component_from_csr(row_ptr, col_idx, n);
}

// [[Rcpp::export]]
IntegerVector cpp_largest_cc_dense(NumericMatrix adj) {
  int n = adj.nrow();
  int m = adj.ncol();

  if (n != m) {
    stop("Input adjacency matrix must be square.");
  }

  if (n == 0) {
    return IntegerVector(0);
  }

  std::vector<unsigned char> visited(n, 0);
  std::vector<int> best;
  std::vector<int> comp;
  std::vector<int> queue;

  best.reserve(n);
  comp.reserve(n);
  queue.reserve(n);

  for (int v = 0; v < n; ++v) {
    if (visited[v]) continue;

    comp.clear();
    queue.clear();

    visited[v] = 1;
    queue.push_back(v);
    comp.push_back(v);

    std::size_t head = 0;

    while (head < queue.size()) {
      int u = queue[head++];

      // Scan row u in increasing column order, matching which(adj[u, ] != 0)
      for (int w = 0; w < n; ++w) {
        if (!visited[w] && is_nonzero(adj(u, w))) {
          visited[w] = 1;
          queue.push_back(w);
          comp.push_back(w);
        }
      }
    }

    // Match R which.max behavior: first largest component wins
    if (comp.size() > best.size()) {
      best = comp;
    }
  }

  IntegerVector out(best.size());
  for (std::size_t i = 0; i < best.size(); ++i) {
    out[i] = best[i] + 1;
  }

  return out;
}
')


find_largest_cc <- function(adj) {
  n <- nrow(adj)

  if (is.null(n)) {
    stop("Input must be a matrix-like object.")
  }

  if (n != ncol(adj)) {
    stop("Input adjacency matrix must be square.")
  }

  # Fast path for sparse Matrix objects
  if (inherits(adj, "sparseMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package Matrix is required for sparse matrices.")
    }

    # Convert to dgCMatrix for fast C++ access.
    # For symmetric Matrix classes, this materializes the implied full matrix.
    adj_c <- if (inherits(adj, "dgCMatrix")) {
      adj
    } else {
      methods::as(adj, "dgCMatrix")
    }

    largest_comp <- cpp_largest_cc_dgC(adj_c)

  } else {
    # Dense Matrix or base matrix path
    largest_comp <- cpp_largest_cc_dense(as.matrix(adj))
  }

  # Preserve the original subsetting behavior/class as much as possible
  largest_submatrix <- adj[largest_comp, largest_comp, drop = FALSE]

  rownames(largest_submatrix) <- largest_comp
  colnames(largest_submatrix) <- largest_comp

  list(
    nodes = largest_comp,
    submatrix = largest_submatrix
  )
}

out.cc <- find_largest_cc(A)
node.id <- twitch_nodes[out.cc$nodes]
twitch.adj.lcc <- A[sort(out.cc$nodes), sort(out.cc$nodes)]

twitch.el.lcc <- as.data.frame(summary(twitch.adj.lcc), what = "edges")

twitch.lang <- twitch_sub |>
  filter(node %in% as.integer(node.id)) |>
  pull(language)

twitch.comm <- factor(twitch.lang) |>
  as.integer()

twitch.lang.family <- case_when(
  twitch.lang %in% c('DA', 'DE', 'FI', 'NL', 'NO', 'SV') ~ 'Germanic',
  twitch.lang %in% c('FR', 'IT') ~ 'French',
  twitch.lang %in% c('ES', 'PT') ~ 'Spanish',
  twitch.lang %in% c('JA', 'KO', 'TH', 'TR', 'ZH') ~ 'Asian',
  .default = 'E.EU'
)

twitch.comm.family <- factor(twitch.lang.family) |>
  as.integer()
