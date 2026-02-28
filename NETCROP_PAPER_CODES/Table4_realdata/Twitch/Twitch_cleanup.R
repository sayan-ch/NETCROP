setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table4_realdata", "Twitch")))

library(tidyverse)
library(Matrix)
# library(igraph)

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

find_largest_cc <- function(adj) {
  # if (!inherits(adj, "dMatrix")) stop("Input must be a dMatrix")
  n <- nrow(adj)
  
  # Convert to list of neighbors for fast lookup
  neighbors <- vector("list", n)
  for (i in 1:n) {
    neighbors[[i]] <- which(adj[i, ] != 0)
  }
  
  visited <- rep(FALSE, n)
  components <- list()
  
  for (v in 1:n) {
    if (!visited[v]) {
      # BFS
      queue <- v
      comp <- v
      visited[v] <- TRUE
      while (length(queue) > 0) {
        u <- queue[1]
        queue <- queue[-1]
        for (w in neighbors[[u]]) {
          if (!visited[w]) {
            visited[w] <- TRUE
            queue <- c(queue, w)
            comp <- c(comp, w)
          }
        }
      }
      components <- append(components, list(unique(comp)))
    }
  }
  
  # Find largest component
  sizes <- sapply(components, length)
  largest_idx <- which.max(sizes)
  largest_comp <- components[[largest_idx]]
  
  # Create submatrix
  largest_submatrix <- adj[largest_comp, largest_comp]
  rownames(largest_submatrix) <- largest_comp
  colnames(largest_submatrix) <- largest_comp
  
  return(list(
    nodes = largest_comp,
    submatrix = largest_submatrix
  ))
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
