# make_graph.R
# Graph construction and ICAR intrinsic precision H

suppressWarnings({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package Matrix is required.", call. = FALSE)
  }
})

make_grid_edges <- function(n1) {
  # Create 4 neighbor grid edges for n1 nodes arranged in near square grid
  nr <- floor(sqrt(n1))
  nc <- ceiling(n1 / nr)
  idx <- matrix(NA_integer_, nrow = nr, ncol = nc)
  idx[seq_len(n1)] <- seq_len(n1)

  edges <- list()
  add_edge <- function(a, b) {
    if (!is.na(a) && !is.na(b)) edges[[length(edges) + 1]] <<- c(a, b)
  }

  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      i <- idx[r, c]
      if (is.na(i)) next
      if (r < nr) add_edge(i, idx[r + 1, c])
      if (c < nc) add_edge(i, idx[r, c + 1])
    }
  }

  if (length(edges) == 0) return(matrix(integer(0), ncol = 2))
  do.call(rbind, edges)
}

make_ring_edges <- function(n1) {
  cbind(seq_len(n1), c(seq_len(n1 - 1) + 1, 1))
}

make_graph <- function(n1, graph_type = c("grid", "ring")) {
  graph_type <- match.arg(graph_type)
  edges <- switch(
    graph_type,
    grid = make_grid_edges(n1),
    ring = make_ring_edges(n1)
  )

  # Make symmetric edges
  if (nrow(edges) > 0) {
    edges_sym <- rbind(edges, edges[, c(2, 1), drop = FALSE])
  } else {
    edges_sym <- matrix(integer(0), ncol = 2)
  }

  # Sparse adjacency W
  W <- Matrix::sparseMatrix(
    i = edges_sym[, 1],
    j = edges_sym[, 2],
    x = 1,
    dims = c(n1, n1)
  )
  diag(W) <- 0

  # Intrinsic precision H = D - W
  d <- Matrix::rowSums(W)
  H <- Matrix::Diagonal(n1, x = d) - W

  list(n1 = n1, graph_type = graph_type, edges = edges, W = W, H = H)
}

