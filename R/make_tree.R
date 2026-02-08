# make_tree.R
# Build balanced or mixed resolution rooted tree for multiscale counts

make_tree <- function(n1, L, branching, refine_prob = 1.0,
                      numbering = c("colmajor", "colmajor_global", "sequential"),
                      seed = NULL) {
  # n1: number of nodes at level 1
  # L: number of levels
  # branching: integer vector of length L - 1; branching[l] is children per refined node at level l
  # refine_prob: probability a node is refined at each nonterminal level; used for mixed resolution trees
  # numbering: node labeling convention for each level.
  #   sequential: children at each level are numbered consecutively in the order parents are visited.
  #   colmajor: children are numbered consecutively within each parent (so parent j has (j-1)k+1:jk)
  #             and are assigned implied grid coordinates using column-major fill within each parent.
  #             This matches the user's desired convention for debugging and plotting.
  #   colmajor_global: nodes at each level are numbered by global column-major order of their implied
  #                   grid coordinates (this may interleave children from different parents).
  set_seed(seed)

  numbering <- match.arg(numbering)

  assert_true(L >= 1, "L must be at least 1.")
  assert_true(n1 >= 1, "n1 must be at least 1.")
  if (L == 1) {
    nr1 <- floor(sqrt(n1))
    nc1 <- ceiling(n1 / nr1)
    coords1 <- data.frame(
      row = ((seq_len(n1) - 1L) %% nr1) + 1L,
      col = ((seq_len(n1) - 1L) %/% nr1) + 1L
    )
    return(list(
      L = 1L,
      n = c(n1),
      children = list(),
      parents = list(),
      coords = list(coords1),
      dims = list(c(nr = as.integer(nr1), nc = as.integer(nc1))),
      numbering = numbering
    ))
  }

  assert_true(length(branching) == (L - 1), "branching must have length L - 1.")
  assert_true(all(branching >= 1), "branching entries must be positive.")
  assert_true(refine_prob >= 0 && refine_prob <= 1, "refine_prob must be in [0,1].")

  n_by_level <- integer(L)
  n_by_level[1] <- n1
  children <- vector("list", L - 1)  # children[[l]][[j]] gives child indices at level l+1
  parents <- vector("list", L - 1)   # parents[[l]] is integer vector of length n_{l+1} giving parent at level l

  # Optional coordinates for debugging and plotting.
  # coords[[l]] is a data.frame with columns row and col, in node id order at level l.
  coords <- vector("list", L)
  dims <- vector("list", L)

  # Level 1 coordinates match make_graph.R: column major fill.
  nr1 <- floor(sqrt(n1))
  nc1 <- ceiling(n1 / nr1)
  coords[[1]] <- data.frame(
    row = ((seq_len(n1) - 1L) %% nr1) + 1L,
    col = ((seq_len(n1) - 1L) %/% nr1) + 1L
  )
  dims[[1]] <- c(nr = as.integer(nr1), nc = as.integer(nc1))

  if (numbering == "sequential") {
    for (l in seq_len(L - 1)) {
      n_parent <- n_by_level[l]
      children_l <- vector("list", n_parent)

      refine_flag <- rep(TRUE, n_parent)
      if (refine_prob < 1) refine_flag <- rbinom(n_parent, size = 1, prob = refine_prob) == 1
      if (all(!refine_flag)) refine_flag[sample.int(n_parent, 1)] <- TRUE

      next_id <- 0L
      parent_of_child <- integer(0)
      for (j in seq_len(n_parent)) {
        if (!refine_flag[j]) {
          children_l[[j]] <- integer(0)
          next
        }
        k <- branching[l]
        child_ids <- (next_id + 1L):(next_id + k)
        children_l[[j]] <- child_ids
        parent_of_child <- c(parent_of_child, rep.int(j, k))
        next_id <- next_id + k
      }

      n_child <- next_id
      n_by_level[l + 1] <- n_child
      children[[l]] <- children_l
      parents[[l]] <- parent_of_child

      # Coordinates exist but numbering is not colmajor beyond level 1.
      coords[[l + 1]] <- NULL
      dims[[l + 1]] <- NULL
    }

    return(list(
      L = as.integer(L),
      n = n_by_level,
      children = children,
      parents = parents,
      coords = coords,
      dims = dims,
      numbering = numbering
    ))
  }

  # numbering == "colmajor" or "colmajor_global"
  for (l in seq_len(L - 1)) {
    n_parent <- n_by_level[l]
    children_l <- vector("list", n_parent)

    refine_flag <- rep(TRUE, n_parent)
    if (refine_prob < 1) refine_flag <- rbinom(n_parent, size = 1, prob = refine_prob) == 1
    if (all(!refine_flag)) refine_flag[sample.int(n_parent, 1)] <- TRUE

    k <- as.integer(branching[l])
    nr_sub <- floor(sqrt(k))
    nc_sub <- ceiling(k / nr_sub)
    assert_true(nr_sub >= 1 && nc_sub >= 1, "Invalid child grid dimensions.")

    if (numbering == "colmajor") {
      # Desired behavior:
      #   - parent j gets consecutive child ids
      #   - within each parent, local indices 1:k correspond to column-major positions
      # This guarantees, for example, branching=4 yields local layout 1,2 in first column; 3,4 in second.

      max_child <- n_parent * k
      parent_of_child <- integer(max_child)
      row_child <- integer(max_child)
      col_child <- integer(max_child)

      next_id <- 0L
      for (j in seq_len(n_parent)) {
        if (!refine_flag[j]) {
          children_l[[j]] <- integer(0)
          next
        }

        prow <- coords[[l]]$row[j]
        pcol <- coords[[l]]$col[j]

        child_ids <- (next_id + 1L):(next_id + k)
        children_l[[j]] <- child_ids

        # Local ids 1:k correspond to column-major fill of an nr_sub by nc_sub subgrid.
        local_id <- seq_len(k)
        lrow <- ((local_id - 1L) %% nr_sub) + 1L
        lcol <- ((local_id - 1L) %/% nr_sub) + 1L

        crow <- (prow - 1L) * nr_sub + lrow
        ccol <- (pcol - 1L) * nc_sub + lcol

        parent_of_child[child_ids] <- j
        row_child[child_ids] <- crow
        col_child[child_ids] <- ccol

        next_id <- next_id + k
      }

      n_child <- next_id
      n_by_level[l + 1] <- n_child

      parent_of_child <- parent_of_child[seq_len(n_child)]
      row_child <- row_child[seq_len(n_child)]
      col_child <- col_child[seq_len(n_child)]

      coords[[l + 1]] <- data.frame(row = row_child, col = col_child)
      dims[[l + 1]] <- c(nr = as.integer(max(coords[[l + 1]]$row)),
                         nc = as.integer(max(coords[[l + 1]]$col)))

      children[[l]] <- children_l
      parents[[l]] <- parent_of_child
    } else {
      # numbering == "colmajor_global"
      # Build child records with implied coordinates, then renumber globally by (col,row).
      rec_parent <- integer(0)
      rec_local <- integer(0)
      rec_row <- integer(0)
      rec_col <- integer(0)

      for (j in seq_len(n_parent)) {
        if (!refine_flag[j]) {
          children_l[[j]] <- integer(0)
          next
        }

        prow <- coords[[l]]$row[j]
        pcol <- coords[[l]]$col[j]

        local_id <- seq_len(k)
        lrow <- ((local_id - 1L) %% nr_sub) + 1L
        lcol <- ((local_id - 1L) %/% nr_sub) + 1L

        crow <- (prow - 1L) * nr_sub + lrow
        ccol <- (pcol - 1L) * nc_sub + lcol

        rec_parent <- c(rec_parent, rep.int(j, k))
        rec_local <- c(rec_local, local_id)
        rec_row <- c(rec_row, crow)
        rec_col <- c(rec_col, ccol)
      }

      n_child <- length(rec_parent)
      n_by_level[l + 1] <- n_child

      if (n_child == 0) {
        coords[[l + 1]] <- data.frame(row = integer(0), col = integer(0))
        dims[[l + 1]] <- c(nr = 0L, nc = 0L)
        children[[l]] <- children_l
        parents[[l]] <- integer(0)
        next
      }

      ord <- order(rec_col, rec_row)
      child_id <- integer(n_child)
      child_id[ord] <- seq_len(n_child)

      parent_of_child <- integer(n_child)
      parent_of_child[child_id] <- rec_parent

      for (j in seq_len(n_parent)) {
        if (!refine_flag[j]) {
          children_l[[j]] <- integer(0)
          next
        }
        idxj <- which(rec_parent == j)
        idxj <- idxj[order(rec_local[idxj])]
        children_l[[j]] <- child_id[idxj]
      }

      coords[[l + 1]] <- data.frame(row = rec_row, col = rec_col)
      coords[[l + 1]] <- coords[[l + 1]][order(child_id), , drop = FALSE]
      dims[[l + 1]] <- c(nr = as.integer(max(coords[[l + 1]]$row)),
                         nc = as.integer(max(coords[[l + 1]]$col)))

      children[[l]] <- children_l
      parents[[l]] <- parent_of_child
    }
  }

  list(
    L = as.integer(L),
    n = n_by_level,
    children = children,
    parents = parents,
    coords = coords,
    dims = dims,
    numbering = numbering
  )
}

is_leaf_node <- function(tree, l, j) {
  if (l >= tree$L) return(TRUE)
  length(tree$children[[l]][[j]]) == 0
}

