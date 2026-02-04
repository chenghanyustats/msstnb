\
# make_tree.R
# Build balanced or mixed resolution rooted tree for multiscale counts

make_tree <- function(n1, L, branching, refine_prob = 1.0, seed = NULL) {
  # n1: number of nodes at level 1
  # L: number of levels
  # branching: integer vector of length L - 1; branching[l] is children per refined node at level l
  # refine_prob: probability a node is refined at each nonterminal level; used for mixed resolution trees
  set_seed(seed)

  assert_true(L >= 1, "L must be at least 1.")
  assert_true(n1 >= 1, "n1 must be at least 1.")
  if (L == 1) {
    return(list(L = 1L, n = c(n1), children = list(), parents = list()))
  }

  assert_true(length(branching) == (L - 1), "branching must have length L - 1.")
  assert_true(all(branching >= 1), "branching entries must be positive.")
  assert_true(refine_prob >= 0 && refine_prob <= 1, "refine_prob must be in [0,1].")

  n_by_level <- integer(L)
  n_by_level[1] <- n1
  children <- vector("list", L - 1)  # children[[l]][[j]] gives child indices at level l+1
  parents <- vector("list", L - 1)   # parents[[l]] is integer vector of length n_{l+1} giving parent at level l

  for (l in seq_len(L - 1)) {
    n_parent <- n_by_level[l]
    children_l <- vector("list", n_parent)

    refine_flag <- rep(TRUE, n_parent)
    if (refine_prob < 1) refine_flag <- rbinom(n_parent, size = 1, prob = refine_prob) == 1

    # Ensure at least one node is refined unless l is last splitting level
    if (all(!refine_flag)) refine_flag[sample.int(n_parent, 1)] <- TRUE

    # Assign children consecutively
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
  }

  list(
    L = as.integer(L),
    n = n_by_level,
    children = children,
    parents = parents
  )
}

is_leaf_node <- function(tree, l, j) {
  if (l >= tree$L) return(TRUE)
  length(tree$children[[l]][[j]]) == 0
}

