# tree_utils.R
# Helper utilities for the multiscale tree (splits) used in MSSTNB.

flatten_internal_nodes <- function(tree) {
    # Returns a data.frame with columns l, j, K for internal nodes with K > 1 children.
    L <- tree$L
    rows <- list()
    idx <- 1L

    for (l in seq_len(L - 1L)) {
        n_parent <- tree$n[l]
        for (j in seq_len(n_parent)) {
            ch <- tree$children[[l]][[j]]
            K <- length(ch)
            if (K > 1L) {
                rows[[idx]] <- data.frame(
                    l = as.integer(l),
                    j = as.integer(j),
                    K = as.integer(K),
                    stringsAsFactors = FALSE
                )
                idx <- idx + 1L
            }
        }
    }

    if (length(rows) == 0L) {
        return(data.frame(l = integer(), j = integer(), K = integer()))
    }

    do.call(rbind, rows)
}

make_delta_map <- function(tree) {
    tab <- flatten_internal_nodes(tree)
    list(
        tab = tab,
        l = as.integer(tab$l),
        j = as.integer(tab$j)
    )
}

delta_to_vec <- function(delta, delta_map) {
    # delta: nested list structure
    # delta_map: output of make_delta_map()
    l <- delta_map$l
    j <- delta_map$j
    out <- numeric(length(l))

    for (i in seq_along(l)) {
        out[i] <- delta[[l[i]]][[j[i]]]
    }

    out
}

vec_to_delta <- function(v, delta_template, delta_map) {
    l <- delta_map$l
    j <- delta_map$j

    assert_true(length(v) == length(l), "vec_to_delta: length mismatch.")

    out <- delta_template
    for (i in seq_along(l)) {
        out[[l[i]]][[j[i]]] <- v[i]
    }

    out
}
