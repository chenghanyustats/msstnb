
# plot_multiscale_regions.R
# Utilities for plotting multiscale regions using the `tree` object returned by make_tree()
# and simulated outputs saved by the data generation pipeline.
#
# Expected inputs
# - tree: list returned by make_tree(), containing at least:
#     tree$L
#     tree$coords[[l]] with columns row and col for each level l
#     tree$children[[l]][[j]] giving children ids of parent j at level l
# - y_levels (optional): list of matrices, y_levels[[l]] is TT by n_l, counts or any scalar values
#
# This script avoids facet_wrap free scales with fixed aspect ratio conflicts by:
# - using facet_wrap without free scales when coord_equal is requested, or
# - using per level plots combined with patchwork when both are requested.

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required. Please install it with install.packages('ggplot2').")
  }
})

# ---------- Small helpers ----------

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

is_scalar_int <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x) && abs(x - round(x)) < 1e-8

# Extract tree from an object that may wrap it in various fields
extract_tree <- function(obj) {
  if (is.list(obj) && !is.null(obj$tree)) return(obj$tree)
  if (is.list(obj) && !is.null(obj$data) && is.list(obj$data) && !is.null(obj$data$tree)) return(obj$data$tree)
  if (is.list(obj) && !is.null(obj$truth) && is.list(obj$truth) && !is.null(obj$truth$tree)) return(obj$truth$tree)
  stop("Could not find `tree` in the provided object. Expected obj$tree or obj$data$tree.")
}

# Extract y_levels from an object that may store it in various fields
extract_y_levels <- function(obj) {
  if (is.list(obj) && !is.null(obj$y_levels)) return(obj$y_levels)
  if (is.list(obj) && !is.null(obj$data) && is.list(obj$data) && !is.null(obj$data$y_levels)) return(obj$data$y_levels)
  if (is.list(obj) && !is.null(obj$y) && is.list(obj$y) && !is.null(obj$y$y_levels)) return(obj$y$y_levels)
  return(NULL)
}

# Read a simulated replicate RDS (data.rds or rep####.rds)
read_sim_rds <- function(path) {
  assert_true(file.exists(path), paste0("File not found: ", path))
  readRDS(path)
}

# Validate tree coords format
validate_tree_coords <- function(tree) {
  assert_true(!is.null(tree$L) && is_scalar_int(tree$L) && tree$L >= 1, "tree$L must be a positive integer.")
  assert_true(!is.null(tree$coords) && length(tree$coords) >= tree$L, "tree$coords must exist for each level.")
  for (l in seq_len(tree$L)) {
    cd <- tree$coords[[l]]
    assert_true(is.data.frame(cd) || is.matrix(cd), paste0("tree$coords[[", l, "]] must be a data frame or matrix."))
    cdn <- colnames(cd)
    assert_true(!is.null(cdn) && all(c("row","col") %in% cdn), paste0("tree$coords[[", l, "]] must have columns row and col."))
  }
  invisible(TRUE)
}

# ---------- Data frame builders ----------

build_level_df <- function(tree, level, values = NULL, t0 = 1L) {
  validate_tree_coords(tree)
  assert_true(level >= 1 && level <= tree$L, "level must be between 1 and tree$L.")

  cd <- as.data.frame(tree$coords[[level]])
  n_l <- nrow(cd)
  out <- data.frame(
    level = factor(level),
    id = seq_len(n_l),
    row = cd$row,
    col = cd$col
  )

  if (!is.null(values)) {
    if (is.list(values)) {
      # values is y_levels
      assert_true(length(values) >= level, "values list does not contain requested level.")
      mat <- values[[level]]
      assert_true(is.matrix(mat) || is.data.frame(mat), "values[[level]] must be a matrix or data frame.")
      mat <- as.matrix(mat)
      assert_true(ncol(mat) == n_l, "values[[level]] has incompatible number of columns for this level.")
      assert_true(t0 >= 1 && t0 <= nrow(mat), "t0 out of range for values[[level]].")
      out$value <- as.numeric(mat[t0, ])
    } else {
      # values is a numeric vector for this level
      assert_true(length(values) == n_l, "values vector must have length equal to number of regions at this level.")
      out$value <- as.numeric(values)
    }
  }

  out
}

build_multiscale_df <- function(tree, y_levels = NULL, t0 = 1L) {
  validate_tree_coords(tree)
  out <- vector("list", tree$L)
  for (l in seq_len(tree$L)) {
    out[[l]] <- build_level_df(tree, level = l, values = y_levels, t0 = t0)
  }
  do.call(rbind, out)
}

# ---------- Plotters ----------

# Plot numbering only (no values needed)
plot_tree_numbering <- function(tree,
                                free_scales = TRUE,
                                square_tiles = TRUE,
                                label_size = 4) {
  df <- build_multiscale_df(tree, y_levels = NULL, t0 = 1L)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row)) +
    ggplot2::geom_tile(fill = "white", color = "grey70", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = id), size = label_size) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "col", y = "row", title = "Region numbering at each level") +
    ggplot2::theme_minimal()

  # Handling the ggplot2 constraint: free scales cannot be combined with fixed ratio
  if (square_tiles && free_scales) {
    # build separate plots per level, then combine if patchwork exists
    plots <- lapply(seq_len(tree$L), function(l) {
      dfl <- subset(df, level == factor(l))
      ggplot2::ggplot(dfl, ggplot2::aes(x = col, y = row)) +
        ggplot2::geom_tile(fill = "white", color = "grey70", linewidth = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = id), size = label_size) +
        ggplot2::scale_y_reverse() +
        ggplot2::coord_equal() +
        ggplot2::labs(title = paste("Level", l), x = "col", y = "row") +
        ggplot2::theme_minimal()
    })

    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(Reduce(`+`, plots) + patchwork::plot_layout(guides = "collect"))
    } else {
      message("patchwork not installed. Returning a list of ggplot objects, one per level.")
      return(plots)
    }
  }

  if (free_scales) {
    p <- p + ggplot2::facet_wrap(~ level, scales = "free") + ggplot2::coord_cartesian()
  } else {
    p <- p + ggplot2::facet_wrap(~ level)
  }

  if (square_tiles && !free_scales) p <- p + ggplot2::coord_equal()
  p
}

# Plot values by level (requires y_levels or a list of level matrices)
plot_multiscale_facets <- function(tree,
                                   y_levels,
                                   t0 = 1L,
                                   free_scales = TRUE,
                                   square_tiles = TRUE,
                                   show_id = TRUE,
                                   id_size = 3) {
  assert_true(!is.null(y_levels), "y_levels must be provided.")
  df <- build_multiscale_df(tree, y_levels = y_levels, t0 = t0)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_tile(color = "grey70", linewidth = 0.3) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "col", y = "row", fill = "value",
                  title = paste("Multiscale values by level at time", t0)) +
    ggplot2::theme_minimal()

  if (show_id) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = id), size = id_size)
  }

  # Same free scales vs fixed ratio handling as numbering plot
  if (square_tiles && free_scales) {
    plots <- lapply(seq_len(tree$L), function(l) {
      dfl <- subset(df, level == factor(l))
      pl <- ggplot2::ggplot(dfl, ggplot2::aes(x = col, y = row, fill = value)) +
        ggplot2::geom_tile(color = "grey70", linewidth = 0.3) +
        ggplot2::scale_y_reverse() +
        ggplot2::coord_equal() +
        ggplot2::labs(title = paste("Level", l, "time", t0), x = "col", y = "row", fill = "value") +
        ggplot2::theme_minimal()
      if (show_id) pl <- pl + ggplot2::geom_text(ggplot2::aes(label = id), size = id_size)
      pl
    })

    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(Reduce(`+`, plots) + patchwork::plot_layout(guides = "collect"))
    } else {
      message("patchwork not installed. Returning a list of ggplot objects, one per level.")
      return(plots)
    }
  }

  if (free_scales) {
    p <- p + ggplot2::facet_wrap(~ level, scales = "free") + ggplot2::coord_cartesian()
  } else {
    p <- p + ggplot2::facet_wrap(~ level)
  }

  if (square_tiles && !free_scales) p <- p + ggplot2::coord_equal()
  p
}

# Nested plot for L = 2 only
# Draws children rectangles inside parents using column major placement within each parent.
plot_nested_L2 <- function(tree, y_levels, t0 = 1L, show_parent_border = TRUE) {
  assert_true(tree$L == 2, "plot_nested_L2 requires tree$L == 2.")
  validate_tree_coords(tree)
  assert_true(!is.null(tree$children) && length(tree$children) >= 1, "tree$children[[1]] must exist.")
  assert_true(!is.null(y_levels) && length(y_levels) >= 2, "y_levels must include levels 1 and 2.")

  cd1 <- as.data.frame(tree$coords[[1]])
  n1 <- nrow(cd1)
  nrow1 <- max(cd1$row)
  ncol1 <- max(cd1$col)

  parent <- data.frame(
    id = seq_len(n1),
    row = cd1$row,
    col = cd1$col
  )
  parent$xmin <- (parent$col - 1) / ncol1
  parent$xmax <- parent$col / ncol1
  parent$ymin <- (parent$row - 1) / nrow1
  parent$ymax <- parent$row / nrow1

  children1 <- tree$children[[1]]
  n2 <- max(unlist(children1))
  mat2 <- as.matrix(y_levels[[2]])
  assert_true(t0 >= 1 && t0 <= nrow(mat2), "t0 out of range for y_levels[[2]].")
  child_vals <- as.numeric(mat2[t0, ])
  assert_true(length(child_vals) == n2, "y_levels[[2]] has incompatible number of columns.")

  rows <- list()
  for (j in seq_len(n1)) {
    ch <- children1[[j]]
    k <- length(ch)
    if (k == 0) next

    # Build a local grid that is filled column major by default in R.
    nr <- floor(sqrt(k))
    nc <- ceiling(k / nr)
    loc <- matrix(NA_integer_, nrow = nr, ncol = nc)
    loc[seq_len(k)] <- seq_len(k)

    pos <- which(!is.na(loc), arr.ind = TRUE)
    local_id <- loc[!is.na(loc)]

    pb <- parent[j, ]
    w <- (pb$xmax - pb$xmin) / nc
    h <- (pb$ymax - pb$ymin) / nr

    for (m in seq_along(local_id)) {
      rr <- pos[m, "row"]
      cc <- pos[m, "col"]
      cid <- ch[local_id[m]]
      rows[[length(rows) + 1L]] <- data.frame(
        id = cid,
        parent = j,
        xmin = pb$xmin + (cc - 1) * w,
        xmax = pb$xmin + cc * w,
        ymin = pb$ymin + (rr - 1) * h,
        ymax = pb$ymin + rr * h,
        value = child_vals[cid]
      )
    }
  }
  child <- do.call(rbind, rows)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = child,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = value),
      color = "grey70",
      linewidth = 0.2
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_equal() +
    ggplot2::labs(title = paste("Nested level 2 inside level 1 at time", t0), fill = "value") +
    ggplot2::theme_void()

  if (show_parent_border) {
    p <- p + ggplot2::geom_rect(
      data = parent,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA,
      color = "black",
      linewidth = 0.8
    )
  }
  p
}

# ---------- Minimal example runner ----------
# Run from terminal or RStudio:
# Rscript plot_multiscale_regions.R path/to/data.rds 6
#
# The script will attempt to read tree and y_levels from the RDS.
# If y_levels is missing, it will plot only numbering.

run_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) return(invisible(NULL))

  path <- args[[1]]
  t0 <- if (length(args) >= 2) as.integer(args[[2]]) else 1L

  obj <- read_sim_rds(path)
  tree <- extract_tree(obj)
  validate_tree_coords(tree)

  y_levels <- extract_y_levels(obj)

  print(plot_tree_numbering(tree, free_scales = TRUE, square_tiles = TRUE, label_size = 4))

  if (!is.null(y_levels)) {
    print(plot_multiscale_facets(tree, y_levels, t0 = t0, free_scales = TRUE, square_tiles = TRUE, show_id = TRUE))
    if (tree$L == 2) print(plot_nested_L2(tree, y_levels, t0 = t0, show_parent_border = TRUE))
  } else {
    message("No y_levels found in the RDS. Only numbering plot was produced.")
  }
  invisible(NULL)
}

if (!interactive()) {
  run_cli()
}
