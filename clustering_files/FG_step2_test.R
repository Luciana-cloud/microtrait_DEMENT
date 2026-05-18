#!/usr/bin/env Rscript

# =============================================================================
# STEP 1 TEST: Single Repetition Clustering Run
# =============================================================================
# One seed, 5000 MAGs, k_start = 160, k_max = 350
# If this completes successfully, the output is directly usable in Step 2
# =============================================================================

library(parallelDist)
library(vegan)
library(dplyr)

# =============================================================================
# SETTINGS
# =============================================================================
input_file      <- "/lustre/nobackup/WUR/ESG/chave013/clustering_YAS/Global_dataset_microtrait.csv"
output_dir      <- "/lustre/nobackup/WUR/ESG/chave013/clustering_YAS/Output_Data/Clustering"
subset_size     <- 5000
seed            <- 1
n_threads       <- 1
perm            <- 999
p_thresh        <- 0.05
k_start         <- 160
k_max           <- 350
gene_cols_start <- 3

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# CUSTOM PAIRWISE ADONIS FUNCTION (replaces pairwiseAdonis package)
# Uses vegan::adonis2 directly — identical logic, no extra dependency
# =============================================================================
pairwise_adonis <- function(dist_matrix, grouping, perm = 999, p_adjust = "BH") {

  groups    <- levels(factor(grouping))
  all_pairs <- combn(groups, 2, simplify = FALSE)
  results   <- list()

  for (pair in all_pairs) {
    idx           <- which(grouping %in% pair)
    sub_dist      <- as.dist(as.matrix(dist_matrix)[idx, idx])
    sub_grouping  <- factor(grouping[idx])

    res <- tryCatch({
      vegan::adonis2(sub_dist ~ sub_grouping, permutations = perm)
    }, error = function(e) NULL)

    if (!is.null(res)) {
      results[[paste(pair, collapse = "_vs_")]] <- data.frame(
        pair    = paste(pair, collapse = " vs "),
        F       = res$F[1],
        R2      = res$R2[1],
        p.value = res$`Pr(>F)`[1]
      )
    }
  }

  out         <- do.call(rbind, results)
  out$p.adj   <- p.adjust(out$p.value, method = p_adjust)
  return(out)
}

# =============================================================================
# AUTOMATED k SELECTION FUNCTION
# =============================================================================
find_max_k <- function(distance_matrix, hclust_obj,
                       k_start = 160, k_max = 350,
                       perm = 999, p_thresh = 0.05) {

  cat("Finding maximum k starting at k =", k_start, "...\n")
  k_accepted <- k_start - 1

  for (k in k_start:k_max) {

    t_k          <- Sys.time()
    assignments  <- cutree(hclust_obj, k = k)
    guild_factor <- factor(assignments)

    # Check all groups have at least 2 members
    group_sizes <- table(guild_factor)
    if (any(group_sizes < 2)) {
      cat("k =", k, "— some groups have < 2 members, stopping.\n")
      break
    }

    # Run pairwise adonis using our custom function
    adonis_result <- tryCatch({
      pairwise_adonis(distance_matrix, guild_factor,
                      perm = perm, p_adjust = "BH")
    }, error = function(e) {
      cat("k =", k, "— adonis error:", conditionMessage(e), "\n")
      return(NULL)
    })

    if (is.null(adonis_result)) break

    all_sig  <- all(adonis_result$p.adj < p_thresh)
    elapsed  <- round(difftime(Sys.time(), t_k, units = "mins"), 2)

    if (all_sig) {
      k_accepted <- k
      cat("k =", k, "— all pairs significant |", elapsed, "mins\n")
    } else {
      cat("k =", k, "— non-significant pair found | stopping |",
          elapsed, "mins\n")
      cat("Maximum valid k =", k_accepted, "\n\n")
      break
    }

    if (k == k_max) {
      cat("Reached k_max =", k_max,
          "— consider increasing. Accepted k =", k_accepted, "\n\n")
    }
  }

  return(k_accepted)
}

# =============================================================================
# SINGLE REPETITION RUN
# =============================================================================
cat("=================================================\n")
cat("TEST RUN | Seed:", seed, "| n_MAGs:", subset_size, "\n")
cat("=================================================\n")

cat("Loading data...\n")
global_dataset <- read.csv(input_file, dec = ".")
gene_cols      <- gene_cols_start:ncol(global_dataset)
genome_ids     <- global_dataset[, 1]
cat("Dataset dimensions:", nrow(global_dataset), "MAGs x",
    length(gene_cols), "genes\n\n")

set.seed(seed)
subset_idx  <- sample(nrow(global_dataset), size = subset_size)
subset_data <- global_dataset[subset_idx, gene_cols]

# Distance matrix
cat("Computing Jaccard distance matrix...\n")
t0       <- Sys.time()
dist_mat <- parDist(
  x       = as.matrix(subset_data),
  method  = "fJaccard",
  threads = n_threads
)
cat("Distance matrix done in",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "mins\n\n")

# Clustering
cat("Clustering with Ward linkage...\n")
t1         <- Sys.time()
hclust_obj <- hclust(dist_mat, method = "ward.D2")
cat("Clustering done in",
    round(difftime(Sys.time(), t1, units = "mins"), 1), "mins\n\n")

# Find maximum valid k
k_final <- find_max_k(
  distance_matrix = dist_mat,
  hclust_obj      = hclust_obj,
  k_start         = k_start,
  k_max           = k_max,
  perm            = perm,
  p_thresh        = p_thresh
)

cat("Final k selected:", k_final, "\n\n")

# =============================================================================
# EXTRACT CENTROIDS
# =============================================================================
cat("Computing centroids for", k_final, "guilds...\n")
assignments <- cutree(hclust_obj, k = k_final)
subset_mat  <- as.matrix(subset_data)

centroids <- matrix(0, nrow = k_final, ncol = ncol(subset_mat))
for (i in 1:k_final) {
  members <- which(assignments == i)
  if (length(members) == 1) {
    centroids[i, ] <- subset_mat[members, ]
  } else {
    centroids[i, ] <- colMeans(subset_mat[members, ])
  }
}
cat("Centroids computed.\n\n")

# =============================================================================
# SAVE OUTPUT
# =============================================================================
out_file <- file.path(output_dir, "FINAL_clustering_for_step2.rds")

saveRDS(list(
  k_final        = k_final,
  centroids      = centroids,
  assignments    = assignments,
  subset_idx     = subset_idx,
  genome_ids_sub = genome_ids[subset_idx],
  gene_cols      = gene_cols,
  seed           = seed,
  n_MAGs         = subset_size
), file = out_file)

cat("Saved:", out_file, "\n")
cat("This file feeds directly into FG_step2_assign_all_MAGs.R\n\n")

total_time <- round(difftime(Sys.time(), t0, units = "hours"), 2)
cat("Total run time:", total_time, "hours\n")
cat("Done.\n")
