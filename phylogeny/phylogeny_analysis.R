# Isolate Phylogenetic Tree Visualization

# Produces a circular phylogenetic tree with three color strips:
#   1. Genome size
#   2. CUE (Yield)
#   3. Functional guild
# And runs phylogenetic signal tests (Blomberg's K and Pagel's lambda)
# on genome size and CUE

# Call packages ----

library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(phytools)
library(dplyr)
library(ggplot2)

# Setting up data ----

tree_file <- "C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/phylogeny/isolate_gtdbtk.unrooted.tree"
iso_file  <- "C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/total_genes.guild.940_ISO.csv"
cue_file  <- "C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset/Input_Data/SOIL_ISOLATES/dement_isolates_CUE.csv"
out_file  <- "C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset/Output_Data/Figures/isolate_phylogeny.pdf"
meta_file_1  <- "C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset/Input_Data/SOIL_ISOLATES/metadata_2.tsv"
meta_file_2  <- "C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset/Input_Data/SOIL_ISOLATES/metadata_1.tsv"

# STEP 1: Load data

cat("Loading data...\n")
tree <- read.tree(tree_file)
iso  <- read.csv(iso_file)
cue  <- read.csv(cue_file)
meta_file_1  <- read_delim(meta_file_1)
meta_file_2  <- read_delim(meta_file_2)

shared_cols <- intersect(colnames(meta_file_1), colnames(meta_file_2))

meta_data <- bind_rows(
  meta_file_1 %>% select(all_of(shared_cols)),
  meta_file_2 %>% select(all_of(shared_cols))
)

# Standardize ID columns to character
tree$tip.label <- as.character(tree$tip.label)
iso$id         <- as.character(iso$id)
cue$id         <- as.character(cue$id)
meta_data$taxon_oid <- as.character(meta_data$taxon_oid)
colnames(meta_data)[1] = "id"

cat("Tree tips:", length(tree$tip.label), "\n")
cat("Isolates in metadata:", nrow(iso), "\n")
cat("Isolates in CUE file:", nrow(cue), "\n")
cat("Isolates in meta_data file:", nrow(meta_data), "\n")

# STEP 2: Join metadata ----

cat("Joining metadata...\n")
meta <- iso %>%
  select(id, genome.size, guild) %>%
  inner_join(
    cue %>% select(id, CUE),
    by = "id"
  ) %>%
  filter(id %in% tree$tip.label)
cat("Working dataset:", nrow(meta), "isolates\n")

cat("Joining additional metadata...\n")
meta.2 <- iso %>%
  select(id, genome.size, guild) %>%
  inner_join(
    cue %>% select(id, CUE),
    by = "id"
  ) %>%
  left_join(
    meta_data %>% select(id, `Specific Ecosystem`, `Geographic Location`),   # choose the columns you want
    by = "id"
  ) %>%
  filter(id %in% tree$tip.label)
cat("Working dataset:", nrow(meta), "isolates\n")

# Add id as explicit column for ggtree mapping
meta_df      <- meta %>%
  mutate(id = as.character(id)) %>%
  tibble::column_to_rownames("id")
meta_df$id   <- rownames(meta_df)

meta_df_2      <- meta.2 %>%
  mutate(id = as.character(id)) %>%
  tibble::column_to_rownames("id")
meta_df_2$id   <- rownames(meta_df_2)

# STEP 3: Root and prune tree ----

cat("Rooting and pruning tree...\n")
tree_rooted <- midpoint.root(tree)
tree_pruned <- keep.tip(tree_rooted, meta_df$id)
cat("Tips after pruning:", length(tree_pruned$tip.label), "\n")
cat("Tree is rooted:", is.rooted(tree_pruned), "\n\n")

cat("Rooting and pruning tree...\n")
tree_rooted.2 <- midpoint.root(tree)
tree_pruned.2 <- keep.tip(tree_rooted.2, meta_df_2$id)
cat("Tips after pruning:", length(tree_pruned.2$tip.label), "\n")
cat("Tree is rooted:", is.rooted(tree_pruned.2), "\n\n")

# STEP 4: Build phylogenetic tree figure ----

cat("Building tree figure...\n")

# Base tree
p <- ggtree(tree_pruned, layout = "circular", size = 0.1)

# Ring 1 — genome size
p1 <- p +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df,
    geom    = geom_tile,
    mapping = aes(y = id, fill = genome.size),
    width   = 0.1,
    offset  = 0.05
  ) +
  scale_fill_gradient(
    low  = "lightblue",
    high = "darkblue",
    name = "Genome size (bp)"
  )

# Ring 2 — CUE (Yield)
p2 <- p1 +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df,
    geom    = geom_tile,
    mapping = aes(y = id, fill = CUE),
    width   = 0.1,
    offset  = 0.05
  ) +
  scale_fill_gradient(
    low  = "lightyellow",
    high = "darkorange",
    name = "CUE (Yield)"
  )

# Ring 3 — functional guild
p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df,
    geom    = geom_tile,
    mapping = aes(y = id, fill = factor(guild)),
    width   = 0.1,
    offset  = 0.05
  ) +
  scale_fill_manual(
    values = colorRampPalette(
      c("red", "blue", "green", "purple", "orange",
        "pink", "brown", "grey", "cyan", "magenta")
    )(length(unique(meta_df$guild))),
    name  = "Functional guild",
    guide = "none"   # too many guilds for a readable legend
  ) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 9)
  )

# Save figure
cat("Saving figure to", out_file, "...\n")
ggsave(out_file, p3, width = 14, height = 14)
cat("Figure saved.\n\n")

# STEP 4.A: Build phylogenetic tree figure ----

cat("Building tree figure...\n")

# Base tree
p.2 <- ggtree(tree_pruned.2, layout = "circular", size = 0.1)

# Ring 1 — Specific Ecosystem
p1.2 <- p.2 +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df_2,
    geom    = geom_tile,
    mapping = aes(y = id, fill = `Specific Ecosystem`),
    width   = 0.1,
    offset  = 0.05
  )

# Ring 2 — functional guild
p3.2 <- p1.2 +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df_2,
    geom    = geom_tile,
    mapping = aes(y = id, fill = factor(guild)),
    width   = 0.1,
    offset  = 0.05
  ) +
  scale_fill_manual(
    values = colorRampPalette(
      c("red", "blue", "green", "purple", "orange",
        "pink", "brown", "grey", "cyan", "magenta")
    )(length(unique(meta_df$guild))),
    name  = "Functional guild",
    guide = "none"   # too many guilds for a readable legend
  ) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 9)
  )
# Save figure
cat("Saving figure to", out_file, "...\n")
ggsave(out_file, p3.2, width = 14, height = 14)
cat("Figure saved.\n\n")

# STEP 5: Phylogenetic signal tests ----

cat("=================================================\n")
cat("PHYLOGENETIC SIGNAL TESTS\n")
cat("=================================================\n")

# Prepare named vectors — required by phytools
genome_size_vec        <- setNames(meta_df$genome.size, rownames(meta_df))
cue_vec                <- setNames(meta_df$CUE,         rownames(meta_df))

# Keep only tips present in pruned tree
genome_size_vec <- genome_size_vec[tree_pruned$tip.label]
cue_vec         <- cue_vec[tree_pruned$tip.label]

# --- Blomberg's K ---
# Fix zero length branches
tree_pruned$edge.length[tree_pruned$edge.length == 0] <- 1e-6

# Blomberg's K — Genome size
cat("\nBlomberg's K — Genome size:\n")
K_genome <- phylosig(tree_pruned, genome_size_vec,
                     method = "K", test = TRUE, nsim = 999)
print(K_genome)

# Blomberg's K — CUE
cat("\nBlomberg's K — CUE (Yield):\n")
K_cue <- phylosig(tree_pruned, cue_vec,
                  method = "K", test = TRUE, nsim = 999)
print(K_cue)

# Pagel's lambda — Genome size
#cat("\nPagel's lambda — Genome size:\n")
#lambda_genome <- phylosig(tree_pruned, genome_size_vec,
#method = "lambda", test = TRUE)
#print(lambda_genome)

# Pagel's lambda — CUE
#cat("\nPagel's lambda — CUE (Yield):\n")
#lambda_cue <- phylosig(tree_pruned, cue_vec,
#method = "lambda", test = TRUE)
#print(lambda_cue)

# Save results
results <- data.frame(
  trait          = c("Genome size", "CUE (Yield)"),
  K              = c(K_genome$K,           K_cue$K),
  K_p_value      = c(K_genome$P,           K_cue$P)
  #  lambda         = c(lambda_genome$lambda,  lambda_cue$lambda),
  #  lambda_p_value = c(lambda_genome$P,       lambda_cue$P)
)

write.csv(results, "C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset/Output_Data/phylogenetic_signal_results.csv", row.names = FALSE)
print(results)

cat("\n=================================================\n")
cat("SUMMARY\n")
cat("=================================================\n")
print(results)
cat("\nSaved: phylogenetic_signal_results.csv\n")
cat("Done.\n")
