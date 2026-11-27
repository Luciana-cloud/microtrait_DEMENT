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

tree_file <- "C:/luciana_datos/project_1_Microtrait_DEMENT/phylogeny/isolate_gtdbtk.unrooted.tree"
iso_file  <- "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/total_genes.guild.940_ISO.csv"
cue_file  <- "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Input_Data/SOIL_ISOLATES/dement_isolates_CUE.csv"
out_file  <- "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Output_Data/Figures/isolate_phylogeny.pdf"
meta_file_1  <- "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Input_Data/SOIL_ISOLATES/metadata_2.tsv"
meta_file_2  <- "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Input_Data/SOIL_ISOLATES/metadata_1.tsv"

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
    meta_data %>% select(id, `Specific Ecosystem`, `Geographic Location`),   
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
p = ggtree(tree_pruned, layout = "circular", size = 0.1)

# Ring 1 — genome size
p1 = p +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df,
    geom    = geom_tile,
    mapping = aes(y = id, fill = genome.size),
    width   = 0.1,
    offset  = 0.05
  ) +
  scale_fill_viridis_c(option = "D", name = "Genome size (bp)", na.value = "grey70")

# Ring 2 — CUE
p2 = p1 +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df,
    geom    = geom_tile,
    mapping = aes(y = id, fill = CUE),
    width   = 0.1,
    offset  = 0.05
  ) +
  scale_fill_viridis_c(option = "C", name = "CUE", na.value = "grey70")   

# Ring 3 — functional guild
p3 = p2 +
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

# --- ecosystem classes (check the grouping before using it) ---
x = trimws(as.character(meta_df_2$`Specific Ecosystem`))
meta_df_2$eco_group = dplyr::case_when(
  is.na(x) | tolower(x) %in% c("", "unclassified")          ~ "Unclassified",
  x %in% c("Agricultural land", "Agricultural soil", "Farm") ~ "Agricultural",
  x %in% c("Contaminated", "Oil-contaminated", "Mine drainage") ~ "Contaminated",
  x == "Desert"                                              ~ "Desert",
  x %in% c("Bulk soil", "Mineral soil")                      ~ "Bulk / mineral soil",
  x == "Permafrost"                                          ~ "Permafrost",
  x %in% c("Forest", "Forest soil", "Grasslands")            ~ "Forest / grassland",
  TRUE                                                       ~ "Other"
)

# check: counts per class and what ended up in "Other"
print(table(meta_df_2$eco_group))
print(table(x[meta_df_2$eco_group == "Other"]))

# legend with n
eco_levels = c("Agricultural", "Contaminated", "Desert", "Bulk / mineral soil",
               "Permafrost", "Forest / grassland", "Other", "Unclassified")
meta_df_2$eco_group = factor(meta_df_2$eco_group, levels = eco_levels)
eco_n   = table(meta_df_2$eco_group)
eco_lab = paste0(eco_levels, " (n = ", eco_n[eco_levels], ")")

eco_cols = c("Agricultural"        = "#E69F00",
             "Contaminated"        = "#D55E00",
             "Desert"              = "#F0E442",
             "Bulk / mineral soil" = "#009E73",
             "Permafrost"          = "#56B4E9",
             "Forest / grassland"  = "#0072B2",
             "Other"               = "#CC79A7",
             "Unclassified"        = "grey85")

# Base tree
p.2 <- ggtree(tree_pruned.2, layout = "circular", size = 0.1)

# Ring 1: source ecosystem
p1.2 <- p.2 +
  new_scale_fill() +
  geom_fruit(
    data    = meta_df_2,
    geom    = geom_tile,
    mapping = aes(y = id, fill = eco_group),
    width   = 0.1,
    offset  = 0.05
  ) +
  scale_fill_manual(values = eco_cols, breaks = eco_levels, labels = eco_lab,
                    name = "Source ecosystem", drop = FALSE)

# Ring 2: functional guild (unchanged)
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
    guide = "none"
  ) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 9)
  )

ggsave(out_file, p3.2, width = 14, height = 14)
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

# Phylogeny for MAGs ----

library(dplyr)

dir_phy = "C:/luciana_datos/project_1_Microtrait_DEMENT/phylogeny/gtdb/gtdb"

cue_MAG = sim_10000 %>%
  transmute(id = as.character(MAG_id),
            CUE = ifelse(BGE <= 0.848, BGE, NA)) %>%
  distinct(id, .keep_all = TRUE)

MAG_trait = MAG_trait %>%
  select(-any_of(c("CUE", "CUE.x", "CUE.y"))) %>%
  left_join(cue_MAG, by = "id")

sum(!is.na(MAG_trait$CUE))                              # expect 14,075
length(unique(MAG_trait$guild[!is.na(MAG_trait$CUE)]))  # expect 1,139

write.csv(MAG_trait[, c("id", "guild", "genome.size", "CUE")],
          file.path(dir_phy, "MAG_trait_for_phylosig.csv"), row.names = FALSE)

library(ape)
library(phytools)

dir_phy  = "C:/luciana_datos/project_1_Microtrait_DEMENT/phylogeny/gtdb/gtdb"
tree_mag = read.tree(file.path(dir_phy, "img_gtdbtk.unrooted.tree"))

MAG_trait = read.csv(file.path(dir_phy, "MAG_trait_for_phylosig.csv"))
MAG_trait$id    = as.character(MAG_trait$id)
MAG_trait$guild = as.character(MAG_trait$guild)

mt = MAG_trait[MAG_trait$id %in% tree_mag$tip.label, ]
cat("MAGs in tree:", nrow(mt), "| with CUE:", sum(!is.na(mt$CUE)),
    "| guilds:", length(unique(mt$guild)), "\n")   # expect 21,077 | 11,200 | 1,256

# FG means from ALL MAGs of the guild (as in the FG tables)
fg_gs  = c(tapply(MAG_trait$genome.size, MAG_trait$guild, mean))
fg_cue = c(tapply(MAG_trait$CUE, MAG_trait$guild, mean, na.rm = TRUE))   # NaN if no CUE

tips_by_guild = split(mt$id, mt$guild)
gs_mag  = setNames(mt$genome.size, mt$id)
cue_mag = setNames(mt$CUE, mt$id)

# Blomberg's K for one draw
sig = function(ids, x) {
  tr = midpoint.root(keep.tip(tree_mag, ids))
  tr$edge.length[tr$edge.length == 0] = 1e-6
  x  = x[tr$tip.label]
  c(K = as.numeric(phylosig(tr, x, method = "K")))
}

# individual MAGs vs one MAG per guild carrying the guild mean
run_trait = function(trait_mag, trait_fg, guilds, R = 20) {
  pool = names(trait_mag)[!is.na(trait_mag)]
  n    = length(guilds)
  res  = lapply(1:R, function(i) {
    ids_b = sample(pool, n)
    ids_a = sapply(tips_by_guild[guilds], function(v) v[sample.int(length(v), 1)])
    rbind(before = sig(ids_b, trait_mag[ids_b]),
          after  = sig(ids_a, setNames(trait_fg[guilds], ids_a)))
  })
  list(n = n, res = res)
}

summarise_run = function(out) {
  K = sapply(out$res, function(m) m[, "K"])
  data.frame(level = c("MAGs", "FGs"), n = out$n,
             K_median = apply(K, 1, median),
             K_min = apply(K, 1, min), K_max = apply(K, 1, max))
}

g_size = names(tips_by_guild)
g_cue  = names(tips_by_guild)[!is.nan(fg_cue[names(tips_by_guild)])]

# test run
set.seed(1)
t0 = Sys.time()
out_test = run_trait(gs_mag, fg_gs, g_size, R = 3)
Sys.time() - t0
summarise_run(out_test)

# full run
set.seed(1)
out_gs  = run_trait(gs_mag,  fg_gs,  g_size, R = 20)
out_cue = run_trait(cue_mag, fg_cue, g_cue,  R = 20)

res_K = rbind(cbind(trait = "Genome size", summarise_run(out_gs)),
              cbind(trait = "CUE",         summarise_run(out_cue)))
tab = res_K
tab[, c("K_median", "K_min", "K_max")] =
  lapply(tab[, c("K_median", "K_min", "K_max")], formatC, format = "e", digits = 1)
names(tab) = c("Trait", "Level", "n tips", "K (median)", "K (min)", "K (max)")
write.csv(tab, "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Output_Data/Table_S_phylosig_K_MAGs.csv", row.names = FALSE)

