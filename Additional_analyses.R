
# Whether A and S are special compared to random traits ----

library(readxl)
library(ggplot2)

total_genes.full = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/total_genes.full.csv")
total_genes.guild.940_ISO = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/total_genes.guild.940_ISO.csv")
total_genes.guild.940_MAG = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/total_genes.guild.940_MAG.csv")
GH_rule = read_excel("Input_Data/microTrait_Rules/microtrait_GH.xlsx")
GH_rule = as.data.frame(rbind("id","guild","genome.size",GH_rule))
PR_rule = read_excel("Input_Data/microTrait_Rules/microtrait_proteins.xlsx")
PR_rule = as.data.frame(rbind("id","guild","genome.size",PR_rule))
transp_rule = read_excel("Input_Data/microTrait_Rules/microtrait_transporters.xlsx")
osmo_rule   = read_excel("Input_Data/microTrait_Rules/microtrait_osmolytes.xlsx")
osmo_rule   = as.data.frame(rbind("id","guild","genome.size",osmo_rule))
biofilm_rule = read_excel("Input_Data/microTrait_Rules/microtrait_biofilm.xlsx")
biofilm_rule = as.data.frame(rbind("id","guild","genome.size",biofilm_rule))
high.T_rule  = read_excel("Input_Data/microTrait_Rules/microtrait_high_T.xlsx")
high.T_rule  = as.data.frame(rbind("id","guild","genome.size",high.T_rule))
low.T_rule   = read_excel("Input_Data/microTrait_Rules/microtrait_low_T.xlsx")          # CHANGED: new
pH_rule      = read_excel("Input_Data/microTrait_Rules/microtrait_pH_stress.xlsx")
pH_rule      = as.data.frame(rbind("id","guild","genome.size",pH_rule))

MAG_gen_trait = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/MAG_gen_trait.csv")

# Add a column to distinguish MAGs and isolates ----
mag_ids = total_genes.guild.940_MAG$id
iso_ids = total_genes.guild.940_ISO$id

total_genes.full$dataset = ifelse(total_genes.full$id %in% mag_ids, "MAG",
                                  ifelse(total_genes.full$id %in% iso_ids, "ISO", NA))

## 0. A and S gene sets from the key files -----------------------------------
key = function(df) setdiff(unique(trimws(as.character(df[[1]]))),
                           c("id", "guild", "genome.size"))                           # CHANGED: drop header strings
A_all = unique(c(key(transp_rule), key(PR_rule), key(GH_rule)))
S_all = unique(c(key(high.T_rule), key(low.T_rule), key(pH_rule),                     # CHANGED: low-T added
                 key(osmo_rule), key(biofilm_rule)))
shared = intersect(A_all, S_all)                       # genes in both A and S (expect 6)

gene_cols      = setdiff(names(total_genes.full), c("X","id","guild","genome.size","dataset","guild_key"))
selected_genes = setdiff(names(total_genes.guild.940_MAG), c("X","id","genome.size","guild"))

A_genes = setdiff(intersect(A_all, selected_genes), shared)
S_genes = setdiff(intersect(S_all, selected_genes), shared)
nA = length(A_genes)
nS = length(S_genes)

cat("A: key", length(A_all), "-> used", nA, "| S: key", length(S_all), "-> used", nS,
    "| shared removed:", length(shared), "\n")
cat("Key genes not found in the matrix:", length(setdiff(c(A_all, S_all), gene_cols)), "\n")
cat("Selected genes in neither A nor S:", length(setdiff(selected_genes, c(A_genes, S_genes, shared))), "\n")

## 1. Find the extra guild (in the matrix but not in your table) ---------------
d_mag = subset(total_genes.full, dataset == "MAG")
extra = setdiff(unique(as.character(d_mag$guild)), as.character(MAG_gen_trait$guild))
extra
subset(d_mag, as.character(guild) %in% extra)[, c("id", "guild", "genome.size")]

## 2. Function: observed A-S correlation + random-gene-set null (weighted) -----
wcor = function(x, y, w) cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]            # CHANGED: new

run_control = function(ds, n_iter = 1000, trait_tab = NULL, drop_guilds = NULL) {
  
  d = subset(total_genes.full, dataset == ds)                 # ds = "MAG" or "ISO"
  if (!is.null(drop_guilds)) d = d[!(as.character(d$guild) %in% as.character(drop_guilds)), ]
  
  sums = rowsum(as.matrix(d[, gene_cols]), group = d$guild)
  n    = tabulate(match(d$guild, rownames(sums)))             # CHANGED: MAGs (or isolates) per FG = weight
  G    = sums / n                                             # guild x gene mean counts
  gs   = as.vector(tapply(d$genome.size, d$guild, mean)[rownames(G)])
  cat(ds, ": guilds =", nrow(G), "\n")
  
  stat = function(Ag, Sg) {
    a = rowSums(G[, Ag, drop = FALSE])
    s = rowSums(G[, Sg, drop = FALSE])
    c(r = wcor(a, s, n),                                                                 # CHANGED: weighted
      partial = wcor(resid(lm(a ~ gs, weights = n)), resid(lm(s ~ gs, weights = n)), n)) # CHANGED: weighted
  }
  observed = stat(A_genes, S_genes)
  
  # check against your existing guild table
  if (!is.null(trait_tab)) {
    m = match(rownames(G), as.character(trait_tab$guild))
    cat("Check vs your table: cor(A) =",
        round(cor(rowSums(G[, A_genes]), trait_tab$A_traits[m], use = "complete.obs"), 3),
        "| cor(S) =",
        round(cor(rowSums(G[, S_genes]), trait_tab$S_traits[m], use = "complete.obs"), 3), "\n")
  }
  
  null_for = function(pool) {
    set.seed(1)
    nl = t(replicate(n_iter, {
      g = sample(pool, nA + nS)
      stat(g[1:nA], g[(nA + 1):(nA + nS)])
    }))
    p = sapply(c("r", "partial"), function(k) (sum(nl[, k] >= observed[k]) + 1) / (n_iter + 1))
    list(quantiles = apply(nl, 2, quantile, c(.025, .5, .975)), p = p, draws = nl)
  }
  
  neither     = setdiff(gene_cols, c(A_genes, S_genes, shared))   # all genes in neither A nor S
  neither_sel = setdiff(selected_genes, shared)                   # random split of the selected genes
  
  list(observed = observed,
       null_all_genes      = null_for(neither),
       null_selected_genes = null_for(neither_sel))
}

## 3. Run (MAG without the extra guild; isolates unchanged) -------------------
res_MAG_w = run_control("MAG", trait_tab = MAG_gen_trait, drop_guilds = extra)
res_ISO_w = run_control("ISO")
res_MAG_w
res_ISO_w

# Additional Figure

res_MAG_w$null_all_genes$quantiles
res_MAG_w$null_selected_genes$quantiles
res_ISO_w$null_all_genes$quantiles
res_ISO_w$null_selected_genes$quantiles

get_df = function(res, ds) {
  rbind(
    data.frame(dataset = ds, null = "Random genes from outside A/S",
               value = res$null_all_genes$draws[, "partial"],
               real = unname(res$observed["partial"])),
    data.frame(dataset = ds, null = "Random splits of selected genes",
               value = res$null_selected_genes$draws[, "partial"],
               real = unname(res$observed["partial"])))
}

plot_df = rbind(get_df(res_MAG_w, "MAGs"), get_df(res_ISO_w, "Isolates"))
plot_df$null = factor(plot_df$null,
                      levels = c("Random genes from outside A/S", "Random splits of selected genes"))
plot_df$dataset = factor(plot_df$dataset, levels = c("MAGs", "Isolates"))

p_ctrl = ggplot(plot_df, aes(x = value)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  geom_vline(aes(xintercept = real), colour = "red", linewidth = 0.9) +
  facet_grid(dataset ~ null) +
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Correlation between two random gene sets (genome size removed)",
       y = "Number of random gene-set pairs") +
  theme_bw(base_size = 12) +
  theme(panel.spacing = unit(1, "lines"))

p_ctrl
ggsave("Output_Data/Figures/Figure_SX_random_control.pdf", p_ctrl, width = 8, height = 6)

# Whether A and S are special compared to RU traits ----

# MAGs

mt = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/data_paper/Data_Base/output/FUNCTIONAL ANNOTATION/mag_functional_traits_clean.csv")
total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")
MAG_gen_trait = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/MAG_gen_trait.csv")

wcor = function(x, y, w) cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
wpartial = function(x, y, z, w) wcor(resid(lm(x ~ z, weights = w)), resid(lm(y ~ z, weights = w)), w)

ru_cols = grep("^Resource.Use", names(mt), value = TRUE)
length(ru_cols)                                                    # how many RU traits
summary(colMeans(mt[, ru_cols] > 0, na.rm = TRUE))                 # share of MAGs with each RU trait

mt_ru = data.frame(id = mt$MAG_id, RU_n = rowSums(mt[, ru_cols] > 0, na.rm = TRUE))
ids   = total_genes.guild.940_MAG[, c("id", "guild")]
cat("MAG ids found in the trait table:", mean(ids$id %in% mt_ru$id), "\n")

tmp      = merge(ids, mt_ru, by = "id")
ru_guild = aggregate(RU_n ~ guild, data = tmp, FUN = mean)
tab      = merge(MAG_gen_trait, ru_guild, by = "guild")
nrow(tab)

w = tab$n_mags
data.frame(pair = c("A-S", "A-RU", "S-RU"),
           r = c(wcor(tab$A_traits, tab$S_traits, w),
                 wcor(tab$A_traits, tab$RU_n, w),
                 wcor(tab$S_traits, tab$RU_n, w)),
           partial = c(wpartial(tab$A_traits, tab$S_traits, tab$genome.size, w),
                       wpartial(tab$A_traits, tab$RU_n, tab$genome.size, w),
                       wpartial(tab$S_traits, tab$RU_n, tab$genome.size, w)))


# Isolates

total_genes.guild.940_ISO = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/total_genes.guild.940_ISO.csv")
mt_iso = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/data_paper/Data_Base/output/FUNCTIONAL ANNOTATION/isolate_functional_traits_clean.csv")
ISO_gen_TOTAL_ab    = read.csv("Intermediate_Results/Isolates_gen_trait.csv",dec=".")

wcor = function(x, y, w) cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
wpartial = function(x, y, z, w) wcor(resid(lm(x ~ z, weights = w)), resid(lm(y ~ z, weights = w)), w)

id_col = names(mt_iso)[1]                  # isolate id column (assumed to be the first column)
id_col
head(mt_iso[[id_col]], 3)                  # compare with: head(total_genes.guild.940_ISO$id, 3)

ru_cols = grep("^Resource.Use", names(mt_iso), value = TRUE)
length(ru_cols)
summary(colMeans(mt_iso[, ru_cols] > 0, na.rm = TRUE))

mt_ru = data.frame(id = mt_iso[[id_col]], RU_n = rowSums(mt_iso[, ru_cols] > 0, na.rm = TRUE))
ids   = total_genes.guild.940_ISO[, c("id", "guild")]
cat("Isolate ids found in the trait table:", mean(ids$id %in% mt_ru$id), "\n")

tmp      = merge(ids, mt_ru, by = "id")
ru_guild = aggregate(RU_n ~ guild, data = tmp, FUN = mean)
tab_iso  = merge(ISO_gen_TOTAL_ab, ru_guild, by = "guild")
nrow(tab_iso)

w = tab_iso$n_mags
data.frame(pair = c("A-S", "A-RU", "S-RU"),
           r = c(wcor(tab_iso$A_traits, tab_iso$S_traits, w),
                 wcor(tab_iso$A_traits, tab_iso$RU_n, w),
                 wcor(tab_iso$S_traits, tab_iso$RU_n, w)),
           partial = c(wpartial(tab_iso$A_traits, tab_iso$S_traits, tab_iso$genome.size, w),
                       wpartial(tab_iso$A_traits, tab_iso$RU_n, tab_iso$genome.size, w),
                       wpartial(tab_iso$S_traits, tab_iso$RU_n, tab_iso$genome.size, w)))

