
library(readxl)

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
    list(quantiles = apply(nl, 2, quantile, c(.025, .5, .975)), p = p)
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
