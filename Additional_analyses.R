
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
A_all = unique(c(key(transp_rule[transp_rule$`function` == "transporter", ]),
                 key(PR_rule), key(GH_rule)))
S_all = unique(c(key(high.T_rule), key(pH_rule),                     # CHANGED: low-T added
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

fmt = function(q) sprintf("%.2f (%.2f-%.2f)", q[2], q[1], q[3])   # median (2.5%-97.5%)
tab_S20 = data.frame(
  Dataset = c("MAGs", "Isolates"),
  `Real A-S (genome size removed)` = sprintf("%.2f", c(res_MAG_w$observed["partial"],
                                                       res_ISO_w$observed["partial"])),
  `Random genes outside A/S` = c(fmt(res_MAG_w$null_all_genes$quantiles[, "partial"]),
                                 fmt(res_ISO_w$null_all_genes$quantiles[, "partial"])),
  `Random splits of selected genes` = c(fmt(res_MAG_w$null_selected_genes$quantiles[, "partial"]),
                                        fmt(res_ISO_w$null_selected_genes$quantiles[, "partial"])),
  check.names = FALSE)
tab_S20
write.csv(tab_S20, "Output_Data/Table_S20.csv", row.names = FALSE)

make_ru = function(tb, label) {
  w = tb$n_mags
  data.frame(
    Dataset = label,
    FGs     = nrow(tb),
    Pair    = c("A-S", "A-RU", "S-RU"),
    r       = c(wcor(tb$A_traits, tb$S_traits, w),
                wcor(tb$A_traits, tb$RU_n, w),
                wcor(tb$S_traits, tb$RU_n, w)),
    partial = c(wpartial(tb$A_traits, tb$S_traits, tb$genome.size, w),
                wpartial(tb$A_traits, tb$RU_n, tb$genome.size, w),
                wpartial(tb$S_traits, tb$RU_n, tb$genome.size, w)))
}
tab_S21 = rbind(make_ru(tab, "MAGs"), make_ru(tab_iso, "Isolates"))
tab_S21[, c("r", "partial")] = round(tab_S21[, c("r", "partial")], 2)
tab_S21
write.csv(tab_S21, "Output_Data/Table_S21.csv", row.names = FALSE)

# ===================== Genome Level Analysis (reviewer figure) =====================
library(readxl); library(readr); library(dplyr); library(ggplot2)

# ===================== Open relevant data =====================
total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv", dec = ".")
GH_rule      = as.data.frame(rbind("id","guild","genome.size", read_excel("Input_Data/microTrait_Rules/microtrait_GH.xlsx")))
PR_rule      = as.data.frame(rbind("id","guild","genome.size", read_excel("Input_Data/microTrait_Rules/microtrait_proteins.xlsx")))
osmo_rule    = as.data.frame(rbind("id","guild","genome.size", read_excel("Input_Data/microTrait_Rules/microtrait_osmolytes.xlsx")))
biofilm_rule = as.data.frame(rbind("id","guild","genome.size", read_excel("Input_Data/microTrait_Rules/microtrait_biofilm.xlsx")))
high.T_rule  = as.data.frame(rbind("id","guild","genome.size", read_excel("Input_Data/microTrait_Rules/microtrait_high_T.xlsx")))
pH_rule      = as.data.frame(rbind("id","guild","genome.size", read_excel("Input_Data/microTrait_Rules/microtrait_pH_stress.xlsx")))
transp_rule  = read_excel("Input_Data/microTrait_Rules/microtrait_transporters.xlsx")

MAG_gen_trait    = read.csv("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Intermediate_Results/MAG_gen_trait.csv")
TOTAL_MAGs_fixed = read_csv("Intermediate_Results/TOTAL_MAGs_fixed.csv")
sim_10000        = read.csv("Output_Data/my_MAGs_BGE_glucose_10000.csv")
microtrait       = read_csv("Input_Data/IMG_JGI_MAGs/total.granularity.3_datasets.csv")

# ===================== Create MAG_trait (one row per MAG) =====================
dat = total_genes.guild.940_MAG
dat$genome.size = as.numeric(as.character(dat$genome.size))

genes_in = function(rule) {
  x = intersect(rule$`microtrait_hmm-name`, names(dat))
  setdiff(x, c("id", "guild", "genome.size"))
}

MAG_trait = data.frame(id = dat$id, guild = dat$guild, genome.size = dat$genome.size)

rules = list(CAZy = GH_rule, protein = PR_rule, osmolyte = osmo_rule,
             biofilm = biofilm_rule, temp = high.T_rule, pH = pH_rule)
for (nm in names(rules)) {
  MAG_trait[[nm]] = rowSums(dat[, genes_in(rules[[nm]]), drop = FALSE])
}

# transporters (same filters as in the FG analysis; AA = aminoacid + peptide, with %in%)
tr0 = transp_rule %>% filter(`function` == "transporter")
g_total = intersect(tr0$`microtrait_hmm-name`, names(dat))
g_ami   = intersect((tr0 %>% filter(class %in% c("aminoacid", "peptide")))$`microtrait_hmm-name`, names(dat))
g_car   = intersect((tr0 %>% filter(class == "carbohydrate"))$`microtrait_hmm-name`, names(dat))
c(total = length(g_total), amino_peptide = length(g_ami), carbohydrate = length(g_car))   # genes per trait

MAG_trait$total = rowSums(dat[, g_total, drop = FALSE])
MAG_trait$AA    = rowSums(dat[, g_ami,   drop = FALSE])
MAG_trait$GH    = rowSums(dat[, g_car,   drop = FALSE])    # "GH" transporters = carbohydrate transporters

# same exclusion as in the FG analysis (guild 247, 2 MAGs)
MAG_trait = MAG_trait[MAG_trait$guild != 247, ]
nrow(MAG_trait)                                             # expect 27212

# ===================== CUE =====================
cue_MAG = sim_10000 %>% transmute(id = MAG_id, CUE = ifelse(BGE <= 0.848, BGE, NA_real_))
MAG_trait = MAG_trait %>% select(-any_of(c("CUE", "CUE.x", "CUE.y"))) %>% left_join(cue_MAG, by = "id")
sum(!is.na(MAG_trait$CUE))                                  # expect 14075

cue_fg = sim_10000 %>%
  inner_join(dat %>% select(id, guild), by = c("MAG_id" = "id")) %>%
  filter(BGE <= 0.848) %>% group_by(guild) %>% summarise(CUE_fg = mean(BGE), n_cue_mags = n())
chk2 = MAG_trait %>% filter(!is.na(CUE)) %>% group_by(guild) %>% summarise(CUE = mean(CUE), n = n())
chk2 = merge(chk2, cue_fg, by = "guild")
cor(chk2$CUE, chk2$CUE_fg); max(abs(chk2$CUE - chk2$CUE_fg)); sum(chk2$n != chk2$n_cue_mags)   # expect 1, 0, 0

# ===================== ogt, mgt, mgr =====================
growth = microtrait %>% filter(id %in% MAG_trait$id) %>% transmute(id, OGT = ogt, MGT = mgt, MGR = 1 / mgt)
MAG_trait = MAG_trait %>%
  select(-any_of(c("OGT", "MGT", "MGR", "OGT.x", "OGT.y", "MGT.x", "MGT.y", "MGR.x", "MGR.y"))) %>%
  left_join(growth, by = "id")

chk3 = MAG_trait %>% group_by(guild) %>% summarise(OGT = mean(OGT, na.rm = TRUE), MGT = mean(MGT, na.rm = TRUE))
chk3 = merge(chk3, MAG_gen_trait[, c("guild", "ogt", "mgt")], by = "guild")
ok = complete.cases(chk3)
cor(chk3$OGT[ok], chk3$ogt[ok]); max(abs(chk3$OGT[ok] - chk3$ogt[ok]))   # expect 1, ~0
cor(chk3$MGT[ok], chk3$mgt[ok]); max(abs(chk3$MGT[ok] - chk3$mgt[ok]))   # expect 1, ~0

# ===================== completeness, contamination, GC =====================
qual = TOTAL_MAGs_fixed %>% filter(id %in% MAG_trait$id) %>%
  transmute(id, completeness, contamination, GC = GC_count)
MAG_trait = MAG_trait %>%
  select(-any_of(c("completeness", "contamination", "GC", "completeness.x", "completeness.y",
                   "contamination.x", "contamination.y", "GC.x", "GC.y"))) %>%
  left_join(qual, by = "id")

# ===================== Check every trait against your FG table =====================
fgmap = c(CAZy = "CAZy", protein = "Protein", osmolyte = "osmolyte", biofilm = "biofilm",
          temp = "temp", pH = "pH", total = "tranport.total", GH = "GH.trasnport",
          completeness = "completeness_mean", contamination = "contamination_mean", GC = "GC_count_mean")
for (tr in names(fgmap)) {
  col = unname(fgmap[tr])
  m = aggregate(MAG_trait[[tr]], list(guild = MAG_trait$guild), mean, na.rm = TRUE)
  m = merge(m, MAG_gen_trait[, c("guild", col)], by = "guild")
  m = m[complete.cases(m), ]
  cat(sprintf("%-14s r = %.4f | max abs diff = %.4f\n", tr, cor(m$x, m[[col]]), max(abs(m$x - m[[col]]))))
}
# expect r = 1.0000 and diff = 0.0000 for all lines (AA is not listed: it has the NEW definition)

# ===================== Reviewer figure: FG-level vs MAG-level correlations =====================
traits = c("OGT", "MGR", "CUE", "AA", "GH", "total", "protein", "CAZy", "pH", "temp",
           "biofilm", "osmolyte", "completeness", "contamination", "GC")

# FG level: mean per guild, number of MAGs, number of MAGs with CUE
fg = MAG_trait %>% group_by(guild) %>%
  summarise(n_mags = n(), n_cue_mags = sum(!is.na(CUE)),
            across(c(OGT, MGT, CUE, AA, GH, total, protein, CAZy, pH, temp, biofilm,
                     osmolyte, completeness, contamination, GC),
                   ~ mean(.x, na.rm = TRUE)))
fg$MGR = 1 / fg$MGT                       # as in the FG analysis: 1 / mean(mgt)

# weighted Pearson correlation (as in S3: weights = number of MAGs per guild)
wcor = function(a, b, w) {
  ok = complete.cases(a, b, w) & w > 0
  cov.wt(cbind(a[ok], b[ok]), wt = w[ok], cor = TRUE)$cor[1, 2]
}

res = as.data.frame(t(combn(traits, 2)), stringsAsFactors = FALSE)
names(res) = c("x", "y")
res$r_FG = unname(mapply(function(a, b) {
  w = if ("CUE" %in% c(a, b)) fg$n_cue_mags else fg$n_mags
  wcor(fg[[a]], fg[[b]], w)
}, res$x, res$y))
res$r_MAG = unname(mapply(function(a, b) cor(MAG_trait[[a]], MAG_trait[[b]], use = "pairwise.complete.obs"),
                          res$x, res$y))
res$n_MAG = unname(mapply(function(a, b) sum(complete.cases(MAG_trait[, c(a, b)])), res$x, res$y))

# CHECK against S3 (CUE pairs; AA is expected to differ):
cue_rows = subset(res, x == "CUE" | y == "CUE")[, c("x", "y", "r_FG")]
cue_rows$r_FG = round(cue_rows$r_FG, 3)
print(cue_rows)

# summary for the answer
c(n_pairs = nrow(res), r_between = cor(res$r_FG, res$r_MAG),
  median_abs_diff = median(abs(res$r_FG - res$r_MAG)),
  same_sign = sum(sign(res$r_FG) == sign(res$r_MAG)))

write.csv(res, "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Output_Data/Table_S_FG_vs_MAG_correlations.csv", row.names = FALSE)

res$type = ifelse(res$x %in% c("OGT", "MGR", "CUE") | res$y %in% c("OGT", "MGR", "CUE"),
                  "Pair includes a Y-like trait", "Other pairs")

p = ggplot(res, aes(r_FG, r_MAG, colour = type)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(x = "Correlation at functional group level (weighted)",
       y = "Correlation at MAG level (unweighted)", colour = NULL) +
  theme_bw() + theme(legend.position = "bottom")
p
ggsave("C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Output_Data/Figures/Figure_S_FG_vs_MAG.pdf", p, width = 6, height = 6.5)

# Additional plots for Loma and Fire project ----

key = function(df) setdiff(unique(trimws(as.character(df[[1]]))),
                           c("id", "guild", "genome.size", "size"))
rd  = function(f) read_excel(file.path("Input_Data/microTrait_Rules", f))

tr = rd("microtrait_transporters.xlsx")
A_all = unique(c(key(tr[tr$`function` == "transporter", ]),
                 key(rd("microtrait_proteins.xlsx")), key(rd("microtrait_GH.xlsx"))))
S_all = unique(c(key(rd("microtrait_high_T.xlsx")), key(rd("microtrait_pH_stress.xlsx")),
                 key(rd("microtrait_osmolytes.xlsx")), key(rd("microtrait_biofilm.xlsx"))))
shared = intersect(A_all, S_all)

wcor = function(x, y, w) cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]

# ---- Random baseline and genome-size-removed values (local datasets) ----

run_control_local = function(d, gene_cols, selected_genes, n_iter = 1000) {
  sums = rowsum(as.matrix(d[, gene_cols]), group = d$guild)
  n    = tabulate(match(d$guild, rownames(sums)))
  G    = sums / n
  gs   = as.vector(tapply(d$genome.size, d$guild, mean, na.rm = TRUE)[rownames(G)])
  A_g  = setdiff(intersect(A_all, selected_genes), shared)
  S_g  = setdiff(intersect(S_all, selected_genes), shared)
  nA = length(A_g); nS = length(S_g)
  
  stat = function(Ag, Sg) {
    a = rowSums(G[, Ag, drop = FALSE]); s = rowSums(G[, Sg, drop = FALSE])
    c(r = wcor(a, s, n),
      partial = wcor(resid(lm(a ~ gs, weights = n)), resid(lm(s ~ gs, weights = n)), n))
  }
  observed = stat(A_g, S_g)
  
  null_for = function(pool) {
    set.seed(1)
    nl = t(replicate(n_iter, { g = sample(pool, nA + nS); stat(g[1:nA], g[(nA+1):(nA+nS)]) }))
    p = sapply(c("r", "partial"), function(k) (sum(nl[, k] >= observed[k]) + 1) / (n_iter + 1))
    list(quantiles = apply(nl, 2, quantile, c(.025, .5, .975)), p = p, draws = nl)
  }
  
  a = rowSums(G[, A_g, drop = FALSE]); s = rowSums(G[, S_g, drop = FALSE])
  list(observed = observed, sizes = c(A = nA, S = nS),
       null_outside = null_for(setdiff(gene_cols, c(A_g, S_g, shared))),
       null_split   = null_for(setdiff(selected_genes, shared)),
       data = data.frame(guild = rownames(G), n = n, gs = gs, a = a, s = s,
                         a_res = resid(lm(a ~ gs, weights = n)),
                         s_res = resid(lm(s ~ gs, weights = n))))
}

# Inputs: one row per MAG, with guild, genome size and ALL gene columns
d_loma = merge(hmm_loma.global[, c("id","guild","size")], hmm_loma, by = "id")
names(d_loma)[names(d_loma) == "size"] = "genome.size"
d_loma$genome.size = as.numeric(as.character(d_loma$genome.size))
gene_cols_loma = setdiff(names(hmm_loma), c("X","id"))
selected_loma  = intersect(Loma_predictors_gene, names(hmm_loma))

d_fire = merge(hmm_fire.global[, c("id","guild","Genome_size")], hmm_fire, by = "id")
names(d_fire)[names(d_fire) == "Genome_size"] = "genome.size"
d_fire$genome.size = as.numeric(as.character(d_fire$genome.size))
gene_cols_fire = setdiff(names(hmm_fire), c("X","id"))
selected_fire  = intersect(Fire_predictors_gene, names(hmm_fire))

ctrl_loma = run_control_local(d_loma, gene_cols_loma, selected_loma)
ctrl_fire = run_control_local(d_fire, gene_cols_fire, selected_fire)

# Panels C and D: A and S traits divided by genome size (as in Fig. S2)
plot_AS_div = function(d) {
  d$A_n = d$A_traits / d$genome.size
  d$S_n = d$S_traits / d$genome.size
  res   = weighted_cor(d$A_n, d$S_n, d$n_mags)
  ggplot(d, aes(A_n, S_n, color = A_S, weight = n_mags)) +
    geom_point(size = 2.5) + stat_poly_line() +
    theme_classic() + theme(text = element_text(size = 14)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("r = ", round(res$r, 2), ", p ",
                            ifelse(res$p < 0.001, "< 0.001", paste0("= ", round(res$p, 3))))) +
    scale_color_gradient(low = "red", high = "blue", limits = c(0.5, 8.5)) +
    xlab("A traits / genome size") + ylab("S traits / genome size")
}

# Panels E and F: A and S traits with genome size removed (residuals, as in Table S20)
plot_AS_res = function(res) {
  d = res$data; d$A_S = d$a / d$s
  n = nrow(d); r = res$observed["partial"]
  p = 2 * pt(-abs(r * sqrt((n - 3) / (1 - r^2))), n - 3)
  ggplot(d, aes(a_res, s_res, color = A_S, weight = n)) +
    geom_point(size = 2.5) + stat_poly_line() +
    theme_classic() + theme(text = element_text(size = 14)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("r = ", round(r, 2), ", p ",
                            ifelse(p < 0.001, "< 0.001", paste0("= ", round(p, 3))))) +
    scale_color_gradient(low = "red", high = "blue", limits = c(0.5, 8.5)) +
    xlab("A traits (residuals)") + ylab("S traits (residuals)")
}

Figure_S11_ext = ggarrange(Figure_Total_S_SA_LOMA, Figure_Total_S_SA_FIRE,
                           plot_AS_div(MAG_gen_trait.agg_LOMA), plot_AS_div(MAG_gen_trait.agg_FIRE),
                           plot_AS_res(ctrl_loma), plot_AS_res(ctrl_fire),
                           labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)
Figure_S11_ext
pdf("Output_Data/Figures/Figure_S11_extended.pdf", width = 12, height = 14)
print(Figure_S11_ext)
dev.off()

# Table S25: real A-S vs random gene sets
fmt = function(q) sprintf("%.2f (%.2f-%.2f)", q[2], q[1], q[3])
mk = function(ctrl, name) data.frame(
  Dataset = name, FGs = nrow(ctrl$data), `A genes` = ctrl$sizes["A"], `S genes` = ctrl$sizes["S"],
  `Real A-S` = sprintf("%.2f", ctrl$observed["r"]),
  `Real A-S (genome size removed)` = sprintf("%.2f", ctrl$observed["partial"]),
  `Random genes outside A/S` = fmt(ctrl$null_outside$quantiles[, "partial"]),
  `p outside` = sprintf("%.3f", mean(ctrl$null_outside$draws[, "partial"] <= ctrl$observed["partial"])),
  `Random splits of selected genes` = fmt(ctrl$null_split$quantiles[, "partial"]),
  `p splits` = sprintf("%.3f", mean(ctrl$null_split$draws[, "partial"] <= ctrl$observed["partial"])),
  check.names = FALSE, row.names = NULL)
tab_S25 = rbind(mk(ctrl_loma, "Loma"), mk(ctrl_fire, "Fire"))
tab_S25
write.csv(tab_S25, "Output_Data/Table_S25_local_control.csv", row.names = FALSE)

# Prevalence and Bias ----

total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")

g = total_genes.guild.940_MAG
is01  = sapply(g, function(v) is.numeric(v) && all(v %in% c(0, 1, NA)))
genes = names(g)[is01]
length(genes)                                   # expect 936

prev = colSums(g[, genes], na.rm = TRUE)        # number of MAGs carrying each gene
summary(prev)
c(only_1_MAG = sum(prev == 1), under_10 = sum(prev < 10), under_100 = sum(prev < 100))

n_mag = nrow(g)

gene_tab = data.frame(Gene = genes,
                      `MAGs with gene` = as.integer(prev),
                      `% of MAGs` = round(100 * prev / n_mag, 2),
                      check.names = FALSE)

gene_tab$Category = ifelse(gene_tab$Gene %in% A_genes & gene_tab$Gene %in% S_genes, "A and S",
                           ifelse(gene_tab$Gene %in% A_genes, "A",
                                  ifelse(gene_tab$Gene %in% S_genes, "S", "Neither")))

table(gene_tab$Category)
gene_tab = gene_tab[order(-gene_tab$`MAGs with gene`), ]
write.csv(gene_tab, "C:/luciana_datos/project_1_Microtrait_DEMENT/Manuscript/Dataset/Output_Data/Table_S_gene_prevalence.csv", row.names = FALSE)


