# Loma Analysis ----

# Call the MAG-Microtrait datasets
# Calling data and preprocessing
hmm_loma    = read.csv("Input_Data/LOMA_MAGs/hmm_Loma.csv",dec=".")
mag_stat    = read.csv("Input_Data/LOMA_MAGs/mag_stats.csv",dec=".") 
mag_abun    = read.csv("Input_Data/LOMA_MAGs/mag_adundance.csv",dec=".") 
mag_stat    = mag_stat %>% full_join(mag_abun)

# Call Genes and merge datasets - Loma
Loma_predictors_files = list.files(path="Intermediate_Results/LOMA_predictors", pattern=NULL, all.files=FALSE, 
                                   full.names=TRUE)
Loma_predictors_list = lapply(Loma_predictors_files, read.csv)
Loma_predictors_gene = do.call(rbind.data.frame, Loma_predictors_list)
Loma_predictors_gene = unique(Loma_predictors_gene$x)

# Select only the best predictors 
hmm_loma.global      = hmm_loma %>% select(id,tidyselect::any_of(Loma_predictors_gene))

# Clustering

FG_test = function(i,k,seed) {
  set.seed(seed)
  rand_df                = hmm_loma.global[sample(nrow(hmm_loma.global), size=i), ]
  distance.total         = parDist(x = as.matrix(rand_df[2:length(rand_df)]),
                                   method = "fJaccard",
                                   threads = 10) # Adapt the number of threads
  cluster.total          = hclust(distance.total, method="ward.D2")
  v                      = cutree(cluster.total,k=k)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                           perm = 999,p.adjust.m='BH')
  return(as.data.frame(adonis_1))
}

# Number of clusters

adonis = FG_test(nrow(hmm_loma.global),41,1000) # 41 is the final number

# Calling distance matrix

set.seed(1)
distance.total  = parDist(x = as.matrix(hmm_loma.global[,2:182]),
                          method = "fJaccard",
                          threads = 1) # Adapt the number of threads
cluster.total   = hclust(distance.total, method="ward.D2")
v                      = cutree(cluster.total,k=round(41))
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
hmm_loma.global        = as.data.frame(cbind(genome2guild,hmm_loma.global))

# Merge with genome size etc

mag_stat = mag_stat %>%
  mutate(id_join = paste0(sample, "_", bin))

hmm_loma.global = hmm_loma.global %>%
  left_join(mag_stat, by = c("id" = "id_join"))

# Delete guild 35 which is extremely big

hmm_loma.global = hmm_loma.global %>% filter(guild != 35)

loma_guild_n = hmm_loma.global %>% 
  filter(guild != 35) %>% 
  count(guild, name = "n_mags")

# Total Traits

# Call Trait Keys 
GH_rule = read_excel("Input_Data/microTrait_Rules/microtrait_GH.xlsx")
GH_rule = as.data.frame(rbind("id","guild","size",GH_rule))
PR_rule = read_excel("Input_Data/microTrait_Rules/microtrait_proteins.xlsx")
PR_rule = as.data.frame(rbind("id","guild","size",PR_rule))
transp_rule = read_excel("Input_Data/microTrait_Rules/microtrait_transporters.xlsx")
osmo_rule   = read_excel("Input_Data/microTrait_Rules/microtrait_osmolytes.xlsx")
osmo_rule   = as.data.frame(rbind("id","guild","size",osmo_rule))
biofilm_rule = read_excel("Input_Data/microTrait_Rules/microtrait_biofilm.xlsx")
biofilm_rule = as.data.frame(rbind("id","guild","size",biofilm_rule))  
high.T_rule  = read_excel("Input_Data/microTrait_Rules/microtrait_high_T.xlsx")
high.T_rule  = as.data.frame(rbind("id","guild","size",high.T_rule))
pH_rule      = read_excel("Input_Data/microTrait_Rules/microtrait_pH_stress.xlsx")
pH_rule      = as.data.frame(rbind("id","guild","size",pH_rule)) 

# GH - genes (MAG)
GH_TOTAL.MAG = hmm_loma.global %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.MAG$size = as.numeric(as.character(GH_TOTAL.MAG$size))
GH_TOTAL_m   = GH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,3:ncol(GH_TOTAL_m)]))
GH_MAG       = GH_TOTAL_m %>% select(guild,size_mean,GH_total)

# Protein - genes (MAG)
PR_TOTAL.MAG = hmm_loma.global %>% select(any_of(PR_rule$`microtrait_hmm-name`))
PR_TOTAL.MAG$size = as.numeric(as.character(PR_TOTAL.MAG$size))
PR_TOTAL_m   = PR_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PR_TOTAL_m   = PR_TOTAL_m %>% mutate(PR_total = rowSums(PR_TOTAL_m[,3:ncol(PR_TOTAL_m)]))
Protein_MAG  = PR_TOTAL_m %>% select(guild,size_mean,PR_total)

# Transport - genes (MAG)
transp_rule  = transp_rule %>% filter(`function` == c("transporter"))
transp_rule  = as.data.frame(rbind("id","guild","size",transp_rule))
TRANSP_TOTAL = hmm_loma.global %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL$size = as.numeric(as.character(TRANSP_TOTAL$size))
TRANSP_TOTAL_m   = TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m   = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,3:ncol(TRANSP_TOTAL_m)]))
Transporter_MAG  = TRANSP_TOTAL_m %>% select(guild,size_mean,transp_total)

# Transport - Aminoacids (MAG)
transp_rule        = transp_rule %>% filter(`function`==c("transporter"))
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","size",transp_rule_ami))
AMI_TRANSP_TOTAL   = hmm_loma.global %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL$size = as.numeric(as.character(AMI_TRANSP_TOTAL$size))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL_m[,3:ncol(AMI_TRANSP_TOTAL_m)]))
Aminoacids_T_MAG   = AMI_TRANSP_TOTAL_m %>% select(guild,size_mean,AMI_TRANSP_TOTAL)

# Transport - Carbohydrate (MAG)
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","size",transp_rule_car))
CAR_TRANSP_TOTAL   = hmm_loma.global %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL$size = as.numeric(as.character(CAR_TRANSP_TOTAL$size))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL_m[,3:ncol(CAR_TRANSP_TOTAL_m)]))
GH_T_MAG           = CAR_TRANSP_TOTAL_m %>% select(guild,size_mean,CAR_TRANSP_TOTAL)

# Osmolytes - genes (MAG)
OSMO_TOTAL.MAG = hmm_loma.global %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.MAG$size = as.numeric(as.character(OSMO_TOTAL.MAG$size))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.MAG_m[,3:ncol(OSMO_TOTAL.MAG_m)]))
Osmolyte_MAG       = OSMO_TOTAL.MAG_m %>% select(guild,size_mean,OSMO_total)

# Biofilm - genes (MAG)
BIO_TOTAL.MAG = hmm_loma.global %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$size = as.numeric(as.character(BIO_TOTAL.MAG$size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,3:ncol(BIO_TOTAL.MAG_m)]))
Biofilm_MAG       = BIO_TOTAL.MAG_m %>% select(guild,size_mean,BIO_total)

# High Temp - genes (MAG)
TEMP.MAG = hmm_loma.global %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP.MAG$size = as.numeric(as.character(TEMP.MAG$size))
TEMP.MAG_m   = TEMP.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                 list(mean=mean), na.rm=TRUE))
TEMP.MAG_m   = TEMP.MAG_m %>% mutate(TEMP_total = rowSums(TEMP.MAG_m[,3:ncol(TEMP.MAG_m)]))
TEMP_MAG     = TEMP.MAG_m %>% select(guild,size_mean,TEMP_total)

# pH - genes (MAG)
PH.MAG = hmm_loma.global %>% select(any_of(pH_rule$`microtrait_hmm-name`))
PH.MAG$size = as.numeric(as.character(PH.MAG$size))
PH.MAG_m   = PH.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
PH.MAG_m   = PH.MAG_m %>% mutate(PH_total = rowSums(PH.MAG_m[,3:ncol(PH.MAG_m)]))
PH_MAG     = PH.MAG_m %>% select(guild,size_mean,PH_total)

MAG_gen_trait    = as.data.frame(cbind(Aminoacids_T_MAG,PH_MAG$PH_total,
                                       TEMP_MAG$TEMP_total,Biofilm_MAG$BIO_total,
                                       Osmolyte_MAG$OSMO_total,GH_T_MAG$CAR_TRANSP_TOTAL,
                                       Transporter_MAG$transp_total,Protein_MAG$PR_total,
                                       GH_MAG$GH_total))
colnames(MAG_gen_trait) = c("guild","genome.size","amino.transport","pH","temp",
                            "biofilm","osmolyte","GH.trasnport","tranport.total",
                            "Protein","CAZy")

loma_guild_n = loma_guild_n %>% mutate(guild = as.numeric(as.character(guild)))

MAG_gen_trait.agg_LOMA = MAG_gen_trait %>%
  mutate(guild = as.numeric(as.character(guild))) %>%
  left_join(loma_guild_n, by = "guild")

cat("Guilds before:", nrow(MAG_gen_trait.agg_LOMA), "\n")
cat("Guilds with n_mags matched:", sum(!is.na(MAG_gen_trait.agg_LOMA$n_mags)), "\n")

# Aggregated traits

MAG_gen_trait.agg_LOMA   =  MAG_gen_trait.agg_LOMA %>% mutate(S_traits  = rowSums(MAG_gen_trait[,4:7]),
                                                A_traits  = rowSums(MAG_gen_trait[,9:11]),
                                                A_enzymes = rowSums(MAG_gen_trait[,10:11]))
MAG_gen_trait.agg_LOMA   = MAG_gen_trait.agg_LOMA %>% mutate(A_S = A_traits/S_traits)

write.csv(MAG_gen_trait.agg_LOMA, file = "Intermediate_Results/MAG_gen_trait.agg_LOMA.csv")

# Plotting

weighted_cor = function(x, y, w) {
  ok = complete.cases(x, y, w)
  if (sum(ok) < 3) return(list(r = NA, p = NA))
  x = x[ok]; y = y[ok]; w = w[ok]
  r = cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
  p = tryCatch({ summary(lm(y ~ x, weights = w))$coefficients[2, 4] }, error = function(e) NA)
  list(r = r, p = p)
}

res_loma = weighted_cor(MAG_gen_trait.agg_LOMA$A_traits, MAG_gen_trait.agg_LOMA$S_traits,
                        MAG_gen_trait.agg_LOMA$n_mags)
cat("LOMA - weighted r:", round(res_loma$r, 3), "| p:", res_loma$p, "\n")

Figure_Total_S_SA_LOMA = ggplot(data = MAG_gen_trait.agg_LOMA, 
                                aes(x = as.numeric(A_traits), y = as.numeric(S_traits),
                                    color = A_S, weight = n_mags)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = 3, y = 20, hjust = 0,
           label = paste0("r = ", round(res_loma$r, 2), ", p ", 
                          ifelse(res_loma$p < 0.001, "< 0.001", paste0("= ", round(res_loma$p, 3))))) +
  xlab("A Traits") + ylab("S Traits") +
  scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5)) +
  ylim(0,25) + xlim(0,40)
Figure_Total_S_SA_LOMA

# Fire Project Analysis ----

# Call the MAG-Microtrait datasets
# Calling data and preprocessing
hmm_fire    = read.csv("Input_Data/FIRE_MAGs/hmm_Fire.csv",dec=".")
mag_abun    = read.csv("Input_Data/FIRE_MAGs/mag_adundance_fire.csv",dec=".") 
mag_abun    = mag_abun[-c(440,546), ]
fire_meta   = read_excel("Input_Data/FIRE_MAGs/MAG_Dataset_BurnSeverity_ARNelson.xlsx")
fire_size   = read_csv("Input_Data/FIRE_MAGs/fire_metadata.csv")

# Call Genes and merge datasets - Fire
Fire_predictors_files = list.files(path="Intermediate_Results/FIRE_predictors", pattern=NULL, all.files=FALSE, 
                                   full.names=TRUE)
Fire_predictors_list = lapply(Fire_predictors_files, read.csv)
Fire_predictors_gene = do.call(rbind.data.frame, Fire_predictors_list)
Fire_predictors_gene = unique(Fire_predictors_gene$x)

# Select only the best predictors 
hmm_fire.global      = hmm_fire %>% select(id,tidyselect::any_of(Fire_predictors_gene))

# Clustering

FG_test = function(i,k,seed) {
  set.seed(seed)
  rand_df                = hmm_fire.global[sample(nrow(hmm_fire.global), size=i), ]
  distance.total         = parDist(x = as.matrix(rand_df[2:length(rand_df)]),
                                   method = "fJaccard",
                                   threads = 10) # Adapt the number of threads
  cluster.total          = hclust(distance.total, method="ward.D2")
  v                      = cutree(cluster.total,k=k)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                           perm = 999,p.adjust.m='BH')
  return(as.data.frame(adonis_1))
}

# Number of clusters

adonis = FG_test(nrow(hmm_fire.global),18,1) # 18 is the final number

# Calling distance matrix

set.seed(1)
distance.total  = parDist(x = as.matrix(hmm_fire.global[,2:214]),
                          method = "fJaccard",
                          threads = 1) # Adapt the number of threads
cluster.total   = hclust(distance.total, method="ward.D2")
v                      = cutree(cluster.total,k=round(18))
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
hmm_fire.global        = as.data.frame(cbind(genome2guild,hmm_fire.global))

# Merge with genome size etc

hmm_fire.global = hmm_fire.global %>%
  left_join(fire_size, by = c("id"))

fire_guild_n = hmm_fire.global %>% 
  filter(guild != 35) %>% 
  count(guild, name = "n_mags")

# Total Traits

# Call Trait Keys 
GH_rule = read_excel("Input_Data/microTrait_Rules/microtrait_GH.xlsx")
GH_rule = as.data.frame(rbind("id","guild","Genome_size",GH_rule))
PR_rule = read_excel("Input_Data/microTrait_Rules/microtrait_proteins.xlsx")
PR_rule = as.data.frame(rbind("id","guild","Genome_size",PR_rule))
transp_rule = read_excel("Input_Data/microTrait_Rules/microtrait_transporters.xlsx")
osmo_rule   = read_excel("Input_Data/microTrait_Rules/microtrait_osmolytes.xlsx")
osmo_rule   = as.data.frame(rbind("id","guild","Genome_size",osmo_rule))
biofilm_rule = read_excel("Input_Data/microTrait_Rules/microtrait_biofilm.xlsx")
biofilm_rule = as.data.frame(rbind("id","guild","Genome_size",biofilm_rule))  
high.T_rule  = read_excel("Input_Data/microTrait_Rules/microtrait_high_T.xlsx")
high.T_rule  = as.data.frame(rbind("id","guild","Genome_size",high.T_rule))
pH_rule      = read_excel("Input_Data/microTrait_Rules/microtrait_pH_stress.xlsx")
pH_rule      = as.data.frame(rbind("id","guild","Genome_size",pH_rule)) 

# GH - genes (MAG)
GH_TOTAL.MAG = hmm_fire.global %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.MAG$Genome_size = as.numeric(as.character(GH_TOTAL.MAG$Genome_size))
GH_TOTAL_m   = GH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,3:ncol(GH_TOTAL_m)]))
GH_MAG       = GH_TOTAL_m %>% select(guild,Genome_size_mean,GH_total)

# Protein - genes (MAG)
PR_TOTAL.MAG = hmm_fire.global %>% select(any_of(PR_rule$`microtrait_hmm-name`))
PR_TOTAL.MAG$Genome_size = as.numeric(as.character(PR_TOTAL.MAG$Genome_size))
PR_TOTAL_m   = PR_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PR_TOTAL_m   = PR_TOTAL_m %>% mutate(PR_total = rowSums(PR_TOTAL_m[,3:ncol(PR_TOTAL_m)]))
Protein_MAG  = PR_TOTAL_m %>% select(guild,Genome_size_mean,PR_total)

# Transport - genes (MAG)
transp_rule  = transp_rule %>% filter(`function` == c("transporter"))
transp_rule  = as.data.frame(rbind("id","guild","Genome_size",transp_rule))
TRANSP_TOTAL = hmm_fire.global %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL$Genome_size = as.numeric(as.character(TRANSP_TOTAL$Genome_size))
TRANSP_TOTAL_m   = TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m   = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,3:ncol(TRANSP_TOTAL_m)]))
Transporter_MAG  = TRANSP_TOTAL_m %>% select(guild,Genome_size_mean,transp_total)

# Transport - Aminoacids (MAG)
transp_rule        = transp_rule %>% filter(`function`==c("transporter"))
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","Genome_size",transp_rule_ami))
AMI_TRANSP_TOTAL   = hmm_fire.global %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL$Genome_size = as.numeric(as.character(AMI_TRANSP_TOTAL$Genome_size))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL_m[,3:ncol(AMI_TRANSP_TOTAL_m)]))
Aminoacids_T_MAG   = AMI_TRANSP_TOTAL_m %>% select(guild,Genome_size_mean,AMI_TRANSP_TOTAL)

# Transport - Carbohydrate (MAG)
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","Genome_size",transp_rule_car))
CAR_TRANSP_TOTAL   = hmm_fire.global %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL$Genome_size = as.numeric(as.character(CAR_TRANSP_TOTAL$Genome_size))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL_m[,3:ncol(CAR_TRANSP_TOTAL_m)]))
GH_T_MAG           = CAR_TRANSP_TOTAL_m %>% select(guild,Genome_size_mean,CAR_TRANSP_TOTAL)

# Osmolytes - genes (MAG)
OSMO_TOTAL.MAG = hmm_fire.global %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.MAG$Genome_size = as.numeric(as.character(OSMO_TOTAL.MAG$Genome_size))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.MAG_m[,3:ncol(OSMO_TOTAL.MAG_m)]))
Osmolyte_MAG       = OSMO_TOTAL.MAG_m %>% select(guild,Genome_size_mean,OSMO_total)

# Biofilm - genes (MAG)
BIO_TOTAL.MAG = hmm_fire.global %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$Genome_size = as.numeric(as.character(BIO_TOTAL.MAG$Genome_size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,3:ncol(BIO_TOTAL.MAG_m)]))
Biofilm_MAG       = BIO_TOTAL.MAG_m %>% select(guild,Genome_size_mean,BIO_total)

# High Temp - genes (MAG)
TEMP.MAG = hmm_fire.global %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP.MAG$Genome_size = as.numeric(as.character(TEMP.MAG$Genome_size))
TEMP.MAG_m   = TEMP.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                 list(mean=mean), na.rm=TRUE))
TEMP.MAG_m   = TEMP.MAG_m %>% mutate(TEMP_total = rowSums(TEMP.MAG_m[,3:ncol(TEMP.MAG_m)]))
TEMP_MAG     = TEMP.MAG_m %>% select(guild,Genome_size_mean,TEMP_total)

# pH - genes (MAG)
PH.MAG = hmm_fire.global %>% select(any_of(pH_rule$`microtrait_hmm-name`))
PH.MAG$Genome_size = as.numeric(as.character(PH.MAG$Genome_size))
PH.MAG_m   = PH.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
PH.MAG_m   = PH.MAG_m %>% mutate(PH_total = rowSums(PH.MAG_m[,3:ncol(PH.MAG_m)]))
PH_MAG     = PH.MAG_m %>% select(guild,Genome_size_mean,PH_total)

MAG_gen_trait    = as.data.frame(cbind(Aminoacids_T_MAG,PH_MAG$PH_total,
                                       TEMP_MAG$TEMP_total,Biofilm_MAG$BIO_total,
                                       Osmolyte_MAG$OSMO_total,GH_T_MAG$CAR_TRANSP_TOTAL,
                                       Transporter_MAG$transp_total,Protein_MAG$PR_total,
                                       GH_MAG$GH_total))
colnames(MAG_gen_trait) = c("guild","genome.size","amino.transport","pH","temp",
                            "biofilm","osmolyte","GH.trasnport","tranport.total",
                            "Protein","CAZy")

fire_guild_n = fire_guild_n %>% mutate(guild = as.numeric(as.character(guild)))

MAG_gen_trait.agg_FIRE = MAG_gen_trait %>%
  mutate(guild = as.numeric(as.character(guild))) %>%
  left_join(fire_guild_n, by = "guild")

cat("Guilds before:", nrow(MAG_gen_trait.agg_FIRE), "\n")
cat("Guilds with n_mags matched:", sum(!is.na(MAG_gen_trait.agg_FIRE$n_mags)), "\n")

write.csv(MAG_gen_trait.agg_FIRE, file = "Intermediate_Results/MAG_gen_trait.agg_FIRE.csv")

# Aggregated traits

MAG_gen_trait.agg_FIRE   =  MAG_gen_trait.agg_FIRE %>% mutate(S_traits  = rowSums(MAG_gen_trait[,4:7]),
                                                A_traits  = rowSums(MAG_gen_trait[,9:11]),
                                                A_enzymes = rowSums(MAG_gen_trait[,10:11]))
MAG_gen_trait.agg_FIRE   = MAG_gen_trait.agg_FIRE %>% mutate(A_S = A_traits/S_traits)

# Plotting

weighted_cor = function(x, y, w) {
  ok = complete.cases(x, y, w)
  if (sum(ok) < 3) return(list(r = NA, p = NA))
  x = x[ok]; y = y[ok]; w = w[ok]
  r = cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
  p = tryCatch({ summary(lm(y ~ x, weights = w))$coefficients[2, 4] }, error = function(e) NA)
  list(r = r, p = p)
}

res_fire = weighted_cor(MAG_gen_trait.agg_FIRE$A_traits, MAG_gen_trait.agg_FIRE$S_traits,
                        MAG_gen_trait.agg_FIRE$n_mags)
cat("FIRE - weighted r:", round(res_fire$r, 3), "| p:", res_fire$p, "\n")

Figure_Total_S_SA_FIRE = ggplot(data = MAG_gen_trait.agg_FIRE, 
                                aes(x = as.numeric(A_traits), y = as.numeric(S_traits),
                                    color = A_S, weight = n_mags)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = 3, y = 20, hjust = 0,
           label = paste0("r = ", round(res_fire$r, 2), ", p ", 
                          ifelse(res_fire$p < 0.001, "< 0.001", paste0("= ", round(res_fire$p, 3))))) +
  xlab("A Traits") + ylab("S Traits") +
  scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5)) +
  ylim(0,25) + xlim(0,40)
Figure_Total_S_SA_FIRE

# General Plot ----

Figure_S11 = ggarrange(Figure_Total_S_SA_LOMA, Figure_Total_S_SA_FIRE, 
                       labels = c("A","B"), ncol = 2, nrow = 1) + 
  theme(panel.background = element_blank())
Figure_S11

pdf("Output_Data/Figures/Figure_S11.pdf", width=12, height=12*2/5)
print(Figure_S11)
dev.off()

