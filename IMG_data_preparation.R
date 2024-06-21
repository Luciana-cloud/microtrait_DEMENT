# IMG MAG ANALYSIS
# Analysis of 33000 MAGs from IMG + LOMAa and FIRE MAGs

# CALLING PACKAGES----

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(readr)
library(vegan)
library(devtools)
library(ggvegan)
library(pairwiseAdonis)
library(EcolUtils)
library(olsrr)
library(stats)
library(corrr)
library(car)
#install.packages("collinear")
library(collinear)
library(glmnet)

# CALLING DATA----

hmm_img       = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/hmm_MAG_IMG.csv",dec=".")
coverage_img  = read_excel("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/bin_avg_coverage.xlsx")
meta_img.nr   = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/IMG_bindata_withmeta_norestricted.csv",dec=".")

# Change column name to match
colnames(hmm_img)[2]      = "bin_id"
colnames(meta_img.nr)[2]  = "bin_id"
colnames(coverage_img)[1] = "bin_id"

# Merge databases
data_combined = meta_img.nr %>% left_join(hmm_img, by='bin_id') 
data_combined = data_combined %>% left_join(coverage_img, by='bin_id')
data_combined = data_combined %>% filter(Domain == "Bacteria")
write.csv(data_combined, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/IMG_global_datasets.csv")

# Summary MAGs by project to select what project to use for extracting filtering genes
data_combined = data_combined[!(is.na(data_combined$avg_coverage)),]
data_combined = data_combined[!(is.na(data_combined$aaeB)),]
data_combined_IMG.ID  = data_combined %>% group_by(IMG.Genome.ID) %>% count()
data_combined_IMG.ID  = data_combined_IMG.ID %>% filter(n > 20)

# Get hmm files from data_combined i.e. after data filtering
hmm_img.filter = data_combined %>% select(names(hmm_img))

# Select and save individual projects with more than 20 MAGs for further analysis
a         = unique(data_combined_IMG.ID$IMG.Genome.ID)
data_project_100 = c()
for(i in a){
  project = data_combined %>% filter(IMG.Genome.ID == i)
  temp    = project %>% summarise(sum(avg_coverage, na.rm = TRUE))
  colnames(temp) = "Total"
  project = project %>% mutate(RelAbund = avg_coverage/temp$Total)
  data_project_100 = rbind(data_project_100, project) 
}
write.csv(data_project_100, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/IMG_global_datasets_relative_abundance.csv")

# Reduce NAs
data_project_20 = data_project_100[!(is.na(data_project_100$RelAbund)),]
# Reduce columns with zero
temp.1            = as.data.frame(cbind((colSums(data_project_20[,c(69:1789,1791)]))))
erase.1           = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
data_project_20   = data_project_20[ , !(names(data_project_20) %in% erase.1$row_names)]
# SELECTION OF GENES PER PROJECT - IMG---- 
a = unique(data_project_20$IMG.Genome.ID)
b = as.data.frame(colnames(data_project_20))
for(i in a){
  project  = data_project_20 %>% filter(IMG.Genome.ID == i)
  mat_file = as.data.frame(project[,c(69:1684,1686)])
  # Reduce columns with zero
  temp.1            = as.data.frame(cbind((colSums(mat_file[,c(1:1616)]))))
  erase.1           = temp.1 %>% filter(V1==0)
  erase.1$row_names = row.names(erase.1)
  mat_file          = mat_file[,!(names(mat_file) %in% erase.1$row_names)]
  # Best predictors statistics - step 1 Reduce perfect correlation
  d    = colnames(mat_file[,1:length(mat_file)-1])
  selected.predictors.1 = cor_select(df = as.data.frame(mat_file),response = "RelAbund",
                                     predictors = d,cor_method = "pearson",max_cor = 0.99)
  mat_file.p            = mat_file[,(names(mat_file) %in% selected.predictors.1)]
  # Best predictors statistics - step 1 Reduce multicollinearity 
  d                     = colnames(mat_file.p[,1:length(mat_file.p)-1])
  selected.predictors   = collinear(df = as.data.frame(mat_file.p),
                                    predictors = d,cor_method = "pearson",max_cor = 1,
                                    max_vif = 5)
  mat_file.p            = mat_file.p[,(names(mat_file.p) %in% selected.predictors)]
  # Best predictors statistics - step 2
  m.all             = lm(mat_file[,length(mat_file)]~.,data = mat_file.p)
  temp.all          = ols_step_forward_p(m.all, details = FALSE)
  best.predictors   = temp.all[["metrics"]][["variable"]]
  assign(paste("best.predictors_", i, sep = ""), best.predictors) 
}

# save each new data frame as an individual .csv file based on its name
list            = lapply(ls(pattern="best.predictors_"), get)
data_list_names = ls(pattern="best.predictors_")
names(list) <- data_list_names

for(i in names(list)){
  setwd("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/best_predictors")
  write.csv(list[[i]], paste0(i,".csv"))
}

# Complete series because site 3300056972 did not work

for(i in a[183:length(a)]){
  project  = data_project_20 %>% filter(IMG.Genome.ID == i)
  mat_file = as.data.frame(project[,c(69:1684,1686)])
  # Reduce columns with zero
  temp.1            = as.data.frame(cbind((colSums(mat_file[,c(1:1616)]))))
  erase.1           = temp.1 %>% filter(V1==0)
  erase.1$row_names = row.names(erase.1)
  mat_file          = mat_file[,!(names(mat_file) %in% erase.1$row_names)]
  # Best predictors statistics - step 1 Reduce perfect correlation
  d    = colnames(mat_file[,1:length(mat_file)-1])
  selected.predictors.1 = cor_select(df = as.data.frame(mat_file),response = "RelAbund",
                                     predictors = d,cor_method = "pearson",max_cor = 0.99)
  mat_file.p            = mat_file[,(names(mat_file) %in% selected.predictors.1)]
  # Best predictors statistics - step 1 Reduce multicollinearity 
  d                     = colnames(mat_file.p[,1:length(mat_file.p)-1])
  selected.predictors   = collinear(df = as.data.frame(mat_file.p),
                                    predictors = d,cor_method = "pearson",max_cor = 1,
                                    max_vif = 5)
  mat_file.p            = mat_file.p[,(names(mat_file.p) %in% selected.predictors)]
  # Best predictors statistics - step 2
  m.all             = lm(mat_file[,length(mat_file)]~.,data = mat_file.p)
  temp.all          = ols_step_forward_p(m.all, details = FALSE)
  best.predictors   = temp.all[["metrics"]][["variable"]]
  assign(paste("best.predictors_", i, sep = ""), best.predictors) 
}

# save each new data frame as an individual .csv file based on its name
list            = lapply(ls(pattern="best.predictors33_"), get)
data_list_names = ls(pattern="best.predictors33_")
names(list) <- data_list_names

for(i in names(list)){
  setwd("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/best_predictors")
  write.csv(list[[i]], paste0(i,".csv"))
}

# SELECTION OF GENES PER PROJECT - LOMA----

# Calling data and preprocessing
hmm_loma    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/hmm_Loma.csv",dec=".")
gene_loma   = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/litter_mags_metadata.txt",dec=".")
mag_stat    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_stats.txt") 
mag_abun    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_adundance.txt") 
mag_stat    = mag_stat %>% full_join(mag_abun)

# Dataframes from each treatment

grassland_ambient  = hmm_loma %>% mutate(Rel.Abund = mag_stat$Average)  
grassland_ambient  = grassland_ambient %>% mutate(treatment = rep("grassland_ambient",each=nrow(grassland_ambient)))
shrubland_ambient  = hmm_loma %>% mutate(Rel.Abund = mag_stat$Average.2)  
shrubland_ambient  = shrubland_ambient %>%  mutate(treatment = rep("shrubland_ambient",each=nrow(shrubland_ambient)))
grassland_drought  = hmm_loma %>% mutate(Rel.Abund = mag_stat$Average.1)  
grassland_drought  = grassland_drought %>%  mutate(treatment = rep("grassland_drought",each=nrow(grassland_drought)))
shrubland_drought  = hmm_loma %>% mutate(Rel.Abund = mag_stat$Average.3)  
shrubland_drought  = shrubland_drought %>%  mutate(treatment = rep("shrubland_drought",each=nrow(shrubland_drought)))
data.loma          = as.data.frame(rbind(grassland_ambient,shrubland_ambient,
                                         grassland_drought,shrubland_drought))

a = unique(data.loma$treatment)
for(i in a[1]){
  project  = data.loma %>% filter(treatment == i)
  mat_file = as.data.frame(project[,c(3:1298,1299)])
  # Reduce columns with zero
  temp.1            = as.data.frame(cbind((colSums(mat_file[,c(1:1296)]))))
  erase.1           = temp.1 %>% filter(V1==0)
  erase.1$row_names = row.names(erase.1)
  mat_file          = mat_file[,!(names(mat_file) %in% erase.1$row_names)]
  # Best predictors statistics - step 1 Reduce perfect correlation
  d    = colnames(mat_file[,1:length(mat_file)-1])
  selected.predictors.1 = cor_select(df = as.data.frame(mat_file),response = "Rel.Abund",
                                     predictors = d,cor_method = "pearson",max_cor = 0.99)
  mat_file.p            = mat_file[,(names(mat_file) %in% selected.predictors.1)]
  # Best predictors statistics - step 1 Reduce multicollinearity 
  d                     = colnames(mat_file.p[,1:length(mat_file.p)-1])
  selected.predictors   = collinear(df = as.data.frame(mat_file.p),
                                    predictors = d,cor_method = "pearson",max_cor = 1,
                                    max_vif = 5)
  mat_file.p            = mat_file.p[,(names(mat_file.p) %in% selected.predictors)]
  # Best predictors statistics - step 2
  m.all             = lm(mat_file[,length(mat_file)]~.,data = mat_file.p)
  temp.all          = ols_step_forward_p(m.all, details = FALSE)
  best.predictors   = temp.all[["metrics"]][["variable"]]
  assign(paste("best.predictors_", i, sep = ""), best.predictors) 
}

# save each new data frame as an individual .csv file based on its name
list            = lapply(ls(pattern="best.predictors_"), get)
data_list_names = ls(pattern="best.predictors_")
names(list) <- data_list_names

for(i in names(list)){
  setwd("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/best_predictors")
  write.csv(list[[i]], paste0(i,".csv"))
}

# SELECTION OF GENES PER PROJECT - FIRE----

# Calling data and preprocessing
hmm_fire    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/hmm_Fire.csv",dec=".")
mag_abun    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/mag_adundance_fire.txt") 
mag_abun    = mag_abun[-c(440,546), ]

# Dataframes from each treatment
Low_shallow   = hmm_fire %>% mutate(Rel.Abund = mag_abun$Avg_low_shallow)
Low_shallow   = Low_shallow %>% mutate(treatment = rep("Low_shallow",each=nrow(Low_shallow)))
Low_deep      = hmm_fire %>% mutate(Rel.Abund = mag_abun$Avg_low_deep)
Low_deep      = Low_deep %>% mutate(treatment = rep("Low_deep",each=nrow(Low_deep)))
High_shallow  = hmm_fire %>% mutate(Rel.Abund = mag_abun$Avg_high_shallow)
High_shallow  = High_shallow %>% mutate(treatment = rep("High_shallow",each=nrow(High_shallow)))
High_deep     = hmm_fire %>% mutate(Rel.Abund = mag_abun$Avg_high_deep)
High_deep     = High_deep %>% mutate(treatment = rep("High_deep",each=nrow(High_deep)))
data.fire     = as.data.frame(rbind(Low_shallow,Low_deep,High_shallow,High_deep))

a = unique(data.fire$treatment)
for(i in a[c(3,4)]){
  project  = data.fire %>% filter(treatment == i)
  mat_file = as.data.frame(project[,c(3:1214,1215)])
  # Reduce columns with zero
  temp.1            = as.data.frame(cbind((colSums(mat_file[,c(1:1212)]))))
  erase.1           = temp.1 %>% filter(V1==0)
  erase.1$row_names = row.names(erase.1)
  mat_file          = mat_file[,!(names(mat_file) %in% erase.1$row_names)]
  # Best predictors statistics - step 1 Reduce perfect correlation
  d    = colnames(mat_file[,1:length(mat_file)-1])
  selected.predictors.1 = cor_select(df = as.data.frame(mat_file),response = "Rel.Abund",
                                     predictors = d,cor_method = "pearson",max_cor = 0.99)
  mat_file.p            = mat_file[,(names(mat_file) %in% selected.predictors.1)]
  # Best predictors statistics - step 1 Reduce multicollinearity 
  d                     = colnames(mat_file.p[,1:length(mat_file.p)-1])
  selected.predictors   = collinear(df = as.data.frame(mat_file.p),
                                    predictors = d,cor_method = "pearson",max_cor = 1,
                                    max_vif = 5)
  mat_file.p            = mat_file.p[,(names(mat_file.p) %in% selected.predictors)]
  # Best predictors statistics - step 2
  m.all             = lm(mat_file[,length(mat_file)]~.,data = mat_file.p)
  temp.all          = ols_step_forward_p(m.all, details = FALSE)
  best.predictors   = temp.all[["metrics"]][["variable"]]
  assign(paste("best.predictors_", i, sep = ""), best.predictors) 
}

# save each new data frame as an individual .csv file based on its name
list            = lapply(ls(pattern="best.predictors_"), get)
data_list_names = ls(pattern="best.predictors_")
names(list) <- data_list_names

for(i in names(list)){
  setwd("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/best_predictors")
  write.csv(list[[i]], paste0(i,".csv"))
}

# FUNCTIONAL GROUP ANALYSIS----

# Call Genes and merge datasets - IMG
IMG_predictors_files = list.files(path="C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/best_predictors", pattern=".csv", all.files=FALSE, 
                            full.names=TRUE)
IMG_predictors_list = lapply(IMG_predictors_files, read.csv)
IMG_predictors_gene = do.call(rbind.data.frame, IMG_predictors_list)
IMG_predictors_gene = unique(IMG_predictors_gene$x)

# Call Genes and merge datasets - Loma
Loma_predictors_files = list.files(path="C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/best_predictors", pattern=NULL, all.files=FALSE, 
                                  full.names=TRUE)
Loma_predictors_list = lapply(Loma_predictors_files, read.csv)
Loma_predictors_gene = do.call(rbind.data.frame, Loma_predictors_list)
Loma_predictors_gene = unique(Loma_predictors_gene$x)

# Call Genes and merge datasets - Fire
Fire_predictors_files = list.files(path="C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/best_predictors", pattern=NULL, all.files=FALSE, 
                                  full.names=TRUE)
Fire_predictors_list = lapply(Fire_predictors_files, read.csv)
Fire_predictors_gene = do.call(rbind.data.frame, Fire_predictors_list)
Fire_predictors_gene = unique(Fire_predictors_gene$x)

# Call Genes and merge datasets - Co-assembly of Harvard Forest metagenomes from Barre Woods
Harvard.Forest       = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/best_predictors/best.predictors_3300050821.csv",dec=".")
Harvard.Forest_gene  = unique(Harvard.Forest$x)

# Call Genes and merge datasets - EAA2017 plots 3/5 metagenome combined assembly
EAA2017.plot       = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/best_predictors/best.predictors_3300043331.csv",dec=".")
EAA2017.plot_gene  = unique(EAA2017.plot$x)

# GLOBAL ANALYSIS----

# Call data (MAGs)
total_genes            = as.data.frame(unique(c(IMG_predictors_gene,Loma_predictors_gene,Fire_predictors_gene)))
colnames(total_genes)  = c("gene")
colnames(hmm_img.filter)[2]   = "id"
total_gene_names       = c("id",total_genes$gene)
hmm_fire.global        = hmm_fire %>% select(tidyselect::any_of(total_gene_names))
hmm_loma.global        = hmm_loma %>% select(tidyselect::any_of(total_gene_names))
hmm_img.global         = hmm_img.filter %>% select(tidyselect::any_of(total_gene_names))

# Call data (Isolates)
isolates        = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/hmm_isolates.csv",dec=".")
df_1            = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils_unpublished.tsv",sep="\t")
df_1            = subset(df_1, select = -c(23))
df_2            = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils.tsv",sep="\t")
df_2            = subset(df_2, select = -c(14,16,18))
iso.metadata    = as.data.frame(rbind(df_1,df_2))
iso.metadata    = filter(iso.metadata, High.Quality %in% c("Yes"))
isolates        = isolates[isolates$id %in% iso.metadata$taxon_oid, ]
isolates.global = isolates %>% select(tidyselect::any_of(total_gene_names))

# Complete missing columns with zeros
hmm_fire.global[setdiff(names(hmm_img.global), names(hmm_fire.global))] = 0
hmm_loma.global[setdiff(names(hmm_img.global), names(hmm_loma.global))] = 0

# Merge data global
MAG_global = as.data.frame(rbind(hmm_fire.global,hmm_loma.global,hmm_img.global,
                                 isolates.global))
write.csv(MAG_global, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_global.csv")

# ANALYSIS WITH LOMA GENES----

# Filter data (MAGs)
LOMA_gene_names        = c("id",Loma_predictors_gene)
hmm_fire.LOMA          = hmm_fire %>% select(tidyselect::any_of(LOMA_gene_names))
hmm_loma.LOMA          = hmm_loma %>% select(tidyselect::any_of(LOMA_gene_names))
hmm_img.LOMA           = hmm_img.filter %>% select(tidyselect::any_of(LOMA_gene_names))
isolates.LOMA          = isolates %>% select(tidyselect::any_of(LOMA_gene_names))

# Complete missing columns with zeros
hmm_fire.LOMA[setdiff(names(hmm_loma.LOMA), names(hmm_fire.LOMA))] = 0

# Merge data global
MAG_LOMA = as.data.frame(rbind(hmm_fire.LOMA,hmm_loma.LOMA,hmm_img.LOMA,
                               isolates.LOMA))
write.csv(MAG_LOMA, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_LOMA.csv")

# ANALYSIS WITH FIRE GENES----

# Filter data (MAGs)
FIRE_gene_names        = c("id",Fire_predictors_gene)
hmm_fire.FIRE          = hmm_fire %>% select(tidyselect::any_of(FIRE_gene_names))
hmm_loma.FIRE          = hmm_loma %>% select(tidyselect::any_of(FIRE_gene_names))
hmm_img.FIRE           = hmm_img.filter %>% select(tidyselect::any_of(FIRE_gene_names))
isolates.FIRE          = isolates %>% select(tidyselect::any_of(FIRE_gene_names))

# Complete missing columns with zeros
hmm_loma.FIRE[setdiff(names(hmm_fire.FIRE), names(hmm_loma.FIRE))] = 0

# Merge data global
MAG_FIRE = as.data.frame(rbind(hmm_fire.FIRE,hmm_loma.FIRE,hmm_img.FIRE,
                               isolates.FIRE))
write.csv(MAG_FIRE, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_FIRE.csv")

# ANALYSIS WITH HARVARD.FOREST GENES----

# Filter data (MAGs)
HARVARD.FOREST_gene_names        = c("id",Harvard.Forest_gene)
hmm_fire.HARVARD.FOREST          = hmm_fire %>% select(tidyselect::any_of(HARVARD.FOREST_gene_names))
hmm_loma.HARVARD.FOREST          = hmm_loma %>% select(tidyselect::any_of(HARVARD.FOREST_gene_names))
hmm_img.HARVARD.FOREST           = hmm_img.filter %>% select(tidyselect::any_of(HARVARD.FOREST_gene_names))
isolates.HARVARD.FOREST          = isolates %>% select(tidyselect::any_of(HARVARD.FOREST_gene_names))

# Complete missing columns with zeros
hmm_fire.HARVARD.FOREST[setdiff(names(hmm_img.HARVARD.FOREST), names(hmm_fire.HARVARD.FOREST))] = 0
hmm_loma.HARVARD.FOREST[setdiff(names(hmm_img.HARVARD.FOREST), names(hmm_loma.HARVARD.FOREST))] = 0

# Merge data global
MAG_HARVARD.FOREST = as.data.frame(rbind(hmm_fire.HARVARD.FOREST,hmm_loma.HARVARD.FOREST,
                                         hmm_img.HARVARD.FOREST,isolates.HARVARD.FOREST))
write.csv(MAG_HARVARD.FOREST, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_HARVARD.FOREST.csv")

# ANALYSIS WITH EAA2017.PLOT GENES----

# Filter data (MAGs)
EAA2017.PLOT_gene_names        = c("id",EAA2017.plot_gene)
hmm_fire.EAA2017.PLOT          = hmm_fire %>% select(tidyselect::any_of(EAA2017.PLOT_gene_names))
hmm_loma.EAA2017.PLOT          = hmm_loma %>% select(tidyselect::any_of(EAA2017.PLOT_gene_names))
hmm_img.EAA2017.PLOT           = hmm_img.filter %>% select(tidyselect::any_of(EAA2017.PLOT_gene_names))
isolates.EAA2017.PLOT          = isolates %>% select(tidyselect::any_of(EAA2017.PLOT_gene_names))

# Complete missing columns with zeros
hmm_fire.EAA2017.PLOT[setdiff(names(hmm_img.EAA2017.PLOT), names(hmm_fire.EAA2017.PLOT))] = 0
hmm_loma.EAA2017.PLOT[setdiff(names(hmm_img.EAA2017.PLOT), names(hmm_loma.EAA2017.PLOT))] = 0

# Merge data global
MAG_EAA2017.PLOT = as.data.frame(rbind(hmm_fire.EAA2017.PLOT,hmm_loma.EAA2017.PLOT,
                                       hmm_img.EAA2017.PLOT,isolates.EAA2017.PLOT))
write.csv(MAG_EAA2017.PLOT, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_EAA2017.PLOT.csv")

# HPC ANALYSIS----

# See file test.file.R


