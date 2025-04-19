# IMG MAG ANALYSIS
# Analysis of 33000 MAGs from IMG + LOMAa and FIRE MAGs + 5000 Soil Isolates

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
library(parallelDist)
library(googledrive)
library(googlesheets4)
library(readxl)
library(ggpmisc)
library(GGally)

# CALLING DATA----

hmm_img       = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/hmm_MAG_IMG.csv",dec=".")
coverage_img  = read_sheet("https://docs.google.com/spreadsheets/d/1KkgbeRIVynmgl_tulWIJ4U4zvNCzVbA-trdRSYVPZCA/edit?gid=709301844#gid=709301844")
meta_img.nr   = read_sheet("https://docs.google.com/spreadsheets/d/1pjlOXDIDhWv8bujJ4P9BdYC-sSA_kueDAi9T8FLgpm8/edit?gid=67047753#gid=67047753")

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
data_combined   = data_combined[!(is.na(data_combined$avg_coverage)),]
data_combined   = data_combined[!(is.na(data_combined$aaeB)),]
data_combined   = apply(data_combined,2,as.character)
write.csv(data_combined, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/IMG_global_dataset.csv")
data_combined.1 = as.data.frame(data_combined) %>% select(2,44,69:1789)
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
hmm_loma    = read_sheet("https://docs.google.com/spreadsheets/d/1pjlOXDIDhWv8bujJ4P9BdYC-sSA_kueDAi9T8FLgpm8/edit?gid=67047753#gid=67047753")
gene_loma   = read_sheet("https://docs.google.com/spreadsheets/d/1j_GH0_JJANgAVytYBc6jPXaQXbiYUkh4SCfJswLJk7Y/edit?gid=349813226#gid=349813226")
mag_stat    = read_sheet("https://docs.google.com/spreadsheets/d/1uwpo3aUPodFvRzg_hPHH_mPNVq-gKUJOAEOakohL7ng/edit?gid=1062499715#gid=1062499715") 
mag_abun    = read_sheet("https://docs.google.com/spreadsheets/d/1IIXeOK2qYqIgU2AFi0Yedx-qkLMisexrLyUFxUuncrg/edit?gid=264725796#gid=264725796") 
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
hmm_fire    = read_sheet("https://docs.google.com/spreadsheets/d/1hudmGyDbKOVfpNNZ1mRTRK0edNPm9henzCvpZA0HJNc/edit?gid=1400821528#gid=1400821528")
mag_abun    = read_sheet("https://docs.google.com/spreadsheets/d/1F6dv4zx0vK83IzfEU02VC9_GJYVERmIV84uwhBV4RL4/edit?gid=1006296082#gid=1006296082") 
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

# FUNCTIONAL GROUP ANALYSIS----

# Call the MAG-Microtrait datasets
global_dataset                = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_global.csv",dec=".")
#global_LOMA                  = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_LOMA.csv",dec=".")
#global_FIRE                  = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_FIRE.csv",dec=".")
#global_HARVARD               = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_HARVARD.FOREST.csv",dec=".")
#global_MAG_EAA2017           = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_EAA2017.PLOT.csv",dec=".")

# Curve for functional groups (1. Random selection of matrix size;
# 2. Distance matrix and pairwise adonis; 3. Store N° of MAGs and N° of functional groups;
# 4. Plot and fit the curve)

MAG_number = c(100,200,500,600,750,800,900,1000,1500,2000,3500,4000,5000,
               6000,7000,8000,9000,10000)

FG_test = function(i,k,seed) {
  set.seed(seed)
  rand_df                = global_dataset[sample(nrow(global_dataset), size=i), ]
  distance.total         = parDist(x = as.matrix(rand_df[3:length(rand_df)]),
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

# i = 100 ----
adonis = FG_test(MAG_number[1],9,1)        # 9 FG
adonis = FG_test(MAG_number[1],11,11)      # 11 FG
adonis = FG_test(MAG_number[1],10,111)     # 10 FG
adonis = FG_test(MAG_number[1],14,2)       # 14 FG
adonis = FG_test(MAG_number[1],13,22)      # 13 FG
adonis = FG_test(MAG_number[1],14,222)     # 14 FG
adonis = FG_test(MAG_number[1],10,3)       # 10 FG
adonis = FG_test(MAG_number[1],10,33)      # 10 FG
adonis = FG_test(MAG_number[1],12,333)     # 12 FG
adonis = FG_test(MAG_number[1],11,27)      # 11 FG
adonis = FG_test(MAG_number[1],5,2707)     # 5 FG
adonis = FG_test(MAG_number[1],9,27071990) # 9 FG
adonis_100 = c(9,11,10,14,13,14,10,10,12,11,5,9)

# i = 200 ----
adonis = FG_test(MAG_number[2],18,1)        # 18 FG
adonis = FG_test(MAG_number[2],15,11)       # 15 FG
adonis = FG_test(MAG_number[2],21,111)      # 21 FG
adonis = FG_test(MAG_number[2],22,2)        # 22 FG
adonis = FG_test(MAG_number[2],22,22)       # 22 FG
adonis = FG_test(MAG_number[2],14,222)      # 14 FG
adonis = FG_test(MAG_number[2],17,3)        # 17 FG
adonis = FG_test(MAG_number[2],14,33)       # 14 FG
adonis = FG_test(MAG_number[2],23,333)      # 23 FG
adonis = FG_test(MAG_number[2],21,27)       # 21 FG
adonis = FG_test(MAG_number[2],17,2707)     # 17 FG
adonis = FG_test(MAG_number[2],21,27071990) # 21 FG
adonis_200 = c(18,15,21,22,22,14,17,14,23,21,14,21)

# i = 1000 ----
adonis = FG_test(MAG_number[8],60,1)        # 60 FG
adonis = FG_test(MAG_number[8],57,11)       # 57 FG
adonis = FG_test(MAG_number[8],65,111)      # 65 FG
adonis = FG_test(MAG_number[8],63,2)        # 63 FG
adonis = FG_test(MAG_number[8],72,22)       # 72 FG
adonis = FG_test(MAG_number[8],64,222)      # 64 FG
adonis = FG_test(MAG_number[8],60,3)        # 60 FG
adonis = FG_test(MAG_number[8],65,33)       # 65 FG
adonis = FG_test(MAG_number[8],52,333)      # 52 FG
adonis = FG_test(MAG_number[8],67,27)       # 67 FG
adonis = FG_test(MAG_number[8],95,2707)     # 95 FG
adonis = FG_test(MAG_number[8],56,27071990) # 56 FG
adonis_1000 = c(60,57,65,63,72,64,60,65,52,67,95,56)

# i = 750 ----
adonis = FG_test(MAG_number[5],48,1)        # 48 FG
adonis = FG_test(MAG_number[5],41,11)       # 41 FG
adonis = FG_test(MAG_number[5],57,111)      # 57 FG
adonis = FG_test(MAG_number[5],45,2)        # 45 FG
adonis = FG_test(MAG_number[5],65,22)       # 65 FG
adonis = FG_test(MAG_number[5],45,222)      # 45 FG
adonis = FG_test(MAG_number[5],42,3)        # 42 FG
adonis = FG_test(MAG_number[5],59,33)       # 59 FG
adonis = FG_test(MAG_number[5],42,333)      # 42 FG
adonis = FG_test(MAG_number[5],51,27)       # 51 FG
adonis = FG_test(MAG_number[5],63,2707)     # 63 FG
adonis = FG_test(MAG_number[5],55,27071990) # 55 FG
adonis_750 = c(48,41,57,45,65,45,42,59,42,41,63,55)

# i = 500 ----
adonis = FG_test(MAG_number[3],33,1)        # 33 FG
adonis = FG_test(MAG_number[3],33,11)       # 33 FG
adonis = FG_test(MAG_number[3],45,111)      # 45 FG
adonis = FG_test(MAG_number[3],41,2)        # 41 FG
adonis = FG_test(MAG_number[3],40,22)       # 40 FG
adonis = FG_test(MAG_number[3],29,222)      # 29 FG
adonis = FG_test(MAG_number[3],35,3)        # 35 FG
adonis = FG_test(MAG_number[3],47,33)       # 47 FG
adonis = FG_test(MAG_number[3],28,333)      # 28 FG
adonis = FG_test(MAG_number[3],36,27)       # 36 FG
adonis = FG_test(MAG_number[3],42,2707)     # 42 FG
adonis = FG_test(MAG_number[3],33,27071990) # 33 FG
adonis_500 = c(33,33,45,41,40,29,35,47,28,36,42,33)

# i = 1500 ----
adonis = FG_test(MAG_number[9],107,1)        # 107 FG
adonis = FG_test(MAG_number[9],86,11)        # 86 FG
adonis = FG_test(MAG_number[9],85,111)       # 85 FG
adonis = FG_test(MAG_number[9],87,2)         # 87 FG
adonis = FG_test(MAG_number[9],87,22)        # 87 FG
adonis = FG_test(MAG_number[9],89,222)       # 89 FG
adonis = FG_test(MAG_number[9],73,3)         # 73 FG
adonis = FG_test(MAG_number[9],79,33)        # 79 FG
adonis = FG_test(MAG_number[9],69,333)       # 69 FG
adonis = FG_test(MAG_number[9],98,27)        # 98 FG
adonis = FG_test(MAG_number[9],82,2707)      # 82 FG
adonis = FG_test(MAG_number[9],76,27071990)  # 76 FG
adonis_1500 = c(107,86,85,87,87,89,73,79,69,98,82,76)

# i = 2000 ----
adonis = FG_test(MAG_number[10],107,1)       # 107 FG
adonis = FG_test(MAG_number[10],119,11)      # 119 FG
adonis = FG_test(MAG_number[10],90,111)      # 90 FG
adonis = FG_test(MAG_number[10],117,2)       # 117 FG
adonis = FG_test(MAG_number[10],99,22)       # 99 FG
adonis = FG_test(MAG_number[10],101,222)     # 101 FG
adonis = FG_test(MAG_number[10],90,3)        # 90 FG
adonis = FG_test(MAG_number[10],123,33)      # 123 FG
adonis = FG_test(MAG_number[10],93,333)      # 93 FG
adonis = FG_test(MAG_number[10],115,27)      # 115 FG
adonis = FG_test(MAG_number[10],126,2707)    # 126 FG
adonis = FG_test(MAG_number[10],96,27071990) # 96 FG
adonis_2000 = c(107,119,90,117,99,101,90,123,93,115,126,96)

# i = 3500 ----
adonis = FG_test(MAG_number[11],181,1)        # 181 FG
adonis = FG_test(MAG_number[11],175,2707)     # 175 FG
adonis = FG_test(MAG_number[11],167,2)        # 167 FG
adonis = FG_test(MAG_number[11],154,27071990) # 154 FG
adonis = FG_test(MAG_number[11],141,3)        # 141 FG 
adonis = FG_test(MAG_number[11],148,222)      # 148 FG 
adonis = FG_test(MAG_number[11],159,22)       # 159 FG 
adonis = FG_test(MAG_number[11],157,11)       # 157 FG 
adonis = FG_test(MAG_number[11],164,33)       # 164 FG 
adonis = FG_test(MAG_number[11],185,27)       # 185 FG 
adonis = FG_test(MAG_number[11],192,111)      # 192 FG 
adonis = FG_test(MAG_number[11],140,333)      # 140 FG 
adonis_3500 = c(181,175,167,154,141,148,159,157,164,185,192,140)

# i = 5000 ----
adonis = FG_test(MAG_number[13],209,11)       # 209 FG
adonis = FG_test(MAG_number[13],229,2707)     # 229 FG
adonis = FG_test(MAG_number[13],205,27071990) # 205 FG
adonis = FG_test(MAG_number[13],205,33)       # 205 FG
adonis = FG_test(MAG_number[13],224,27)       # 224 FG
adonis = FG_test(MAG_number[13],190,222)      # 190 FG
adonis = FG_test(MAG_number[13],199,1)        # 199 FG
adonis = FG_test(MAG_number[13],177,3)        # 177 FG
adonis = FG_test(MAG_number[13],269,2)        # 269 FG
adonis = FG_test(MAG_number[13],179,333)      # 179 FG
adonis = FG_test(MAG_number[13],239,111)      # 239 FG
adonis = FG_test(MAG_number[13],245,22)       # 245 FG

adonis_5000 = c(209,229,205,205,224,190,199,177,269,239,245)

# Extrapolating the number of functional groups for the whole datasets ----

plot_file = as.data.frame(cbind(adonis_100,adonis_200,adonis_500,adonis_750,
                            adonis_1000,adonis_1500,adonis_2000,adonis_3500,
                            adonis_5000)) %>% 
  summarise_each(funs( mean( .,na.rm = TRUE)))

MAG_number.1 = as.data.frame(c(100,200,500,750,1000,1500,2000,3500,5000))
colnames(MAG_number.1) = "mag"
plot_file = as.data.frame(t(plot_file))

data = as.data.frame(cbind(MAG_number.1$mag,plot_file$V1))

model = lm(plot_file$V1 ~ MAG_number.1$mag)

Figure_S1 = ggplot(data,aes(V1,V2)) +
  geom_point() +
  geom_smooth(method='lm') + 
  stat_cor(label.x = 30, label.y = 130, size = 4) +
  stat_regline_equation(label.x = 30, label.y = 150, size = 4) + 
  xlab("# of MAGs") + ylab("# of Functional Groups") + theme_classic() + 
  theme(text = element_text(size=20))

png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_S1.png",
    width=3500,height=1969,res=300)
print(Figure_S1)
dev.off()

# Final functional groups for the total dataset----
FG = 16.597382 + 0.041440*32515 # Number of Functional groups for whole dataset

# Calling distance matrix
load("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/FG_parallel/distance.total.RData")
cluster.total   = hclust(distance.total, method="ward.D2")
v                      = cutree(cluster.total,k=round(FG))
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
global_dataset.1       = as.data.frame(cbind(genome2guild,global_dataset))
write.csv(global_dataset.1, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/global_datasets_guild.csv")

# Gene overlapping analysis----

# Genome size analysis----

# Call data again
# IMG-MAG
hmm_img     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/IMG_global_dataset.csv",dec=".")
hmm_img     = as.data.frame(hmm_img) %>% select(3,26,70:1790)
colnames(hmm_img)[1]  = "id"
colnames(hmm_img)[2]  = "genome.size"

# Loma-MAG
hmm_loma    = read_sheet("https://docs.google.com/spreadsheets/d/1pjlOXDIDhWv8bujJ4P9BdYC-sSA_kueDAi9T8FLgpm8/edit?gid=67047753#gid=67047753")
gene_loma   = read_sheet("https://docs.google.com/spreadsheets/d/1j_GH0_JJANgAVytYBc6jPXaQXbiYUkh4SCfJswLJk7Y/edit?gid=349813226#gid=349813226")
hmm_loma    = gene_loma %>% full_join(hmm_loma)
hmm_loma[3] = NULL
colnames(hmm_loma)[2]  = "genome.size"

# Wildfire-MAG
hmm_fire    = read_sheet("https://docs.google.com/spreadsheets/d/1hudmGyDbKOVfpNNZ1mRTRK0edNPm9henzCvpZA0HJNc/edit?gid=1400821528#gid=1400821528")
fire_stat   = read_sheet("https://docs.google.com/spreadsheets/d/1Ck_WNZbSZHJ5aVSSH7rdq9uNUYRqyQ5b-0RC7MtvdVk/edit?gid=963628746#gid=963628746")
fire_stat   = fire_stat[-c(440,546), ]
hmm_fire    = fire_stat %>% full_join(hmm_fire)
hmm_fire[3] = NULL
colnames(hmm_fire)[2]  = "genome.size"

# IMG-Isolates
isolates      = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/hmm_isolates.csv",dec=".")
df_1          = read_sheet("https://docs.google.com/spreadsheets/d/1QFSw5U0tf_G6LYYUmkCIkVAnCKaPehnIbYDvu0x4Bzc/edit?gid=1504418445#gid=1504418445")
df_1          = subset(df_1, select = -c(23))
df_2          = read_sheet("https://docs.google.com/spreadsheets/d/12sKqTlsMWOSDus66dfz3ao6M9IBK_iPnpfGPMCnCt1E/edit?gid=2065174436#gid=2065174436")
df_2          = subset(df_2, select = -c(14,16,18))
iso.metadata  = as.data.frame(rbind(df_1,df_2))
iso.metadata  = filter(iso.metadata, High.Quality %in% c("Yes"))
test          = as.data.frame(isolates[3:2300])
sums          = as.data.frame(rowSums(test, na.rm = FALSE, dims = 1))
temp          = cbind(sums,isolates$id,1:nrow(isolates))
isolates      = isolates[isolates$id %in% iso.metadata$taxon_oid, ]
iso.metadata  = iso.metadata[iso.metadata$taxon_oid %in% isolates$id, ]
colnames(iso.metadata)[1]  = "id"
isolates      = iso.metadata %>% full_join(isolates)
isolates      = as.data.frame(isolates) %>% select(1,23,26:2323)
colnames(isolates)[2]  = "genome.size"

# Erase temporary files
rm(gene_loma,fire_stat,df_1,df_2,iso.metadata,test,sums,temp)

# Bind the big matrix
hmm_img     = mutate(hmm_img, across(everything(), as.factor))
hmm_loma    = mutate(hmm_loma, across(everything(), as.factor))
hmm_fire    = mutate(hmm_fire, across(everything(), as.factor))
isolates    = mutate(isolates, across(everything(), as.factor))
total_genes = bind_rows(hmm_img,hmm_loma,hmm_fire,isolates)
total_genes[is.na(total_genes)] = 0

# Guild number
total_genes.guild     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/global_datasets_guild.csv",dec=".")
total.guild           = as.data.frame(total_genes.guild) %>% select(4,2)
total.g.size          = as.data.frame(total_genes) %>% select(1,2)

# Total guilds for the 940 selected genes
total_genes.guild     = left_join(total.g.size,total_genes.guild, by=c('id'))
total_genes.guild.940 = as.data.frame(total_genes.guild) %>% select(1,2,4,6:941)
write.csv(total_genes.guild.940, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.940.csv")

# Total guilds for the complete MAG-genes matrix
total_genes           = left_join(total.guild, total_genes, by=c('id'))
write.csv(total_genes, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.full.csv")

# 940 MAG-Gene matrix----
# Call Trait Keys 
GH_rule = read_sheet("https://docs.google.com/spreadsheets/d/1U6qWJJHossiK3kIwV2XYNKgsXkVrUc42QHiigMhX8Wg/edit?gid=1910000517#gid=1910000517")
GH_rule = as.data.frame(rbind("id","guild","genome.size",GH_rule))
PR_rule = read_sheet("https://docs.google.com/spreadsheets/d/1y8kqtT9mBdf-34wDUO9LN11kN0Ugz7vk9G2ZlZSc31Y/edit?gid=1394473486#gid=1394473486")
PR_rule = as.data.frame(rbind("id","guild","genome.size",PR_rule))
transp_rule = read_sheet("https://docs.google.com/spreadsheets/d/1NJGOHKHM8IpKEZNAs69R_XlftkZ_GrhzOhQvTZfPz_w/edit?gid=1971062019#gid=1971062019")
osmo_rule   = read_sheet("https://docs.google.com/spreadsheets/d/1WQ27I2Hd9jtCOv3z3cZnSB_A3xTt9bP5tZlOjiU4wVg/edit?gid=379499472#gid=379499472")
osmo_rule   = as.data.frame(rbind("id","guild","genome.size",osmo_rule))
biofilm_rule = read_sheet("https://docs.google.com/spreadsheets/d/1-FR1s9-txuPZWg8uJ21wCn8S0hYP3JDAZpF3Qmx93QU/edit?gid=1459756284#gid=1459756284")
biofilm_rule = as.data.frame(rbind("id","guild","genome.size",biofilm_rule))  
high.T_rule  = read_sheet("https://docs.google.com/spreadsheets/d/1-PeZ-F2RnXVFMa0Hkmg8Q-zlp38bCf4lxC-gpI9jEu4/edit?gid=123098806#gid=123098806")
high.T_rule  = as.data.frame(rbind("id","guild","genome.size",high.T_rule))
pH_rule      = read_sheet("https://docs.google.com/spreadsheets/d/1kANfYGvbb8tDiYEkJ_Whb9kocfhv8oDWpbNITMdZFQg/edit?gid=1661588961#gid=1661588961")
pH_rule      = as.data.frame(rbind("id","guild","genome.size",pH_rule)) 
  
# Genome-size for MAGs----
key_MAGs = bind_rows(hmm_img,hmm_loma,hmm_fire)
key_MAGs = as.data.frame(key_MAGs$id)
colnames(key_MAGs) = "key"
total_genes.guild.940_MAG = total_genes.guild.940[total_genes.guild.940$id %in% key_MAGs$key, ]
write.csv(total_genes.guild.940_MAG, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.940_MAG.csv")

total_genes.guild.940_MAG = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.940_MAG.csv",dec=".")

# Genome-size for Isolates----
key_isolates = as.data.frame(isolates$id)
colnames(key_isolates) = "key"
total_genes.guild.940_ISO = total_genes.guild.940[total_genes.guild.940$id %in% key_isolates$key, ]
write.csv(total_genes.guild.940_ISO, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.940_ISO.csv")

total_genes.guild.940_ISO = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.940_ISO.csv",dec=".")

# LINEAR REGRESSIONS----

# GH - genes (MAG)----
GH_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.MAG$genome.size = as.numeric(as.character(GH_TOTAL.MAG$genome.size))
GH_TOTAL_m   = GH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                               list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,3:ncol(GH_TOTAL_m)]))

ggplot(data = GH_TOTAL_m, aes(x = genome.size_mean, y = GH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("CAZy") +
  geom_point() + theme(text = element_text(size=20))

GH_TOTAL_m.1 = GH_TOTAL_m[-c(247), ] # Erasing super big Functional Group

sheet_write(GH_TOTAL_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1BvCCY4WmRk9Gaz01XGvhUgBpd6ZE8gHLzQSZsNv_J-4/edit?gid=0#gid=0",
            sheet = "GH_MAG")

GH_TOTAL_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1BvCCY4WmRk9Gaz01XGvhUgBpd6ZE8gHLzQSZsNv_J-4/edit?gid=0#gid=0")
MAG_gen_trait.2  = read_sheet("https://docs.google.com/spreadsheets/d/1Gtb4oLsNibF-yPwA9f1y78axGuKCO70Y1cZKNw0KRtk/edit?gid=0#gid=0")
MAG_gen_trait.3  = read_sheet("https://docs.google.com/spreadsheets/d/1Y6yMpQygXJof8EdDAE967ufLbNxmJcBqSZzUgMBnpYY/edit?gid=0#gid=0")

colnames(MAG_gen_trait.2) = c("Guild", "Genome-Size", "Amino-transporter", "pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","OGT",
                              "MGT")
colnames(MAG_gen_trait.3) = c("Guild", "Genome-Size", "Amino-transporter", "pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","OGT",
                              "MGT")
Total = as.data.frame(rbind(MAG_gen_trait.2,MAG_gen_trait.3))
Total = Total %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                           MGT  > 5 ~ "slow"))

Figure_test_1 = ggplot(data = GH_TOTAL_m.1, aes(x = genome.size_mean/1e6, y = (GH_total))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("CAZy") +
  geom_point() + theme_classic() + theme(text = element_text(size=16)) + 
  xlim(0,1.25e7/1e6) + ylim(0,30)
Figure_test_1

Figure_test_1.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (CAZy),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("CAZy (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,30) + theme(legend.position="none")
Figure_test_1.new

# Protein - genes (MAG)----
PR_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(PR_rule$`microtrait_rule-name`))
PR_TOTAL.MAG$genome.size = as.numeric(as.character(PR_TOTAL.MAG$genome.size))
PR_TOTAL_m   = PR_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PR_TOTAL_m   = PR_TOTAL_m %>% mutate(PR_total = rowSums(PR_TOTAL_m[,3:ncol(PR_TOTAL_m)]))

ggplot(data = PR_TOTAL_m, aes(x = genome.size_mean, y = PR_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total protein costs per FG") +
  geom_point() + theme(text = element_text(size=20))

PR_TOTAL_m.1 = PR_TOTAL_m[-c(247), ] # Erasing super big Functional Group

sheet_write(PR_TOTAL_m.1,
            ss = "https://docs.google.com/spreadsheets/d/18x3e2MPXYvhWrQ5HoGj_MmAXJfhHko-bP7TfDp9vCQM/edit?gid=0#gid=0",
            sheet = "Protein_MAG")

PR_TOTAL_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/18x3e2MPXYvhWrQ5HoGj_MmAXJfhHko-bP7TfDp9vCQM/edit?gid=0#gid=0")

Figure_test_2 = ggplot(data = PR_TOTAL_m.1, aes(x = genome.size_mean/1e6, y = PR_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Protein Enzymes") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,50)
Figure_test_2

Figure_test_2.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Protein-enzyme`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Protein Enzymes (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,50) + theme(legend.position="none")
Figure_test_2.new

# Transport - genes (MAG)----
transp_rule  = transp_rule %>% filter(function_gene == c("transporter"))
transp_rule  = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL = total_genes.guild.940_MAG %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL$genome.size = as.numeric(as.character(TRANSP_TOTAL$genome.size))
TRANSP_TOTAL_m   = TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m   = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,3:ncol(TRANSP_TOTAL_m)]))

ggplot(data = TRANSP_TOTAL_m, aes(x = genome.size_mean, y = transp_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total transporter costs per FG") +
  geom_point() + theme(text = element_text(size=20))

TRANSP_TOTAL_m.1 = TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group

sheet_write(TRANSP_TOTAL_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1SHPSWR5bB-31cVpd98uaJ0225PYDiPEXh1jeCeez-Dw/edit?gid=0#gid=0",
            sheet = "Transporter_MAG")

TRANSP_TOTAL_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1SHPSWR5bB-31cVpd98uaJ0225PYDiPEXh1jeCeez-Dw/edit?gid=0#gid=0")

Figure_test_3 = ggplot(data = TRANSP_TOTAL_m.1, aes(x = genome.size_mean/1e6, y = transp_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Total Transporter") +
  geom_point() + theme_classic()+ theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,200)
Figure_test_3

Figure_test_3.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Total-transporter`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Total Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,200) + theme(legend.position="none")
Figure_test_3.new

# Transport - Aminoacids (MAG)----
transp_rule        = transp_rule %>% filter(function_gene==c("transporter"))
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL   = total_genes.guild.940_MAG %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL$genome.size))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL_m[,3:ncol(AMI_TRANSP_TOTAL_m)]))

ggplot(data = AMI_TRANSP_TOTAL_m, aes(x = genome.size_mean, y = AMI_TRANSP_TOTAL)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total aminoacid transporter costs per FG") +
  geom_point() + theme(text = element_text(size=20))

AMI_TRANSP_TOTAL_m.1 = AMI_TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group

sheet_write(AMI_TRANSP_TOTAL_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1roWnESdGSQ3DiCGBogrfLBFzjwrLrPUHYEYqVNpNEsI/edit?gid=0#gid=0",
            sheet = "Aminoacids_T_MAG")

AMI_TRANSP_TOTAL_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1roWnESdGSQ3DiCGBogrfLBFzjwrLrPUHYEYqVNpNEsI/edit?gid=0#gid=0")

Figure_test_4 = ggplot(data = AMI_TRANSP_TOTAL_m.1, aes(x = genome.size_mean/1e6, y = (AMI_TRANSP_TOTAL))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Aminoacid Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20)
Figure_test_4

Figure_test_4.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Amino-transporter`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Aminoacid Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_4.new

# Transport - Carbohydrate (MAG)----
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL   = total_genes.guild.940_MAG %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL$genome.size))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL_m[,3:ncol(CAR_TRANSP_TOTAL_m)]))

ggplot(data = CAR_TRANSP_TOTAL_m, aes(x = genome.size_mean, y = CAR_TRANSP_TOTAL)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Carbohydrate Transport") +
  geom_point() + theme(text = element_text(size=14)) + xlim(0,1.25e7) + ylim(0,25)

CAR_TRANSP_TOTAL_m.1 = CAR_TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group

sheet_write(CAR_TRANSP_TOTAL_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1c08rBTbHN6G8mVN8ybqMWvsJQgZzQ8Z1B72I_AlF06M/edit?gid=0#gid=0",
            sheet = "GH_T_MAG")

CAR_TRANSP_TOTAL_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1c08rBTbHN6G8mVN8ybqMWvsJQgZzQ8Z1B72I_AlF06M/edit?gid=0#gid=0")

Figure_test_5 = ggplot(data = CAR_TRANSP_TOTAL_m.1, aes(x = genome.size_mean/1e6, y = (CAR_TRANSP_TOTAL))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Carbohydrate Transport") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,55)
Figure_test_5

Figure_test_5.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`GH-transporter`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Carbohydrate Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,55) + theme(legend.position="none")
Figure_test_5.new

# Osmolytes - genes (MAG)----
OSMO_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.MAG$genome.size = as.numeric(as.character(OSMO_TOTAL.MAG$genome.size))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.MAG_m[,3:ncol(OSMO_TOTAL.MAG_m)]))

ggplot(data = OSMO_TOTAL.MAG_m, aes(x = genome.size_mean, y = OSMO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

OSMO_TOTAL.MAG_m.1 = OSMO_TOTAL.MAG_m[-c(247), ] # Erasing super big Functional Group

sheet_write(OSMO_TOTAL.MAG_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1dzeG2XKS0fhznOMWXEa521XzC4YO5IUGfc4ShiUJ5yM/edit?gid=0#gid=0",
            sheet = "Osmolyte_MAG")

OSMO_TOTAL.MAG_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1dzeG2XKS0fhznOMWXEa521XzC4YO5IUGfc4ShiUJ5yM/edit?gid=0#gid=0")

Figure_test_6 = ggplot(data = OSMO_TOTAL.MAG_m.1, aes(x = genome.size_mean/1e6, y = OSMO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Osmolytes") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32)
Figure_test_6

Figure_test_6.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (Osmolyte),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Osmolytes (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32) + theme(legend.position="none")
Figure_test_6.new

# Biofilm - genes (MAG)----
BIO_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,3:ncol(BIO_TOTAL.MAG_m)]))

ggplot(data = BIO_TOTAL.MAG_m, aes(x = genome.size_mean, y = BIO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

BIO_TOTAL.MAG_m.1 = BIO_TOTAL.MAG_m[-c(247), ] # Erasing super big Functional Group

sheet_write(BIO_TOTAL.MAG_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1k32wnkrBHdG31a79sFvxVboZ-DwhWe28n9qk7H9lnOI/edit?gid=0#gid=0",
            sheet = "Biofilm_MAG")

BIO_TOTAL.MAG_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1k32wnkrBHdG31a79sFvxVboZ-DwhWe28n9qk7H9lnOI/edit?gid=0#gid=0")

Figure_test_7 = ggplot(data = BIO_TOTAL.MAG_m.1, aes(x = genome.size_mean/1e6, y = BIO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Biofilm") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,16)
Figure_test_7

Figure_test_7.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (Biofilm),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Biofilm (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,16) + theme(legend.position="none")
Figure_test_7.new

# High Temp - genes (MAG)----
TEMP.MAG = total_genes.guild.940_MAG %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP.MAG$genome.size = as.numeric(as.character(TEMP.MAG$genome.size))
TEMP.MAG_m   = TEMP.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
TEMP.MAG_m   = TEMP.MAG_m %>% mutate(TEMP_total = rowSums(TEMP.MAG_m[,3:ncol(TEMP.MAG_m)]))

ggplot(data = TEMP.MAG_m, aes(x = genome.size_mean, y = TEMP_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

TEMP.MAG_m.1 = TEMP.MAG_m[-c(247), ] # Erasing super big Functional Group

sheet_write(TEMP.MAG_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1WYvJIxHLaksDbD0thH3DWOVViQ7gC8KU9OEb9cVqp50/edit?gid=0#gid=0",
            sheet = "TEMP_MAG")

TEMP.MAG_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1WYvJIxHLaksDbD0thH3DWOVViQ7gC8KU9OEb9cVqp50/edit?gid=0#gid=0")

Figure_test_8 = ggplot(data = TEMP.MAG_m.1, aes(x = genome.size_mean/1e6, y = TEMP_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Heat Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,15)
Figure_test_8

Figure_test_8.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Temp-Tol`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Heat Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,15) + theme(legend.position="none")
Figure_test_8.new

# pH - genes (MAG)----
PH.MAG = total_genes.guild.940_MAG %>% select(any_of(pH_rule$`microtrait_hmm-name`))
PH.MAG$genome.size = as.numeric(as.character(PH.MAG$genome.size))
PH.MAG_m   = PH.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                 list(mean=mean), na.rm=TRUE))
PH.MAG_m   = PH.MAG_m %>% mutate(PH_total = rowSums(PH.MAG_m[,3:ncol(PH.MAG_m)]))

ggplot(data = PH.MAG_m, aes(x = genome.size_mean, y = PH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

PH.MAG_m.1 = PH.MAG_m[-c(247), ] # Erasing super big Functional Group

sheet_write(PH.MAG_m.1,
            ss = "https://docs.google.com/spreadsheets/d/1za5oqXIRLFYkg6eV-BhL6MVlftTkDv78jQCZoYOnYBQ/edit?gid=0#gid=0",
            sheet = "PH_MAG")

PH.MAG_m.1     = read_sheet("https://docs.google.com/spreadsheets/d/1za5oqXIRLFYkg6eV-BhL6MVlftTkDv78jQCZoYOnYBQ/edit?gid=0#gid=0")

Figure_test_9 = ggplot(data = PH.MAG_m.1, aes(x = genome.size_mean/1e6, y = PH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("pH Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20)
Figure_test_9

Figure_test_9.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`pH-Tol`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("pH Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_9.new

# GH - genes (Isolates)----
GH_TOTAL.ISO = total_genes.guild.940_ISO %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.ISO$genome.size = as.numeric(as.character(GH_TOTAL.ISO$genome.size))
GH_TOTAL_m   = GH_TOTAL.ISO %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,4:ncol(GH_TOTAL_m)]))

sheet_write(GH_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/110N6SnbJWrh3bLcS3rwL8wQsntuTUHAeqYqfB2_R7d0/edit?gid=0#gid=0",
            sheet = "GH_Isolates")

ISO_gen_trait.1  = read_sheet("https://docs.google.com/spreadsheets/d/1xT3vB2K2tupkeHZwUlVzh4ZYWhmUpOmMQtESNH2Mst0/edit?gid=0#gid=0",
                              sheet = "ISO_gen_trait.1")
GH_TOTAL_m       = read_sheet("https://docs.google.com/spreadsheets/d/110N6SnbJWrh3bLcS3rwL8wQsntuTUHAeqYqfB2_R7d0/edit?gid=0#gid=0")
Total.iso        = ISO_gen_trait.1 %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                           MGT  > 5 ~ "slow"))
# Erase outlier isolates that has a 0.99 yield value
Total.iso        = Total.iso %>% filter(yield <= 0.9)

Figure_test_1.ISO = ggplot(data = GH_TOTAL_m, aes(x = genome.size_mean/1e6, y = GH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CAZy") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,30)
Figure_test_1.ISO

Figure_test_1.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (CAZy),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("CAZy (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,30) + theme(legend.position="none")
Figure_test_1.ISO.new

# Protein - genes (Isolates)----
PH_TOTAL.ISO = total_genes.guild.940_ISO %>% select(any_of(PR_rule$`microtrait_rule-name`))
PH_TOTAL.ISO$genome.size = as.numeric(as.character(PH_TOTAL.ISO$genome.size))
PH_TOTAL_m   = PH_TOTAL.ISO %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                       list(mean=mean), na.rm=TRUE))
PH_TOTAL_m   = PH_TOTAL_m %>% mutate(PR_total = rowSums(PH_TOTAL_m[,4:ncol(PH_TOTAL_m)]))

sheet_write(PH_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1G3xtk7sHCy15e_0n2vLkbGc57fQFYyqCUghLDPQTuy0/edit?gid=0#gid=0",
            sheet = "Protein_Isolates")

PH_TOTAL_m       = read_sheet("https://docs.google.com/spreadsheets/d/1G3xtk7sHCy15e_0n2vLkbGc57fQFYyqCUghLDPQTuy0/edit?gid=0#gid=0")

Figure_test_2.ISO = ggplot(data = PH_TOTAL_m, aes(x = genome.size_mean/1e6, y = PR_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Protein Enzyme") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,50)
Figure_test_2.ISO

Figure_test_2.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (Protein),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Protein Enzyme (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,50) + theme(legend.position="none")
Figure_test_2.ISO.new

# Transport - genes (Isolates)----
transp_rule    = transp_rule %>% filter(function_gene == c("transporter"))
transp_rule    = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL.i$genome.size = as.numeric(as.character(TRANSP_TOTAL.i$genome.size))
TRANSP_TOTAL_m = TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,4:ncol(TRANSP_TOTAL_m)]))

sheet_write(TRANSP_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1pSIRE1LZknqQz5RSn_jZvnid2mZsns4_X6AiK7AHl8A/edit?gid=0#gid=0",
            sheet = "Transporter_Isolates")

TRANSP_TOTAL_m       = read_sheet("https://docs.google.com/spreadsheets/d/1pSIRE1LZknqQz5RSn_jZvnid2mZsns4_X6AiK7AHl8A/edit?gid=0#gid=0")

Figure_test_3.ISO = ggplot(data = TRANSP_TOTAL_m, aes(x = genome.size_mean/1e6, y = transp_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Total Transporter") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,200)
Figure_test_3.ISO

Figure_test_3.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (tranport.total),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Total Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,200) + theme(legend.position="none")
Figure_test_3.ISO.new

# Transport - Aminoacids (Isolates)----
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL.i$genome.size))
AMI_TRANSP_TOTAL.i_m = AMI_TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL.i_m = AMI_TRANSP_TOTAL.i_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL.i_m[,4:ncol(AMI_TRANSP_TOTAL.i_m)]))

sheet_write(AMI_TRANSP_TOTAL.i_m,
            ss = "https://docs.google.com/spreadsheets/d/1_-aWScfsDb-TNSFQ9bzbWVgEjdleCD_woBNRKHGy8o4/edit?gid=0#gid=0",
            sheet = "Aminoacids_T_Isolates")

AMI_TRANSP_TOTAL.i_m       = read_sheet("https://docs.google.com/spreadsheets/d/1_-aWScfsDb-TNSFQ9bzbWVgEjdleCD_woBNRKHGy8o4/edit?gid=0#gid=0")

Figure_test_4.ISO = ggplot(data = AMI_TRANSP_TOTAL.i_m, aes(x = genome.size_mean/1e6, y = AMI_TRANSP_TOTAL)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Aminoacid Transporter") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20)
Figure_test_4.ISO

Figure_test_4.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (amino.transport),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Aminoacid Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_4.ISO.new

# Transport - Carbohydrate (Isolates)----
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL.i$genome.size))
CAR_TRANSP_TOTAL.i_m = CAR_TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL.i_m = CAR_TRANSP_TOTAL.i_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL.i_m[,4:ncol(CAR_TRANSP_TOTAL.i_m)]))

sheet_write(CAR_TRANSP_TOTAL.i_m,
            ss = "https://docs.google.com/spreadsheets/d/14RSIx5VpsXPvxV3Awrj_XRwAc1uYx5tdpjK1vjJGJAo/edit?gid=0#gid=0",
            sheet = "GH_T_Isolates")

CAR_TRANSP_TOTAL.i_m       = read_sheet("https://docs.google.com/spreadsheets/d/14RSIx5VpsXPvxV3Awrj_XRwAc1uYx5tdpjK1vjJGJAo/edit?gid=0#gid=0")

Figure_test_5.ISO = ggplot(data = CAR_TRANSP_TOTAL.i_m, aes(x = genome.size_mean/1e6, y = CAR_TRANSP_TOTAL)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Carbohydrate Transporter") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,55)
Figure_test_5.ISO

Figure_test_5.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (GH.trasnport),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Carbohydrate Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,55) + theme(legend.position="none")
Figure_test_5.ISO.new

# Osmolytes - genes (Isolates)----
OSMO_TOTAL.TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.TOTAL.i$genome.size = as.numeric(as.character(OSMO_TOTAL.TOTAL.i$genome.size))
OSMO_TOTAL.TOTAL.i_m   = OSMO_TOTAL.TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
OSMO_TOTAL_m   = OSMO_TOTAL.TOTAL.i_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.TOTAL.i_m[,4:ncol(OSMO_TOTAL.TOTAL.i_m)]))

sheet_write(OSMO_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1FsE_lt-gHphdH1Ndcce7IODTXyBVwu6uckVY97ipovo/edit?gid=0#gid=0",
            sheet = "Osmolyte_Isolates")

OSMO_TOTAL_m       = read_sheet("https://docs.google.com/spreadsheets/d/1FsE_lt-gHphdH1Ndcce7IODTXyBVwu6uckVY97ipovo/edit?gid=0#gid=0")

Figure_test_6.ISO = ggplot(data = OSMO_TOTAL_m, aes(x = genome.size_mean/1e6, y = OSMO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Osmolytes") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32)
Figure_test_6.ISO

Figure_test_6.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (osmolyte),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Osmolytes (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32) + theme(legend.position="none")
Figure_test_6.ISO.new

# Biofilm - genes (Isolates)----
BIO_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,4:ncol(BIO_TOTAL.MAG_m)]))

sheet_write(BIO_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1m5BCkSa5CoWY7sLeJl0K16Ef1P8lPuli7HG6KpSxZPE/edit?gid=0#gid=0",
            sheet = "Biofilm_Isolates")

BIO_TOTAL.MAG_m       = read_sheet("https://docs.google.com/spreadsheets/d/1m5BCkSa5CoWY7sLeJl0K16Ef1P8lPuli7HG6KpSxZPE/edit?gid=0#gid=0")

Figure_test_7.ISO = ggplot(data = BIO_TOTAL.MAG_m, aes(x = genome.size_mean/1e6, y = BIO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Biofilm") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,16)
Figure_test_7.ISO

Figure_test_7.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (biofilm),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Biofilm (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32) + theme(legend.position="none")
Figure_test_7.ISO.new

# High Temperature - genes (Isolates)----
TEMP_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP_TOTAL.MAG$genome.size = as.numeric(as.character(TEMP_TOTAL.MAG$genome.size))
TEMP_TOTAL.MAG_m   = TEMP_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
TEMP_TOTAL.MAG_m   = TEMP_TOTAL.MAG_m %>% mutate(TEMP_total = rowSums(TEMP_TOTAL.MAG_m[,4:ncol(TEMP_TOTAL.MAG_m)]))

sheet_write(TEMP_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1GXAwzFiMHrZAkIxYqDxGbj0VHzB8SHWxF8mnGqU3jGo/edit?gid=0#gid=0",
            sheet = "TEMP_Isolates")

TEMP_TOTAL.MAG_m       = read_sheet("https://docs.google.com/spreadsheets/d/1GXAwzFiMHrZAkIxYqDxGbj0VHzB8SHWxF8mnGqU3jGo/edit?gid=0#gid=0")

Figure_test_8.ISO = ggplot(data = TEMP_TOTAL.MAG_m, aes(x = genome.size_mean/1e6, y = TEMP_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Heat Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,15)
Figure_test_8.ISO

Figure_test_8.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (temp),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Heat Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32) + theme(legend.position="none")
Figure_test_8.ISO.new

# pH - genes (Isolates)----
pH_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(pH_rule$`microtrait_hmm-name`))
pH_TOTAL.MAG$genome.size = as.numeric(as.character(pH_TOTAL.MAG$genome.size))
pH_TOTAL.MAG_m   = pH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
pH_TOTAL.MAG_m   = pH_TOTAL.MAG_m %>% mutate(PH_total = rowSums(pH_TOTAL.MAG_m[,4:ncol(pH_TOTAL.MAG_m)]))

sheet_write(pH_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1p58fzYUw9I-jA4FkWQz6u8-qtNVi7GMVVcGCtYPdxHk/edit?gid=0#gid=0",
            sheet = "PH_Isolates")

pH_TOTAL.MAG_m       = read_sheet("https://docs.google.com/spreadsheets/d/1p58fzYUw9I-jA4FkWQz6u8-qtNVi7GMVVcGCtYPdxHk/edit?gid=0#gid=0")

Figure_test_9.ISO = ggplot(data = pH_TOTAL.MAG_m, aes(x = genome.size_mean/1e6, y = PH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("pH Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20)
Figure_test_9.ISO

Figure_test_9.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (pH),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("pH Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_9.ISO.new

# Data Preparation - Isolates
temp    = drive_ls("https://drive.google.com/drive/u/0/folders/1C9gqH5mSadOGdZv05jKPu8tELph2isd4")
a       = as.data.frame(temp$name)
PH_Isolates       = read_sheet(temp$id[1])
PH_Isolates       = PH_Isolates %>% select(guild,genome.size_mean,PH_total)
TEMP_Isolates     = read_sheet(temp$id[2])
TEMP_Isolates     = TEMP_Isolates %>% select(guild,genome.size_mean,TEMP_total)
Biofilm_Isolates  = read_sheet(temp$id[3])
Biofilm_Isolates  = Biofilm_Isolates %>% select(guild,genome.size_mean,BIO_total)
Osmolyte_Isolates = read_sheet(temp$id[4])
Osmolyte_Isolates = Osmolyte_Isolates %>% select(guild,genome.size_mean,OSMO_total)
GH_T_Isolates     = read_sheet(temp$id[5])
GH_T_Isolates     = GH_T_Isolates %>% select(guild,genome.size_mean,CAR_TRANSP_TOTAL)
Aminoacids_T_Isolates = read_sheet(temp$id[6])
Aminoacids_T_Isolates = Aminoacids_T_Isolates %>% select(guild,genome.size_mean,AMI_TRANSP_TOTAL)
Transporter_Isolates  = read_sheet(temp$id[7])
Transporter_Isolates  = Transporter_Isolates %>% select(guild,genome.size_mean,transp_total)
Protein_Isolates      = read_sheet(temp$id[8])
Protein_Isolates      = Protein_Isolates %>% select(guild,genome.size_mean,PR_total)
GH_Isolates           = read_sheet(temp$id[9])
GH_Isolates           = GH_Isolates %>% select(guild,genome.size_mean,GH_total)
Isolates_gen_trait    = as.data.frame(cbind(Aminoacids_T_Isolates,PH_Isolates$PH_total,
                                            TEMP_Isolates$TEMP_total,Biofilm_Isolates$BIO_total,
                                            Osmolyte_Isolates$OSMO_total,GH_T_Isolates$CAR_TRANSP_TOTAL,
                                            Transporter_Isolates$transp_total,Protein_Isolates$PR_total,
                                            GH_Isolates$GH_total))
colnames(Isolates_gen_trait) = c("guild","genome.size","amino.transport","pH","temp",
                            "biofilm","osmolyte","GH.trasnport","tranport.total",
                            "Protein","CAZy")
rm(temp,a,Aminoacids_T_Isolates,PH_Isolates,TEMP_Isolates,Biofilm_Isolates,
   Osmolyte_Isolates,GH_T_Isolates,Transporter_Isolates,Protein_Isolates,
   GH_Isolates)

sheet_write(Isolates_gen_trait,
            ss = "https://docs.google.com/spreadsheets/d/1xT3vB2K2tupkeHZwUlVzh4ZYWhmUpOmMQtESNH2Mst0/edit?gid=0#gid=0",
            sheet = "Isolates_gen_trait")

# LOMA GENES ANALYSIS----

global_LOMA                  = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_LOMA.csv",dec=".")

# Curve for functional groups (1. Random selection of matrix size;
# 2. Distance matrix and pairwise adonis; 3. Store N° of MAGs and N° of functional groups;
# 4. Plot and fit the curve)

MAG_number = c(100,200,500,600,750,800,900,1000,1500,2000,3500,4000,5000,
               6000,7000,8000,9000,10000)

FG_test = function(i,k,seed) {
  set.seed(seed)
  rand_df                = global_LOMA[sample(nrow(global_LOMA), size=i), ]
  distance.total         = parDist(x = as.matrix(rand_df[3:length(rand_df)]),
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

# i = 100 ----
adonis = FG_test(MAG_number[1],7,1)        # 7 FG
adonis = FG_test(MAG_number[1],8,11)       # 8 FG
adonis = FG_test(MAG_number[1],10,111)     # 10 FG
adonis = FG_test(MAG_number[1],12,2)       # 12 FG
adonis = FG_test(MAG_number[1],10,22)      # 10 FG
adonis = FG_test(MAG_number[1],11,222)     # 11 FG
adonis = FG_test(MAG_number[1],10,3)       # 10 FG
adonis = FG_test(MAG_number[1],12,33)      # 12 FG
adonis = FG_test(MAG_number[1],15,333)     # 15 FG
adonis = FG_test(MAG_number[1],11,27)      # 11 FG
adonis = FG_test(MAG_number[1],7,2707)     # 7 FG
adonis = FG_test(MAG_number[1],8,27071990) # 8 FG
adonis_100 = c(7,8,10,12,10,11,10,12,15,11,7,8)

# i = 1000 ----
adonis = FG_test(MAG_number[8],63,1)        # 63 FG
adonis = FG_test(MAG_number[8],52,11)       # 52 FG
adonis = FG_test(MAG_number[8],59,111)      # 59 FG
adonis = FG_test(MAG_number[8],69,2)        # 69 FG
adonis = FG_test(MAG_number[8],59,22)       # 59 FG
adonis = FG_test(MAG_number[8],55,222)      # 55 FG
adonis = FG_test(MAG_number[8],60,3)        # 60 FG
adonis = FG_test(MAG_number[8],59,33)       # 59 FG
adonis = FG_test(MAG_number[8],45,333)      # 45 FG
adonis = FG_test(MAG_number[8],62,27)       # 62 FG
adonis = FG_test(MAG_number[8],44,2707)     # 44 FG
adonis = FG_test(MAG_number[8],70,27071990) # 70 FG
adonis_1000 = c(63,52,59,69,59,55,60,59,45,62,44,70)

# i = 2000 ----
adonis = FG_test(MAG_number[10],101,1)       # 101 FG
adonis = FG_test(MAG_number[10],123,11)      # 123 FG 
adonis = FG_test(MAG_number[10],83,111)      # 83 FG
adonis = FG_test(MAG_number[10],107,2)       # 107 FG
adonis = FG_test(MAG_number[10],90,22)       # 90 FG
adonis = FG_test(MAG_number[10],120,222)     # 120 FG
adonis = FG_test(MAG_number[10],96,3)        # 96 FG
adonis = FG_test(MAG_number[10],131,33)      # 131 FG
adonis = FG_test(MAG_number[10],118,333)     # 118 FG
adonis = FG_test(MAG_number[10],112,27)      # 112 FG
adonis = FG_test(MAG_number[10],138,2707)    # 138 FG
adonis = FG_test(MAG_number[10],82,27071990) # 82 FG
adonis_2000 = c(101,123,83,107,90,120,96,131,118,112,138,82)

# i = 500 ----
adonis = FG_test(MAG_number[3],37,1)        # 37 FG
adonis = FG_test(MAG_number[3],24,11)       # 24 FG
adonis = FG_test(MAG_number[3],34,111)      # 34 FG
adonis = FG_test(MAG_number[3],28,2)        # 28 FG
adonis = FG_test(MAG_number[3],36,22)       # 36 FG
adonis = FG_test(MAG_number[3],28,222)      # 28 FG
adonis = FG_test(MAG_number[3],33,3)        # 33 FG
adonis = FG_test(MAG_number[3],43,33)       # 43 FG
adonis = FG_test(MAG_number[3],25,333)      # 28 FG
adonis = FG_test(MAG_number[3],34,27)       # 34 FG
adonis = FG_test(MAG_number[3],37,2707)     # 37 FG
adonis = FG_test(MAG_number[3],43,27071990) # 43 FG
adonis_500 = c(37,24,34,28,36,28,33,43,28,34,37,43)

# i = 5000 ----
adonis = FG_test(MAG_number[13],207,27071990) # 207 FG
adonis = FG_test(MAG_number[13],172,33)       # 172 FG
adonis = FG_test(MAG_number[13],190,333)      # 190 FG
adonis = FG_test(MAG_number[13],223,22)       # 223 FG
adonis = FG_test(MAG_number[13],212,27)       # 212 FG
adonis = FG_test(MAG_number[13],229,11)       # 229 FG
adonis = FG_test(MAG_number[13],222,2)        # 222 FG
adonis = FG_test(MAG_number[13],228,1)        # 228 FG
adonis = FG_test(MAG_number[13],231,222)      # 231 FG
adonis = FG_test(MAG_number[13],245,3)        # 245 FG
adonis = FG_test(MAG_number[13],259,2707)     # 259 FG
adonis = FG_test(MAG_number[13],276,111)      # 276 FG
adonis_5000 = c(207,172,190,223,212,229,222,228,231,245,259,276)

# Extrapolating the number of functional groups for the whole datasets ----

plot_file = as.data.frame(cbind(adonis_100,adonis_500,adonis_1000,adonis_2000,
                                adonis_5000)) %>% 
  summarise_each(funs( mean( .,na.rm = TRUE)))
rm(adonis_100,adonis_500,adonis_1000,adonis_2000,
   adonis_5000)
MAG_number.1 = as.data.frame(c(100,500,1000,2000,5000))
colnames(MAG_number.1) = "mag"
plot_file = as.data.frame(t(plot_file))

data = as.data.frame(cbind(MAG_number.1$mag,plot_file$V1))

model = lm(plot_file$V1 ~ MAG_number.1$mag)

Figure_S2 = ggplot(data,aes(V1,V2)) +
  geom_point() +
  geom_smooth(method='lm') + 
  stat_cor(label.x = 30, label.y = 130, size = 4) +
  stat_regline_equation(label.x = 30, label.y = 150, size = 4) + 
  xlab("# of MAGs") + ylab("# of Functional Groups") + theme_classic() + 
  theme(text = element_text(size=20))

png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_S2.png",
    width=3500,height=1969,res=300)
print(Figure_S2)
dev.off()

# Final functional groups for the total dataset----
FG = 13 + 0.043*32515 # Number of Functional groups for whole dataset

# Distance matrix
distance.total  = parDist(x = as.matrix(global_LOMA[,3:183]),
                          method = "fJaccard",
                          threads = 1) # Adapt the number of threads
cluster.total   = hclust(distance.total, method="ward.D2")
v                      = cutree(cluster.total,k=round(FG))
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
global_dataset.1       = as.data.frame(cbind(genome2guild,global_LOMA))
write.csv(global_dataset.1, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/LOMA_datasets_guild.csv")

# Genome size analysis----

# Call data again
# IMG-MAG
hmm_img     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/IMG_global_dataset.csv",dec=".")
hmm_img     = as.data.frame(hmm_img) %>% select(3,26,70:1790)
colnames(hmm_img)[1]  = "id"
colnames(hmm_img)[2]  = "genome.size"

# Loma-MAG
hmm_loma    = read_sheet("https://docs.google.com/spreadsheets/d/1ixRY5OwvfrPk5MIc1LdZPHCKj5YYSqkBSA-OOo9pnoo/edit?gid=601691898#gid=601691898")
gene_loma   = read_sheet("https://docs.google.com/spreadsheets/d/1j_GH0_JJANgAVytYBc6jPXaQXbiYUkh4SCfJswLJk7Y/edit?gid=349813226#gid=349813226")
hmm_loma    = gene_loma %>% full_join(hmm_loma)
hmm_loma[3] = NULL
colnames(hmm_loma)[2]  = "genome.size"

# Wildfire-MAG
hmm_fire    = read_sheet("https://docs.google.com/spreadsheets/d/1hudmGyDbKOVfpNNZ1mRTRK0edNPm9henzCvpZA0HJNc/edit?gid=1400821528#gid=1400821528")
fire_stat   = read_sheet("https://docs.google.com/spreadsheets/d/1Ck_WNZbSZHJ5aVSSH7rdq9uNUYRqyQ5b-0RC7MtvdVk/edit?gid=963628746#gid=963628746")
fire_stat   = fire_stat[-c(440,546), ]
hmm_fire    = fire_stat %>% full_join(hmm_fire)
hmm_fire[3] = NULL
colnames(hmm_fire)[2]  = "genome.size"

# IMG-Isolates
isolates      = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/hmm_isolates.csv",dec=".")
df_1          = read_sheet("https://docs.google.com/spreadsheets/d/1QFSw5U0tf_G6LYYUmkCIkVAnCKaPehnIbYDvu0x4Bzc/edit?gid=1504418445#gid=1504418445")
df_1          = subset(df_1, select = -c(23))
df_2          = read_sheet("https://docs.google.com/spreadsheets/d/12sKqTlsMWOSDus66dfz3ao6M9IBK_iPnpfGPMCnCt1E/edit?gid=2065174436#gid=2065174436")
df_2          = subset(df_2, select = -c(14,16,18))
iso.metadata  = as.data.frame(rbind(df_1,df_2))
iso.metadata  = filter(iso.metadata, High.Quality %in% c("Yes"))
test          = as.data.frame(isolates[3:2300])
sums          = as.data.frame(rowSums(test, na.rm = FALSE, dims = 1))
temp          = cbind(sums,isolates$id,1:nrow(isolates))
isolates      = isolates[isolates$id %in% iso.metadata$taxon_oid, ]
iso.metadata  = iso.metadata[iso.metadata$taxon_oid %in% isolates$id, ]
colnames(iso.metadata)[1]  = "id"
isolates      = iso.metadata %>% full_join(isolates)
isolates      = as.data.frame(isolates) %>% select(1,23,26:2323)
colnames(isolates)[2]  = "genome.size"

# Erase temporary files
rm(gene_loma,fire_stat,df_1,df_2,iso.metadata,test,sums,temp)

# Bind the big matrix
hmm_img     = mutate(hmm_img, across(everything(), as.factor))
hmm_loma    = mutate(hmm_loma, across(everything(), as.factor))
hmm_fire    = mutate(hmm_fire, across(everything(), as.factor))
isolates    = mutate(isolates, across(everything(), as.factor))
total_genes = bind_rows(hmm_img,hmm_loma,hmm_fire,isolates)
total_genes[is.na(total_genes)] = 0

# Guild number
total_genes.guild     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/LOMA_datasets_guild.csv",dec=".")
total.guild           = as.data.frame(total_genes.guild) %>% select(4,2)
total.g.size          = as.data.frame(total_genes) %>% select(1,2)

# Total guilds for the 940 selected genes
total_genes.guild     = left_join(total.g.size,total_genes.guild, by=c('id'))
total_genes.guild.186 = as.data.frame(total_genes.guild) %>% select(1,2,4,6:186)
write.csv(total_genes.guild.186, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.LOMA.csv")

# Total guilds for the complete MAG-genes matrix
total_genes           = left_join(total.guild, total_genes, by=c('id'))
write.csv(total_genes, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.LOMA.csv")

# 186 MAG-Gene matrix----
# Call Trait Keys 
GH_rule = read_sheet("https://docs.google.com/spreadsheets/d/1U6qWJJHossiK3kIwV2XYNKgsXkVrUc42QHiigMhX8Wg/edit?gid=1910000517#gid=1910000517")
GH_rule = as.data.frame(rbind("id","guild","genome.size",GH_rule))
PR_rule = read_sheet("https://docs.google.com/spreadsheets/d/1y8kqtT9mBdf-34wDUO9LN11kN0Ugz7vk9G2ZlZSc31Y/edit?gid=1394473486#gid=1394473486")
PR_rule = as.data.frame(rbind("id","guild","genome.size",PR_rule))
transp_rule = read_sheet("https://docs.google.com/spreadsheets/d/1NJGOHKHM8IpKEZNAs69R_XlftkZ_GrhzOhQvTZfPz_w/edit?gid=1971062019#gid=1971062019")
osmo_rule   = read_sheet("https://docs.google.com/spreadsheets/d/1WQ27I2Hd9jtCOv3z3cZnSB_A3xTt9bP5tZlOjiU4wVg/edit?gid=379499472#gid=379499472")
osmo_rule   = as.data.frame(rbind("id","guild","genome.size",osmo_rule))
biofilm_rule = read_sheet("https://docs.google.com/spreadsheets/d/1-FR1s9-txuPZWg8uJ21wCn8S0hYP3JDAZpF3Qmx93QU/edit?gid=1459756284#gid=1459756284")
biofilm_rule = as.data.frame(rbind("id","guild","genome.size",biofilm_rule))  
high.T_rule  = read_sheet("https://docs.google.com/spreadsheets/d/1-PeZ-F2RnXVFMa0Hkmg8Q-zlp38bCf4lxC-gpI9jEu4/edit?gid=123098806#gid=123098806")
high.T_rule  = as.data.frame(rbind("id","guild","genome.size",high.T_rule))
pH_rule      = read_sheet("https://docs.google.com/spreadsheets/d/1kANfYGvbb8tDiYEkJ_Whb9kocfhv8oDWpbNITMdZFQg/edit?gid=1661588961#gid=1661588961")
pH_rule      = as.data.frame(rbind("id","guild","genome.size",pH_rule)) 

# Genome-size for MAGs----
key_MAGs = bind_rows(hmm_img,hmm_loma,hmm_fire)
key_MAGs = as.data.frame(key_MAGs$id)
colnames(key_MAGs) = "key"
total_genes.guild.186_MAG = total_genes.guild.186[total_genes.guild.186$id %in% key_MAGs$key, ]
write.csv(total_genes.guild.186, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.186_MAG.csv")

total_genes.guild.186_MAG = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.186_MAG.csv",dec=".")

# Genome-size for Isolates----
key_isolates = as.data.frame(isolates$id)
colnames(key_isolates) = "key"
total_genes.guild.186_ISO = total_genes.guild.186[total_genes.guild.186$id %in% key_isolates$key, ]
write.csv(total_genes.guild.186_ISO, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.186_ISO.csv")

total_genes.guild.186_ISO = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.186_ISO.csv",dec=".")

# LINEAR REGRESSIONS----

# GH - genes (MAG)----
GH_TOTAL.MAG = total_genes.guild.186_MAG %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.MAG$genome.size = as.numeric(as.character(GH_TOTAL.MAG$genome.size))
GH_TOTAL_m   = GH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,3:ncol(GH_TOTAL_m)]))

sheet_write(GH_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1jUgB0CcVwmFekcEL5K2Ebr-jHAI-hm7fqiXi6vBmXec/edit?gid=0#gid=0",
            sheet = "GH_MAG")

GH_TOTAL_m       = read_sheet("https://docs.google.com/spreadsheets/d/1jUgB0CcVwmFekcEL5K2Ebr-jHAI-hm7fqiXi6vBmXec/edit?gid=0#gid=0")
MAG_gen_trait.2  = read_sheet("https://docs.google.com/spreadsheets/d/1Gtb4oLsNibF-yPwA9f1y78axGuKCO70Y1cZKNw0KRtk/edit?gid=0#gid=0")
MAG_gen_trait.3  = read_sheet("https://docs.google.com/spreadsheets/d/1Y6yMpQygXJof8EdDAE967ufLbNxmJcBqSZzUgMBnpYY/edit?gid=0#gid=0")

colnames(MAG_gen_trait.2) = c("Guild", "Genome-Size", "Amino-transporter", "pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","OGT",
                              "MGT")
colnames(MAG_gen_trait.3) = c("Guild", "Genome-Size", "Amino-transporter", "pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","OGT",
                              "MGT")
Total = as.data.frame(rbind(MAG_gen_trait.2,MAG_gen_trait.3))
Total = Total %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                           MGT  > 5 ~ "slow")) 

Figure_test_1 = ggplot(data = GH_TOTAL_m, aes(x = genome.size_mean/1e6, y = (GH_total))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("CAZy") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20)
Figure_test_1

Figure_test_1.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (CAZy),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("CAZy") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.0e7/1e6) + ylim(0,30) + theme(legend.position="none")
Figure_test_1.new

# Protein - genes (MAG)----
PR_TOTAL.MAG = total_genes.guild.186_MAG %>% select(any_of(PR_rule$`microtrait_rule-name`))
PR_TOTAL.MAG$genome.size = as.numeric(as.character(PR_TOTAL.MAG$genome.size))
PR_TOTAL_m   = PR_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PR_TOTAL_m   = PR_TOTAL_m %>% mutate(PR_total = rowSums(PR_TOTAL_m[,3:ncol(PR_TOTAL_m)]))

sheet_write(PR_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1PEHIyhkOi-hxLbeMDbtD_12Qq2hVmFUyNLF4GASYVRQ/edit?gid=0#gid=0",
            sheet = "Protein_MAG")

PR_TOTAL_m    = read_sheet("https://docs.google.com/spreadsheets/d/1PEHIyhkOi-hxLbeMDbtD_12Qq2hVmFUyNLF4GASYVRQ/edit?gid=0#gid=0")
  
Figure_test_2 = ggplot(data = PR_TOTAL_m, aes(x = genome.size_mean/1e6, y = PR_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Protein Enzymes") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20)
Figure_test_2

Figure_test_2.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Protein-enzyme`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("CAZy") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.0e7/1e6) + ylim(0,50) + theme(legend.position="none")
Figure_test_2.new

# Transport - genes (MAG)----
transp_rule  = transp_rule %>% filter(function_gene == c("transporter"))
transp_rule  = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL = total_genes.guild.186_MAG %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL$genome.size = as.numeric(as.character(TRANSP_TOTAL$genome.size))
TRANSP_TOTAL_m   = TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m   = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,3:ncol(TRANSP_TOTAL_m)]))

sheet_write(TRANSP_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1frAIOBIQewDGxpUFJpdcP_5pAcyYfK8idKr2eHk3-ns/edit?gid=0#gid=0",
            sheet = "Transporter_MAG")

TRANSP_TOTAL_m    = read_sheet("https://docs.google.com/spreadsheets/d/1frAIOBIQewDGxpUFJpdcP_5pAcyYfK8idKr2eHk3-ns/edit?gid=0#gid=0")

Figure_test_3 = ggplot(data = TRANSP_TOTAL_m, aes(x = genome.size_mean/1e6, y = transp_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Total Transporter") +
  geom_point() + theme_classic()+ theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,30)
Figure_test_3

Figure_test_3.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Total-transporter`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("CAZy") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.0e7/1e6) + ylim(0,100) + theme(legend.position="none")
Figure_test_3.new

# Transport - Aminoacids (MAG)----
transp_rule        = transp_rule %>% filter(function_gene==c("transporter"))
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL   = total_genes.guild.186_MAG %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL$genome.size))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL_m[,3:ncol(AMI_TRANSP_TOTAL_m)]))

sheet_write(AMI_TRANSP_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1IKW-R9_T9H81P9BdBXsO1tymK1CYmI_vj7CDptw4LCc/edit?gid=0#gid=0",
            sheet = "AMI_TRANSP_TOTAL_m")

Figure_test_4 = ggplot(data = AMI_TRANSP_TOTAL_m, aes(x = genome.size_mean, y = (AMI_TRANSP_TOTAL))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Aminoacid Transporter") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,5)

# Transport - Carbohydrate (MAG)----
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL   = total_genes.guild.186_MAG %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL$genome.size))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL_m[,3:ncol(CAR_TRANSP_TOTAL_m)]))

ggplot(data = CAR_TRANSP_TOTAL_m, aes(x = genome.size_mean, y = CAR_TRANSP_TOTAL)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Carbohydrate Transport") +
  geom_point() + theme(text = element_text(size=14)) 

sheet_write(CAR_TRANSP_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1YwP4g8R3mi4TkTJI6p8DUm8QgwB4dVbX3i_AL1wNB5w/edit?gid=0#gid=0",
            sheet = "GH_T_MAG")

Figure_test_5 = ggplot(data = CAR_TRANSP_TOTAL_m, aes(x = genome.size_mean, y = (CAR_TRANSP_TOTAL))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Carbohydrate Transport") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,15)

# Osmolytes - genes (MAG)----
OSMO_TOTAL.MAG = total_genes.guild.186_MAG %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.MAG$genome.size = as.numeric(as.character(OSMO_TOTAL.MAG$genome.size))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.MAG_m[,3:ncol(OSMO_TOTAL.MAG_m)]))

ggplot(data = OSMO_TOTAL.MAG_m, aes(x = genome.size_mean, y = OSMO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

sheet_write(OSMO_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/16P4PkdWsxXiL1KUXQ7SPfV0-IOxkncvgmVfL1APTgJE/edit?gid=0#gid=0",
            sheet = "Osmolyte_MAG")

Figure_test_6 = ggplot(data = OSMO_TOTAL.MAG_m, aes(x = genome.size_mean, y = OSMO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Osmolytes") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,15)

# Biofilm - genes (MAG)----
BIO_TOTAL.MAG = total_genes.guild.186_MAG %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,3:ncol(BIO_TOTAL.MAG_m)]))

ggplot(data = BIO_TOTAL.MAG_m, aes(x = genome.size_mean, y = BIO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

sheet_write(BIO_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1OtixkEN_UAPIU7kNJqJjEkvOePwk1e0cNo7tidm9dMY/edit?gid=0#gid=0",
            sheet = "Biofilm_MAG")

Figure_test_7 = ggplot(data = BIO_TOTAL.MAG_m, aes(x = genome.size_mean, y = BIO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Biofilm") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,10)

# High Temp - genes (MAG)----
TEMP.MAG = total_genes.guild.186_MAG %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP.MAG$genome.size = as.numeric(as.character(TEMP.MAG$genome.size))
TEMP.MAG_m   = TEMP.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                 list(mean=mean), na.rm=TRUE))
TEMP.MAG_m   = TEMP.MAG_m %>% mutate(TEMP_total = rowSums(TEMP.MAG_m[,3:ncol(TEMP.MAG_m)]))

ggplot(data = TEMP.MAG_m, aes(x = genome.size_mean, y = TEMP_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

sheet_write(TEMP.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1spmhor9psOs_OXbMWbque6hldXX-zSIaLIKJE77b84Q/edit?gid=0#gid=0",
            sheet = "TEMP_MAG")

Figure_test_8 = ggplot(data = TEMP.MAG_m, aes(x = genome.size_mean, y = TEMP_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Heat Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,5)

# pH - genes (MAG)----
PH.MAG = total_genes.guild.186_MAG %>% select(any_of(pH_rule$`microtrait_hmm-name`))
PH.MAG$genome.size = as.numeric(as.character(PH.MAG$genome.size))
PH.MAG_m   = PH.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
PH.MAG_m   = PH.MAG_m %>% mutate(PH_total = rowSums(PH.MAG_m[,3:ncol(PH.MAG_m)]))

ggplot(data = PH.MAG_m, aes(x = genome.size_mean, y = PH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mean genome size per FG") + 
  ylab("Mean total osmolytes costs per FG") +
  geom_point() + theme(text = element_text(size=20))

sheet_write(PH.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/16mtVRN2sLusJC2_YU-To51oSmB_2EX_XxUwSkBn7YRI/edit?gid=0#gid=0",
            sheet = "PH_MAG")

Figure_test_9 = ggplot(data = PH.MAG_m, aes(x = genome.size_mean, y = PH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("pH Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,10)

# GH - genes (Isolates)----
GH_TOTAL.ISO = total_genes.guild.186_ISO %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.ISO$genome.size = as.numeric(as.character(GH_TOTAL.ISO$genome.size))
GH_TOTAL_m   = GH_TOTAL.ISO %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,4:ncol(GH_TOTAL_m)]))

sheet_write(GH_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1vxaR5T4dreJO5q_vVCCesBnBQvpcuD4kVYXyTqpcbXg/edit?gid=0#gid=0",
            sheet = "GH_Isolates")

Figure_test_1.ISO = ggplot(data = GH_TOTAL_m, aes(x = genome.size_mean, y = GH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CAZy") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,20)

# Protein - genes (Isolates)----
PH_TOTAL.ISO = total_genes.guild.186_ISO %>% select(any_of(PR_rule$`microtrait_rule-name`))
PH_TOTAL.ISO$genome.size = as.numeric(as.character(PH_TOTAL.ISO$genome.size))
PH_TOTAL_m   = PH_TOTAL.ISO %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PH_TOTAL_m   = PH_TOTAL_m %>% mutate(PR_total = rowSums(PH_TOTAL_m[,4:ncol(PH_TOTAL_m)]))

sheet_write(PH_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1wugrm2Jj77S3Q8jao0o_ILntnMGDsJdlfdZpvRYMr28/edit?gid=0#gid=0",
            sheet = "Protein_Isolates")

Figure_test_2.ISO = ggplot(data = PH_TOTAL_m, aes(x = genome.size_mean, y = PR_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Protein Enzyme") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,20)

# Transport - genes (Isolates)----
transp_rule    = transp_rule %>% filter(function_gene == c("transporter"))
transp_rule    = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL.i = total_genes.guild.186_ISO %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL.i$genome.size = as.numeric(as.character(TRANSP_TOTAL.i$genome.size))
TRANSP_TOTAL_m = TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,4:ncol(TRANSP_TOTAL_m)]))

sheet_write(TRANSP_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1ov3RjzE9NoiFGsT4B01atfJWvwYuZAqtXR3fOF8zl0Y/edit?gid=0#gid=0",
            sheet = "Transporter_Isolates")

Figure_test_3.ISO = ggplot(data = TRANSP_TOTAL_m, aes(x = genome.size_mean, y = transp_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Total Transporter") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,30)

# Transport - Aminoacids (Isolates)----
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL.i = total_genes.guild.186_ISO %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL.i$genome.size))
AMI_TRANSP_TOTAL.i_m = AMI_TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL.i_m = AMI_TRANSP_TOTAL.i_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL.i_m[,4:ncol(AMI_TRANSP_TOTAL.i_m)]))

sheet_write(AMI_TRANSP_TOTAL.i_m,
            ss = "https://docs.google.com/spreadsheets/d/14n8ydCBAwO_CtGuCjP0hMBcfZ_3qI8Wth7A2azK9xy0/edit?gid=0#gid=0",
            sheet = "Aminoacids_T_Isolates")

Figure_test_4.ISO = ggplot(data = AMI_TRANSP_TOTAL.i_m, aes(x = genome.size_mean, y = AMI_TRANSP_TOTAL)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Aminoacid Transporter") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,5)

# Transport - Carbohydrate (Isolates)----
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL.i = total_genes.guild.186_ISO %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL.i$genome.size))
CAR_TRANSP_TOTAL.i_m = CAR_TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL.i_m = CAR_TRANSP_TOTAL.i_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL.i_m[,4:ncol(CAR_TRANSP_TOTAL.i_m)]))

sheet_write(CAR_TRANSP_TOTAL.i_m,
            ss = "https://docs.google.com/spreadsheets/d/1ZKHklNdpi9SZjOizTMHsL2ScZqNIL1cFEnhS74e43C0/edit?gid=0#gid=0",
            sheet = "GH_T_Isolates")

Figure_test_5.ISO = ggplot(data = CAR_TRANSP_TOTAL.i_m, aes(x = genome.size_mean, y = CAR_TRANSP_TOTAL)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Carbohydrate Transporter") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,14)

# Osmolytes - genes (Isolates)----
OSMO_TOTAL.TOTAL.i = total_genes.guild.186_ISO %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.TOTAL.i$genome.size = as.numeric(as.character(OSMO_TOTAL.TOTAL.i$genome.size))
OSMO_TOTAL.TOTAL.i_m   = OSMO_TOTAL.TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                     list(mean=mean), na.rm=TRUE))
OSMO_TOTAL_m   = OSMO_TOTAL.TOTAL.i_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.TOTAL.i_m[,4:ncol(OSMO_TOTAL.TOTAL.i_m)]))

sheet_write(OSMO_TOTAL_m,
            ss = "https://docs.google.com/spreadsheets/d/1PpQdrLCFbzkanDn_R-PTmZAxTGQFbo3nk_77Tcxy0Yo/edit?gid=0#gid=0",
            sheet = "Osmolyte_Isolates")

Figure_test_6.ISO = ggplot(data = OSMO_TOTAL_m, aes(x = genome.size_mean, y = OSMO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Osmolytes") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,15)

# Biofilm - genes (Isolates)----
BIO_TOTAL.MAG = total_genes.guild.186_ISO %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,4:ncol(BIO_TOTAL.MAG_m)]))

sheet_write(BIO_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1fZD2srLcGZX6axd05_hlRZyL23izL4U7PIxeoL0HuuE/edit?gid=0#gid=0",
            sheet = "Biofilm_Isolates")

Figure_test_7.ISO = ggplot(data = BIO_TOTAL.MAG_m, aes(x = genome.size_mean, y = BIO_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Biofilm") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,10)

# High Temperature - genes (Isolates)----
TEMP_TOTAL.MAG = total_genes.guild.186_ISO %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP_TOTAL.MAG$genome.size = as.numeric(as.character(TEMP_TOTAL.MAG$genome.size))
TEMP_TOTAL.MAG_m   = TEMP_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
TEMP_TOTAL.MAG_m   = TEMP_TOTAL.MAG_m %>% mutate(TEMP_total = rowSums(TEMP_TOTAL.MAG_m[,4:ncol(TEMP_TOTAL.MAG_m)]))

sheet_write(TEMP_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1aBiAnokl323Lw8dj6k_B_uV4BIWNJvwVhoLuQeVRBU4/edit?gid=0#gid=0",
            sheet = "TEMP_Isolates")

Figure_test_8.ISO = ggplot(data = TEMP_TOTAL.MAG_m, aes(x = genome.size_mean, y = TEMP_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Heat Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,5)

# pH - genes (Isolates)----
pH_TOTAL.MAG = total_genes.guild.186_ISO %>% select(any_of(pH_rule$`microtrait_hmm-name`))
pH_TOTAL.MAG$genome.size = as.numeric(as.character(pH_TOTAL.MAG$genome.size))
pH_TOTAL.MAG_m   = pH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
pH_TOTAL.MAG_m   = pH_TOTAL.MAG_m %>% mutate(PH_total = rowSums(pH_TOTAL.MAG_m[,4:ncol(pH_TOTAL.MAG_m)]))

sheet_write(pH_TOTAL.MAG_m,
            ss = "https://docs.google.com/spreadsheets/d/1PPcpcweUZkPBwH8LpxBS7ZEuEkjOy3MTA4KI8QBfFkg/edit?gid=0#gid=0",
            sheet = "PH_Isolates")

Figure_test_9.ISO = ggplot(data = pH_TOTAL.MAG_m, aes(x = genome.size_mean, y = PH_total)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("pH Resistance") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,10)

# Approximation for Y strategy (for MAG)

total.granularity.3 = read.csv("Input_Data/IMG_JGI_MAGs/total.granularity.3_datasets.csv",dec=".") 
MAG_id              = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")
MAG_id              = as.data.frame(MAG_id$id)
colnames(MAG_id)    = "id"
Y_MAG               = merge(total.granularity.3,MAG_id,by = "id")
Y_MAG               = Y_MAG %>% select(c("guild","genome.size","mgt","ogt"))
Y_MAG.A             = Y_MAG %>% select(c("guild","genome.size","mgt"))

# Minimum generation time > 5
Y_MAG.A     = Y_MAG.A %>% filter(mgt > 5)
Y_MAG.A.1   = Y_MAG.A %>% group_by(guild) %>% summarise(across(where(is.numeric),list(mean=mean), na.rm=TRUE))

Figure_test_10a = ggplot(data = Y_MAG.A.1, aes(x = genome.size_mean/1e6, y = mgt_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Minimum generation time (hrs)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.5e7/1e6) + ylim(0,120)
Figure_test_10a

# Minimum generation time <= 5
Y_MAG.B   = Y_MAG %>% select(c("guild","genome.size","mgt"))
Y_MAG.B   = Y_MAG.B %>% filter(mgt <= 5)
Y_MAG.B.1 = Y_MAG.B %>% group_by(guild) %>% summarise(across(where(is.numeric),
                                                             list(mean=mean), na.rm=TRUE))

Figure_test_10b = ggplot(data = Y_MAG.B.1, aes(x = genome.size_mean/1e6, y = mgt_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Maximum growth rate (hrs)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,8)
Figure_test_10b

# Optimal Growth Temperature

Y_MAG.C     = Y_MAG %>% select(c("guild","genome.size","ogt"))
Y_MAG.C     = na.omit(Y_MAG.C)
Y_MAG.C.1   = Y_MAG.C %>% group_by(guild) %>% summarise(across(where(is.numeric),
                                                             list(mean=mean), na.rm=TRUE)) 

Figure_test_11 = ggplot(data = Y_MAG.C.1, aes(x = genome.size_mean/1e6, y = ogt_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("Optimal Growth Temperature (°C)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,80)
Figure_test_11

# Cross correlations----

# Data Preparation - MAGs
temp    = drive_ls("https://drive.google.com/drive/u/0/folders/1OlIc_JaS5a9s5jz7m1hTduTtYV7wAcG3")
a       = as.data.frame(temp$name)
Aminoacids_T_MAG = read_sheet(temp$id[1])
Aminoacids_T_MAG = Aminoacids_T_MAG %>% select(guild,genome.size_mean,AMI_TRANSP_TOTAL)
PH_MAG           = read_sheet(temp$id[2])
PH_MAG           = PH_MAG %>% select(guild,genome.size_mean,PH_total)
TEMP_MAG         = read_sheet(temp$id[3])
TEMP_MAG         = TEMP_MAG %>% select(guild,genome.size_mean,TEMP_total)
Biofilm_MAG      = read_sheet(temp$id[4])
Biofilm_MAG      = Biofilm_MAG %>% select(guild,genome.size_mean,BIO_total)
Osmolyte_MAG     = read_sheet(temp$id[5])
Osmolyte_MAG     = Osmolyte_MAG %>% select(guild,genome.size_mean,OSMO_total)
GH_T_MAG         = read_sheet(temp$id[6])
GH_T_MAG         = GH_T_MAG %>% select(guild,genome.size_mean,CAR_TRANSP_TOTAL)
Transporter_MAG  = read_sheet(temp$id[7])
Transporter_MAG  = Transporter_MAG %>% select(guild,genome.size_mean,transp_total)
Protein_MAG      = read_sheet(temp$id[8])
Protein_MAG      = Protein_MAG %>% select(guild,genome.size_mean,PR_total)
GH_MAG           = read_sheet(temp$id[9])
GH_MAG           = GH_MAG %>% select(guild,genome.size_mean,GH_total)
MAG_gen_trait    = as.data.frame(cbind(Aminoacids_T_MAG,PH_MAG$PH_total,
                                       TEMP_MAG$TEMP_total,Biofilm_MAG$BIO_total,
                                       Osmolyte_MAG$OSMO_total,GH_T_MAG$CAR_TRANSP_TOTAL,
                                       Transporter_MAG$transp_total,Protein_MAG$PR_total,
                                       GH_MAG$GH_total))
colnames(MAG_gen_trait) = c("guild","genome.size","amino.transport","pH","temp",
                            "biofilm","osmolyte","GH.trasnport","tranport.total",
                            "Protein","CAZy")
rm(temp,a,Aminoacids_T_MAG,PH_MAG,TEMP_MAG,Biofilm_MAG,Osmolyte_MAG,GH_T_MAG,
   Transporter_MAG,Protein_MAG,GH_MAG)

# Call data

MAG_gen_trait       = read_sheet("https://docs.google.com/spreadsheets/d/1Y6yMpQygXJof8EdDAE967ufLbNxmJcBqSZzUgMBnpYY/edit?gid=0#gid=0")
colnames(Y_MAG.B.1) = c("guild","genome.size","mgr")
colnames(Y_MAG.A.1) = c("guild","genome.size","mgr")
colnames(Y_MAG.C.1) = c("guild","genome.size","ogt")
MAG_gen_trait.1     = merge(MAG_gen_trait,Y_MAG.C.1,by = c("guild","genome.size"))
MAG_gen_trait.2     = merge(MAG_gen_trait.1,Y_MAG.B.1,by = c("guild","genome.size"))
MAG_gen_trait.3     = merge(MAG_gen_trait.1,Y_MAG.A.1,by = c("guild","genome.size"))

sheet_write(MAG_gen_trait.2,
            ss = "https://docs.google.com/spreadsheets/d/1Gtb4oLsNibF-yPwA9f1y78axGuKCO70Y1cZKNw0KRtk/edit?gid=0#gid=0",
            sheet = "MAG_gen_trait")
sheet_write(MAG_gen_trait.3,
            ss = "https://docs.google.com/spreadsheets/d/1Y6yMpQygXJof8EdDAE967ufLbNxmJcBqSZzUgMBnpYY/edit?gid=0#gid=0",
            sheet = "MAG_gen_trait.5")

# Plotting cross-correlations

MAG_gen_trait.2  = read_sheet("https://docs.google.com/spreadsheets/d/1Gtb4oLsNibF-yPwA9f1y78axGuKCO70Y1cZKNw0KRtk/edit?gid=0#gid=0")
MAG_gen_trait.3  = read_sheet("https://docs.google.com/spreadsheets/d/1Y6yMpQygXJof8EdDAE967ufLbNxmJcBqSZzUgMBnpYY/edit?gid=0#gid=0")

colnames(MAG_gen_trait.2) = c("Guild", "Genome Size", "Amino-transporter", "pH-Tol",
                               "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                               "Total-transporter","Protein-enzyme","CAZy","OGT",
                               "MGT")
colnames(MAG_gen_trait.3) = c("Guild", "Genome Size", "Amino-transporter", "pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","OGT",
                              "MGT")

# Reordering columns

col_order = c("Guild", "Genome Size","OGT","MGT","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte")

MAG_gen_trait.2 = MAG_gen_trait.2[, col_order]
MAG_gen_trait.3 = MAG_gen_trait.3[, col_order]

# Function for merging colors and correlation plot
my_fn <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("red", "white", "blue"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

# Actual figure
Figure_5 = ggpairs(MAG_gen_trait.2[,3:13], 
              upper = list(continuous = my_fn),
              lower = list(continuous = "smooth"))  
Figure_5

Figure_5.a = ggpairs(MAG_gen_trait.3[,3:13], 
                     upper = list(continuous = my_fn),
                     lower = list(continuous = "smooth")) 
Figure_5.a

# CUE ----

Isolates_gen_trait  = read_sheet("https://docs.google.com/spreadsheets/d/1xT3vB2K2tupkeHZwUlVzh4ZYWhmUpOmMQtESNH2Mst0/edit?gid=0#gid=0")
isolates     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/dement_isolates_CUE.csv",dec=".")
isolatest    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/total.granularity.3_datasets.csv",dec=".")
isolates.1   = subset(isolates, isolates$CUE !='NaN')
total_genes.guild.940_ISO = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.940_ISO.csv",dec=".")

# Merge

isolates.2   = isolates.1 %>% left_join(total_genes.guild.940_ISO, by='id')
isolates.2.m = isolates.2 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
isolates.2.m = isolates.2.m %>% select(c("guild","CUE_mean","genome.size_mean",
                                         "mingentime_mean","optimumT_mean"))
colnames(isolates.2.m) = c("guild","yield","genome.size","MGT","OGT")
ISO_gen_trait.1        = merge(Isolates_gen_trait,isolates.2.m,by = c("guild","genome.size"))

sheet_write(ISO_gen_trait.1,
            ss = "https://docs.google.com/spreadsheets/d/1xT3vB2K2tupkeHZwUlVzh4ZYWhmUpOmMQtESNH2Mst0/edit?gid=0#gid=0v",
            sheet = "ISO_gen_trait.1")

# Plotting cross correlations

ISO_gen_trait.1  = read_sheet("https://docs.google.com/spreadsheets/d/1xT3vB2K2tupkeHZwUlVzh4ZYWhmUpOmMQtESNH2Mst0/edit?gid=0#gid=0",
                              sheet = "ISO_gen_trait.1")

colnames(ISO_gen_trait.1) = c("Guild","Genome Size","Amino-transporter","pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","Yield","MGT","OGT")
# Reordering columns

col_order = c("Guild", "Genome Size","Yield","OGT","MGT","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte")
ISO_gen_trait.1   = ISO_gen_trait.1[, col_order]
ISO_gen_trait.1   = ISO_gen_trait.1 %>% filter(Yield <= 0.9)

# Minimum generation time <= 5
ISO_gen_trait.1.f = ISO_gen_trait.1 %>% filter(MGT <= 5)

Figure_6 = ggpairs(ISO_gen_trait.1.f[,3:14], 
                   upper = list(continuous = my_fn),
                   lower = list(continuous = "smooth"))  
Figure_6

# Minimum generation time > 5
ISO_gen_trait.1.s = ISO_gen_trait.1 %>% filter(MGT > 5)

Figure_6.a = ggpairs(ISO_gen_trait.1.s[,3:14], 
                   upper = list(continuous = my_fn),
                   lower = list(continuous = "smooth"))  
Figure_6.a

# CUE for different taxonomic levels----

# Metadata

meta_iso   = readr::read_tsv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils_unpublished.tsv")
meta_iso   = as.data.frame(cbind(meta_iso$`IMG Genome ID`,meta_iso$Phylum,
                                 meta_iso$Class,meta_iso$Order,meta_iso$Family,
                                 meta_iso$Genus))
meta_iso.1 = readr::read_tsv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils.tsv")
meta_iso.1 = as.data.frame(cbind(meta_iso.1$`IMG Genome ID`,meta_iso.1$Phylum,
                                 meta_iso.1$Class,meta_iso.1$Order,meta_iso.1$Family,
                                 meta_iso.1$Genus))

# final meta
final_meta   = as.data.frame(rbind(meta_iso,meta_iso.1))
colnames(final_meta) = c("id","Phylum","Class","Order","Family","Genus")

merged_df    = merge(isolates.1,final_meta,by = "id")
merged_df    = merged_df %>% filter(CUE <= 0.9)

# Phylum

phylum = merged_df %>% group_by(Phylum) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
phylum = phylum %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                             mingentime_mean  > 5 ~ "slow"))

Figure_test_1.phylum = ggplot(data = phylum, aes(x = genome_length_mean/1e6, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Phylum")
Figure_test_1.phylum

Figure_test_1.phylum.new = ggplot(data = phylum, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                     color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Phylum") + theme(legend.position="none")
Figure_test_1.phylum.new

# Class

Class = merged_df %>% group_by(Class) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Class = Class %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                             mingentime_mean  > 5 ~ "slow"))
Figure_test_1.Class = ggplot(data = Class, aes(x = genome_length_mean/1e6, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Class")
Figure_test_1.Class

Figure_test_1.Class.new = ggplot(data = Class, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                     color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Class") + theme(legend.position="none")
Figure_test_1.Class.new

# Order

Order = merged_df %>% group_by(Order) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Order = Order %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Order = ggplot(data = Order, aes(x = genome_length_mean/1e6, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Order")
Figure_test_1.Order

Figure_test_1.Order.new = ggplot(data = Order, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                   color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Order") + theme(legend.position="none")
Figure_test_1.Order.new

# Family

Family = merged_df %>% group_by(Family) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
Family = Family %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Family = ggplot(data = Family, aes(x = genome_length_mean/1e6, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Family")
Figure_test_1.Family

Figure_test_1.Family.new = ggplot(data = Family, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                   color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Family") + theme(legend.position="none")
Figure_test_1.Family.new

# Genus

Genus = merged_df %>% group_by(Genus) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Genus = Genus %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                             mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Genus = ggplot(data = Genus, aes(x = genome_length_mean/1e6, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Genus")
Figure_test_1.Genus

Figure_test_1.Genus.new = ggplot(data = Genus, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                     color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Genus") + theme(legend.position="none")
Figure_test_1.Genus.new


