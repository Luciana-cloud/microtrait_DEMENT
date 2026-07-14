# IMG MAG ANALYSIS
# Analysis of 33000 MAGs from IMG + LOMAa and FIRE MAGs + 5000 Soil Isolates

# Calling Packages ----

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(readr)
library(vegan)
library(devtools)
#library(ggvegan)
library(pairwiseAdonis)
#library(EcolUtils)
library(olsrr)
#library(stats)
library(corrr)
library(car)
#install.packages("collinear")
library(collinear)
library(glmnet)
library(parallelDist)
#library(googledrive)
#library(googlesheets4)
#library(readxl)
library(ggpmisc)
library(GGally)
library(reshape2)
library(factoextra)
library(domir)

# Set directory----

setwd("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset")

# Calling data IMG-JGI----

hmm_img       = read.csv("Input_Data/IMG_JGI_MAGs/hmm_MAG_IMG.csv",dec=".")
coverage_img  = read.csv("Input_Data/IMG_JGI_MAGs/bin_coverage-IMG_MAGs.csv",dec=".")
meta_img.nr   = read.csv("Input_Data/IMG_JGI_MAGs/IMG_bindata_withmeta_norestricted.csv",dec=".")

# Change column name to match
colnames(hmm_img)[2]      = "bin_id"
colnames(meta_img.nr)[2]  = "bin_id"
colnames(coverage_img)[1] = "bin_id"

# Merge databases
data_combined = meta_img.nr %>% left_join(hmm_img, by='bin_id') 
data_combined = data_combined %>% left_join(coverage_img, by='bin_id')
data_combined = data_combined %>% filter(Domain == "Bacteria")
#write.csv(data_combined, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/IMG_global_datasets.csv")

# Data preparation ----

# Summary MAGs by project to select what project to use for extracting filtering genes
data_combined   = data_combined[!(is.na(data_combined$avg_coverage)),]
data_combined   = data_combined[!(is.na(data_combined$aaeB)),]
data_combined   = as.data.frame(apply(data_combined,2,as.character))
# write.csv(data_combined, file = "Intermediate_Results/Intermediate_results/IMG_global_dataset.csv")
data_combined.1   = as.data.frame(data_combined) %>% select(1,2,44,69:1789)
data_combined_IMG.ID  = data_combined %>% group_by(IMG.Genome.ID) %>% count()
data_combined_IMG.ID  = data_combined_IMG.ID %>% filter(n > 20)

# Get hmm files from data_combined i.e. after data filtering
hmm_img        = hmm_img[,-1]
hmm_img.filter = data_combined %>% select(names(hmm_img))

# Select and save individual projects with more than 20 MAGs for further analysis
a         = unique(data_combined_IMG.ID$IMG.Genome.ID)
data_project_100 = c()
for(i in a){
  project = data_combined %>% filter(IMG.Genome.ID == i)
  temp    = project %>% summarise(sum(as.numeric(avg_coverage), na.rm = TRUE))
  colnames(temp) = "Total"
  project = project %>% mutate(RelAbund = as.numeric(avg_coverage)/temp$Total)
  data_project_100 = rbind(data_project_100, project) 
}
# write.csv(data_project_100, file = "Intermediate_Results/IMG_global_datasets_relative_abundance.csv")

# Reduce NAs
data_project_20 = data_project_100[!(is.na(data_project_100$RelAbund)),]
data_project_20[,1790] = as.numeric(data_project_20[,1790])
data_project_20[,69:1789] = as.numeric(unlist(data_project_20[,69:1789]))

# Reduce columns with zero
temp.1            = as.data.frame(cbind((colSums(data_project_20[,69:1789]))))
erase.1           = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
data_project_20   = data_project_20[ , !(names(data_project_20) %in% erase.1$row_names)]

# Selection of genes per project (best predictors of MAG abundances) for IMG projects ----

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
  setwd("Intermediate_Results/IMG_MAG_predictors")
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
  setwd("Intermediate_Results/IMG_MAG_predictors")
  write.csv(list[[i]], paste0(i,".csv"))
}

# Selection of genes per project (best predictors of MAG abundances) for LOMA project ----

# Calling data and preprocessing
hmm_loma    = read.csv("Input_Data/LOMA_MAGs/hmm_Loma.csv",dec=".")
mag_stat    = read.csv("Input_Data/LOMA_MAGs/mag_stats.csv",dec=".") 
mag_abun    = read.csv("Input_Data/LOMA_MAGs/mag_adundance.csv",dec=".") 
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
for(i in a){
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
  setwd("Intermediate_Results/LOMA_predictors")
  write.csv(list[[i]], paste0(i,".csv"))
}

# Selection of genes per project (best predictors of MAG abundances) for FIRE project ----

# Calling data and preprocessing
hmm_fire    = read.csv("Input_Data/FIRE_MAGs/hmm_Fire.csv",dec=".")
mag_abun    = read.csv("Input_Data/FIRE_MAGs/mag_adundance_fire.csv",dec=".") 
mag_abun    = mag_abun[-c(440,546), ]
fire_meta   = read_excel("Input_Data/FIRE_MAGs/MAG_Dataset_BurnSeverity_ARNelson.xlsx")

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
for(i in a){
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
  setwd("Intermediate_Results/FIRE_predictors")
  write.csv(list[[i]], paste0(i,".csv"))
}

# Functional group analysis ----

# Call Genes and merge datasets - IMG
IMG_predictors_files = list.files(path="Intermediate_Results/IMG_MAG_predictors", pattern=".csv", all.files=FALSE, 
                                  full.names=TRUE)
IMG_predictors_list = lapply(IMG_predictors_files, read.csv)
IMG_predictors_gene = do.call(rbind.data.frame, IMG_predictors_list)
IMG_predictors_gene = unique(IMG_predictors_gene$x)

# Call Genes and merge datasets - Loma
Loma_predictors_files = list.files(path="Intermediate_Results/LOMA_predictors", pattern=NULL, all.files=FALSE, 
                                   full.names=TRUE)
Loma_predictors_list = lapply(Loma_predictors_files, read.csv)
Loma_predictors_gene = do.call(rbind.data.frame, Loma_predictors_list)
Loma_predictors_gene = unique(Loma_predictors_gene$x)

# Call Genes and merge datasets - Fire
Fire_predictors_files = list.files(path="Intermediate_Results/FIRE_predictors", pattern=NULL, all.files=FALSE, 
                                   full.names=TRUE)
Fire_predictors_list = lapply(Fire_predictors_files, read.csv)
Fire_predictors_gene = do.call(rbind.data.frame, Fire_predictors_list)
Fire_predictors_gene = unique(Fire_predictors_gene$x)

# Call data (MAGs)
total_genes            = as.data.frame(unique(c(IMG_predictors_gene,Loma_predictors_gene,Fire_predictors_gene)))
colnames(total_genes)  = c("gene")
write.csv(total_genes, file = "Output_Data/Best_predictors/Selected_genes.csv")
colnames(hmm_img.filter)[2]   = "id"
total_gene_names       = c("id",total_genes$gene)
hmm_fire.global        = hmm_fire %>% select(tidyselect::any_of(total_gene_names))
hmm_loma.global        = hmm_loma %>% select(tidyselect::any_of(total_gene_names))
hmm_img.global         = hmm_img.filter %>% select(tidyselect::any_of(total_gene_names))

# Call data (Isolates)
isolates        = read.csv("Input_Data/SOIL_ISOLATES/hmm_isolates.csv",dec=".")
df_1            = read.delim("Input_Data/SOIL_ISOLATES/metadata_2.tsv",sep="\t")
df_1            = subset(df_1, select = -c(23))
df_2            = read.delim("Input_Data/SOIL_ISOLATES/metadata_1.tsv",sep="\t")
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
write.csv(MAG_global, file = "Output_Data/Global_dataset_microtrait.csv")

# Call the MAG-Microtrait datasets
global_dataset              = read.csv("Output_Data/Global_dataset_microtrait.csv",dec=".")

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

# Fig_S1 ----

Figure_S1 = ggplot(data,aes(V1,V2)) +
  geom_point() +
  geom_smooth(method='lm') + 
  stat_cor(label.x = 30, label.y = 130, size = 4) +
  stat_regline_equation(label.x = 30, label.y = 150, size = 4) + 
  xlab("# of MAGs") + ylab("# of Functional Groups") + theme_classic() + 
  theme(text = element_text(size=20))
Figure_S1

pdf("Output_Data/Figures/Figure_S1.pdf",
    width=5,height=5)
print(Figure_S1)
dev.off()

# Final functional groups for the total dataset----
FG = 16.597382 + 0.041440*32515 # Number of Functional groups for whole dataset

# Calling distance matrix
#load("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/FG_parallel/distance.total.RData")
set.seed(1)
distance.total  = parDist(x = as.matrix(global_dataset[,3:938]),
                          method = "fJaccard",
                          threads = 1) # Adapt the number of threads
cluster.total   = hclust(distance.total, method="ward.D2")
v                      = cutree(cluster.total,k=round(FG))
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
global_dataset.1       = as.data.frame(cbind(genome2guild,global_dataset))
write.csv(global_dataset.1, file = "Intermediate_Results/global_datasets_guild.csv")

# Genome size analysis----

# Call data again
# IMG-MAG
hmm_img     = read.csv("Intermediate_Results/IMG_global_dataset.csv",dec=".")
hmm_img     = as.data.frame(hmm_img) %>% select(3,26,70:1790)
colnames(hmm_img)[1]  = "id"
colnames(hmm_img)[2]  = "genome.size"

# Loma-MAG
hmm_loma    = read.csv("Input_Data/LOMA_MAGs/hmm_Loma.csv",dec=".")
gene_loma   = read.csv("Input_Data/LOMA_MAGs/litter_mags_metadata.csv",dec=".")
hmm_loma    = gene_loma %>% full_join(hmm_loma)
hmm_loma[3] = NULL
colnames(hmm_loma)[2]  = "genome.size"

# Wildfire-MAG
hmm_fire    = read.csv("Input_Data/FIRE_MAGs/hmm_Fire.csv",dec=".")
fire_stat   = read.csv("Input_Data/FIRE_MAGs/fire_metadata.csv",dec=".")
fire_stat   = fire_stat[-c(440,546), ]
hmm_fire    = fire_stat %>% full_join(hmm_fire)
hmm_fire[3] = NULL
colnames(hmm_fire)[2]  = "genome.size"

# IMG-Isolates
isolates      = read.csv("Input_Data/SOIL_ISOLATES/hmm_isolates.csv",dec=".")
df_1          = read.delim("Input_Data/SOIL_ISOLATES/metadata_2.tsv",sep="\t")
df_1          = subset(df_1, select = -c(23))
df_2          = read.delim("Input_Data/SOIL_ISOLATES/metadata_1.tsv",sep="\t")
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
total_genes.guild     = read.csv("Intermediate_Results/global_datasets_guild.csv",dec=".")
total.guild           = as.data.frame(total_genes.guild) %>% select(4,2)
total.g.size          = as.data.frame(total_genes) %>% select(1,2)

# Total guilds for the 940 selected genes
total_genes.guild     = left_join(total.g.size,total_genes.guild, by=c('id'))
total_genes.guild.940 = as.data.frame(total_genes.guild) %>% select(1,2,4,6:941)
write.csv(total_genes.guild.940, file = "Intermediate_Results/total_genes.guild.940.csv")

# Total guilds for the complete MAG-genes matrix
total_genes           = left_join(total.guild, total_genes, by=c('id'))
write.csv(total_genes, file = "Intermediate_Results/total_genes.full.csv")

# 940 MAG-Gene matrix----
# Call Trait Keys 
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
pH_rule      = read_excel("Input_Data/microTrait_Rules/microtrait_pH_stress.xlsx")
pH_rule      = as.data.frame(rbind("id","guild","genome.size",pH_rule)) 

# Genome-size for MAGs----
key_MAGs = bind_rows(hmm_img,hmm_loma,hmm_fire)
key_MAGs = as.data.frame(key_MAGs$id)
colnames(key_MAGs) = "key"
total_genes.guild.940_MAG = total_genes.guild.940[total_genes.guild.940$id %in% key_MAGs$key, ]
write.csv(total_genes.guild.940_MAG, file = "Intermediate_Results/total_genes.guild.940_MAG.csv")

# Genome-size for Isolates----
key_isolates = as.data.frame(isolates$id)
colnames(key_isolates) = "key"
total_genes.guild.940_ISO = total_genes.guild.940[total_genes.guild.940$id %in% key_isolates$key, ]
write.csv(total_genes.guild.940_ISO, file = "Intermediate_Results/total_genes.guild.940_ISO.csv")

# Trait Data Preparation (MAGs) ----

total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")

# GH - genes (MAG)
GH_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.MAG$genome.size = as.numeric(as.character(GH_TOTAL.MAG$genome.size))
GH_TOTAL_m   = GH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,3:ncol(GH_TOTAL_m)]))
GH_TOTAL_m.1 = GH_TOTAL_m[-c(247), ] # Erasing super big Functional Group
GH_MAG       = GH_TOTAL_m.1 %>% select(guild,genome.size_mean,GH_total)

# Protein - genes (MAG)
PR_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(PR_rule$`microtrait_hmm-name`))
PR_TOTAL.MAG$genome.size = as.numeric(as.character(PR_TOTAL.MAG$genome.size))
PR_TOTAL_m   = PR_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PR_TOTAL_m   = PR_TOTAL_m %>% mutate(PR_total = rowSums(PR_TOTAL_m[,3:ncol(PR_TOTAL_m)]))
PR_TOTAL_m.1 = PR_TOTAL_m[-c(247), ] # Erasing super big Functional Group
Protein_MAG  = PR_TOTAL_m.1 %>% select(guild,genome.size_mean,PR_total)

# Transport - genes (MAG)
transp_rule  = transp_rule %>% filter(`function` == c("transporter"))
transp_rule  = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL = total_genes.guild.940_MAG %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL$genome.size = as.numeric(as.character(TRANSP_TOTAL$genome.size))
TRANSP_TOTAL_m   = TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m   = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,3:ncol(TRANSP_TOTAL_m)]))
TRANSP_TOTAL_m.1 = TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group
Transporter_MAG  = TRANSP_TOTAL_m.1 %>% select(guild,genome.size_mean,transp_total)

# Transport - Aminoacids (MAG)
transp_rule        = transp_rule %>% filter(`function`==c("transporter"))
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL   = total_genes.guild.940_MAG %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL$genome.size))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL_m[,3:ncol(AMI_TRANSP_TOTAL_m)]))
AMI_TRANSP_TOTAL_m.1 = AMI_TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group
Aminoacids_T_MAG   = AMI_TRANSP_TOTAL_m.1 %>% select(guild,genome.size_mean,AMI_TRANSP_TOTAL)

# Transport - Carbohydrate (MAG)
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL   = total_genes.guild.940_MAG %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL$genome.size))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL_m[,3:ncol(CAR_TRANSP_TOTAL_m)]))
CAR_TRANSP_TOTAL_m.1 = CAR_TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group
GH_T_MAG           = CAR_TRANSP_TOTAL_m.1 %>% select(guild,genome.size_mean,CAR_TRANSP_TOTAL)

# Osmolytes - genes (MAG)
OSMO_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.MAG$genome.size = as.numeric(as.character(OSMO_TOTAL.MAG$genome.size))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.MAG_m[,3:ncol(OSMO_TOTAL.MAG_m)]))
OSMO_TOTAL.MAG_m.1 = OSMO_TOTAL.MAG_m[-c(247), ] # Erasing super big Functional Group
Osmolyte_MAG       = OSMO_TOTAL.MAG_m.1 %>% select(guild,genome.size_mean,OSMO_total)

# Biofilm - genes (MAG)
BIO_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,3:ncol(BIO_TOTAL.MAG_m)]))
BIO_TOTAL.MAG_m.1 = BIO_TOTAL.MAG_m[-c(247), ] # Erasing super big Functional Group
Biofilm_MAG       = BIO_TOTAL.MAG_m.1 %>% select(guild,genome.size_mean,BIO_total)

# High Temp - genes (MAG)
TEMP.MAG = total_genes.guild.940_MAG %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP.MAG$genome.size = as.numeric(as.character(TEMP.MAG$genome.size))
TEMP.MAG_m   = TEMP.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                 list(mean=mean), na.rm=TRUE))
TEMP.MAG_m   = TEMP.MAG_m %>% mutate(TEMP_total = rowSums(TEMP.MAG_m[,3:ncol(TEMP.MAG_m)]))
TEMP.MAG_m.1 = TEMP.MAG_m[-c(247), ] # Erasing super big Functional Group
TEMP_MAG     = TEMP.MAG_m.1 %>% select(guild,genome.size_mean,TEMP_total)

# pH - genes (MAG)
PH.MAG = total_genes.guild.940_MAG %>% select(any_of(pH_rule$`microtrait_hmm-name`))
PH.MAG$genome.size = as.numeric(as.character(PH.MAG$genome.size))
PH.MAG_m   = PH.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
PH.MAG_m   = PH.MAG_m %>% mutate(PH_total = rowSums(PH.MAG_m[,3:ncol(PH.MAG_m)]))
PH.MAG_m.1 = PH.MAG_m[-c(247), ] # Erasing super big Functional Group
PH_MAG     = PH.MAG_m.1 %>% select(guild,genome.size_mean,PH_total)

MAG_gen_trait    = as.data.frame(cbind(Aminoacids_T_MAG,PH_MAG$PH_total,
                                       TEMP_MAG$TEMP_total,Biofilm_MAG$BIO_total,
                                       Osmolyte_MAG$OSMO_total,GH_T_MAG$CAR_TRANSP_TOTAL,
                                       Transporter_MAG$transp_total,Protein_MAG$PR_total,
                                       GH_MAG$GH_total))
colnames(MAG_gen_trait) = c("guild","genome.size","amino.transport","pH","temp",
                            "biofilm","osmolyte","GH.trasnport","tranport.total",
                            "Protein","CAZy")

write.csv(MAG_gen_trait, file = "Intermediate_Results/MAG_gen_trait.csv")

# Trait Data Preparation (Isolates) ----

total_genes.guild.940_ISO = read.csv("Intermediate_Results/total_genes.guild.940_ISO.csv",dec=".")

# GH - genes (Isolates)
GH_TOTAL.ISO = total_genes.guild.940_ISO %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.ISO$genome.size = as.numeric(as.character(GH_TOTAL.ISO$genome.size))
GH_TOTAL_m   = GH_TOTAL.ISO %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
GH_TOTAL_m   = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,4:ncol(GH_TOTAL_m)]))
GH_Isolates  = GH_TOTAL_m %>% select(guild,genome.size_mean,GH_total)

# Protein - genes (Isolates)
PH_TOTAL.ISO = total_genes.guild.940_ISO %>% select(any_of(PR_rule$`microtrait_hmm-name`))
PH_TOTAL.ISO$genome.size = as.numeric(as.character(PH_TOTAL.ISO$genome.size))
PH_TOTAL_m   = PH_TOTAL.ISO %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PH_TOTAL_m   = PH_TOTAL_m %>% mutate(PR_total = rowSums(PH_TOTAL_m[,4:ncol(PH_TOTAL_m)]))
Protein_Isolates      = PH_TOTAL_m %>% select(guild,genome.size_mean,PR_total)

# Transport - genes (Isolates)
transp_rule    = transp_rule %>% filter(`function` == c("transporter"))
transp_rule    = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL.i$genome.size = as.numeric(as.character(TRANSP_TOTAL.i$genome.size))
TRANSP_TOTAL_m = TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m = TRANSP_TOTAL_m %>% mutate(transp_total = rowSums(TRANSP_TOTAL_m[,4:ncol(TRANSP_TOTAL_m)]))
Transporter_Isolates  = TRANSP_TOTAL_m %>% select(guild,genome.size_mean,transp_total)

# Transport - Aminoacids (Isolates)
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL.i$genome.size))
AMI_TRANSP_TOTAL.i_m = AMI_TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL.i_m = AMI_TRANSP_TOTAL.i_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL.i_m[,4:ncol(AMI_TRANSP_TOTAL.i_m)]))
Aminoacids_T_Isolates = AMI_TRANSP_TOTAL.i_m %>% select(guild,genome.size_mean,AMI_TRANSP_TOTAL)

# Transport - Carbohydrate (Isolates)
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL.i$genome.size))
CAR_TRANSP_TOTAL.i_m = CAR_TRANSP_TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL.i_m = CAR_TRANSP_TOTAL.i_m %>% mutate(CAR_TRANSP_TOTAL = rowSums(CAR_TRANSP_TOTAL.i_m[,4:ncol(CAR_TRANSP_TOTAL.i_m)]))
GH_T_Isolates     = CAR_TRANSP_TOTAL.i_m %>% select(guild,genome.size_mean,CAR_TRANSP_TOTAL)

# Osmolytes - genes (Isolates)
OSMO_TOTAL.TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.TOTAL.i$genome.size = as.numeric(as.character(OSMO_TOTAL.TOTAL.i$genome.size))
OSMO_TOTAL.TOTAL.i_m   = OSMO_TOTAL.TOTAL.i %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                     list(mean=mean), na.rm=TRUE))
OSMO_TOTAL_m   = OSMO_TOTAL.TOTAL.i_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.TOTAL.i_m[,4:ncol(OSMO_TOTAL.TOTAL.i_m)]))
Osmolyte_Isolates = OSMO_TOTAL_m %>% select(guild,genome.size_mean,OSMO_total)

# Biofilm - genes (Isolates)
BIO_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_m %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG_m[,4:ncol(BIO_TOTAL.MAG_m)]))
Biofilm_Isolates  = BIO_TOTAL.MAG_m %>% select(guild,genome.size_mean,BIO_total)

# High Temperature - genes (Isolates)
TEMP_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP_TOTAL.MAG$genome.size = as.numeric(as.character(TEMP_TOTAL.MAG$genome.size))
TEMP_TOTAL.MAG_m   = TEMP_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
TEMP_TOTAL.MAG_m   = TEMP_TOTAL.MAG_m %>% mutate(TEMP_total = rowSums(TEMP_TOTAL.MAG_m[,4:ncol(TEMP_TOTAL.MAG_m)]))
TEMP_Isolates      = TEMP_TOTAL.MAG_m %>% select(guild,genome.size_mean,TEMP_total)

# pH - genes (Isolates)
pH_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(pH_rule$`microtrait_hmm-name`))
pH_TOTAL.MAG$genome.size = as.numeric(as.character(pH_TOTAL.MAG$genome.size))
pH_TOTAL.MAG_m   = pH_TOTAL.MAG %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
pH_TOTAL.MAG_m   = pH_TOTAL.MAG_m %>% mutate(PH_total = rowSums(pH_TOTAL.MAG_m[,4:ncol(pH_TOTAL.MAG_m)]))
PH_Isolates      = pH_TOTAL.MAG_m %>% select(guild,genome.size_mean,PH_total)

Isolates_gen_trait    = as.data.frame(cbind(Aminoacids_T_Isolates,PH_Isolates$PH_total,
                                            TEMP_Isolates$TEMP_total,Biofilm_Isolates$BIO_total,
                                            Osmolyte_Isolates$OSMO_total,GH_T_Isolates$CAR_TRANSP_TOTAL,
                                            Transporter_Isolates$transp_total,Protein_Isolates$PR_total,
                                            GH_Isolates$GH_total))
colnames(Isolates_gen_trait) = c("guild","genome.size","amino.transport","pH","temp",
                                 "biofilm","osmolyte","GH.trasnport","tranport.total",
                                 "Protein","CAZy")

write.csv(Isolates_gen_trait, file = "Intermediate_Results/Isolates_gen_trait.csv")

# Additional preparations for MAGs ----

# Adding ogt and mgt ----

microtrait = read_csv("Input_Data/IMG_JGI_MAGs/total.granularity.3_datasets.csv")
total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")
guild_n = total_genes.guild.940_MAG %>% count(guild, name = "n_mags")

# One aggregation step for ogt + mgt + genome.size together ----

MAG_growth = microtrait %>% select(guild, genome.size, ogt, mgt) %>%
  group_by(guild) %>% summarise(across(everything(), mean, na.rm = TRUE))

# Calculating maximum growth rate ----
MAG_growth = MAG_growth %>% mutate(mgr = 1/mgt)   

# Updated MAG_gen_trait file ----

MAG_gen_trait_v2 = merge(MAG_gen_trait, MAG_growth %>% select(guild, ogt, mgt, mgr),
                         by = "guild")
MAG_gen_trait_v2 = merge(MAG_gen_trait_v2, guild_n, by = "guild")
MAG_gen_trait_v2 = MAG_gen_trait_v2 %>% select(-`...1`)

write.csv(MAG_gen_trait_v2, file = "Intermediate_Results/MAG_gen_trait.csv")

# Adding contamination, GC count and completeness ----

# Calling Data

# IMG MAGs
metadata = read.csv("Input_Data/IMG_JGI_MAGs/IMG_bindata_withmeta_norestricted.csv", dec=".")  # or wherever `metadata` came from originally

# Fire
fire_meta   = readxl::read_excel("Input_Data/FIRE_MAGs/MAG_Dataset_BurnSeverity_ARNelson.xlsx")
colnames(fire_meta)[1] = "id"
fire_stat   = read.csv("Input_Data/FIRE_MAGs/fire_metadata.csv",dec=".")

# Loma
mag_stat = read.csv("Input_Data/LOMA_MAGs/mag_stats.csv", dec=".")
mag_abun = read.csv("Input_Data/LOMA_MAGs/mag_adundance.csv", dec=".")
mag_stat = mag_stat %>% full_join(mag_abun)

# Aligning the IDs

# Split the guild table's ids by source
guild_ids = total_genes.guild.940_MAG %>% select(id, guild, genome.size)

guild_ids_fire = guild_ids %>% filter(grepl("^RYN", id))
guild_ids_loma = guild_ids %>% filter(grepl("shrub|grass", id, ignore.case = TRUE)) %>%
  mutate(id_stripped = sub("\\.orig$|\\.strict$", "", id))
guild_ids_img  = guild_ids %>% filter(!grepl("^RYN", id) & !grepl("shrub|grass", id, ignore.case = TRUE))

cat("Routing check - should sum to", nrow(guild_ids), "total rows:\n")
cat("  IMG: ", nrow(guild_ids_img),
    "| Fire:", nrow(guild_ids_fire),
    "| Loma:", nrow(guild_ids_loma),
    "| Sum:", nrow(guild_ids_img) + nrow(guild_ids_fire) + nrow(guild_ids_loma), "\n\n")

# IMG: join by Bin.ID
IMG_meta = metadata %>% select(Bin.ID, Bin.Completeness, Bin.Contamination, GC.....assembled) %>%
  rename(id = Bin.ID, completeness = Bin.Completeness,
         contamination = Bin.Contamination, GC_count = GC.....assembled)
IMG_joined = left_join(guild_ids_img, IMG_meta, by = "id")
cat("IMG - rows before:", nrow(guild_ids_img), "| after join:", nrow(IMG_joined),
    "| unmatched (NA completeness):", sum(is.na(IMG_joined$completeness)), "\n\n")

# Fire: join by id (RYN_n) directly
FIRE_meta = fire_meta %>%
  merge(fire_stat, by = "id") %>%
  select(id, Completeness, Contamination, GC) %>%
  rename(completeness = Completeness, contamination = Contamination, GC_count = GC)
FIRE_joined = left_join(guild_ids_fire, FIRE_meta, by = "id")
cat("Fire - rows before:", nrow(guild_ids_fire), "| after join:", nrow(FIRE_joined),
    "| unmatched (NA completeness):", sum(is.na(FIRE_joined$completeness)), "\n\n")

# Loma: join by stripped id
LOMA_meta = mag_stat %>% select(id, completeness, contamination) %>% mutate(GC_count = NA)
LOMA_joined = left_join(guild_ids_loma, LOMA_meta, by = c("id_stripped" = "id"))
cat("Loma - rows before:", nrow(guild_ids_loma), "| after join:", nrow(LOMA_joined),
    "| unmatched (NA completeness):", sum(is.na(LOMA_joined$completeness)), "\n\n")
LOMA_joined = LOMA_joined %>% select(-id_stripped)

# Recombine
TOTAL_MAGs_fixed = bind_rows(IMG_joined, FIRE_joined, LOMA_joined)
cat("=== TOTAL (all sources) ===\n")
cat("Rows:", nrow(TOTAL_MAGs_fixed), "| unique ids:", n_distinct(TOTAL_MAGs_fixed$id),
    "| unique guilds:", n_distinct(TOTAL_MAGs_fixed$guild), "\n")
cat("Rows with any NA completeness/contamination:",
    sum(is.na(TOTAL_MAGs_fixed$completeness) | is.na(TOTAL_MAGs_fixed$contamination)), "\n\n")

write.csv(TOTAL_MAGs_fixed, file = "Intermediate_Results/TOTAL_MAGs_fixed.csv")

# Guild-level aggregation
TOTAL_MAGs_guild = TOTAL_MAGs_fixed %>%
  group_by(guild) %>%
  summarise(completeness_mean = mean(completeness, na.rm = TRUE),
            contamination_mean = mean(contamination, na.rm = TRUE),
            GC_count_mean = mean(GC_count, na.rm = TRUE))
cat("Guild-level completeness/contamination/GC table - guilds:", nrow(TOTAL_MAGs_guild), "\n")
cat("(compare to MAG_gen_trait_v2's 1342 guilds)\n")

# Merge with MAG_gen_trait_v2

MAG_gen_trait_v3 = merge(MAG_gen_trait_v2, TOTAL_MAGs_guild, by = "guild")

write.csv(MAG_gen_trait_v3, file = "Intermediate_Results/MAG_gen_trait.csv")

# Adding CUE ----

sim_10000 = read.csv("Output_Data/my_MAGs_BGE_glucose_10000.csv")

# Map MAG-level CUE to guild via the real id
guild_lookup = total_genes.guild.940_MAG %>% select(id, guild)

cue_with_guild = sim_10000 %>%
  left_join(guild_lookup, by = c("MAG_id" = "id"))

# Guild-level mean CUE (BGE)
CUE_guild = cue_with_guild %>%
  filter(!is.na(guild)) %>%
  group_by(guild) %>%
  summarise(CUE_mean = mean(BGE, na.rm = TRUE),
            CUE_n_mags = n())

# Merge into the running MAG base
MAG_gen_trait_v4 = merge(MAG_gen_trait_v3, CUE_guild, by = "guild", all.x = TRUE)

write.csv(MAG_gen_trait_v4, file = "Intermediate_Results/MAG_gen_trait.csv")

# Adding A and S traits ----

MAG_gen_trait_v5 = MAG_gen_trait_v4 %>%
  mutate(
    S_traits  = rowSums(across(c(pH, temp, biofilm, osmolyte)), na.rm = TRUE),
    A_traits  = rowSums(across(c(tranport.total, Protein, CAZy)), na.rm = TRUE),
    A_enzymes = rowSums(across(c(Protein, CAZy)), na.rm = TRUE),
    A_S       = A_traits / S_traits
  )

write.csv(MAG_gen_trait_v5, file = "Intermediate_Results/MAG_gen_trait.csv")

# Normalizing by genome size ----

MAG_gen_trait_v6 = MAG_gen_trait_v5 %>%
  mutate(across(
    c(amino.transport, pH, temp, biofilm, osmolyte,
      GH.trasnport, tranport.total, Protein, CAZy,
      S_traits, A_traits, A_enzymes),
    ~ .x / genome.size,          # uses the piped df's own genome.size - safe against future filters
    .names = "{.col}_nor"
  ))

write.csv(MAG_gen_trait_v6, file = "Intermediate_Results/MAG_gen_trait_normalized.csv")

# Additional preparations for ISOLATES ----

# Adding ogt and mgt ----

isolates_microtrait = read_csv("Input_Data/SOIL_ISOLATES/dement_isolates_CUE.csv")
my_isolates_BGE_glucose = read_csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/DEBmicroTrait.jl-main/DEBmicroTrait.jl-main/files/my_isolates_BGE_glucose.csv")

guild_lookup_iso = total_genes.guild.940_ISO %>% select(id, guild, genome.size)

guild_n_iso = total_genes.guild.940_ISO %>% count(guild, name = "n_mags")

#ogt/mgt from isolates_microtrait 
iso_with_guild = left_join(guild_lookup_iso,
                           isolates_microtrait %>% select(id, mingentime, optimumT),
                           by = "id")

# CUE joined separately
iso_with_cue = left_join(guild_lookup_iso %>% select(id, guild),
                         my_isolates_BGE_glucose %>% select(isolate_id, BGE),
                         by = c("id" = "isolate_id"))

# Aggregate growth traits (genome.size + ogt + mgt together)
ISO_growth = iso_with_guild %>%
  group_by(guild) %>%
  summarise(
    genome.size_check = mean(genome.size, na.rm = TRUE),
    mgt = mean(mingentime, na.rm = TRUE),
    ogt = mean(optimumT, na.rm = TRUE)
  ) %>%
  mutate(mgr = 1/mgt)
cat("Guild-level growth table - guilds:", nrow(ISO_growth), "\n\n")

# Aggregate CUE separately (own simulation)
ISO_CUE = iso_with_cue %>%
  group_by(guild) %>%
  summarise(
    CUE_mean = mean(BGE, na.rm = TRUE),
    CUE_n_isolates = sum(!is.na(BGE))
  )

# Merge both into checkpoint-1 base
Isolates_gen_trait_v2 = merge(Isolates_gen_trait,
                              ISO_growth %>% select(guild, mgt, ogt, mgr),
                              by = "guild")
Isolates_gen_trait_v2 = merge(Isolates_gen_trait_v2,
                              ISO_CUE %>% select(guild, CUE_mean, CUE_n_isolates),
                              by = "guild", all.x = TRUE)

Isolates_gen_trait_v2 = merge(Isolates_gen_trait_v2, guild_n_iso, by = "guild")

cat("Guilds before:", n_distinct(Isolates_gen_trait$guild), "\n")
cat("Guilds after: ", n_distinct(Isolates_gen_trait_v2$guild), "\n")
cat("Any NA n_mags:", sum(is.na(Isolates_gen_trait_v2$n_mags)), "\n")

write.csv(Isolates_gen_trait_v2, file = "Intermediate_Results/Isolates_gen_trait.csv")

# Add S and A traits

Isolates_gen_trait_v3 = Isolates_gen_trait_v2 %>%
  mutate(
    S_traits  = rowSums(across(c(pH, temp, biofilm, osmolyte)), na.rm = TRUE),
    A_traits  = rowSums(across(c(tranport.total, Protein, CAZy)), na.rm = TRUE),
    A_enzymes = rowSums(across(c(Protein, CAZy)), na.rm = TRUE),
    A_S       = A_traits / S_traits
  )

write.csv(Isolates_gen_trait_v3, file = "Intermediate_Results/Isolates_gen_trait.csv")

# Normalization

Isolates_gen_trait_v4 = Isolates_gen_trait_v3 %>%
  mutate(across(
    c(amino.transport, pH, temp, biofilm, osmolyte,
      GH.trasnport, tranport.total, Protein, CAZy,
      S_traits, A_traits, A_enzymes),
    ~ .x / genome.size,
    .names = "{.col}_nor"
  ))

write.csv(Isolates_gen_trait_v4, file = "Intermediate_Results/Isolates_gen_trait_normalized.csv")

# FIGURES ----

# Call data

MAG_gen_trait_total = read.csv("Intermediate_Results/MAG_gen_trait.csv",dec=".")
ISO_gen_TOTAL_ab    = read.csv("Intermediate_Results/Isolates_gen_trait.csv",dec=".")

MAG_gen_trait_total_nor = read.csv("Intermediate_Results/MAG_gen_trait_normalized.csv",dec=".")
ISO_gen_TOTAL_ab_nor    = read.csv("Intermediate_Results/Isolates_gen_trait_normalized.csv",dec=".")

# Fig_2 ----

lm_weighted   = lm(S_traits ~ A_traits, data = MAG_gen_trait_total, weights = n_mags)
wcor = cov.wt(cbind(MAG_gen_trait_total$A_traits, MAG_gen_trait_total$S_traits),
              wt = MAG_gen_trait_total$n_mags, cor = TRUE)
r_weighted = wcor$cor[1, 2]

cat("Weighted Pearson r:", round(r_weighted, 3), "\n")
cat("Unweighted Pearson r (for comparison):", 
    round(cor(MAG_gen_trait_total$A_traits, MAG_gen_trait_total$S_traits), 3), "\n")
p_weighted = summary(lm_weighted)$coefficients[2, 4]
cat("p-value:", p_weighted, "\n")

Figure_Total_S_SA_MAGs = ggplot(data = MAG_gen_trait_total, 
                                aes(x = as.numeric(A_traits), 
                                    y = as.numeric(S_traits),
                                    color = A_S,
                                    weight = n_mags)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = 3, y = 80, hjust = 0,
           label = paste0("r = ", round(r_weighted, 2), ", p ", 
                          ifelse(p_weighted < 0.001, "< 0.001", paste0("= ", round(p_weighted, 3))))) +
  xlab("A Traits") + 
  ylab("S Traits") +
  scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5)) + 
  ylim(0,80) + xlim(0,250)
Figure_Total_S_SA_MAGs

lm_weighted   = lm(S_traits ~ A_traits, data = ISO_gen_TOTAL_ab, weights = n_mags)
wcor = cov.wt(cbind(ISO_gen_TOTAL_ab$A_traits, ISO_gen_TOTAL_ab$S_traits),
              wt = ISO_gen_TOTAL_ab$n_mags, cor = TRUE)
r_weighted = wcor$cor[1, 2]

cat("Weighted Pearson r:", round(r_weighted, 3), "\n")
cat("Unweighted Pearson r (for comparison):", 
    round(cor(ISO_gen_TOTAL_ab$A_traits, ISO_gen_TOTAL_ab$S_traits), 3), "\n")
p_weighted = summary(lm_weighted)$coefficients[2, 4]
cat("p-value:", p_weighted, "\n")

Figure_Total_S_SA_ISO   = ggplot(data = ISO_gen_TOTAL_ab, 
                                 aes(x = as.numeric(A_traits), 
                                     y = as.numeric(S_traits),
                                     color = A_S,
                                     weight = n_mags)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = 3, y = 80, hjust = 0,
           label = paste0("r = ", round(r_weighted, 2), ", p ", 
                          ifelse(p_weighted < 0.001, "< 0.001", paste0("= ", round(p_weighted, 3))))) +
  xlab("A Traits") + 
  ylab("S Traits") +
  scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5)) + 
  ylim(0,80) + xlim(0,250) 
Figure_Total_S_SA_ISO

Fig_2 = ggarrange(Figure_Total_S_SA_MAGs, Figure_Total_S_SA_ISO, 
                  labels = c("A","B"),
                  ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Fig_2
pdf("Output_Data/Figures/Fig_2.pdf",
    width=12,height=12*2/5)
print(Fig_2)
dev.off()

# Fig_S2 ----

lm_weighted   = lm(S_traits_nor ~ A_traits_nor, data = MAG_gen_trait_total_nor, weights = n_mags)
wcor = cov.wt(cbind(MAG_gen_trait_total_nor$A_traits_nor, MAG_gen_trait_total_nor$S_traits_nor),
              wt = MAG_gen_trait_total_nor$n_mags, cor = TRUE)
r_weighted = wcor$cor[1, 2]

cat("Weighted Pearson r:", round(r_weighted, 3), "\n")
cat("Unweighted Pearson r (for comparison):", 
    round(cor(MAG_gen_trait_total_nor$A_traits_nor, MAG_gen_trait_total_nor$S_traits_nor), 3), "\n")
p_weighted = summary(lm_weighted)$coefficients[2, 4]
cat("p-value:", p_weighted, "\n")

Figure_Total_S_SA_MAG_nor = ggplot(data = MAG_gen_trait_total_nor, 
                                   aes(x = as.numeric(A_traits_nor),
                                       y = as.numeric(S_traits_nor),
                                       color = A_S,
                                       weight = n_mags)) + 
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = min(MAG_gen_trait_total_nor$A_traits_nor, na.rm=TRUE), 
           y = max(MAG_gen_trait_total_nor$S_traits_nor, na.rm=TRUE), 
           hjust = 0, vjust = 1,
           label = paste0("r = ", round(r_weighted, 2), ", p ", 
                          ifelse(p_weighted < 0.001, "< 0.001", paste0("= ", round(p_weighted, 3))))) +
  xlab("A Traits") + 
  ylab("S Traits") +
  theme() + scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5))
Figure_Total_S_SA_MAG_nor

lm_weighted   = lm(S_traits_nor ~ A_traits_nor, data = ISO_gen_TOTAL_ab_nor, weights = n_mags)
wcor = cov.wt(cbind(ISO_gen_TOTAL_ab_nor$A_traits_nor, ISO_gen_TOTAL_ab_nor$S_traits_nor),
              wt = ISO_gen_TOTAL_ab_nor$n_mags, cor = TRUE)
r_weighted = wcor$cor[1, 2]

cat("Weighted Pearson r:", round(r_weighted, 3), "\n")
cat("Unweighted Pearson r (for comparison):", 
    round(cor(ISO_gen_TOTAL_ab_nor$A_traits_nor, ISO_gen_TOTAL_ab_nor$S_traits_nor), 3), "\n")
p_weighted = summary(lm_weighted)$coefficients[2, 4]
cat("p-value:", p_weighted, "\n")

Figure_Total_S_SA_ISO_nor   = ggplot(data = ISO_gen_TOTAL_ab_nor, 
                                     aes(x = as.numeric(A_traits_nor),
                                         y = as.numeric(S_traits_nor),
                                         color = A_S,
                                         weight = n_mags)) + 
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = min(ISO_gen_TOTAL_ab_nor$A_traits_nor, na.rm=TRUE), 
           y = max(ISO_gen_TOTAL_ab_nor$S_traits_nor, na.rm=TRUE), 
           hjust = 0, vjust = 1,
           label = paste0("r = ", round(r_weighted, 2), ", p ", 
                          ifelse(p_weighted < 0.001, "< 0.001", paste0("= ", round(p_weighted, 3))))) +
  xlab("A Traits") + 
  ylab("S Traits") +
  theme() + scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5))
Figure_Total_S_SA_ISO_nor

Figure_S2 = ggarrange(Figure_Total_S_SA_MAG_nor, Figure_Total_S_SA_ISO_nor, 
                      labels = c("A","B"),
                      ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Figure_S2
pdf("Output_Data/Figures/Figure_S2.pdf",
    width=12,height=12*2/5)
print(Figure_S2)
dev.off()

# MAGs ----

# CAZy----

Figure_CAZy = ggplot(data = MAG_gen_trait_total, 
                     aes(x = as.numeric(log(genome.size)), 
                         y = as.numeric(log(CAZy + 1)),
                         weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(CAZy)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.8,3.5)
Figure_CAZy

# Protein enzyme ----

Figure_Protein   = ggplot(data = MAG_gen_trait_total, 
                          aes(x = as.numeric(log(genome.size)), 
                              y = as.numeric(log(Protein + 1)),
                              weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Protein)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.8,4)
Figure_Protein

# Transport total Transporters ----

Figure_tranport.total   = ggplot(data = MAG_gen_trait_total, 
                                 aes(x = as.numeric(log(genome.size)), 
                                     y = as.numeric(log(tranport.total + 1)),
                                     weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.5,5.5)
Figure_tranport.total

# GH total Transporters ----

Figure_GH.trasnport   = ggplot(data = MAG_gen_trait_total, 
                               aes(x = as.numeric(log(genome.size)), 
                                   y = as.numeric(log(GH.trasnport + 1)),
                                   weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(GH transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.5,4.0)
Figure_GH.trasnport

# Amino total Transporters ----

Figure_amino.transport   = ggplot(data = MAG_gen_trait_total, 
                                  aes(x = as.numeric(log(genome.size)), 
                                      y = as.numeric(log(amino.transport + 1)),
                                      weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Amino transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0,3.5)
Figure_amino.transport

# Osmolytes ----

Figure_osmolyte  = ggplot(data = MAG_gen_trait_total, 
                          aes(x = as.numeric(log(genome.size)), 
                              y = as.numeric(log(osmolyte + 1)),
                              weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Osmolytes)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.7,3.5)
Figure_osmolyte

# Biofilm ----

Figure_biofilm  = ggplot(data = MAG_gen_trait_total, 
                         aes(x = as.numeric(log(genome.size)), 
                             y = as.numeric(log(biofilm + 1)),
                             weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Biofilm)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.2,2.8)
Figure_biofilm

# Heat Tolerance ----

Figure_temp  = ggplot(data = MAG_gen_trait_total, 
                      aes(x = as.numeric(log(genome.size)), 
                          y = as.numeric(log(temp + 1)),
                          weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Temperature Tolerance)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.5,2.6)
Figure_temp

# pH Tolerance Transporters ----

Figure_pH  = ggplot(data = MAG_gen_trait_total, 
                    aes(x = as.numeric(log(genome.size)), 
                        y = as.numeric(log(pH + 1)),
                        weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(pH Tolerance)") +
  theme(legend.position="none") + xlim(13,16.5) # + ylim(-2.7,2.85)
Figure_pH

# Maximum Growth Rate ----

Figure_mgr = ggplot(data = MAG_gen_trait_total, 
                    aes(x = as.numeric(log(genome.size)), 
                        y = as.numeric(log(mgr)),
                        weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(mgr)") +
  theme(legend.position="none") + xlim(13,16.5)
Figure_mgr

# Optimum Growth Temperature ----

Figure_ogt = ggplot(data = MAG_gen_trait_total %>% 
                      filter(is.finite(ogt)),  # 9 guilds excluded: entirely/mostly Loma MAGs, 
                    # which lack ogt data (consistent with other 
                    # documented Loma metadata gaps)
                    aes(x = as.numeric(log(genome.size)), 
                        y = as.numeric(log(ogt)),
                        weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + 
  xlab("Log (Genome Size)") + 
  ylab("Log (ogt)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(2.9,4.3)
Figure_ogt

# Yield ----

guild_lookup = total_genes.guild.940_MAG %>% select(id, guild)
cue_with_guild = sim_10000 %>%
  left_join(guild_lookup, by = c("MAG_id" = "id")) %>%
  filter(!is.na(guild))

guild_cue_without = cue_with_guild %>%
  filter(BGE <= 0.848) %>%
  group_by(guild) %>%
  summarise(CUE_mean = mean(BGE, na.rm = TRUE),
            genome.size = mean(genomesize, na.rm = TRUE),
            n_cue_mags = n())          # <-- how many MAGs survived the filter for THIS guild

guild_cue_without = guild_cue_without %>%
  mutate(CUE_logit = log(CUE_mean / (1 - CUE_mean)))

cat("Any non-finite logit values:", sum(!is.finite(guild_cue_without$CUE_logit)), "\n\n")
cat("Any single-MAG guilds still present:", sum(guild_cue_without$n_cue_mags == 1), "\n\n")

Figure_yield_logit = ggplot(data = guild_cue_without,
                            aes(x = log(genome.size), y = CUE_logit,
                                weight = n_cue_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") +
  xlab("Log (Genome Size)") +
  ylab("Logit (Yield)") +
  theme(legend.position="none") +
  coord_cartesian(xlim = c(13, 16.5))
Figure_yield_logit

# ISOLATES ----

# CAZy----

Figure_CAZy.total.ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                               aes(x = as.numeric(log(genome.size)), 
                                   y = as.numeric(log(CAZy + 1)),
                                   weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(CAZy)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(-0.8,3.5))
Figure_CAZy.total.ISO

# Protein enzyme ----

Figure_Protein.enzyme.total.ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                                         aes(x = as.numeric(log(genome.size)), 
                                             y = as.numeric(log(Protein + 1)),
                                             weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Protein)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(0.8,4))
Figure_Protein.enzyme.total.ISO

# Transport total Transporters ----

Figure_tranport.total_ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                                   aes(x = as.numeric(log(genome.size)), 
                                       y = as.numeric(log(tranport.total + 1)),
                                       weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Transporters)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(-0.5,5.5))
Figure_tranport.total_ISO

# GH transporter ----

Figure_GH.transporter = ggplot(data = ISO_gen_TOTAL_ab, 
                               aes(x = as.numeric(log(genome.size)), 
                                   y = as.numeric(log(GH.trasnport + 1)),
                                   weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5, 
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(GH transporters)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(-3.5,4.0))
Figure_GH.transporter

# Amino transporter ----

Figure_Amino.transporter = ggplot(data = ISO_gen_TOTAL_ab, 
                                  aes(x = as.numeric(log(genome.size)), 
                                      y = as.numeric(log(amino.transport + 1)),
                                      weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Amino transporters)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(0,3.5))
Figure_Amino.transporter

# Osmolyte ----

Figure_osmolyte_ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                             aes(x = as.numeric(log(genome.size)), 
                                 y = as.numeric(log(osmolyte + 1)),
                                 weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Osmolytes)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(-2.7,3.5))
Figure_osmolyte_ISO

# Biofilm ----

Figure_Biofilm_ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                            aes(x = as.numeric(log(genome.size)), 
                                y = as.numeric(log(biofilm + 1)),
                                weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Biofilm)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(-3.2,2.8))
Figure_Biofilm_ISO

# Heat Tolerance ----

Figure_Temp.Tol_ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                             aes(x = as.numeric(log(genome.size)), 
                                 y = as.numeric(log(temp + 1)),
                                 weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(Temperature Tolerance)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(0.5,2.6))
Figure_Temp.Tol_ISO

# pH Tolerance ----

Figure_pH.Tol_ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                           aes(x = as.numeric(log(genome.size)), 
                               y = as.numeric(log(pH + 1)),
                               weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4,
               label.x = "left", label.y = "bottom") + 
  xlab("Log(Genome Size)") + 
  ylab("Log(pH Tolerance)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(-2.7,2.85))
Figure_pH.Tol_ISO

# Minimum Generation Time ----

Figure_mgr_ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                        aes(x = as.numeric(log(genome.size)), 
                            y = as.numeric(log(mgt)),
                            weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(mgr)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(-2.55,3.45))
Figure_mgr_ISO

# Optimum Growth Temperature ----

Figure_OGT_ISO = ggplot(data = ISO_gen_TOTAL_ab, 
                        aes(x = as.numeric(log(genome.size)), 
                            y = as.numeric(log(ogt)),
                            weight = n_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4) + 
  xlab("Log(Genome Size)") + 
  ylab("Log(ogt)") +
  theme(legend.position="none") + coord_cartesian(xlim = c(13,16.5), ylim = c(2.9,4.))
Figure_OGT_ISO

# Yield ----

guild_lookup_iso = total_genes.guild.940_ISO %>% select(id, guild)

iso_cue_with_guild = my_isolates_BGE_glucose %>%
  left_join(guild_lookup_iso, by = c("isolate_id" = "id")) %>%
  filter(!is.na(guild))

# Genome size skew check for the CORRECTED threshold
iso_cue_with_guild %>%
  mutate(near_ceiling = BGE > 0.90) %>%
  filter(!is.na(near_ceiling)) %>%
  group_by(near_ceiling) %>%
  summarise(mean_genome_size = mean(genomesize, na.rm=TRUE), n = n())

guild_cue_iso = iso_cue_with_guild %>%
  filter(BGE <= 0.90) %>%
  group_by(guild) %>%
  summarise(CUE_mean = mean(BGE, na.rm = TRUE),
            genome.size = mean(genomesize, na.rm = TRUE),
            n_cue_mags = n())

guild_cue_iso = guild_cue_iso %>%
  mutate(CUE_logit = log(CUE_mean / (1 - CUE_mean)))

cat("Any non-finite logit values:", sum(!is.finite(guild_cue_iso$CUE_logit)), "\n")
cat("Any single-isolate guilds:", sum(guild_cue_iso$n_cue_mags == 1), "\n")
cat("Guilds remaining:", nrow(guild_cue_iso), "\n\n")

Figure_yield_logit_ISO = ggplot(data = guild_cue_iso,
                                aes(x = log(genome.size), y = CUE_logit,
                                    weight = n_cue_mags)) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") +
  xlab("Log (Genome Size)") +
  ylab("Logit (Yield)") +
  theme(legend.position="none") +
  coord_cartesian(xlim = c(13, 16.5))
Figure_yield_logit_ISO

# Fig_3 ----

Fig_3 = ggarrange(Figure_CAZy, Figure_Protein,Figure_CAZy.total.ISO,
                  Figure_Protein.enzyme.total.ISO,
                  labels = c("A","B","C","D"),ncol = 2, nrow = 2) + 
  theme(panel.background = element_blank())

Fig_3
pdf("Output_Data/Figures/Fig_3.pdf",
    width=14,height=14*5/5)
print(Fig_3)
dev.off()

# Fig_4 ----

Fig_4 = ggarrange(Figure_tranport.total, Figure_GH.trasnport,
                  Figure_amino.transport,Figure_tranport.total_ISO,
                  Figure_GH.transporter,Figure_Amino.transporter,
                  labels = c("A","B","C","D","E","F"),ncol = 3, nrow = 2) + 
  theme(panel.background = element_blank())

Fig_4
pdf("Output_Data/Figures/Fig_4.pdf",
    width=15,height=15*4/5)
print(Fig_4)
dev.off()

# Fig_5 ----

Fig_5 = ggarrange(Figure_osmolyte,Figure_biofilm,Figure_temp,
                  Figure_pH,Figure_osmolyte_ISO,Figure_Biofilm_ISO,
                  Figure_Temp.Tol_ISO,Figure_pH.Tol_ISO,
                  labels = c("A","B","C","D","E","F","G","H"),ncol = 4, 
                  nrow = 2) + theme(panel.background = element_blank())
Fig_5

pdf("Output_Data/Figures/Fig_5.pdf",
    width=16,height=16*3/5)
print(Fig_5)
dev.off()

# Fig_6 ----

Fig_6 = ggarrange(Figure_mgr, Figure_mgr_ISO, Figure_yield_logit,
                  Figure_yield_logit_ISO,
                  labels = c("A","B","C","D"),ncol = 2, nrow = 2) + 
  theme(panel.background = element_blank())

Fig_6

pdf("Output_Data/Figures/Fig_6.pdf",
    width=12,height=12*5/5)
print(Fig_6)
dev.off()

# Fig_7 ----

Fig_7 = ggarrange(Figure_ogt, Figure_OGT_ISO,labels = c("A","B"),
                  ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Fig_7

pdf("Output_Data/Figures/Fig_7.pdf",
    width=12,height=12*1/2)
print(Fig_7)
dev.off()

# Fig_S3 ----

# Weighted correlation matrix

total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv", dec=".")
sim_10000 = read.csv("Output_Data/my_MAGs_BGE_glucose_10000.csv")

guild_lookup = total_genes.guild.940_MAG %>% select(id, guild)

cue_with_guild = sim_10000 %>%
  left_join(guild_lookup, by = c("MAG_id" = "id")) %>%
  filter(!is.na(guild))

guild_cue_without = cue_with_guild %>%
  filter(BGE <= 0.848) %>%
  group_by(guild) %>%
  summarise(CUE_mean = mean(BGE, na.rm = TRUE),
            genome.size = mean(genomesize, na.rm = TRUE),
            n_cue_mags = n())

# Step 1: merge full-trait guild table with the CORRECTED CUE-only guild table.
# IMPORTANT: MAG_gen_trait_total already has CUE_mean/CUE_n_mags columns, but
# these were computed BEFORE the BGE-artifact fix (BR≈0 rows excluded at
# BGE > 0.848). Confirmed 744/1342 guilds differ meaningfully between the two
# versions (some by up to 0.85!) - so we DROP the stale columns and replace
# them with the corrected ones from guild_cue_without, rather than adding
# alongside (which would invite accidentally using the wrong one later).

MAG_corr_data = MAG_gen_trait_total %>%
  select(-CUE_mean, -CUE_n_mags) %>%
  left_join(guild_cue_without %>% select(guild, CUE_mean, n_cue_mags), by = "guild")

cat("Guilds before merge:", n_distinct(MAG_gen_trait_total$guild), "\n")
cat("Guilds after merge: ", n_distinct(MAG_corr_data$guild), "(should be unchanged)\n")
cat("Guilds with NA CUE:  ", sum(is.na(MAG_corr_data$CUE_mean)), "\n\n")

# Readable column names

MAG_corr_data = MAG_corr_data %>%
  rename(
    "Guild" = guild,
    "Genome Size" = genome.size,
    "Amino-transporter" = amino.transport,
    "pH-Tol" = pH,
    "Temp-Tol" = temp,
    "Biofilm" = biofilm,
    "Osmolyte" = osmolyte,
    "GH-transporter" = GH.trasnport,
    "Total-transporter" = tranport.total,
    "Protein-enzyme" = Protein,
    "CAZy" = CAZy,
    "OGT" = ogt,
    "MGR" = mgr,
    "Completeness" = completeness_mean,
    "Contamination" = contamination_mean,
    "GC-count" = GC_count_mean,
    "Yield" = CUE_mean
  )

# Step 2: weighted correlation function

# WEIGHTING NOTE: MAG_corr_data has TWO possible weight columns:
#   - n_mags: total MAGs in the guild (used for all non-CUE trait pairs)
#   - n_cue_mags: MAGs that survived the BGE artifact filter for THIS guild's
#     CUE (only meaningful/available for pairs involving CUE_mean)
# The function below picks the right one automatically based on whether
# CUE_mean is one of the two variables being correlated.

weighted_cor = function(x, y, w) {
  ok = complete.cases(x, y, w)
  if (sum(ok) < 3) return(list(r = NA, p = NA))
  x = x[ok]; y = y[ok]; w = w[ok]
  r = cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
  # p-value via weighted lm - mathematically equivalent to testing correlation significance
  p = tryCatch({
    summary(lm(y ~ x, weights = w))$coefficients[2, 4]
  }, error = function(e) NA)
  list(r = r, p = p)
}

sig_stars = function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

# Step 3: custom ggpairs panel function - weighted correlation,

my_fn_weighted = function(data, mapping, ...) {
  x_col = rlang::as_label(mapping$x)
  y_col = rlang::as_label(mapping$y)
  x = eval_data_col(data, mapping$x)
  y = eval_data_col(data, mapping$y)
  
  # Use n_cue_mags if Yield is involved in this pair, otherwise n_mags
  # (checking against "Yield" - the renamed display column, not "CUE_mean")
  w = if (x_col == "Yield" | y_col == "Yield") data$n_cue_mags else data$n_mags
  
  result = weighted_cor(x, y, w)
  corr = result$r
  stars = sig_stars(result$p)
  
  colFn = colorRampPalette(c("red", "white", "blue"), interpolate = 'spline')
  fill = if (is.na(corr)) "grey90" else colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  
  label_text = if (is.na(corr)) "NA" else sprintf("Corr: %.3f%s", corr, stars)
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, size = 3.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = fill)) +
    xlim(0,1) + ylim(0,1)
}
# Step 4: build the plot
# IMPORTANT: n_mags AND n_cue_mags must both be present in MAG_corr_data
# (not listed in trait_cols) so my_fn_weighted can access them via
# data$n_mags / data$n_cue_mags even though they aren't plotted themselves.

trait_cols = c("Yield","OGT","MGR","Amino-transporter","GH-transporter",
               "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
               "Biofilm", "Osmolyte", "Completeness","Contamination","GC-count")

Figure_S3_weighted = ggpairs(MAG_corr_data, 
                             columns = trait_cols,
                             upper = list(continuous = my_fn_weighted),
                             lower = list(continuous = "smooth"))
Figure_S3_weighted

pdf("Output_Data/Figures/Figure_S3.pdf",
    width=4*4,height=4*4)
print(Figure_S3_weighted)
dev.off()

# Fig_S7 ----

# Step 1: merge full-trait guild table with the CORRECTED CUE-only guild table.
MAG_corr_data = MAG_gen_trait_total_nor %>%
  select(-CUE_mean, -CUE_n_mags) %>%
  left_join(guild_cue_without %>% select(guild, CUE_mean, n_cue_mags), by = "guild")

cat("Guilds before merge:", n_distinct(MAG_gen_trait_total_nor$guild), "\n")
cat("Guilds after merge: ", n_distinct(MAG_corr_data$guild), "(should be unchanged)\n")
cat("Guilds with NA CUE:  ", sum(is.na(MAG_corr_data$CUE_mean)), "\n\n")

# Readable column names

MAG_corr_data = MAG_corr_data %>%
  rename(
    "Guild" = guild,
    "Genome Size" = genome.size,
    "Amino-transporter-n" = amino.transport_nor,
    "pH-Tol-n" = pH_nor,
    "Temp-Tol-n" = temp_nor,
    "Biofilm-n" = biofilm_nor,
    "Osmolyte-n" = osmolyte_nor,
    "GH-transporter-n" = GH.trasnport_nor,
    "Total-transporter-n" = tranport.total_nor,
    "Protein-enzyme-n" = Protein_nor,
    "CAZy-n" = CAZy_nor,
    "OGT" = ogt,
    "MGR" = mgr,
    "Completeness" = completeness_mean,
    "Contamination" = contamination_mean,
    "GC-count" = GC_count_mean,
    "Yield" = CUE_mean
  )

# Step 2: weighted correlation function

weighted_cor = function(x, y, w) {
  ok = complete.cases(x, y, w)
  if (sum(ok) < 3) return(list(r = NA, p = NA))
  x = x[ok]; y = y[ok]; w = w[ok]
  r = cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
  # p-value via weighted lm - mathematically equivalent to testing correlation significance
  p = tryCatch({
    summary(lm(y ~ x, weights = w))$coefficients[2, 4]
  }, error = function(e) NA)
  list(r = r, p = p)
}

sig_stars = function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

# Step 3: custom ggpairs panel function - weighted correlation,

my_fn_weighted = function(data, mapping, ...) {
  x_col = rlang::as_label(mapping$x)
  y_col = rlang::as_label(mapping$y)
  x = eval_data_col(data, mapping$x)
  y = eval_data_col(data, mapping$y)
  
  # Use n_cue_mags if Yield is involved in this pair, otherwise n_mags
  # (checking against "Yield" - the renamed display column, not "CUE_mean")
  w = if (x_col == "Yield" | y_col == "Yield") data$n_cue_mags else data$n_mags
  
  result = weighted_cor(x, y, w)
  corr = result$r
  stars = sig_stars(result$p)
  
  colFn = colorRampPalette(c("red", "white", "blue"), interpolate = 'spline')
  fill = if (is.na(corr)) "grey90" else colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  
  label_text = if (is.na(corr)) "NA" else sprintf("Corr: %.3f%s", corr, stars)
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, size = 3.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = fill)) +
    xlim(0,1) + ylim(0,1)
}
# Step 4: build the plot

trait_cols = c("Yield","OGT","MGR","Amino-transporter-n","GH-transporter-n",
               "Total-transporter-n","Protein-enzyme-n","CAZy-n","pH-Tol-n","Temp-Tol-n", 
               "Biofilm-n", "Osmolyte-n", "Completeness","Contamination","GC-count")

Figure_S3A = ggpairs(MAG_corr_data, 
                     columns = trait_cols,
                     upper = list(continuous = my_fn_weighted),
                     lower = list(continuous = "smooth"))
Figure_S3A

pdf("Output_Data/Figures/Figure_S3A.pdf",
    width=4*4,height=4*4)
print(Figure_S3A)
dev.off()

# Fig_S4 ----

# Step 1: check for the same stale-CUE problem we found in MAGs
cat("Does ISO_gen_TOTAL_ab have old CUE columns?\n")
cat("  CUE_mean present:", "CUE_mean" %in% colnames(ISO_gen_TOTAL_ab), "\n")
cat("  n_mags present:  ", "n_mags" %in% colnames(ISO_gen_TOTAL_ab), "\n\n")

# Step 2: merge in the CORRECTED isolate Yield (built with 0.90 threshold)

ISO_corr_data = ISO_gen_TOTAL_ab %>%
  select(-any_of(c("CUE_mean", "CUE_n_isolates"))) %>%   # drop stale columns if present
  left_join(guild_cue_iso %>% select(guild, CUE_mean, n_cue_mags), by = "guild")

cat("Guilds before merge:", n_distinct(ISO_gen_TOTAL_ab$guild), "\n")
cat("Guilds after merge: ", n_distinct(ISO_corr_data$guild), "(should be unchanged)\n")
cat("Guilds with NA Yield:", sum(is.na(ISO_corr_data$CUE_mean)), "\n\n")

# Step 3: weighted correlation + significance (same functions as MAGs)
weighted_cor = function(x, y, w) {
  ok = complete.cases(x, y, w)
  if (sum(ok) < 3) return(list(r = NA, p = NA))
  x = x[ok]; y = y[ok]; w = w[ok]
  r = cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
  p = tryCatch({
    summary(lm(y ~ x, weights = w))$coefficients[2, 4]
  }, error = function(e) NA)
  list(r = r, p = p)
}

sig_stars = function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

my_fn_weighted_iso = function(data, mapping, ...) {
  x_col = rlang::as_label(mapping$x)
  y_col = rlang::as_label(mapping$y)
  x = eval_data_col(data, mapping$x)
  y = eval_data_col(data, mapping$y)
  
  # Use n_cue_mags if Yield is involved, otherwise n_mags
  w = if (x_col == "Yield" | y_col == "Yield") data$n_cue_mags else data$n_mags
  
  result = weighted_cor(x, y, w)
  corr = result$r
  stars = sig_stars(result$p)
  
  colFn = colorRampPalette(c("red", "white", "blue"), interpolate = 'spline')
  fill = if (is.na(corr)) "grey90" else colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  label_text = if (is.na(corr)) "NA" else sprintf("Corr: %.3f%s", corr, stars)
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, size = 3.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = fill)) +
    xlim(0,1) + ylim(0,1)
}

# Step 4: readable column names
# NOTE: n_mags and n_cue_mags NOT renamed - accessed programmatically.
# NOTE: no Completeness/Contamination/GC-count - see header note above.

ISO_corr_data = ISO_corr_data %>%
  rename(
    "Guild" = guild,
    "Genome Size" = genome.size,
    "Amino-transporter" = amino.transport,
    "pH-Tol" = pH,
    "Temp-Tol" = temp,
    "Biofilm" = biofilm,
    "Osmolyte" = osmolyte,
    "GH-transporter" = GH.trasnport,
    "Total-transporter" = tranport.total,
    "Protein-enzyme" = Protein,
    "CAZy" = CAZy,
    "OGT" = ogt,
    "MGR" = mgr,
    "Yield" = CUE_mean
  )

col_order = c("Yield","OGT","MGR","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte")

# Step 5: build the plot

Figure_S3_ISO_weighted = ggpairs(ISO_corr_data, 
                                 columns = col_order,
                                 upper = list(continuous = my_fn_weighted_iso),
                                 lower = list(continuous = "smooth"))
Figure_S3_ISO_weighted

pdf("Output_Data/Figures/Figure_S4.pdf",
    width=4*4,height=4*4)
print(Figure_S3_ISO_weighted)
dev.off()

# Fig_S8 ----

# Step 1: check for the same stale-CUE problem we found in MAGs
cat("Does ISO_gen_TOTAL_ab_nor have old CUE columns?\n")
cat("  CUE_mean present:", "CUE_mean" %in% colnames(ISO_gen_TOTAL_ab_nor), "\n")
cat("  n_mags present:  ", "n_mags" %in% colnames(ISO_gen_TOTAL_ab_nor), "\n\n")

# Step 2: merge in the CORRECTED isolate Yield (built with 0.90 threshold)

ISO_corr_data = ISO_gen_TOTAL_ab_nor %>%
  select(-any_of(c("CUE_mean", "CUE_n_isolates"))) %>%   # drop stale columns if present
  left_join(guild_cue_iso %>% select(guild, CUE_mean, n_cue_mags), by = "guild")

cat("Guilds before merge:", n_distinct(ISO_gen_TOTAL_ab_nor$guild), "\n")
cat("Guilds after merge: ", n_distinct(ISO_corr_data$guild), "(should be unchanged)\n")
cat("Guilds with NA Yield:", sum(is.na(ISO_corr_data$CUE_mean)), "\n\n")

# Step 3: weighted correlation + significance (same functions as MAGs)
weighted_cor = function(x, y, w) {
  ok = complete.cases(x, y, w)
  if (sum(ok) < 3) return(list(r = NA, p = NA))
  x = x[ok]; y = y[ok]; w = w[ok]
  r = cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
  p = tryCatch({
    summary(lm(y ~ x, weights = w))$coefficients[2, 4]
  }, error = function(e) NA)
  list(r = r, p = p)
}

sig_stars = function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

my_fn_weighted_iso = function(data, mapping, ...) {
  x_col = rlang::as_label(mapping$x)
  y_col = rlang::as_label(mapping$y)
  x = eval_data_col(data, mapping$x)
  y = eval_data_col(data, mapping$y)
  
  # Use n_cue_mags if Yield is involved, otherwise n_mags
  w = if (x_col == "Yield" | y_col == "Yield") data$n_cue_mags else data$n_mags
  
  result = weighted_cor(x, y, w)
  corr = result$r
  stars = sig_stars(result$p)
  
  colFn = colorRampPalette(c("red", "white", "blue"), interpolate = 'spline')
  fill = if (is.na(corr)) "grey90" else colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]
  label_text = if (is.na(corr)) "NA" else sprintf("Corr: %.3f%s", corr, stars)
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, size = 3.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = fill)) +
    xlim(0,1) + ylim(0,1)
}

# Step 4: readable column names
# NOTE: n_mags and n_cue_mags NOT renamed - accessed programmatically.
# NOTE: no Completeness/Contamination/GC-count - see header note above.

ISO_corr_data = ISO_corr_data %>%
  rename(
    "Guild" = guild,
    "Genome Size" = genome.size,
    "Amino-transporter-n" = amino.transport_nor,
    "pH-Tol-n" = pH_nor,
    "Temp-Tol-n" = temp_nor,
    "Biofilm-n" = biofilm_nor,
    "Osmolyte-n" = osmolyte_nor,
    "GH-transporter-n" = GH.trasnport_nor,
    "Total-transporter-n" = tranport.total_nor,
    "Protein-enzyme-n" = Protein_nor,
    "CAZy-n" = CAZy_nor,
    "OGT" = ogt,
    "MGR" = mgr,
    "Yield" = CUE_mean
  )

col_order = c("Yield","OGT","MGR","Amino-transporter-n","GH-transporter-n",
              "Total-transporter-n","Protein-enzyme-n","CAZy-n","pH-Tol-n","Temp-Tol-n", 
              "Biofilm-n", "Osmolyte-n")

# Step 5: build the plot

Figure_S4A = ggpairs(ISO_corr_data, 
                     columns = col_order,
                     upper = list(continuous = my_fn_weighted_iso),
                     lower = list(continuous = "smooth"))
Figure_S4A

pdf("Output_Data/Figures/Figure_S4A4.pdf",
    width=4*4,height=4*4)
print(Figure_S4A)
dev.off()

# Fig_S5 -----

# Load raw files

my_isolates_BGE_glucose = read.csv("Output_Data/my_isolates_BGE_glucose.csv")  # CONFIRM PATH
total_genes.guild.940_ISO = read.csv("Intermediate_Results/total_genes.guild.940_ISO.csv", dec=".")

meta_iso   = read.delim("Input_Data/SOIL_ISOLATES/metadata_2.tsv", sep="\t")
meta_iso   = as.data.frame(cbind(meta_iso$IMG.Genome.ID, meta_iso$Phylum,
                                 meta_iso$Class, meta_iso$Order, meta_iso$Family,
                                 meta_iso$Genus))
meta_iso.1 = read.delim("Input_Data/SOIL_ISOLATES/metadata_1.tsv", sep="\t")
meta_iso.1 = as.data.frame(cbind(meta_iso.1$IMG.Genome.ID, meta_iso.1$Phylum,
                                 meta_iso.1$Class, meta_iso.1$Order, meta_iso.1$Family,
                                 meta_iso.1$Genus))
final_meta = as.data.frame(rbind(meta_iso, meta_iso.1))
colnames(final_meta) = c("id","Phylum","Class","Order","Family","Genus")

# Isolate-level data: OUR OWN BGE + genome size, correct artifact filter
iso_cue_clean = my_isolates_BGE_glucose %>%
  filter(BGE <= 0.90) %>%              # correct, evidence-based threshold
  select(isolate_id, genomesize, BGE)

cat("Isolates before artifact filter:", nrow(my_isolates_BGE_glucose), "\n")
cat("Isolates after filter (BGE<=0.90):", nrow(iso_cue_clean), "\n\n")

# Attach taxonomy
final_meta = final_meta %>% mutate(id = as.numeric(id))

iso_cue_tax = iso_cue_clean %>%
  left_join(final_meta, by = c("isolate_id" = "id"))

cat("Isolates with taxonomy match:", sum(!is.na(iso_cue_tax$Phylum)), 
    "out of", nrow(iso_cue_tax), "\n\n")

# Attach guild (for Figure S9, which is guild-level, not taxonomy-level)
guild_lookup_iso = total_genes.guild.940_ISO %>% select(id, guild)

iso_cue_guild = iso_cue_clean %>%
  left_join(guild_lookup_iso, by = c("isolate_id" = "id")) %>%
  filter(!is.na(guild))

cat("Isolates with guild match:", nrow(iso_cue_guild), "\n\n")

# Figure

build_tax_panel = function(data, rank_col, title) {
  agg = data %>%
    group_by(.data[[rank_col]]) %>%
    summarise(genome_length_mean = mean(genomesize, na.rm = TRUE),
              CUE_mean = mean(BGE, na.rm = TRUE),
              n = n()) %>%
    filter(!is.na(.data[[rank_col]]), .data[[rank_col]] != "") %>%
    mutate(CUE_logit = log(CUE_mean / (1 - CUE_mean)))
  
  cat(title, "- non-finite logit values:", sum(!is.finite(agg$CUE_logit)), "\n")
  
  ggplot(data = agg, aes(x = log(genome_length_mean), y = CUE_logit, weight = n)) +
    geom_point() + theme_classic() + theme(text = element_text(size = 12)) +
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 3, 
                 label.y = "bottom", label.x = "left") +
    xlab("Log(Genome size)") +
    ylab("Logit(Yield)") +
    labs(title = title)
}

Figure_test_1.phylum = build_tax_panel(iso_cue_tax, "Phylum", "Phylum")
Figure_test_1.Class  = build_tax_panel(iso_cue_tax, "Class", "Class")
Figure_test_1.Order  = build_tax_panel(iso_cue_tax, "Order", "Order")
Figure_test_1.Family = build_tax_panel(iso_cue_tax, "Family", "Family")
Figure_test_1.Genus  = build_tax_panel(iso_cue_tax, "Genus", "Genus")

Figure_test_1.phylum
Figure_test_1.Class
Figure_test_1.Order
Figure_test_1.Family
Figure_test_1.Genus

# Combined figure
Figure_S5 = ggarrange(Figure_test_1.phylum, Figure_test_1.Class,
                      Figure_test_1.Order, Figure_test_1.Family, Figure_test_1.Genus,
                      labels = c("A","B","C","D","E"),
                      ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_S5

pdf("Output_Data/Figures/Figure_S5.pdf",
    width=14,height=14*3/5)
print(Figure_S5)
dev.off()

# Fig_S9 -----

# MAGs

# Load raw files
MAG_gen_trait         = read.csv("Intermediate_Results/MAG_gen_trait.csv",dec=".")
total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv", dec=".")
sim_10000 = read.csv("Output_Data/my_MAGs_BGE_glucose_10000.csv")

# S/A traits (from MAG_gen_trait, already confirmed 1342 guilds)
MAG_traits = MAG_gen_trait %>%
  mutate(
    S_traits = rowSums(across(c(pH, temp, biofilm, osmolyte)), na.rm = TRUE),
    A_traits = rowSums(across(c(tranport.total, Protein, CAZy)), na.rm = TRUE)
  ) %>%
  select(guild, S_traits, A_traits)

# n_mags (total MAG count per guild)
guild_n_mag = total_genes.guild.940_MAG %>% count(guild, name = "n_mags")

# Corrected Yield (BGE <= 0.848, real artifact threshold)
guild_lookup_mag = total_genes.guild.940_MAG %>% select(id, guild)
mag_cue_guild = sim_10000 %>%
  left_join(guild_lookup_mag, by = c("MAG_id" = "id")) %>%
  filter(!is.na(guild), BGE <= 0.848) %>%
  group_by(guild) %>%
  summarise(Yield = mean(BGE, na.rm = TRUE), n_cue_mags = n())

# Merge everything
MAG_S9_data = MAG_traits %>%
  left_join(guild_n_mag, by = "guild") %>%
  left_join(mag_cue_guild, by = "guild")

# Isolates

Isolates_gen_trait    = read.csv("Intermediate_Results/Isolates_gen_trait.csv",dec=".")
total_genes.guild.940_ISO = read.csv("Intermediate_Results/total_genes.guild.940_ISO.csv", dec=".")
my_isolates_BGE_glucose = read.csv("Output_Data/my_isolates_BGE_glucose.csv")  # CONFIRM PATH

# S/A traits
ISO_traits = Isolates_gen_trait %>%
  mutate(
    S_traits = rowSums(across(c(pH, temp, biofilm, osmolyte)), na.rm = TRUE),
    A_traits = rowSums(across(c(tranport.total, Protein, CAZy)), na.rm = TRUE)
  ) %>%
  select(guild, S_traits, A_traits)

# n_mags (total isolate count per guild)
guild_n_iso = total_genes.guild.940_ISO %>% count(guild, name = "n_mags")

# Corrected Yield (BGE <= 0.90, isolate-specific artifact threshold)
guild_lookup_iso = total_genes.guild.940_ISO %>% select(id, guild)
iso_cue_guild = my_isolates_BGE_glucose %>%
  left_join(guild_lookup_iso, by = c("isolate_id" = "id")) %>%
  filter(!is.na(guild), BGE <= 0.90) %>%
  group_by(guild) %>%
  summarise(Yield = mean(BGE, na.rm = TRUE), n_cue_mags = n())

# Merge everything
ISO_S9_data = ISO_traits %>%
  left_join(guild_n_iso, by = "guild") %>%
  left_join(iso_cue_guild, by = "guild")

# MAGs
MAG_S9_data = MAG_S9_data %>%
  mutate(AS_total = A_traits + S_traits)

MAG_S9_complete = MAG_S9_data %>% filter(!is.na(Yield))

wcor_MAG = cov.wt(cbind(MAG_S9_complete$AS_total, MAG_S9_complete$Yield),
                  wt = MAG_S9_complete$n_cue_mags, cor = TRUE)
r_weighted_MAG = wcor_MAG$cor[1, 2]

lm_weighted_MAG = lm(Yield ~ AS_total, data = MAG_S9_complete, weights = n_cue_mags)
p_weighted_MAG = summary(lm_weighted_MAG)$coefficients[2, 4]

Figure_S9_MAG = ggplot(data = MAG_S9_complete,
                       aes(x = AS_total, y = Yield, weight = n_cue_mags)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = min(MAG_S9_complete$AS_total, na.rm=TRUE), 
           y = min(MAG_S9_complete$Yield, na.rm=TRUE), hjust = 0, vjust = 0,
           label = paste0("r = ", round(r_weighted_MAG, 3), ", p ", 
                          ifelse(p_weighted_MAG < 0.001, "< 0.001", 
                                 paste0("= ", round(p_weighted_MAG, 3))))) +
  xlab("Total S+A gene counts per functional group") +
  ylab("Yield (CUE)") 
Figure_S9_MAG

# Isolates
ISO_S9_data = ISO_S9_data %>%
  mutate(AS_total = A_traits + S_traits)

ISO_S9_complete = ISO_S9_data %>% filter(!is.na(Yield))

wcor_ISO = cov.wt(cbind(ISO_S9_complete$AS_total, ISO_S9_complete$Yield),
                  wt = ISO_S9_complete$n_cue_mags, cor = TRUE)
r_weighted_ISO = wcor_ISO$cor[1, 2]

lm_weighted_ISO = lm(Yield ~ AS_total, data = ISO_S9_complete, weights = n_cue_mags)
p_weighted_ISO = summary(lm_weighted_ISO)$coefficients[2, 4]

Figure_S9_ISO = ggplot(data = ISO_S9_complete,
                       aes(x = AS_total, y = Yield, weight = n_cue_mags)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  annotate("text", x = min(ISO_S9_complete$AS_total, na.rm=TRUE), 
           y = min(ISO_S9_complete$Yield, na.rm=TRUE), hjust = 0, vjust = 0,
           label = paste0("r = ", round(r_weighted_ISO, 3), ", p ", 
                          ifelse(p_weighted_ISO < 0.001, "< 0.001", 
                                 paste0("= ", round(p_weighted_ISO, 3))))) +
  xlab("Total S+A gene counts per functional group") +
  ylab("Yield (CUE)")
Figure_S9_ISO

# --- Combined ---
Figure_S9 = ggarrange(Figure_S9_MAG, Figure_S9_ISO,
                      labels = c("A","B"), ncol = 2, nrow = 1) +
  theme(panel.background = element_blank())
Figure_S9

pdf("Output_Data/Figures/Figure_S9.pdf", width=10, height=5)
print(Figure_S9)
dev.off()

# Fig_S6 ----

total_genes.guild     = read.csv("Intermediate_Results/global_datasets_guild.csv",dec=".")
total_genes.guild     = total_genes.guild %>% select(id,guild)
total_genes.guild.940_ISO = read.csv("Intermediate_Results/total_genes.guild.940_ISO.csv",dec=".")
total_genes.guild.940_ISO = total_genes.guild.940_ISO %>% select(id,guild)

# MAGs
meta_img.nr  = read.csv("Input_Data/IMG_JGI_MAGs/IMG_bindata_withmeta_norestricted.csv",dec=".")
meta_img.nr  = meta_img.nr %>% select(Bin.ID, Bin.Completeness, Bin.Contamination)
colnames(meta_img.nr) = c("id","completeness","contamination")

loma_stat    = read.csv("Input_Data/LOMA_MAGs/mag_stats.csv",dec=".") 
loma_stat    = loma_stat %>% select(id,completeness,contamination)

fire_meta    = read_excel("Input_Data/FIRE_MAGs/MAG_Dataset_BurnSeverity_ARNelson.xlsx")
fire_meta    = fire_meta %>% select("MAG Id #",Completeness,Contamination)
colnames(fire_meta) = c("id","completeness","contamination")

MAGs_combined    = as.data.frame(rbind(meta_img.nr,loma_stat,fire_meta))
MAGs_combined    = merge(MAGs_combined,total_genes.guild,by = "id")

# Adding Loma guilds
total.genes.loma = total_genes.guild[636:1168,]
loma_stat1       = cbind(loma_stat,total.genes.loma)
loma_stat1       = loma_stat1[,c(1:3,5)]
MAGs_combined    = as.data.frame(rbind(MAGs_combined,loma_stat1))
high_quality     = MAGs_combined %>% filter(completeness > 90 & contamination < 5)
(nrow(high_quality)+100)*100/(27214)# 22.85221

write.csv(MAGs_combined, file = "Output_Data/MAGs_FG_table.csv")

MAGs_combined.FG = MAGs_combined %>% group_by(guild) %>% summarise(completeness = mean(completeness),
                                                                   contamination = mean(contamination))
MAGs_combined    = MAGs_combined %>% mutate(cat = rep("MAG", nrow(MAGs_combined)))
MAGs_combined    = MAGs_combined %>% select(id,completeness,contamination,cat)
MAGs_combined.FG = MAGs_combined.FG %>% mutate(cat = rep("Funtional Group", nrow(MAGs_combined.FG)))
colnames(MAGs_combined.FG) = c("id","completeness","contamination","cat")
total_data       = as.data.frame(rbind(MAGs_combined,MAGs_combined.FG))

Figure_S6.A      = ggplot(data = total_data, aes(x=cat, y=completeness)) + 
  geom_boxplot() + theme_light() + theme(legend.position = "none") + 
  labs(x = "", y = "Completeness (%)") + ylim(0,100)
Figure_S6.A  

Figure_S6.B      = ggplot(data = total_data, aes(x=cat, y=contamination)) + 
  geom_boxplot() + theme_light() + theme(legend.position = "none") + 
  labs(x = "", y = "Contamination (%)")
Figure_S6.B  

# Final figure

Figure_S6 = ggarrange(Figure_S6.A,Figure_S6.B,
                      labels = c("A","B"),
                      ncol = 2, nrow = 1) + theme(panel.background = element_blank())
Figure_S6

pdf("Output_Data/Figures/Figure_S6.pdf",
    width=4.81*1.5,height=2.21*1.5)
print(Figure_S6)
dev.off()

# Fig_S11 -----

# Loma Analysis

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

# Fire Project Analysis

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

# General Plot

Figure_S11 = ggarrange(Figure_Total_S_SA_LOMA, Figure_Total_S_SA_FIRE, 
                       labels = c("A","B"), ncol = 2, nrow = 1) + 
  theme(panel.background = element_blank())
Figure_S11

pdf("Output_Data/Figures/Figure_S11.pdf", width=12, height=12*2/5)
print(Figure_S11)
dev.off()

# Table S1 ----

# Call data

MAG_gen_trait_total = read.csv("Intermediate_Results/MAG_gen_trait.csv",dec=".")
MAG_gen_test = MAG_gen_trait_total
TOTAL_MAGs_fixed = read.csv("Intermediate_Results/TOTAL_MAGs_fixed.csv",dec=".")

# Loop Function

weighted_r2 = function(actual, predicted, weights) {
  # Defensive: drop any row where actual/predicted/weights is NA or non-finite,
  # rather than letting a single missing guild (e.g. mgr's 1 NA, ogt's 9 NA)
  # silently propagate NA through sum() and poison the whole R2.
  ok = stats::complete.cases(actual, predicted, weights) &
    is.finite(actual) & is.finite(predicted) & is.finite(weights)
  n_dropped <- sum(!ok)
  if (n_dropped > 0) {
    message(sprintf("weighted_r2: dropping %d row(s) with NA/non-finite actual, predicted, or weight", n_dropped))
  }
  actual    = actual[ok]
  predicted = predicted[ok]
  weights   = weights[ok]
  
  wmean    = weighted.mean(actual, weights)
  ss_res   = sum(weights * (actual - predicted)^2)
  ss_tot   = sum(weights * (actual - wmean)^2)
  1 - ss_res / ss_tot
}

# 1. Build n_GC — the per-guild count of MAGs

guild_n_GC = TOTAL_MAGs_fixed %>%
  filter(!is.na(GC_count)) %>%
  count(guild, name = "n_GC")

MAG_gen_test_complete = MAG_gen_test %>%
  filter(!is.na(GC_count_mean)) %>%
  left_join(guild_n_GC, by = "guild")

# 2. Four-model fitter for one trait (MAGs)

run_four_models = function(trait, data, weight_mags = "n_mags", weight_gc = "n_GC",
                           transform = c("log1p", "log", "logit")) {
  transform = match.arg(transform)
  df = data
  df$.trait_raw  = df[[trait]]
  
  if (transform == "log1p") {
    df$.log_trait = log(df[[trait]] + 1)   # gene-count traits: handles zero counts
  } else if (transform == "log") {
    n_nonpos = sum(df[[trait]] <= 0, na.rm = TRUE)
    if (n_nonpos > 0) {
      warning(sprintf(
        "%s: %d rows have trait <= 0 -> log() will produce -Inf/NaN. Filter these before fitting.",
        trait, n_nonpos))
    }
    df$.log_trait = log(df[[trait]])       # continuous positive variable (mgr, ogt): no pseudocount
  } else if (transform == "logit") {
    n_bad = sum(df[[trait]] <= 0 | df[[trait]] >= 1, na.rm = TRUE)
    if (n_bad > 0) {
      warning(sprintf(
        "%s: %d rows have trait <= 0 or >= 1 -> qlogis() will produce -Inf/Inf. Filter these before fitting.",
        trait, n_bad))
    }
    df$.log_trait = qlogis(df[[trait]])    # logit transform for bounded proportions (CUE)
  }
  
  df$.log_genome = log(df$genome.size)
  df$.w_mags     = df[[weight_mags]]
  df$.w_gc       = df[[weight_gc]]
  
  # Model 1: linear, raw scale
  m1 = lm(.trait_raw ~ genome.size, weights = .w_mags, data = df)
  
  # Model 2: power law / logit-linear, R2 back-transformed to raw scale for
  # direct comparison with Model 1
  m2 = lm(.log_trait ~ .log_genome, weights = .w_mags, data = df)
  # newdata=df forces evaluation on every row of df, aligned positionally,
  # regardless of the session's na.action setting (na.omit would otherwise
  # silently return a SHORTER vector, misaligning row-by-row with df$.trait_raw)
  pred_link_m2 = predict(m2, newdata = df)
  pred_raw_m2  = switch(transform,
                        log1p = exp(pred_link_m2) - 1,
                        log   = exp(pred_link_m2),
                        logit = plogis(pred_link_m2))
  m2_r2_bt    = weighted_r2(df$.trait_raw, pred_raw_m2, df$.w_mags)
  
  # Model 3: confound check (completeness/contamination), same scale as M2/M4
  m3 = lm(.log_trait ~ .log_genome + completeness_mean + contamination_mean,
          weights = .w_mags, data = df)
  
  # Model 4: GC alternative — uses its OWN weight (n_GC), not n_mags
  m4 = lm(.log_trait ~ GC_count_mean, weights = .w_gc, data = df)
  
  list(trait = trait, n = nrow(df),
       model1 = m1, model2 = m2, model2_R2_backtransformed = m2_r2_bt,
       model3 = m3, model4 = m4)
}

# Extract 

extract_summary = function(model, r2_override = NULL) {
  s  = summary(model)
  co = coef(model)
  n  = length(residuals(model))
  p  = length(co) - 1L  # number of predictors, excluding intercept
  
  if (!is.null(r2_override)) {
    R2     = r2_override
    R2_adj = 1 - (1 - r2_override) * (n - 1) / (n - p - 1)
  } else {
    R2     = s$r.squared
    R2_adj = s$adj.r.squared
  }
  
  data.frame(
    formula   = deparse(formula(model)),
    n         = n,
    R2        = R2,
    R2_adj    = R2_adj,
    intercept = unname(co[1]),
    slope     = unname(co[2]),
    p_slope   = coef(s)[2, 4],
    AIC       = AIC(model)
  )
}

# 4. Generalize into a loop across all remaining MAG traits

mag_traits = c("pH", "temp", "biofilm", "osmolyte",
               "CAZy", "Protein", "tranport.total",
               "GH.trasnport", "amino.transport")

mag_results_list = lapply(mag_traits, function(tr) {
  fit = run_four_models(tr, MAG_gen_test_complete)
  bind_rows(
    extract_summary(fit$model1) |> mutate(trait = tr, model = "1_linear"),
    extract_summary(fit$model2, fit$model2_R2_backtransformed) |>
      mutate(trait = tr, model = "2_powerlaw_backtransformed"),
    extract_summary(fit$model3) |> mutate(trait = tr, model = "3_confound_check"),
    extract_summary(fit$model4) |> mutate(trait = tr, model = "4_GC_alternative")
  )
})
MAG_TableS1 = bind_rows(mag_results_list) |>
  select(trait, model, n, formula, intercept, slope, p_slope, R2, R2_adj, AIC)

# 5. mgr, ogt, CUE

mgr_fit = run_four_models("mgr", MAG_gen_test_complete, transform = "log")
ogt_fit = run_four_models("ogt", MAG_gen_test_complete, transform = "log")

mgr_ogt_table = bind_rows(
  bind_rows(
    extract_summary(mgr_fit$model1) |> mutate(model = "1_linear"),
    extract_summary(mgr_fit$model2, mgr_fit$model2_R2_backtransformed) |>
      mutate(model = "2_powerlaw_backtransformed"),
    extract_summary(mgr_fit$model3) |> mutate(model = "3_confound_check"),
    extract_summary(mgr_fit$model4) |> mutate(model = "4_GC_alternative")
  ) |> mutate(trait = "mgr"),
  bind_rows(
    extract_summary(ogt_fit$model1) |> mutate(model = "1_linear"),
    extract_summary(ogt_fit$model2, ogt_fit$model2_R2_backtransformed) |>
      mutate(model = "2_powerlaw_backtransformed"),
    extract_summary(ogt_fit$model3) |> mutate(model = "3_confound_check"),
    extract_summary(ogt_fit$model4) |> mutate(model = "4_GC_alternative")
  ) |> mutate(trait = "ogt")
)

# CUE — proportion, log-transformed 

MAG_gen_trait_total = read.csv("Intermediate_Results/MAG_gen_trait.csv", dec = ".")
MAG_gen_test = MAG_gen_trait_total
TOTAL_MAGs_fixed = read.csv("Intermediate_Results/TOTAL_MAGs_fixed.csv", dec = ".")
total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv", dec = ".")
sim_10000 = read.csv("Output_Data/my_MAGs_BGE_glucose_10000.csv")

guild_n_GC = TOTAL_MAGs_fixed %>%
  filter(!is.na(GC_count)) %>%
  count(guild, name = "n_GC")

# Used by all traits EXCEPT CUE (pH, temp, CAZy, Protein, transporters, mgr, ogt, etc.)
MAG_gen_test_complete = MAG_gen_test %>%
  filter(!is.na(GC_count_mean)) %>%
  left_join(guild_n_GC, by = "guild")

# CUE only

guild_lookup = total_genes.guild.940_MAG %>% select(id, guild)

cue_with_guild = sim_10000 %>%
  left_join(guild_lookup, by = c("MAG_id" = "id")) %>%
  filter(!is.na(guild))

guild_cue_without = cue_with_guild %>%
  filter(BGE <= 0.848) %>%
  group_by(guild) %>%
  summarise(CUE_mean = mean(BGE, na.rm = TRUE),
            genome.size = mean(genomesize, na.rm = TRUE),
            n_cue_mags = n())

MAG_gen_test_CUE_complete <- guild_cue_without %>%
  left_join(
    MAG_gen_test %>% select(guild, completeness_mean, contamination_mean, GC_count_mean),
    by = "guild"
  ) %>%
  left_join(guild_n_GC, by = "guild")

cue_fit = run_four_models("CUE_mean", MAG_gen_test_CUE_complete,
                          weight_mags = "n_cue_mags",
                          weight_gc   = "n_GC",
                          transform   = "logit")

cue_table = bind_rows(
  extract_summary(cue_fit$model1) |> mutate(model = "1_linear"),
  extract_summary(cue_fit$model2, cue_fit$model2_R2_backtransformed) |>
    mutate(model = "2_powerlaw_backtransformed"),
  extract_summary(cue_fit$model3) |> mutate(model = "3_confound_check"),
  extract_summary(cue_fit$model4) |> mutate(model = "4_GC_alternative")
) |> mutate(trait = "CUE_mean")

# Merge tables 

Total_Models = as.data.frame(rbind(MAG_TableS1,mgr_ogt_table,cue_table))

write.csv(Total_Models, file = "Output_Data/MAG_model_comparison.csv")

# Table S2 ----

ISO_gen_TOTAL_ab    = read.csv("Intermediate_Results/Isolates_gen_trait.csv",dec=".")

# Traits

run_two_models_isolates = function(trait, data, weight_col = NULL,
                                   transform = c("log1p", "log", "logit")) {
  transform = match.arg(transform)
  df = data
  df$.trait_raw  = df[[trait]]
  
  if (transform == "log1p") {
    df$.log_trait = log(df[[trait]] + 1)
  } else if (transform == "log") {
    n_nonpos <- sum(df[[trait]] <= 0, na.rm = TRUE)
    if (n_nonpos > 0) {
      warning(sprintf("%s: %d rows have trait <= 0 -> log() will produce -Inf/NaN.", trait, n_nonpos))
    }
    df$.log_trait = log(df[[trait]])
  } else if (transform == "logit") {
    n_bad <- sum(df[[trait]] <= 0 | df[[trait]] >= 1, na.rm = TRUE)
    if (n_bad > 0) {
      warning(sprintf("%s: %d rows have trait <= 0 or >= 1 -> qlogis() will produce -Inf/Inf.", trait, n_bad))
    }
    df$.log_trait = qlogis(df[[trait]])
  }
  
  df$.log_genome = log(df$genome.size)
  
  if (!is.null(weight_col) && weight_col %in% names(df)) {
    df$.w = df[[weight_col]]
    m1 = lm(.trait_raw ~ genome.size, weights = .w, data = df)
    m2 = lm(.log_trait ~ .log_genome, weights = .w, data = df)
  } else {
    df$.w = 1
    m1 = lm(.trait_raw ~ genome.size, data = df)
    m2 = lm(.log_trait ~ .log_genome, data = df)
  }
  
  pred_link_m2 = predict(m2, newdata = df)
  pred_raw_m2  = switch(transform,
                        log1p = exp(pred_link_m2) - 1,
                        log   = exp(pred_link_m2),
                        logit = plogis(pred_link_m2))
  m2_r2_bt = weighted_r2(df$.trait_raw, pred_raw_m2, df$.w)
  
  list(trait = trait, n = nrow(df), model1 = m1, model2 = m2,
       model2_R2_backtransformed = m2_r2_bt)
}

iso_traits = c("pH", "temp", "biofilm", "osmolyte",
               "CAZy", "Protein", "tranport.total",
               "GH.trasnport", "amino.transport")  

iso_results_list = lapply(iso_traits, function(tr) {
  fit = run_two_models_isolates(tr, ISO_gen_TOTAL_ab, "n_mags")
  bind_rows(
    extract_summary(fit$model1) |> mutate(trait = tr, model = "1_linear"),
    extract_summary(fit$model2, fit$model2_R2_backtransformed) |>
      mutate(trait = tr, model = "2_powerlaw_backtransformed")
  )
})

ISO_TableS2 = bind_rows(iso_results_list) |>
  select(trait, model, n, formula, intercept, slope, p_slope, R2, R2_adj, AIC)

# MGR, OGT and CUE

mgr_iso_fit = run_two_models_isolates("mgr", ISO_gen_TOTAL_ab, "n_mags", transform = "log")
ogt_iso_fit = run_two_models_isolates("ogt", ISO_gen_TOTAL_ab, "n_mags", transform = "log")

mgr_ogt_iso_table = bind_rows(
  bind_rows(
    extract_summary(mgr_iso_fit$model1) |> mutate(model = "1_linear"),
    extract_summary(mgr_iso_fit$model2, mgr_iso_fit$model2_R2_backtransformed) |>
      mutate(model = "2_powerlaw_backtransformed")
  ) |> mutate(trait = "mgr"),
  bind_rows(
    extract_summary(ogt_iso_fit$model1) |> mutate(model = "1_linear"),
    extract_summary(ogt_iso_fit$model2, ogt_iso_fit$model2_R2_backtransformed) |>
      mutate(model = "2_powerlaw_backtransformed")
  ) |> mutate(trait = "ogt")
)

# CUE

range(ISO_gen_TOTAL_ab$CUE_mean, na.rm = TRUE)
sum(ISO_gen_TOTAL_ab$CUE_mean > 0.848, na.rm = TRUE)

cue_iso_fit = run_two_models_isolates("CUE_mean", ISO_gen_TOTAL_ab, "CUE_n_isolates", transform = "logit")

cue_iso_table = bind_rows(
  extract_summary(cue_iso_fit$model1) |> mutate(model = "1_linear"),
  extract_summary(cue_iso_fit$model2, cue_iso_fit$model2_R2_backtransformed) |>
    mutate(model = "2_powerlaw_backtransformed")
) |> mutate(trait = "CUE_mean")

# Merge tables 

Total_Models_ISO = as.data.frame(rbind(ISO_TableS2,mgr_ogt_iso_table,cue_iso_table))

write.csv(Total_Models_ISO, file = "Output_Data/ISOLATES_model_comparison.csv")