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

# Figure S1 ----

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

# Adding genome size to the trait dataset - MAG

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

# Minimum generation time <= 5
Y_MAG.B   = Y_MAG %>% select(c("guild","genome.size","mgt"))
Y_MAG.B   = Y_MAG.B %>% filter(mgt <= 5)
Y_MAG.B.1 = Y_MAG.B %>% group_by(guild) %>% summarise(across(where(is.numeric),
                                                             list(mean=mean), na.rm=TRUE))

# Optimal Growth Temperature
Y_MAG.C     = Y_MAG %>% select(c("guild","genome.size","ogt"))
Y_MAG.C     = na.omit(Y_MAG.C)
Y_MAG.C.1   = Y_MAG.C %>% group_by(guild) %>% summarise(across(where(is.numeric),
                                                               list(mean=mean), na.rm=TRUE)) 
colnames(Y_MAG.B.1) = c("guild","genome.size","mgr")
colnames(Y_MAG.A.1) = c("guild","genome.size","mgr")
colnames(Y_MAG.C.1) = c("guild","genome.size","ogt")
MAG_gen_trait.1     = merge(MAG_gen_trait,Y_MAG.C.1,by = c("guild","genome.size"))
MAG_gen_trait.2     = merge(MAG_gen_trait.1,Y_MAG.B.1,by = c("guild","genome.size"))
MAG_gen_trait.3     = merge(MAG_gen_trait.1,Y_MAG.A.1,by = c("guild","genome.size"))

write.csv(MAG_gen_trait.2, file = "Intermediate_Results/MAG_gen_trait.csv")
write.csv(MAG_gen_trait.3, file = "Intermediate_Results/MAG_gen_trait.5.csv")

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

# Adding genome size to the trait dataset - Isolates

Isolates_gen_trait    = as.data.frame(cbind(Aminoacids_T_Isolates,PH_Isolates$PH_total,
                                            TEMP_Isolates$TEMP_total,Biofilm_Isolates$BIO_total,
                                            Osmolyte_Isolates$OSMO_total,GH_T_Isolates$CAR_TRANSP_TOTAL,
                                            Transporter_Isolates$transp_total,Protein_Isolates$PR_total,
                                            GH_Isolates$GH_total))
colnames(Isolates_gen_trait) = c("guild","genome.size","amino.transport","pH","temp",
                                 "biofilm","osmolyte","GH.trasnport","tranport.total",
                                 "Protein","CAZy")

# CUE ----

isolates     = read.csv("Input_Data/SOIL_ISOLATES/dement_isolates_CUE.csv",dec=".")
isolatest    = read.csv("Input_Data/SOIL_ISOLATES/total.granularity.3_datasets.csv",dec=".")
isolates.1   = subset(isolates, isolates$CUE !='NaN')

isolates.2   = isolates.1 %>% left_join(total_genes.guild.940_ISO, by='id')
isolates.2.m = isolates.2 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
isolates.2.m = isolates.2.m %>% select(c("guild","CUE_mean","genome.size_mean",
                                         "mingentime_mean","optimumT_mean"))
colnames(isolates.2.m) = c("guild","yield","genome.size","MGT","OGT")
ISO_gen_trait.1        = merge(Isolates_gen_trait,isolates.2.m,by = c("guild","genome.size"))
colnames(ISO_gen_trait.1) = c("Guild","Genome Size","Amino-transporter","pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","Yield","MGT","OGT")
col_order = c("Guild", "Genome Size","Yield","OGT","MGT","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte")
ISO_gen_trait.1   = ISO_gen_trait.1[, col_order]
ISO_gen_trait.1   = ISO_gen_trait.1 %>% filter(Yield <= 0.9)

write.csv(ISO_gen_trait.1, file = "Intermediate_Results/ISO_gen_trait.1.csv")

# Aggregated Traits ----

# MAGs no normalized ----

MAG_gen_trait.fast   = read.csv("Intermediate_Results/MAG_gen_trait.csv",dec=".")
MAG_gen_trait.slow   = read.csv("Intermediate_Results/MAG_gen_trait.5.csv",dec=".")

# A-S-traits ----

MAG_gen_trait.fast   =  MAG_gen_trait.fast %>% mutate(S_traits  = rowSums(MAG_gen_trait.fast[,5:8]),
                                                      A_traits  = rowSums(MAG_gen_trait.fast[,10:12]),
                                                      A_enzymes = rowSums(MAG_gen_trait.fast[,11:12]),
                                                      speed     = rep("fast",nrow(MAG_gen_trait.fast)))

MAG_gen_trait.slow   =  MAG_gen_trait.slow %>% mutate(S_traits  = rowSums(MAG_gen_trait.slow[,5:8]),
                                                      A_traits  = rowSums(MAG_gen_trait.slow[,10:12]),
                                                      A_enzymes = rowSums(MAG_gen_trait.slow[,11:12]),
                                                      speed     = rep("low",nrow(MAG_gen_trait.slow)))


MAG_gen_trait_total  = as.data.frame(rbind(MAG_gen_trait.fast,MAG_gen_trait.slow))
MAG_gen_trait_total  = MAG_gen_trait_total %>% mutate(A_S = A_traits/S_traits)

# Isolates ----

ISO_gen_TOTAL_ab = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_TOTAL_ab = ISO_gen_TOTAL_ab %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                                                 MGT  > 5 ~ "slow"))
# A-S-traits ----

ISO_gen_TOTAL_ab   =  ISO_gen_TOTAL_ab %>% mutate(S_traits      = rowSums(ISO_gen_TOTAL_ab[,12:15]),
                                                  A_traits  = rowSums(ISO_gen_TOTAL_ab[,9:11]),
                                                  A_enzymes = rowSums(ISO_gen_TOTAL_ab[,10:11]))
ISO_gen_TOTAL_ab   = ISO_gen_TOTAL_ab %>% mutate(A_S = A_traits/S_traits)

# Figure 1 ----

Figure_Total_S_SA_MAGs   = ggplot(data = MAG_gen_trait_total, aes(x = as.numeric((A_traits)), 
                                                                  y = as.numeric((S_traits)),color = A_S)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_cor(method = "pearson", label.x = 3, label.y = 80) + 
  xlab("A Traits") + 
  ylab("S Traits") +
  theme() + scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5))  + 
  ylim(0,80)  + xlim(0,250) 
Figure_Total_S_SA_MAGs

Figure_Total_S_SA_ISO   = ggplot(data = ISO_gen_TOTAL_ab, aes(x = as.numeric((A_traits)), 
                                                              y = as.numeric((S_traits)),color = A_S)) +
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_cor(method = "pearson", label.x = 3, label.y = 80) + 
  xlab("A Traits") + 
  ylab("S Traits") +
  theme() + scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5))  + 
  ylim(0,80)  + xlim(0,250) 
Figure_Total_S_SA_ISO

Figure_1 = ggarrange(Figure_Total_S_SA_MAGs, Figure_Total_S_SA_ISO, 
                     labels = c("A","B"),
                     ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Figure_1
pdf("Output_Data/Figures/Figure_1.pdf",
    width=12,height=12*2/5)
print(Figure_1)
dev.off()

# MAGs normalized ----

MAG_gen_trait_total_nor = MAG_gen_trait_total %>% mutate (S_traits_nor  = S_traits/genome.size,
                                                          A_traits_nor  = A_traits/genome.size,
                                                          A_enzymes_nor = A_enzymes/genome.size)

# Isolates normalized ----

ISO_gen_TOTAL_ab_nor       = ISO_gen_TOTAL_ab %>% mutate (S_traits_nor  = S_traits/Genome.Size,
                                                          A_traits_nor  = A_traits/Genome.Size,
                                                          A_enzymes_nor = A_enzymes/Genome.Size)

# Figure S2 ----

Figure_Total_S_SA_MAG_nor   = ggplot(data = MAG_gen_trait_total_nor, 
                                     aes(x = as.numeric((A_traits_nor)),
                                         y = as.numeric((S_traits_nor)),
                                         color = A_S)) + 
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_cor(method = "pearson") + 
  xlab("A Traits") + 
  ylab("S Traits") +
  theme() + scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5))# + ylim(3.5e-6,1.5e-5) + xlim(8e-6,4e-5) 
Figure_Total_S_SA_MAG_nor

Figure_Total_S_SA_ISO_nor   = ggplot(data = ISO_gen_TOTAL_ab_nor, 
                                     aes(x = as.numeric((A_traits_nor)),
                                         y = as.numeric((S_traits_nor)),
                                         color = A_S)) + 
  geom_point(size = 2.5) + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_cor(method = "pearson") + 
  xlab("A Traits") + 
  ylab("S Traits") +
  theme() + scale_color_gradient(low="red", high="blue", limits = c(0.5,8.5))# + ylim(3.5e-6,1.5e-5) + xlim(8e-6,4e-5) 
Figure_Total_S_SA_ISO_nor

Figure_S2 = ggarrange(Figure_Total_S_SA_MAG_nor, Figure_Total_S_SA_ISO_nor, 
                      labels = c("A","B"),
                      ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Figure_S2
pdf("Output_Data/Figures/Figure_S2.pdf",
    width=12,height=12*2/5)
print(Figure_S2)
dev.off()

# Adding GC content to the dataset of traits ----

# MAGs----

# All MAGs
meta_img.nr   = read.csv("Input_Data/IMG_JGI_MAGs/IMG_bindata_withmeta_norestricted.csv",dec=".")
IMG_MAGs      = meta_img.nr %>% select(Total.Number.of.Bases,Bin.Completeness,Bin.Contamination,GC.....assembled)
colnames(IMG_MAGs) = c("genome.size", "completeness",
                       "contamination","GC_count")
hmm_img     = read.csv("Intermediate_Results/IMG_global_dataset.csv",dec=".")
hmm_img     = as.data.frame(hmm_img) %>% select(3,26,70:1790)
colnames(hmm_img)[1]  = "id"
colnames(hmm_img)[2]  = "genome.size"
IMG_MAGs    = left_join(hmm_img, IMG_MAGs, by=c('genome.size'))
IMG_MAGs    = IMG_MAGs %>% select(genome.size,completeness,contamination,GC_count)

# Loma
mag_stat    = read.csv("Input_Data/LOMA_MAGs/mag_stats.csv",dec=".") 
mag_abun    = read.csv("Input_Data/LOMA_MAGs/mag_adundance.csv",dec=".") 
mag_stat    = mag_stat %>% full_join(mag_abun)
LOMA_MAGs   = mag_stat %>% select(size,completeness,contamination)
LOMA_MAGs   = LOMA_MAGs %>% mutate(GC_count=NA)
colnames(LOMA_MAGs) = c("genome.size", "completeness",
                        "contamination","GC_count")
# Fire
fire_meta   = read_excel("Input_Data/FIRE_MAGs/MAG_Dataset_BurnSeverity_ARNelson.xlsx")
colnames(fire_meta)[1] = "id"
fire_stat   = read.csv("Input_Data/FIRE_MAGs/fire_metadata.csv",dec=".")
FIRE_MAGs   = merge(x = fire_meta, y = fire_stat, by = "id")
FIRE_MAGs   = FIRE_MAGs %>% select(Genome_size,Completeness,Contamination,GC)
colnames(FIRE_MAGs) = c("genome.size", "completeness",
                        "contamination","GC_count")
# Merge datasets
TOTAL_MAGs = as.data.frame(rbind(IMG_MAGs,LOMA_MAGs,FIRE_MAGs))

# Guild codes
total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")
# Extracting completeness, GC counts and contamination
TOTAL_MAGs_1 = left_join(total_genes.guild.940_MAG, TOTAL_MAGs, by=c('genome.size'))
# Remove duplicaes
TOTAL_MAGs_1 = TOTAL_MAGs_1[!duplicated(TOTAL_MAGs_1$id), ]
TOTAL_MAGs_1 = TOTAL_MAGs_1 %>% select(id,guild,genome.size,completeness,contamination,GC_count)
# Summarize by guild
TOTAL_MAGS_2 = TOTAL_MAGs_1 %>% group_by(guild) %>% summarise(across(where(is.numeric), list(mean=mean), na.rm=TRUE))
# Merge with trait data (normalized)
# Fast growing
MAG_gen_trait.2  = read.csv("Intermediate_Results/MAG_gen_trait_norm.csv",dec=".")
MAG_gen_trait.2  = MAG_gen_trait.2[, -1]
MAG_gen_trait.fast  = merge(x = MAG_gen_trait.2, y = TOTAL_MAGS_2, by = "guild")
write.csv(MAG_gen_trait.fast, file = "Intermediate_Results/MAG_gen_trait.fast_GC_norm.csv")

# Slow growing
MAG_gen_trait.3  = read.csv("Intermediate_Results/MAG_gen_trait.5_norm.csv",dec=".")
MAG_gen_trait.3  = MAG_gen_trait.3[, -1]
MAG_gen_trait.slow  = merge(x = MAG_gen_trait.3, y = TOTAL_MAGS_2, by = "guild")
write.csv(MAG_gen_trait.slow, file = "Intermediate_Results/MAG_gen_trait.slow_GC_norm.csv")

# Total
MAG_gen_TOTAL = as.data.frame(rbind(MAG_gen_trait.fast,MAG_gen_trait.slow))
write.csv(MAG_gen_TOTAL, file = "Intermediate_Results/MAG_gen_TOTAL_GC_norm.csv")

# Merge with trait data (gene counts)
# Fast growing
MAG_gen_trait.2  = read.csv("Intermediate_Results/MAG_gen_trait.csv",dec=".")
MAG_gen_trait.2  = MAG_gen_trait.2[, -1]
MAG_gen_trait.fast  = merge(x = MAG_gen_trait.2, y = TOTAL_MAGS_2, by = "guild")
write.csv(MAG_gen_trait.fast, file = "Intermediate_Results/MAG_gen_trait.fast_GC.csv")

# Slow growing
MAG_gen_trait.3  = read.csv("Intermediate_Results/MAG_gen_trait.5.csv",dec=".")
MAG_gen_trait.3  = MAG_gen_trait.3[, -1]
MAG_gen_trait.slow  = merge(x = MAG_gen_trait.3, y = TOTAL_MAGS_2, by = "guild")
write.csv(MAG_gen_trait.slow, file = "Intermediate_Results/MAG_gen_trait.slow_GC.csv")

# Total
MAG_gen_TOTAL = as.data.frame(rbind(MAG_gen_trait.fast,MAG_gen_trait.slow))
write.csv(MAG_gen_TOTAL, file = "Intermediate_Results/MAG_gen_TOTAL_GC.csv")

# Correlation tables ----

# MAGs gene counts ----
MAG_gen_TOTAL_ab = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")
MAG_gen_TOTAL_ab = MAG_gen_TOTAL_ab[, -c(1,15)]

# Convert MGT to MGR ----
MAG_gen_TOTAL_ab = MAG_gen_TOTAL_ab %>% mutate(mgr = 1/mgr)

# Colnames ----
colnames(MAG_gen_TOTAL_ab) = c("Guild", "Genome Size", "Amino-transporter", "pH-Tol",
                               "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                               "Total-transporter","Protein-enzyme","CAZy","OGT",
                               "MGR","Completeness","Contamination","GC-count")

# Reordering columns ----

col_order = c("Guild", "Genome Size","OGT","MGR","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte", "Completeness","Contamination","GC-count")

MAG_gen_TOTAL_ab = MAG_gen_TOTAL_ab[, col_order]

# Correlations

cor_1 = as.data.frame(cor(MAG_gen_TOTAL_ab[,3:16]))
write.csv(cor_1, file = "Output_Data/MAG.correlation_complete.csv")

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

# Figure S3 - Cross correlation MAGs ----

Figure_S3 = ggpairs(MAG_gen_TOTAL_ab[,3:16], 
                    upper = list(continuous = my_fn),
                    lower = list(continuous = "smooth"))  
Figure_S3

pdf("Output_Data/Figures/Figure_S3.pdf",
    width=4*4,height=4*4)
print(Figure_S3)
dev.off()

# Isolates ----

ISO_gen_trait.1  = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_trait.1  = ISO_gen_trait.1[, -1]

# Convert MGT to MGR ----
ISO_gen_trait.1 = ISO_gen_trait.1 %>% mutate(MGT = 1/MGT)
colnames(ISO_gen_trait.1)[5] = "MGR"

# Correlation tables----

cor_1 = as.data.frame(cor(ISO_gen_trait.1[,3:14]))
write.csv(cor_1, file = "Output_Data/ISOLATE.correlation_complete.csv")

# Figure S4 - Cross correlation ISOLATES ----

Figure_S4 = ggpairs(ISO_gen_trait.1[,3:14], 
                    upper = list(continuous = my_fn),
                    lower = list(continuous = "smooth"))  
Figure_S4

pdf("Output_Data/Figures/Figure_S4.pdf",
    width=4*4,height=4*4)
print(Figure_S4)
dev.off()

# MAGs normalized 2 ----

MAG_gen_trait_total_nor_2 =
  MAG_gen_trait_total %>%
  select(guild, genome.size, amino.transport,
         pH, temp, biofilm, osmolyte,
         GH.trasnport, tranport.total,
         Protein, CAZy, ogt, mgr) %>%
  mutate(across(
    c(amino.transport, pH, temp, biofilm, osmolyte,
      GH.trasnport, tranport.total, Protein, CAZy),
    ~ .x / MAG_gen_trait_total$genome.size
  ))

# Add contamination, GC content etc ----

MAG_gen_TOTAL_ab = MAG_gen_TOTAL_ab %>% select(guild,completeness_mean,
                                               contamination_mean,GC_count_mean)
MAG_gen_trait_total_nor_2 = merge(MAG_gen_trait_total_nor_2,MAG_gen_TOTAL_ab,by="guild")

# Colnames ----
colnames(MAG_gen_trait_total_nor_2) = c("Guild", "Genome Size", "Amino-transporter", "pH-Tol",
                                        "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                                        "Total-transporter","Protein-enzyme","CAZy","OGT","MGR",
                                        "Completeness","Contamination","GC-count")

# Reordering columns ----

col_order = c("Guild", "Genome Size","OGT","MGR","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte","Completeness","Contamination","GC-count")
MAG_gen_trait_total_nor_2 = MAG_gen_trait_total_nor_2[, col_order]
# Convert MGT to MGR ----
MAG_gen_trait_total_nor_2 = MAG_gen_trait_total_nor_2 %>% mutate(MGR = 1/MGR)

# Figure S3A - Cross correlation MAGs normalized ----

Figure_S3A = ggpairs(MAG_gen_trait_total_nor_2[,3:16], 
                     upper = list(continuous = my_fn),
                     lower = list(continuous = "smooth"))  
Figure_S3A

pdf("Output_Data/Figures/Figure_S3A.pdf",
    width=4*4,height=4*4)
print(Figure_S3A)
dev.off()

# Isolates normalized ----

ISO_gen_TOTAL_ab_nor_2    = ISO_gen_TOTAL_ab %>%
  select(Guild, Genome.Size, Amino.transporter,
         pH.Tol, Temp.Tol, Biofilm, Osmolyte,
         GH.transporter, Total.transporter,
         Protein.enzyme, CAZy, OGT, MGT,Yield) %>%
  mutate(across(
    c(Amino.transporter,
      pH.Tol, Temp.Tol, Biofilm, Osmolyte,
      GH.transporter, Total.transporter,
      Protein.enzyme, CAZy),
    ~ .x / ISO_gen_TOTAL_ab$Genome.Size
  ))

# Colnames ----
colnames(ISO_gen_TOTAL_ab_nor_2) = c("Guild", "Genome Size", "Amino-transporter", "pH-Tol",
                                     "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                                     "Total-transporter","Protein-enzyme","CAZy","OGT","MGR",
                                     "Yield")

# Reordering columns ----

col_order = c("Guild", "Genome Size","Yield","OGT","MGR","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte")
ISO_gen_TOTAL_ab_nor_2 = ISO_gen_TOTAL_ab_nor_2[, col_order]
# Convert MGT to MGR ----
ISO_gen_TOTAL_ab_nor_2 = ISO_gen_TOTAL_ab_nor_2 %>% mutate(MGR = 1/MGR)

# Figure S4A - Cross correlation ISOLATES ----

Figure_S4A = ggpairs(ISO_gen_TOTAL_ab_nor_2[,3:14], 
                     upper = list(continuous = my_fn),
                     lower = list(continuous = "smooth"))  
Figure_S4A

pdf("Output_Data/Figures/Figure_S4A.pdf",
    width=4*4,height=4*4)
print(Figure_S4A)
dev.off()

# Empirical models between Genome Size and functional groups for MAGs (Table S1) ----

# S vs A tradeoffs (MAGs) ----

MAG_gen_TOTAL_ab = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")
MAG_gen_test     = MAG_gen_TOTAL_ab %>% mutate(A_trait = rowSums(MAG_gen_TOTAL_ab[,10:12]),
                                               S_trait = rowSums(MAG_gen_TOTAL_ab[,5:8]))
MAG_gen_test     = MAG_gen_test %>% mutate(A_S = A_trait/S_trait)

# CAZy enzyme ----

# Linear Model 
cazy_model_A       = lm(CAZy ~ genome.size, data = MAG_gen_test)
summary(cazy_model_A)

# Genome Size and Completeness
cazy_model_B       = lm(CAZy ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(cazy_model_B)

domin(CAZy ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
cazy_model_C       = lm(CAZy ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(cazy_model_C)

domin(CAZy ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
cazy_model_D       = lm(CAZy ~ GC_count_mean, data = MAG_gen_test)
summary(cazy_model_D)

# Power Law
CAZy_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                         log(MAG_gen_test$CAZy)))

CAZy_data           = subset(CAZy_data, CAZy_data$V2 != "-Inf") 
colnames(CAZy_data) = c("genome.size","CAZy")

cazy_model_E       = lm(CAZy ~ genome.size, data = CAZy_data)
summary(cazy_model_E)

AIC(cazy_model_A, cazy_model_B, cazy_model_C, cazy_model_D, cazy_model_E)

# Protein enzyme ----

# Linear Model 
Protein_model_A       = lm(Protein ~ genome.size, data = MAG_gen_test)
summary(Protein_model_A)

# Genome Size and Completeness
Protein_model_B       = lm(Protein ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(Protein_model_B)

domin(Protein ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
Protein_model_C       = lm(Protein ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(Protein_model_C)

domin(Protein ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
Protein_model_D       = lm(Protein ~ GC_count_mean, data = MAG_gen_test)
summary(Protein_model_D)

# Power Law
Protein_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                            log(MAG_gen_test$Protein)))

Protein_data           = subset(Protein_data, Protein_data$V2 != "-Inf") 
colnames(Protein_data) = c("genome.size","Protein")

Protein_model_E       = lm(Protein ~ genome.size, data = Protein_data)
summary(Protein_model_E)

AIC(Protein_model_A, Protein_model_B, Protein_model_C, Protein_model_D, Protein_model_E)

# tranport.total enzyme ----

# Linear Model 
tranport.total_model_A       = lm(tranport.total ~ genome.size, data = MAG_gen_test)
summary(tranport.total_model_A)

# Genome Size and Completeness
tranport.total_model_B       = lm(tranport.total ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(tranport.total_model_B)

domin(tranport.total ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
tranport.total_model_C       = lm(tranport.total ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(tranport.total_model_C)

domin(tranport.total ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
tranport.total_model_D       = lm(tranport.total ~ GC_count_mean, data = MAG_gen_test)
summary(tranport.total_model_D)

# Power Law
tranport.total_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                                   log(MAG_gen_test$tranport.total)))

tranport.total_data           = subset(tranport.total_data, tranport.total_data$V2 != "-Inf") 
colnames(tranport.total_data) = c("genome.size","tranport.total")

tranport.total_model_E       = lm(tranport.total ~ genome.size, data = tranport.total_data)
summary(tranport.total_model_E)

AIC(tranport.total_model_A, tranport.total_model_B, tranport.total_model_C, tranport.total_model_D, tranport.total_model_E)

# GH.trasnport enzyme ----

# Linear Model 
GH.trasnport_model_A       = lm(GH.trasnport ~ genome.size, data = MAG_gen_test)
summary(GH.trasnport_model_A)

# Genome Size and Completeness
GH.trasnport_model_B       = lm(GH.trasnport ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(GH.trasnport_model_B)

domin(GH.trasnport ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
GH.trasnport_model_C       = lm(GH.trasnport ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(GH.trasnport_model_C)

domin(GH.trasnport ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
GH.trasnport_model_D       = lm(GH.trasnport ~ GC_count_mean, data = MAG_gen_test)
summary(GH.trasnport_model_D)

# Power Law
GH.trasnport_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                                 log(MAG_gen_test$GH.trasnport)))

GH.trasnport_data           = subset(GH.trasnport_data, GH.trasnport_data$V2 != "-Inf") 
colnames(GH.trasnport_data) = c("genome.size","GH.trasnport")

GH.trasnport_model_E       = lm(GH.trasnport ~ genome.size, data = GH.trasnport_data)
summary(GH.trasnport_model_E)

AIC(GH.trasnport_model_A, GH.trasnport_model_B, GH.trasnport_model_C, GH.trasnport_model_D, GH.trasnport_model_E)

# amino.transport enzyme ----

# Linear Model 
amino.transport_model_A       = lm(amino.transport ~ genome.size, data = MAG_gen_test)
summary(amino.transport_model_A)

# Genome Size and Completeness
amino.transport_model_B       = lm(amino.transport ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(amino.transport_model_B)

domin(amino.transport ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
amino.transport_model_C       = lm(amino.transport ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(amino.transport_model_C)

domin(amino.transport ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
amino.transport_model_D       = lm(amino.transport ~ GC_count_mean, data = MAG_gen_test)
summary(amino.transport_model_D)

# Power Law
amino.transport_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                                    log(MAG_gen_test$amino.transport)))

amino.transport_data           = subset(amino.transport_data, amino.transport_data$V2 != "-Inf") 
colnames(amino.transport_data) = c("genome.size","amino.transport")

amino.transport_model_E       = lm(amino.transport ~ genome.size, data = amino.transport_data)
summary(amino.transport_model_E)

AIC(amino.transport_model_A, amino.transport_model_B, amino.transport_model_C, amino.transport_model_D, amino.transport_model_E)

# osmolyte enzyme ----

# Linear Model 
osmolyte_model_A       = lm(osmolyte ~ genome.size, data = MAG_gen_test)
summary(osmolyte_model_A)

# Genome Size and Completeness
osmolyte_model_B       = lm(osmolyte ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(osmolyte_model_B)

domin(osmolyte ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
osmolyte_model_C       = lm(osmolyte ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(osmolyte_model_C)

domin(osmolyte ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
osmolyte_model_D       = lm(osmolyte ~ GC_count_mean, data = MAG_gen_test)
summary(osmolyte_model_D)

# Power Law
osmolyte_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                             log(MAG_gen_test$osmolyte)))

osmolyte_data           = subset(osmolyte_data, osmolyte_data$V2 != "-Inf") 
colnames(osmolyte_data) = c("genome.size","osmolyte")

osmolyte_model_E       = lm(osmolyte ~ genome.size, data = osmolyte_data)
summary(osmolyte_model_E)

AIC(osmolyte_model_A, osmolyte_model_B, osmolyte_model_C, osmolyte_model_D, osmolyte_model_E)

# biofilm enzyme ----

# Linear Model 
biofilm_model_A       = lm(biofilm ~ genome.size, data = MAG_gen_test)
summary(biofilm_model_A)

# Genome Size and Completeness
biofilm_model_B       = lm(biofilm ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(biofilm_model_B)

domin(biofilm ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
biofilm_model_C       = lm(biofilm ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(biofilm_model_C)

domin(biofilm ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
biofilm_model_D       = lm(biofilm ~ GC_count_mean, data = MAG_gen_test)
summary(biofilm_model_D)

# Power Law
biofilm_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                            log(MAG_gen_test$biofilm)))

biofilm_data           = subset(biofilm_data, biofilm_data$V2 != "-Inf") 
colnames(biofilm_data) = c("genome.size","biofilm")

biofilm_model_E       = lm(biofilm ~ genome.size, data = biofilm_data)
summary(biofilm_model_E)

AIC(biofilm_model_A, biofilm_model_B, biofilm_model_C, biofilm_model_D, biofilm_model_E)

# temp enzyme ----

# Linear Model 
temp_model_A       = lm(temp ~ genome.size, data = MAG_gen_test)
summary(temp_model_A)

# Genome Size and Completeness
temp_model_B       = lm(temp ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(temp_model_B)

domin(temp ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
temp_model_C       = lm(temp ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(temp_model_C)

domin(temp ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
temp_model_D       = lm(temp ~ GC_count_mean, data = MAG_gen_test)
summary(temp_model_D)

# Power Law
temp_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                         log(MAG_gen_test$temp)))

temp_data           = subset(temp_data, temp_data$V2 != "-Inf") 
colnames(temp_data) = c("genome.size","temp")

temp_model_E       = lm(temp ~ genome.size, data = temp_data)
summary(temp_model_E)

AIC(temp_model_A, temp_model_B, temp_model_C, temp_model_D, temp_model_E)

# pH enzyme ----

# Linear Model 
pH_model_A       = lm(pH ~ genome.size, data = MAG_gen_test)
summary(pH_model_A)

# Genome Size and Completeness
pH_model_B       = lm(pH ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(pH_model_B)

domin(pH ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
pH_model_C       = lm(pH ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(pH_model_C)

domin(pH ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
pH_model_D       = lm(pH ~ GC_count_mean, data = MAG_gen_test)
summary(pH_model_D)

# Power Law
pH_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                       log(MAG_gen_test$pH)))

pH_data           = subset(pH_data, pH_data$V2 != "-Inf") 
colnames(pH_data) = c("genome.size","pH")

pH_model_E       = lm(pH ~ genome.size, data = pH_data)
summary(pH_model_E)

AIC(pH_model_A, pH_model_B, pH_model_C, pH_model_D, pH_model_E)

# mgr ----

MAG_gen_test     = MAG_gen_test %>% mutate(mgr = 1/mgr)

# Linear Model 
mgr_model_A       = lm(mgr ~ genome.size, data = MAG_gen_test)
summary(mgr_model_A)

# Genome Size and Completeness
mgr_model_B       = lm(mgr ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(mgr_model_B)

domin(mgr ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
mgr_model_C       = lm(mgr ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(mgr_model_C)

domin(mgr ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
mgr_model_D       = lm(mgr ~ GC_count_mean, data = MAG_gen_test)
summary(mgr_model_D)

# Power Law
mgr_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                        log(MAG_gen_test$mgr)))

mgr_data           = subset(mgr_data, mgr_data$V2 != "-Inf") 
colnames(mgr_data) = c("genome.size","mgr")

mgr_model_E       = lm(mgr ~ genome.size, data = mgr_data)
summary(mgr_model_E)

AIC(mgr_model_A, mgr_model_B, mgr_model_C, mgr_model_D, mgr_model_E)

# ogt ----

# Linear Model 
ogt_model_A       = lm(ogt ~ genome.size, data = MAG_gen_test)
summary(ogt_model_A)

# Genome Size and Completeness
ogt_model_B       = lm(ogt ~ genome.size*completeness_mean, data = MAG_gen_test)
summary(ogt_model_B)

domin(ogt ~ genome.size*completeness_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# Genome size and Contamination
ogt_model_C       = lm(ogt ~ genome.size*contamination_mean, data = MAG_gen_test)
summary(ogt_model_C)

domin(ogt ~ genome.size*contamination_mean, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

# GC count
ogt_model_D       = lm(ogt ~ GC_count_mean, data = MAG_gen_test)
summary(ogt_model_D)

# Power Law
ogt_data          = as.data.frame(cbind(log(MAG_gen_test$genome.size), 
                                        log(MAG_gen_test$ogt)))

ogt_data           = subset(ogt_data, ogt_data$V2 != "-Inf") 
colnames(ogt_data) = c("genome.size","ogt")

ogt_model_E       = lm(ogt ~ genome.size, data = ogt_data)
summary(ogt_model_E)

AIC(ogt_model_A, ogt_model_B, ogt_model_C, ogt_model_D, ogt_model_E)

# S vs A tradeoffs (Isolates) ----

ISO_gen_TOTAL_ab = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_test     = ISO_gen_TOTAL_ab %>% mutate(A_trait = rowSums(ISO_gen_TOTAL_ab[,9:11]),
                                               S_trait = rowSums(ISO_gen_TOTAL_ab[,12:15]))
ISO_gen_test     = ISO_gen_test %>% mutate(A_S = A_trait/S_trait)

ISO_gen_test     = ISO_gen_test %>% mutate(MGR = 1/MGT)

# CAZy enzyme ----

# Linear Model 
cazy_model_A       = lm(CAZy ~ Genome.Size, data = ISO_gen_test)
summary(cazy_model_A)

# Power Law
CAZy_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                         log(ISO_gen_test$CAZy)))

CAZy_data           = subset(CAZy_data, CAZy_data$V2 != "-Inf") 
colnames(CAZy_data) = c("genome.size","CAZy")

cazy_model_E       = lm(CAZy ~ genome.size, data = CAZy_data)
summary(cazy_model_E)

AIC(cazy_model_A, cazy_model_E)

# Protein.enzyme enzyme ----

# Linear Model 
Protein.enzyme_model_A       = lm(Protein.enzyme ~ Genome.Size, data = ISO_gen_test)
summary(Protein.enzyme_model_A)

# Power Law
Protein.enzyme_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                                   log(ISO_gen_test$Protein.enzyme)))

Protein.enzyme_data           = subset(Protein.enzyme_data, Protein.enzyme_data$V2 != "-Inf") 
colnames(Protein.enzyme_data) = c("genome.size","Protein.enzyme")

Protein.enzyme_model_E       = lm(Protein.enzyme ~ genome.size, data = Protein.enzyme_data)
summary(Protein.enzyme_model_E)

AIC(Protein.enzyme_model_A, Protein.enzyme_model_E)

# Total.transporter enzyme ----

# Linear Model 
Total.transporter_model_A       = lm(Total.transporter ~ Genome.Size, data = ISO_gen_test)
summary(Total.transporter_model_A)

# Power Law
Total.transporter_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                                      log(ISO_gen_test$Total.transporter)))

Total.transporter_data           = subset(Total.transporter_data, Total.transporter_data$V2 != "-Inf") 
colnames(Total.transporter_data) = c("genome.size","Total.transporter")

Total.transporter_model_E       = lm(Total.transporter ~ genome.size, data = Total.transporter_data)
summary(Total.transporter_model_E)

AIC(Total.transporter_model_A, Total.transporter_model_E)

# GH.transporter enzyme ----

# Linear Model 
GH.transporter_model_A       = lm(GH.transporter ~ Genome.Size, data = ISO_gen_test)
summary(GH.transporter_model_A)

# Power Law
GH.transporter_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                                   log(ISO_gen_test$GH.transporter)))

GH.transporter_data           = subset(GH.transporter_data, GH.transporter_data$V2 != "-Inf") 
colnames(GH.transporter_data) = c("genome.size","GH.transporter")

GH.transporter_model_E       = lm(GH.transporter ~ genome.size, data = GH.transporter_data)
summary(GH.transporter_model_E)

AIC(GH.transporter_model_A, GH.transporter_model_E)

# Amino.transporter enzyme ----

# Linear Model 
Amino.transporter_model_A       = lm(Amino.transporter ~ Genome.Size, data = ISO_gen_test)
summary(Amino.transporter_model_A)

# Power Law
Amino.transporter_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                                      log(ISO_gen_test$Amino.transporter)))

Amino.transporter_data           = subset(Amino.transporter_data, Amino.transporter_data$V2 != "-Inf") 
colnames(Amino.transporter_data) = c("genome.size","Amino.transporter")

Amino.transporter_model_E       = lm(Amino.transporter ~ genome.size, data = Amino.transporter_data)
summary(Amino.transporter_model_E)

AIC(Amino.transporter_model_A, Amino.transporter_model_E)

# Osmolyte enzyme ----

# Linear Model 
Osmolyte_model_A       = lm(Osmolyte ~ Genome.Size, data = ISO_gen_test)
summary(Osmolyte_model_A)

# Power Law
Osmolyte_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                             log(ISO_gen_test$Osmolyte)))

Osmolyte_data           = subset(Osmolyte_data, Osmolyte_data$V2 != "-Inf") 
colnames(Osmolyte_data) = c("genome.size","Osmolyte")

Osmolyte_model_E       = lm(Osmolyte ~ genome.size, data = Osmolyte_data)
summary(Osmolyte_model_E)

AIC(Osmolyte_model_A, Osmolyte_model_E)

# Biofilm enzyme ----

# Linear Model 
Biofilm_model_A       = lm(Biofilm ~ Genome.Size, data = ISO_gen_test)
summary(Biofilm_model_A)

# Power Law
Biofilm_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                            log(ISO_gen_test$Biofilm)))

Biofilm_data           = subset(Biofilm_data, Biofilm_data$V2 != "-Inf") 
colnames(Biofilm_data) = c("genome.size","Biofilm")

Biofilm_model_E       = lm(Biofilm ~ genome.size, data = Biofilm_data)
summary(Biofilm_model_E)

AIC(Biofilm_model_A, Biofilm_model_E)

# Temp.Tol enzyme ----

# Linear Model 
Temp.Tol_model_A       = lm(Temp.Tol ~ Genome.Size, data = ISO_gen_test)
summary(Temp.Tol_model_A)

# Power Law
Temp.Tol_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                             log(ISO_gen_test$Temp.Tol)))

Temp.Tol_data           = subset(Temp.Tol_data, Temp.Tol_data$V2 != "-Inf") 
colnames(Temp.Tol_data) = c("genome.size","Temp.Tol")

Temp.Tol_model_E       = lm(Temp.Tol ~ genome.size, data = Temp.Tol_data)
summary(Temp.Tol_model_E)

AIC(Temp.Tol_model_A, Temp.Tol_model_E)

# pH.Tol enzyme ----

# Linear Model 
pH.Tol_model_A       = lm(pH.Tol ~ Genome.Size, data = ISO_gen_test)
summary(pH.Tol_model_A)

# Power Law
pH.Tol_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                           log(ISO_gen_test$pH.Tol)))

pH.Tol_data           = subset(pH.Tol_data, pH.Tol_data$V2 != "-Inf") 
colnames(pH.Tol_data) = c("genome.size","pH.Tol")

pH.Tol_model_E       = lm(pH.Tol ~ genome.size, data = pH.Tol_data)
summary(pH.Tol_model_E)

AIC(pH.Tol_model_A, pH.Tol_model_E)

# Yield ----

# Linear Model 
Yield_model_A       = lm(Yield ~ Genome.Size, data = ISO_gen_test)
summary(Yield_model_A)

# Power Law
Yield_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                          log(ISO_gen_test$Yield)))

Yield_data           = subset(Yield_data, Yield_data$V2 != "-Inf") 
colnames(Yield_data) = c("genome.size","Yield")

Yield_model_E       = lm(Yield ~ genome.size, data = Yield_data)
summary(Yield_model_E)

AIC(Yield_model_A, Yield_model_E)

# MGR ----

# Linear Model 
MGR_model_A       = lm(MGR ~ Genome.Size, data = ISO_gen_test)
summary(MGR_model_A)

# Power Law
MGR_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                        log(ISO_gen_test$MGR)))

MGR_data           = subset(MGR_data, MGR_data$V2 != "-Inf") 
colnames(MGR_data) = c("genome.size","MGR")

MGR_model_E       = lm(MGR ~ genome.size, data = MGR_data)
summary(MGR_model_E)

AIC(MGR_model_A, MGR_model_E)

# OGT ----

# Linear Model 
OGT_model_A       = lm(OGT ~ Genome.Size, data = ISO_gen_test)
summary(OGT_model_A)

# Power Law
OGT_data          = as.data.frame(cbind(log(ISO_gen_test$Genome.Size), 
                                        log(ISO_gen_test$OGT)))

OGT_data           = subset(OGT_data, OGT_data$V2 != "-Inf") 
colnames(OGT_data) = c("genome.size","OGT")

OGT_model_E       = lm(OGT ~ genome.size, data = OGT_data)
summary(OGT_model_E)

AIC(OGT_model_A, OGT_model_E)

# Power Law models (Selected) - MAGs ----

# MAGs gene counts ----

MAG_gen_TOTAL_ab = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")
MAG_gen_TOTAL_ab = MAG_gen_TOTAL_ab %>% mutate(speed = case_when(mgr <= 5 ~ "fast",
                                                                 mgr  > 5 ~ "slow"))

MAG_gen_TOTAL_ab_fast = MAG_gen_TOTAL_ab %>% filter(speed == "fast")
MAG_gen_TOTAL_ab_slow = MAG_gen_TOTAL_ab %>% filter(speed == "slow")

# CAZy----

CAZy_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                     log(MAG_gen_TOTAL_ab_fast$CAZy), 
                                     rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(CAZy_data_fast) = c("genome.size","CAZy","speed")
CAZy_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                     log(MAG_gen_TOTAL_ab_slow$CAZy), 
                                     rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(CAZy_data_slow) = c("genome.size","CAZy","speed")

CAZy_data     = as.data.frame(rbind(CAZy_data_fast,CAZy_data_slow))

Figure_CAZy   = ggplot(data = CAZy_data, aes(x = as.numeric(genome.size), 
                                             y = as.numeric(CAZy))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log(Genome Size)") + 
  ylab("Log(CAZy)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.8,3.5)
Figure_CAZy

# Protein enzyme ----

Protein_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                        log(MAG_gen_TOTAL_ab_fast$Protein), 
                                        rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(Protein_data_fast) = c("genome.size","Protein","speed")
Protein_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                        log(MAG_gen_TOTAL_ab_slow$Protein), 
                                        rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(Protein_data_slow) = c("genome.size","Protein","speed")

Protein_data     = as.data.frame(rbind(Protein_data_fast,Protein_data_slow))

Figure_Protein   = ggplot(data = Protein_data, aes(x = as.numeric(genome.size), 
                                                   y = as.numeric(Protein))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log(Genome Size)") + 
  ylab("Log(Protein)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.8,4)
Figure_Protein

# Transport total Transporters ----

tranport.total_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                               log(MAG_gen_TOTAL_ab_fast$tranport.total), 
                                               rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(tranport.total_data_fast) = c("genome.size","tranport.total","speed")
tranport.total_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                               log(MAG_gen_TOTAL_ab_slow$tranport.total), 
                                               rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(tranport.total_data_slow) = c("genome.size","tranport.total","speed")

tranport.total_data.1     = as.data.frame(rbind(tranport.total_data_fast,
                                                tranport.total_data_slow))

Figure_tranport.total   = ggplot(data = tranport.total_data.1, aes(x = as.numeric(genome.size), 
                                                                   y = as.numeric(tranport.total))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.5,5.5)
Figure_tranport.total

# GH total Transporters ----

GH.trasnport_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                             log(MAG_gen_TOTAL_ab_fast$GH.trasnport), 
                                             rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(GH.trasnport_data_fast) = c("genome.size","GH.trasnport","speed")
GH.trasnport_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                             log(MAG_gen_TOTAL_ab_slow$GH.trasnport), 
                                             rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(GH.trasnport_data_slow) = c("genome.size","GH.trasnport","speed")

GH.trasnport_data     = as.data.frame(rbind(GH.trasnport_data_fast,
                                            GH.trasnport_data_slow))

Figure_GH.trasnport   = ggplot(data = GH.trasnport_data, aes(x = as.numeric(genome.size), 
                                                             y = as.numeric(GH.trasnport))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (GH transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.5,4.0)
Figure_GH.trasnport

# Amino total Transporters ----

amino.transport_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                                log(MAG_gen_TOTAL_ab_fast$amino.transport), 
                                                rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(amino.transport_data_fast) = c("genome.size","amino.transport","speed")
amino.transport_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                                log(MAG_gen_TOTAL_ab_slow$amino.transport), 
                                                rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(amino.transport_data_slow) = c("genome.size","amino.transport","speed")

amino.transport_data     = as.data.frame(rbind(amino.transport_data_fast,
                                               amino.transport_data_slow))

Figure_amino.transport   = ggplot(data = amino.transport_data, aes(x = as.numeric(genome.size), 
                                                                   y = as.numeric(amino.transport))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log (Genome Size)") + 
  ylab("Log (Animo transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0,3.5)
Figure_amino.transport

# Osmolytes ----

osmolyte_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                         log(MAG_gen_TOTAL_ab_fast$osmolyte), 
                                         rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(osmolyte_data_fast) = c("genome.size","osmolyte","speed")
osmolyte_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                         log(MAG_gen_TOTAL_ab_slow$osmolyte), 
                                         rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(osmolyte_data_slow) = c("genome.size","osmolyte","speed")

osmolyte_data     = as.data.frame(rbind(osmolyte_data_fast,
                                        osmolyte_data_slow))

Figure_osmolyte  = ggplot(data = osmolyte_data, aes(x = as.numeric(genome.size), 
                                                    y = as.numeric(osmolyte))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Osmolytes)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.7,3.5)
Figure_osmolyte

# Biofilm ----

biofilm_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                        log(MAG_gen_TOTAL_ab_fast$biofilm), 
                                        rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(biofilm_data_fast) = c("genome.size","biofilm","speed")
biofilm_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                        log(MAG_gen_TOTAL_ab_slow$biofilm), 
                                        rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(biofilm_data_slow) = c("genome.size","biofilm","speed")

biofilm_data     = as.data.frame(rbind(biofilm_data_fast,
                                       biofilm_data_slow))

biofilm_data     = subset(biofilm_data, biofilm!="-Inf") 

Figure_biofilm  = ggplot(data = biofilm_data, aes(x = as.numeric(genome.size), 
                                                  y = as.numeric(biofilm))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Biofilm)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.2,2.8)
Figure_biofilm

# Heat Tolerance ----

temp_fast_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                          log(MAG_gen_TOTAL_ab_fast$temp), 
                                          rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(temp_fast_data_fast) = c("genome.size","temp_fast","speed")
temp_fast_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                          log(MAG_gen_TOTAL_ab_slow$temp), 
                                          rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(temp_fast_data_slow) = c("genome.size","temp_fast","speed")

temp_fast_data     = as.data.frame(rbind(temp_fast_data_fast,
                                         temp_fast_data_slow))

# temp_fast_data     = subset(biofilm_data, biofilm!="-Inf") 

Figure_temp_fast  = ggplot(data = temp_fast_data, aes(x = as.numeric(genome.size), 
                                                      y = as.numeric(temp_fast))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Temperature Tolerance)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.5,2.6)
Figure_temp_fast

# pH Tolerance Transporters ----

pH_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                   log(MAG_gen_TOTAL_ab_fast$pH), 
                                   rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(pH_data_fast) = c("genome.size","pH","speed")
pH_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                   log(MAG_gen_TOTAL_ab_slow$pH), 
                                   rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(pH_data_slow) = c("genome.size","pH","speed")

pH_data     = as.data.frame(rbind(pH_data_fast,
                                  pH_data_slow))

pH_data     = subset(pH_data, pH!="-Inf") 

Figure_pH_fast  = ggplot(data = pH_data, aes(x = as.numeric(genome.size), 
                                             y = as.numeric(pH))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (pH Tolerance)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.7,2.85)
Figure_pH_fast

# Minimum Generation Time ----

mgr_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                    log(1/MAG_gen_TOTAL_ab_fast$mgr), 
                                    rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(mgr_data_fast) = c("genome.size","mgr","speed")
mgr_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                    log(1/MAG_gen_TOTAL_ab_slow$mgr), 
                                    rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(mgr_data_slow) = c("genome.size","mgr","speed")

mgr_data     = as.data.frame(rbind(mgr_data_fast,
                                   mgr_data_slow))

Figure_mgr_fast  = ggplot(data = mgr_data, aes(x = as.numeric(genome.size), 
                                               y = as.numeric(mgr))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "top", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (mgr)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.55,3.45)
Figure_mgr_fast

# Optimum Growth Temperature ----

ogt_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                    log(MAG_gen_TOTAL_ab_fast$ogt), 
                                    rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(ogt_data_fast) = c("genome.size","ogt","speed")
ogt_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                    log(MAG_gen_TOTAL_ab_slow$ogt), 
                                    rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(ogt_data_slow) = c("genome.size","ogt","speed")

ogt_data     = as.data.frame(rbind(ogt_data_fast,
                                   ogt_data_slow))

Figure_ogt_fast  = ggplot(data = ogt_data, aes(x = as.numeric(genome.size), 
                                               y = as.numeric(ogt))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (OGT)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(2.9,4.3)
Figure_ogt_fast

# Power Law models - Isolates ----

# Isolate gene counts ----

ISO_gen_TOTAL_ab = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_TOTAL_ab = ISO_gen_TOTAL_ab %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                                                 MGT  > 5 ~ "slow"))

ISO_gen_TOTAL_ab_fast = ISO_gen_TOTAL_ab %>% filter(speed == "fast")
ISO_gen_TOTAL_ab_slow = ISO_gen_TOTAL_ab %>% filter(speed == "slow")

# CAZy----

CAZy.total_data_fast = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                           log(ISO_gen_TOTAL_ab_fast$CAZy), 
                                           rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(CAZy.total_data_fast) = c("genome.size","CAZy","speed")
CAZy.total_data_slow = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                           log(ISO_gen_TOTAL_ab_slow$CAZy), 
                                           rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(CAZy.total_data_slow) = c("genome.size","CAZy","speed")

CAZy.total_data.ISO     = as.data.frame(rbind(CAZy.total_data_fast,
                                              CAZy.total_data_slow))

Figure_CAZy.total.ISO   = ggplot(data = CAZy.total_data.ISO, aes(x = as.numeric(genome.size), 
                                                                 y = as.numeric(CAZy))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log (Genome Size)") + 
  ylab("Log (CAZy)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.8,3.5)
Figure_CAZy.total.ISO

# Protein enzyme ----

Protein.enzyme.total_data_fast = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                                     log(ISO_gen_TOTAL_ab_fast$Protein.enzyme), 
                                                     rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(Protein.enzyme.total_data_fast) = c("genome.size","Protein.enzyme","speed")
Protein.enzyme.total_data_slow = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                                     log(ISO_gen_TOTAL_ab_slow$Protein.enzyme), 
                                                     rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(Protein.enzyme.total_data_slow) = c("genome.size","Protein.enzyme","speed")

Protein.enzyme.total_data.ISO     = as.data.frame(rbind(Protein.enzyme.total_data_fast,
                                                        Protein.enzyme.total_data_slow))

Figure_Protein.enzyme.total.ISO   = ggplot(data = Protein.enzyme.total_data.ISO, aes(x = as.numeric(genome.size), 
                                                                                     y = as.numeric(Protein.enzyme))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log (Genome Size)") + 
  ylab("Log (Protein)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.8,4)
Figure_Protein.enzyme.total.ISO

# Transport total Transporters ----

tranport.total_data_fast = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                               log(ISO_gen_TOTAL_ab_fast$Total.transporter), 
                                               rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(tranport.total_data_fast) = c("genome.size","tranport.total","speed")
tranport.total_data_slow = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                               log(ISO_gen_TOTAL_ab_slow$Total.transporter), 
                                               rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(tranport.total_data_slow) = c("genome.size","tranport.total","speed")

tranport.total_data     = as.data.frame(rbind(tranport.total_data_fast,
                                              tranport.total_data_slow))

Figure_tranport.total_ISO   = ggplot(data = tranport.total_data, aes(x = as.numeric(genome.size), 
                                                                     y = as.numeric(tranport.total))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.5,5.5)
Figure_tranport.total_ISO

# GH transporter ----

GH.transporter_data_fast = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                               log(ISO_gen_TOTAL_ab_fast$GH.transporter), 
                                               rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(GH.transporter_data_fast) = c("genome.size","GH.transporter","speed")
GH.transporter_data_slow = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                               log(ISO_gen_TOTAL_ab_slow$GH.transporter), 
                                               rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(GH.transporter_data_slow) = c("genome.size","GH.transporter","speed")

GH.transporter_data     = as.data.frame(rbind(GH.transporter_data_fast,
                                              GH.transporter_data_slow))

Figure_GH.transporter   = ggplot(data = GH.transporter_data, aes(x = as.numeric(genome.size), 
                                                                 y = as.numeric(GH.transporter))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (GH transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.5,4.0)
Figure_GH.transporter

# Amino transporter ----

Amino.transporter_data_fast = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                                  log(ISO_gen_TOTAL_ab_fast$Amino.transporter), 
                                                  rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(Amino.transporter_data_fast) = c("genome.size","Amino.transporter","speed")
Amino.transporter_data_slow = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                                  log(ISO_gen_TOTAL_ab_slow$Amino.transporter), 
                                                  rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(Amino.transporter_data_slow) = c("genome.size","Amino.transporter","speed")

Amino.transporter_data     = as.data.frame(rbind(Amino.transporter_data_fast,
                                                 Amino.transporter_data_slow))

Figure_Amino.transporter   = ggplot(data = Amino.transporter_data, aes(x = as.numeric(genome.size), 
                                                                       y = as.numeric(Amino.transporter))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log (Genome Size)") + 
  ylab("Log (Amino transporters)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0,3.5)
Figure_Amino.transporter

# Osmolyte ----

osmolyte_data_fast_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                             log(ISO_gen_TOTAL_ab_fast$Osmolyte), 
                                             rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(osmolyte_data_fast_ISO) = c("genome.size","osmolyte","speed")
osmolyte_data_slow_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                             log(ISO_gen_TOTAL_ab_slow$Osmolyte), 
                                             rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(osmolyte_data_slow_ISO) = c("genome.size","osmolyte","speed")

osmolyte_data_ISO     = as.data.frame(rbind(osmolyte_data_fast_ISO,
                                            osmolyte_data_slow_ISO))

Figure_osmolyte_ISO   = ggplot(data = osmolyte_data_ISO, aes(x = as.numeric(genome.size), 
                                                             y = as.numeric(osmolyte))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Osmolytes)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.7,3.5)
Figure_osmolyte_ISO

# Biofilm ----

Biofilm_data_fast_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                            log(ISO_gen_TOTAL_ab_fast$Biofilm), 
                                            rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(Biofilm_data_fast_ISO) = c("genome.size","Biofilm","speed")
Biofilm_data_slow_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                            log(ISO_gen_TOTAL_ab_slow$Biofilm), 
                                            rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(Biofilm_data_slow_ISO) = c("genome.size","Biofilm","speed")

Biofilm_data_ISO     = as.data.frame(rbind(Biofilm_data_fast_ISO,
                                           Biofilm_data_slow_ISO))

Figure_Biofilm_ISO   = ggplot(data = Biofilm_data_ISO, aes(x = as.numeric(genome.size), 
                                                           y = as.numeric(Biofilm))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Biofilm)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.2,2.8)
Figure_Biofilm_ISO

# Heat Tolerance ----

Temp.Tol_data_fast_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                             log(ISO_gen_TOTAL_ab_fast$Temp.Tol), 
                                             rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(Temp.Tol_data_fast_ISO) = c("genome.size","Temp.Tol","speed")
Temp.Tol_data_slow_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                             log(ISO_gen_TOTAL_ab_slow$Temp.Tol), 
                                             rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(Temp.Tol_data_slow_ISO) = c("genome.size","Temp.Tol","speed")

Temp.Tol_data_ISO     = as.data.frame(rbind(Temp.Tol_data_fast_ISO,
                                            Temp.Tol_data_slow_ISO))

Figure_Temp.Tol_ISO   = ggplot(data = Temp.Tol_data_ISO, aes(x = as.numeric(genome.size), 
                                                             y = as.numeric(Temp.Tol))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Temperature Tolerance)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.5,2.6)
Figure_Temp.Tol_ISO

# pH Tolerance ----

pH.Tol_data_fast_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                           log(ISO_gen_TOTAL_ab_fast$pH.Tol), 
                                           rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(pH.Tol_data_fast_ISO) = c("genome.size","pH.Tol","speed")
pH.Tol_data_slow_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                           log(ISO_gen_TOTAL_ab_slow$pH.Tol), 
                                           rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(pH.Tol_data_slow_ISO) = c("genome.size","pH.Tol","speed")

pH.Tol_data_ISO     = as.data.frame(rbind(pH.Tol_data_fast_ISO,
                                          pH.Tol_data_slow_ISO))

Figure_pH.Tol_ISO   = ggplot(data = pH.Tol_data_ISO, aes(x = as.numeric(genome.size), 
                                                         y = as.numeric(pH.Tol))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (pH Tolerance)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.7,2.85)
Figure_pH.Tol_ISO

# Yield ----

Yield_data_fast = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                      log(ISO_gen_TOTAL_ab_fast$Yield), 
                                      rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(Yield_data_fast) = c("genome.size","Yield","speed")
Yield_data_slow = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                      log(ISO_gen_TOTAL_ab_slow$Yield), 
                                      rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(Yield_data_slow) = c("genome.size","Yield","speed")

Yield_data     = as.data.frame(rbind(Yield_data_fast,
                                     Yield_data_slow))

Figure_Yield_ISO   = ggplot(data = Yield_data, aes(x = as.numeric(genome.size), 
                                                   y = as.numeric(Yield))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Yield)") +
  theme(legend.position="none") + xlim(13,16.5)
Figure_Yield_ISO

# Minimum Generation Time ----

mgr_data_fast_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                        log(1/ISO_gen_TOTAL_ab_fast$MGT), 
                                        rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(mgr_data_fast_ISO) = c("genome.size","mgr","speed")
mgr_data_slow_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                        log(1/ISO_gen_TOTAL_ab_slow$MGT), 
                                        rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(mgr_data_slow_ISO) = c("genome.size","mgr","speed")

mgr_data_ISO     = as.data.frame(rbind(mgr_data_fast_ISO,
                                       mgr_data_slow_ISO))

Figure_mgr_ISO   = ggplot(data = mgr_data_ISO, aes(x = as.numeric(genome.size), 
                                                   y = as.numeric(mgr))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (mgr)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.55,3.45)
Figure_mgr_ISO

# Optimum Growth Temperature ----

OGT_data_fast_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                        log(ISO_gen_TOTAL_ab_fast$OGT), 
                                        rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(OGT_data_fast_ISO) = c("genome.size","OGT","speed")
OGT_data_slow_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                        log(ISO_gen_TOTAL_ab_slow$OGT), 
                                        rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(OGT_data_slow_ISO) = c("genome.size","OGT","speed")

OGT_data_ISO     = as.data.frame(rbind(OGT_data_fast_ISO,
                                       OGT_data_slow_ISO))

Figure_OGT_ISO   = ggplot(data = OGT_data_ISO, aes(x = as.numeric(genome.size), 
                                                   y = as.numeric(OGT))) +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (OGT)") +
  theme(legend.position="none") + xlim(13,16.5) + ylim(2.9,4.3)
Figure_OGT_ISO

# Figure 2 (Enzyme encoding genes) ----

Figure_2 = ggarrange(Figure_CAZy, Figure_Protein,Figure_CAZy.total.ISO,
                     Figure_Protein.enzyme.total.ISO,
                     labels = c("A","B","C","D"),ncol = 2, nrow = 2) + 
  theme(panel.background = element_blank())

Figure_2
pdf("Output_Data/Figures/Figure_2.pdf",
    width=14,height=14*5/5)
print(Figure_2)
dev.off()

# Figure 3 (Transporters encoding genes) ----

Figure_3 = ggarrange(Figure_tranport.total, Figure_GH.trasnport,
                     Figure_amino.transport,Figure_tranport.total_ISO,
                     Figure_GH.transporter,Figure_Amino.transporter,
                     labels = c("A","B","C","D","E","F"),ncol = 3, nrow = 2) + 
  theme(panel.background = element_blank())

Figure_3
pdf("Output_Data/Figures/Figure_3.pdf",
    width=15,height=15*4/5)
print(Figure_3)
dev.off()

# Figure 4(Stress Tolerance encoding genes) ----

Figure_4 = ggarrange(Figure_osmolyte,Figure_biofilm,Figure_temp_fast,
                     Figure_pH_fast,Figure_osmolyte_ISO,Figure_Biofilm_ISO,
                     Figure_Temp.Tol_ISO,Figure_pH.Tol_ISO,
                     labels = c("A","B","C","D","E","F","G","H"),ncol = 4, 
                     nrow = 2) + theme(panel.background = element_blank())
Figure_4

pdf("Output_Data/Figures/Figure_4.pdf",
    width=16,height=16*3/5)
print(Figure_4)
dev.off()

# Figure 5 (Life History encoding genes (Yield)) ----

Figure_5 = ggarrange(Figure_mgr_fast, Figure_mgr_ISO, Figure_Yield_ISO,
                     labels = c("A","B","C"),ncol = 3, nrow = 1) + 
  theme(panel.background = element_blank())

Figure_5

pdf("Output_Data/Figures/Figure_5.pdf",
    width=12,height=12*1/3)
print(Figure_5)
dev.off()

# Figure 6 (Life History encoding genes (Temperature)) ----

Figure_6 = ggarrange(Figure_ogt_fast, Figure_OGT_ISO,labels = c("A","B"),
                     ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Figure_6

pdf("Output_Data/Figures/Figure_6.pdf",
    width=12,height=12*1/2)
print(Figure_6)
dev.off()

# CUE for different taxonomic levels ----

isolates     = read.csv("Input_Data/SOIL_ISOLATES/dement_isolates_CUE.csv",dec=".")
isolatest    = read.csv("Input_Data/SOIL_ISOLATES/total.granularity.3_datasets.csv",dec=".")
isolates.1   = subset(isolates, isolates$CUE !='NaN')
meta_iso     = read.delim("Input_Data/SOIL_ISOLATES/metadata_2.tsv",sep="\t")
meta_iso     = as.data.frame(cbind(meta_iso$IMG.Genome.ID,meta_iso$Phylum,
                                 meta_iso$Class,meta_iso$Order,meta_iso$Family,
                                 meta_iso$Genus))
meta_iso.1   = read.delim("Input_Data/SOIL_ISOLATES/metadata_1.tsv",sep="\t")
meta_iso.1   = as.data.frame(cbind(meta_iso.1$IMG.Genome.ID,meta_iso.1$Phylum,
                                 meta_iso.1$Class,meta_iso.1$Order,meta_iso.1$Family,
                                 meta_iso.1$Genus))

# final meta
final_meta   = as.data.frame(rbind(meta_iso,meta_iso.1))
colnames(final_meta) = c("id","Phylum","Class","Order","Family","Genus")

merged_df    = merge(isolates.1,final_meta,by = "id")
merged_df    = merged_df %>% filter(CUE <= 0.9)

# Figure S5A - Phylum ----

phylum = merged_df %>% group_by(Phylum) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
phylum = phylum %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                             mingentime_mean  > 5 ~ "slow"))

Figure_test_1.phylum = ggplot(data = phylum, aes(x = log(genome_length_mean), 
                                                 y = log(CUE_mean))) +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  labs(title = "Phylum")
Figure_test_1.phylum

# Figure S5B - Class ----

Class = merged_df %>% group_by(Class) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Class = Class %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))
Figure_test_1.Class = ggplot(data = Class, aes(x = log(genome_length_mean), 
                                               y = log(CUE_mean))) +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") + 
  labs(title = "Class")
Figure_test_1.Class

# Figure S5C - Order ----

Order = merged_df %>% group_by(Order) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Order = Order %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Order = ggplot(data = Order, aes(x = log(genome_length_mean), 
                                               y = log(CUE_mean))) + 
  geom_point() + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  labs(title = "Order")
Figure_test_1.Order

# Figure S5D - Family ----

Family = merged_df %>% group_by(Family) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
Family = Family %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                             mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Family = ggplot(data = Family, aes(x = log(genome_length_mean), 
                                                 y = log(CUE_mean))) +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  labs(title = "Family")
Figure_test_1.Family

# Figure S5E - Genus ----

Genus = merged_df %>% group_by(Genus) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Genus = Genus %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Genus = ggplot(data = Genus, aes(x = log(genome_length_mean), 
                                               y = log(CUE_mean))) +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  labs(title = "Genus")
Figure_test_1.Genus

# Figure S5 ----

Figure_S5 = ggarrange(Figure_test_1.phylum,Figure_test_1.Class,
                      Figure_test_1.Order,Figure_test_1.Family,Figure_test_1.Genus, 
                      labels = c("A","B","C","D","E"),
                      ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_S5

pdf("Output_Data/Figures/Figure_S5.pdf",
    width=14,height=14*3/5)
print(Figure_S5)
dev.off()

# Summary statistics from the MAGs ----

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

# Figure S6 ----

Figure_S6 = ggarrange(Figure_S6.A,Figure_S6.B,
                      labels = c("A","B"),
                      ncol = 2, nrow = 1) + theme(panel.background = element_blank())
Figure_S6

pdf("Output_Data/Figures/Figure_S6.pdf",
    width=4.81*1.5,height=2.21*1.5)
print(Figure_S6)
dev.off()

# New Figures arised during revision ----

# Scatterplot Yield vs A+S traits ----

pdf("Output_Data/Figures/Figure_S9.pdf",
    width=2.5*2.5,height=2.5*2.5)
plot(ISO_gen_TOTAL_ab$A_traits+ISO_gen_TOTAL_ab$S_traits, ISO_gen_TOTAL_ab$Yield,
     xlab = "Total S+A gene counts per functional group",
     ylab = "Yield (CUE)",
     main = "Yield vs S+A gene counts — isolate functional groups")
abline(lm(ISO_gen_TOTAL_ab$Yield ~ ISO_gen_TOTAL_ab$A_traits+ISO_gen_TOTAL_ab$S_traits), col = "red", lwd = 2)
legend("bottomleft", 
       legend = "r = -0.136, p = 0.042",
       bty = "n")
dev.off()


