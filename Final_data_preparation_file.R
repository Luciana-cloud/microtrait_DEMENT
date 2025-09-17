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

Figure_S1 = ggplot(data,aes(V1,V2)) +
  geom_point() +
  geom_smooth(method='lm') + 
  stat_cor(label.x = 30, label.y = 130, size = 4) +
  stat_regline_equation(label.x = 30, label.y = 150, size = 4) + 
  xlab("# of MAGs") + ylab("# of Functional Groups") + theme_classic() + 
  theme(text = element_text(size=20))

png("Output_Data/Figures/Figure_S1.png",
    width=3500,height=1969,res=300)
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

# Linear regressions with genome size for MAGs ----

MAG_gen_trait.2  = read.csv("Intermediate_Results/MAG_gen_trait.csv",dec=".")
MAG_gen_trait.2  = MAG_gen_trait.2[, -1]
MAG_gen_trait.3  = read.csv("Intermediate_Results/MAG_gen_trait.5.csv",dec=".")
MAG_gen_trait.3  = MAG_gen_trait.3[, -1]

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

# Figure 3A - CAZy enzyme ----

Figure_test_1.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (CAZy),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("CAZy (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,30) + theme(legend.position="none")
Figure_test_1.new

# Figure 3B - Protein enzyme ----

Figure_test_2.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Protein-enzyme`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Protein Enzymes (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,50) + theme(legend.position="none")
Figure_test_2.new

# Figure 4A - Transport ----

Figure_test_3.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Total-transporter`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Total Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,200) + theme(legend.position="none")
Figure_test_3.new

# Figure 4B - Aminoacids Transport ----

Figure_test_4.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Amino-transporter`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Aminoacid Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_4.new

# Figure 4C - Carbohydrate Transport ----

Figure_test_5.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`GH-transporter`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Carbohydrate Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,55) + theme(legend.position="none")
Figure_test_5.new

# Figure 5A - Osmolytes ----

Figure_test_6.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (Osmolyte),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Osmolytes (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32) + theme(legend.position="none")
Figure_test_6.new

# Figure 5B - Biofilm ----

Figure_test_7.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (Biofilm),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Biofilm (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,16) + theme(legend.position="none")
Figure_test_7.new

# Figure 5C - Heat Resistance ----

Figure_test_8.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`Temp-Tol`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("Heat Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,15) + theme(legend.position="none")
Figure_test_8.new

# Figure 5D - pH Resistance ----

Figure_test_9.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (`pH-Tol`),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("") + 
  ylab("pH Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_9.new

# Figure 6A - Minimum generation time ----

Figure_test_10a.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (MGT),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Minimum generation time (hrs)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,40) + theme(legend.position="none")
Figure_test_10a.new

# Figure 7A - Optimal Growth Temperature ----

Figure_test_11.new = ggplot(data = Total, aes(x = `Genome-Size`/1e6, y = (OGT),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Optimal Growth Temperature (°C)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,80) + theme(legend.position="none")
Figure_test_11.new

# Linear regressions with genome size for Isolates ----

ISO_gen_trait.1  = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_trait.1  = ISO_gen_trait.1[, -1]
Total.iso        = ISO_gen_trait.1 %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                                                MGT  > 5 ~ "slow"))
# Erase outlier isolates that has a 0.99 yield value
Total.iso        = Total.iso %>% filter(Yield <= 0.9)

# Figure 3C - CAZy enzyme ----

Figure_test_1.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (CAZy),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("CAZy (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,30) + theme(legend.position="none")
Figure_test_1.ISO.new

# Figure 3D - Protein enzyme ----

Figure_test_2.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (Protein),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Protein Enzyme (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,50) + theme(legend.position="none")
Figure_test_2.ISO.new

# Figure 4D - Transport ----

Figure_test_3.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (tranport.total),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Total Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,200) + theme(legend.position="none")
Figure_test_3.ISO.new

# Figure 4E - Aminoacids Transport ----

Figure_test_4.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (amino.transport),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Aminoacid Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_4.ISO.new

# Figure 4F - Carbohydrate Transport ----

Figure_test_5.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (GH.trasnport),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Carbohydrate Transporter (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,55) + theme(legend.position="none")
Figure_test_5.ISO.new

# Figure 5E - Osmolytes ----

Figure_test_6.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (osmolyte),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Osmolytes (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,32) + theme(legend.position="none")
Figure_test_6.ISO.new

# Figure 5F - Biofilm ----

Figure_test_7.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (biofilm),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Biofilm (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,16) + theme(legend.position="none")
Figure_test_7.ISO.new

# Figure 5G - Heat Resistance ----

Figure_test_8.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (temp),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Heat Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,15) + theme(legend.position="none")
Figure_test_8.ISO.new

# Figure 5H - pH Resistance ----

Figure_test_9.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (pH),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("pH Resistance (genes)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_9.ISO.new

# Figure 6B - Minimum generation time ----

Figure_test_10.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (MGT),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Minimum generation time (hrs)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,20) + theme(legend.position="none")
Figure_test_10.ISO.new

# Figure 6C - Yield ----

Figure_test_11.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (yield),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1) + theme(legend.position="none")
Figure_test_11.ISO.new

# Figure 7B - Optimal Growth Temperature ----

Figure_test_12.ISO.new = ggplot(data = Total.iso, aes(x = genome.size/1e6, y = (OGT),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Optimal Growth Temperature (°C)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,80) + theme(legend.position="none")
Figure_test_12.ISO.new

# Cross correlations MAGs ----

# Call data

MAG_gen_trait.2  = read.csv("Intermediate_Results/MAG_gen_trait.csv",dec=".")
MAG_gen_trait.2  = MAG_gen_trait.2[, -1]
MAG_gen_trait.3  = read.csv("Intermediate_Results/MAG_gen_trait.5.csv",dec=".")
MAG_gen_trait.3  = MAG_gen_trait.3[, -1]

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

# Correlation tables----

# Fast-growing bacteria

cor_1 = as.data.frame(cor(MAG_gen_trait.2[,3:13]))
write.csv(cor_1, file = "Output_Data/Fast.growing_correlation.MAG.csv")

# Slow-growing bacteria

cor_2 = as.data.frame(cor(MAG_gen_trait.3[,3:13]))
write.csv(cor_2, file = "Output_Data/Slow.growing_correlation.MAG.csv")

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

# Figure 1 - Cross correlation MAGs ----

Figure_5 = ggpairs(MAG_gen_trait.2[,3:13], 
                   upper = list(continuous = my_fn),
                   lower = list(continuous = "smooth"))  
Figure_5

# Figure S1 - Cross correlation MAGs ----

Figure_5.a = ggpairs(MAG_gen_trait.3[,3:13], 
                     upper = list(continuous = my_fn),
                     lower = list(continuous = "smooth")) 
Figure_5.a

# Cross correlations Isolates ----

# Plotting cross correlations

ISO_gen_trait.1  = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_trait.1  = ISO_gen_trait.1[, -1]

# Figure 2 - Cross correlation Isolates ----

ISO_gen_trait.1.f = ISO_gen_trait.1 %>% filter(MGT <= 5)

Figure_6 = ggpairs(ISO_gen_trait.1.f[,3:14], 
                   upper = list(continuous = my_fn),
                   lower = list(continuous = "smooth"))  
Figure_6

# Figure S2 - Cross correlation Isolates ----

ISO_gen_trait.1.s = ISO_gen_trait.1 %>% filter(MGT > 5)

Figure_6.a = ggpairs(ISO_gen_trait.1.s[,3:14], 
                     upper = list(continuous = my_fn),
                     lower = list(continuous = "smooth"))  
Figure_6.a

# Correlation tables----

# Fast-growing bacteria

cor_1 = as.data.frame(cor(ISO_gen_trait.1.f[,3:14]))
write.csv(cor_1, file = "Output_Data/Fast.growing_correlation.isolates.csv")

# Slow-growing bacteria

cor_2 = as.data.frame(cor(ISO_gen_trait.1.s[,3:14]))
write.csv(cor_2, file = "Output_Data/Slow.growing_correlation.isolates.csv")

# CUE for different taxonomic levels ----

isolates     = read.csv("Input_Data/SOIL_ISOLATES/dement_isolates_CUE.csv",dec=".")
isolatest    = read.csv("Input_Data/SOIL_ISOLATES/total.granularity.3_datasets.csv",dec=".")
isolates.1   = subset(isolates, isolates$CUE !='NaN')
meta_iso   = read.delim("Input_Data/SOIL_ISOLATES/metadata_2.tsv",sep="\t")
meta_iso   = as.data.frame(cbind(meta_iso$IMG.Genome.ID,meta_iso$Phylum,
                                 meta_iso$Class,meta_iso$Order,meta_iso$Family,
                                 meta_iso$Genus))
meta_iso.1 = read.delim("Input_Data/SOIL_ISOLATES/metadata_1.tsv",sep="\t")
meta_iso.1 = as.data.frame(cbind(meta_iso.1$IMG.Genome.ID,meta_iso.1$Phylum,
                                 meta_iso.1$Class,meta_iso.1$Order,meta_iso.1$Family,
                                 meta_iso.1$Genus))

# final meta
final_meta   = as.data.frame(rbind(meta_iso,meta_iso.1))
colnames(final_meta) = c("id","Phylum","Class","Order","Family","Genus")

merged_df    = merge(isolates.1,final_meta,by = "id")
merged_df    = merged_df %>% filter(CUE <= 0.9)

# Figure S4A - Phylum ----

phylum = merged_df %>% group_by(Phylum) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
phylum = phylum %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                             mingentime_mean  > 5 ~ "slow"))

Figure_test_1.phylum = ggplot(data = phylum, aes(x = log(genome_length_mean), y = log(CUE_mean))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  labs(title = "Phylum")
Figure_test_1.phylum

Figure_test_1.phylum.new = ggplot(data = phylum, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                     color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Phylum") + theme(legend.position="none")
Figure_test_1.phylum.new

# Figure S4B - Class ----

Class = merged_df %>% group_by(Class) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Class = Class %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))
Figure_test_1.Class = ggplot(data = Class, aes(x = log(genome_length_mean), y = log(CUE_mean))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") + 
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  labs(title = "Class")
Figure_test_1.Class

Figure_test_1.Class.new = ggplot(data = Class, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                   color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Class") + theme(legend.position="none")
Figure_test_1.Class.new

# Figure S4C - Order ----

Order = merged_df %>% group_by(Order) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Order = Order %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Order = ggplot(data = Order, aes(x = log(genome_length_mean), y = log(CUE_mean))) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  labs(title = "Order")
Figure_test_1.Order

Figure_test_1.Order.new = ggplot(data = Order, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                   color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Order") + theme(legend.position="none")
Figure_test_1.Order.new

# Figure S4D - Family ----

Family = merged_df %>% group_by(Family) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
Family = Family %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                             mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Family = ggplot(data = Family, aes(x = log(genome_length_mean), y = log(CUE_mean))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  labs(title = "Family")
Figure_test_1.Family

Figure_test_1.Family.new = ggplot(data = Family, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                     color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Family") + theme(legend.position="none")
Figure_test_1.Family.new

# Figure S4E - Genus ----

Genus = merged_df %>% group_by(Genus) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Genus = Genus %>% mutate(speed = case_when(mingentime_mean <= 5 ~ "fast",
                                           mingentime_mean  > 5 ~ "slow"))

Figure_test_1.Genus = ggplot(data = Genus, aes(x = log(genome_length_mean), y = log(CUE_mean))) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 5), label.y = "bottom", label.x = "left") + 
  xlab("Log(Genome size)") + 
  ylab("Log(Yield)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  labs(title = "Genus")
Figure_test_1.Genus

Figure_test_1.Genus.new = ggplot(data = Genus, aes(x = genome_length_mean/1e6, y = CUE_mean,
                                                   color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Mbp") + 
  ylab("Yield (-)") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7/1e6) + ylim(0,1.5) + labs(title = "Genus") + theme(legend.position="none")
Figure_test_1.Genus.new

# Summary statistics from the MAGs and isolates ----

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

Figure_S5.A      = ggplot(data = total_data, aes(x=cat, y=completeness)) + 
  geom_boxplot() + theme_light() + theme(legend.position = "none") + 
  labs(x = "", y = "Completeness (%)") + ylim(0,100)
Figure_S5.A  

Figure_S5.B      = ggplot(data = total_data, aes(x=cat, y=contamination)) + 
  geom_boxplot() + theme_light() + theme(legend.position = "none") + 
  labs(x = "", y = "Contamination (%)")
Figure_S5.B  

# Isolates
meta_iso   = read.delim("Input_Data/SOIL_ISOLATES/metadata_1.tsv",sep="\t")
meta_iso   = meta_iso %>% select(taxon_oid,High.Quality)
meta_iso.1 = read.delim("Input_Data/SOIL_ISOLATES/metadata_2.tsv",sep="\t")
meta_iso.1 = meta_iso.1 %>% select(taxon_oid,High.Quality)
final_meta = as.data.frame(rbind(meta_iso,meta_iso.1))
colnames(final_meta) = c("id","quality")
Iso_combined  = merge(final_meta,total_genes.guild.940_ISO,by = "id")



