# Data Preparation

# Normalization of Traits by Genome size ----

# Call Trait Keys ----

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

# Trait Data Preparation (MAGs) ----

total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")

# GH - genes (MAG)----
GH_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.MAG$genome.size = as.numeric(as.character(GH_TOTAL.MAG$genome.size))
# Normalize by genome size
GH_TOTAL.MAG_n = GH_TOTAL.MAG %>% mutate(GH_total = rowSums(GH_TOTAL.MAG[,4:ncol(GH_TOTAL.MAG)])/GH_TOTAL.MAG$genome.size)
# Aggregation
GH_TOTAL_m     = GH_TOTAL.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
GH_TOTAL_m.1   = GH_TOTAL_m[-c(247), ] # Erasing super big Functional Group
GH_MAG         = GH_TOTAL_m.1 %>% select(guild,genome.size_mean,GH_total_mean)

# Protein - genes (MAG)----
PR_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(PR_rule$`microtrait_hmm-name`))
PR_TOTAL.MAG$genome.size = as.numeric(as.character(PR_TOTAL.MAG$genome.size))
# Normalize by genome size
Protein_MAG_n = PR_TOTAL.MAG %>% mutate(PR_total = rowSums(PR_TOTAL.MAG[,4:ncol(PR_TOTAL.MAG)])/PR_TOTAL.MAG$genome.size)
# Aggregation
PR_TOTAL_m   = Protein_MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                      list(mean=mean), na.rm=TRUE))
PR_TOTAL_m.1 = PR_TOTAL_m[-c(247), ] # Erasing super big Functional Group
Protein_MAG  = PR_TOTAL_m.1 %>% select(guild,genome.size_mean,PR_total_mean)

# Transport - genes (MAG)----
transp_rule  = transp_rule %>% filter(`function` == c("transporter"))
transp_rule  = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL = total_genes.guild.940_MAG %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL$genome.size = as.numeric(as.character(TRANSP_TOTAL$genome.size))
# Normalize by genome size
TRANSP_TOTAL_n = TRANSP_TOTAL %>% mutate(transp_total = rowSums(TRANSP_TOTAL[,4:ncol(TRANSP_TOTAL)])/TRANSP_TOTAL$genome.size)
# Aggregation
TRANSP_TOTAL_m   = TRANSP_TOTAL_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m.1 = TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group
Transporter_MAG  = TRANSP_TOTAL_m.1 %>% select(guild,genome.size_mean,transp_total_mean)

# Transport - Aminoacids (MAG)----
transp_rule        = transp_rule %>% filter(`function`==c("transporter"))
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL   = total_genes.guild.940_MAG %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL$genome.size))
# Normalize by genome size
AMI_TRANSP_TOTAL_n = AMI_TRANSP_TOTAL %>% mutate(transp_total = rowSums(AMI_TRANSP_TOTAL[,4:ncol(AMI_TRANSP_TOTAL)])/AMI_TRANSP_TOTAL$genome.size)
# Aggregation
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                 list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m.1 = AMI_TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group
Aminoacids_T_MAG   = AMI_TRANSP_TOTAL_m.1 %>% select(guild,genome.size_mean,transp_total_mean)

# Transport - Carbohydrate (MAG)----
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL   = total_genes.guild.940_MAG %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL$genome.size))
# Normalize by genome size
CAR_TRANSP_TOTAL_n = CAR_TRANSP_TOTAL %>% mutate(CAR_TRANSP = rowSums(CAR_TRANSP_TOTAL[,4:ncol(CAR_TRANSP_TOTAL)])/CAR_TRANSP_TOTAL$genome.size)
# Aggregation
CAR_TRANSP_TOTAL_m = CAR_TRANSP_TOTAL_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                 list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL_m.1 = CAR_TRANSP_TOTAL_m[-c(247), ] # Erasing super big Functional Group
CAR_T_MAG            = CAR_TRANSP_TOTAL_m.1 %>% select(guild,genome.size_mean,CAR_TRANSP_mean)

# Osmolytes - genes (MAG)----
OSMO_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.MAG$genome.size = as.numeric(as.character(OSMO_TOTAL.MAG$genome.size))
# Normalize by genome size
OSMO_TOTAL.MAG_n = OSMO_TOTAL.MAG %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.MAG[,4:ncol(OSMO_TOTAL.MAG)])/OSMO_TOTAL.MAG$genome.size)
# Aggregation
OSMO_TOTAL.MAG_m   = OSMO_TOTAL.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
OSMO_TOTAL.MAG_m.1 = OSMO_TOTAL.MAG_m[-c(247), ] # Erasing super big Functional Group
Osmolyte_MAG       = OSMO_TOTAL.MAG_m.1 %>% select(guild,genome.size_mean,OSMO_total_mean)

# Biofilm - genes (MAG)----
BIO_TOTAL.MAG = total_genes.guild.940_MAG %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
# Normalize by genome size
BIO_TOTAL.MAG_n = BIO_TOTAL.MAG %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG[,4:ncol(BIO_TOTAL.MAG)])/BIO_TOTAL.MAG$genome.size)
# Aggregation
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
BIO_TOTAL.MAG_m.1 = BIO_TOTAL.MAG_m[-c(247), ] # Erasing super big Functional Group
Biofilm_MAG       = BIO_TOTAL.MAG_m.1 %>% select(guild,genome.size_mean,BIO_total_mean)

# High Temp - genes (MAG)----
TEMP.MAG = total_genes.guild.940_MAG %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP.MAG$genome.size = as.numeric(as.character(TEMP.MAG$genome.size))
# Normalize by genome size
TEMP.MAG_n   = TEMP.MAG %>% mutate(TEMP_total = rowSums(TEMP.MAG[,4:ncol(TEMP.MAG)])/TEMP.MAG$genome.size)
# Aggregation
TEMP.MAG_m   = TEMP.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
TEMP.MAG_m.1 = TEMP.MAG_m[-c(247), ] # Erasing super big Functional Group
TEMP_MAG     = TEMP.MAG_m.1 %>% select(guild,genome.size_mean,TEMP_total_mean)

# pH - genes (MAG)----
PH.MAG = total_genes.guild.940_MAG %>% select(any_of(pH_rule$`microtrait_hmm-name`))
PH.MAG$genome.size = as.numeric(as.character(PH.MAG$genome.size))
# Normalize by genome size
PH.MAG_n   = PH.MAG %>% mutate(PH_total = rowSums(PH.MAG[,4:ncol(PH.MAG)])/PH.MAG$genome.size)
# Aggregation
PH.MAG_m   = PH.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                               list(mean=mean), na.rm=TRUE))
PH.MAG_m.1 = PH.MAG_m[-c(247), ] # Erasing super big Functional Group
PH_MAG     = PH.MAG_m.1 %>% select(guild,genome.size_mean,PH_total_mean)

# Combine files----

MAG_gen_trait_norm    = as.data.frame(cbind(Aminoacids_T_MAG,PH_MAG$PH_total_mean,
                                            TEMP_MAG$TEMP_total_mean,Biofilm_MAG$BIO_total_mean,
                                            Osmolyte_MAG$OSMO_total_mean,CAR_T_MAG$CAR_TRANSP_mean,
                                            Transporter_MAG$transp_total_mean,Protein_MAG$PR_total_mean,
                                            GH_MAG$GH_total_mean))
colnames(MAG_gen_trait_norm) = c("guild","genome.size","amino.transport","pH","temp",
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
MAG_gen_trait.1     = merge(MAG_gen_trait_norm,Y_MAG.C.1,by = c("guild","genome.size"))
MAG_gen_trait.2     = merge(MAG_gen_trait.1,Y_MAG.B.1,by = c("guild","genome.size"))
MAG_gen_trait.3     = merge(MAG_gen_trait.1,Y_MAG.A.1,by = c("guild","genome.size"))

write.csv(MAG_gen_trait.2, file = "Intermediate_Results/MAG_gen_trait_norm.csv")
write.csv(MAG_gen_trait.3, file = "Intermediate_Results/MAG_gen_trait.5_norm.csv")

total_MAGs = as.data.frame(rbind(MAG_gen_trait.2,MAG_gen_trait.3))
write.csv(total_MAGs, file = "Intermediate_Results/total_MAGs_norm.csv")

# Trait Data Preparation (Isolates) ----

total_genes.guild.940_ISO = read.csv("Intermediate_Results/total_genes.guild.940_ISO.csv",dec=".")

# GH - genes (Isolates) ----
GH_TOTAL.ISO = total_genes.guild.940_ISO %>% select(any_of(GH_rule$`microtrait_hmm-name`))
GH_TOTAL.ISO$genome.size = as.numeric(as.character(GH_TOTAL.ISO$genome.size))
# Normalize by genome size
GH_TOTAL.ISO_n = GH_TOTAL.ISO %>% mutate(GH_total = rowSums(GH_TOTAL.ISO[,4:ncol(GH_TOTAL.ISO)])/GH_TOTAL.ISO$genome.size)
# Aggregation
GH_TOTAL_m   = GH_TOTAL.ISO_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                       list(mean=mean), na.rm=TRUE))
GH_Isolates  = GH_TOTAL_m %>% select(guild,genome.size_mean,GH_total_mean)

# Protein - genes (Isolates) ----
PH_TOTAL.ISO = total_genes.guild.940_ISO %>% select(any_of(PR_rule$`microtrait_hmm-name`))
PH_TOTAL.ISO$genome.size = as.numeric(as.character(PH_TOTAL.ISO$genome.size))
# Normalize by genome size
PH_TOTAL.ISO_n = PH_TOTAL.ISO %>% mutate(PR_total = rowSums(PH_TOTAL.ISO[,4:ncol(PH_TOTAL.ISO)])/PH_TOTAL.ISO$genome.size)
# Aggregation
PH_TOTAL_m   = PH_TOTAL.ISO_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                       list(mean=mean), na.rm=TRUE))
Protein_Isolates      = PH_TOTAL_m %>% select(guild,genome.size_mean,PR_total_mean)

# Transport - genes (Isolates)----
transp_rule    = transp_rule %>% filter(`function` == c("transporter"))
transp_rule    = as.data.frame(rbind("id","guild","genome.size",transp_rule))
TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule$`microtrait_hmm-name`))
TRANSP_TOTAL.i$genome.size = as.numeric(as.character(TRANSP_TOTAL.i$genome.size))
# Normalize by genome size
TRANSP_TOTAL_n = TRANSP_TOTAL.i %>% mutate(transp_total = rowSums(TRANSP_TOTAL.i[,4:ncol(TRANSP_TOTAL.i)])/TRANSP_TOTAL.i$genome.size)
# Aggregation
TRANSP_TOTAL_m = TRANSP_TOTAL_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
Transporter_Isolates  = TRANSP_TOTAL_m %>% select(guild,genome.size_mean,transp_total_mean)

# Transport - Aminoacids (Isolates)----
transp_rule_ami    = transp_rule %>% filter(class==c("aminoacid","peptide"))
transp_rule_ami    = as.data.frame(rbind("id","guild","genome.size",transp_rule_ami))
AMI_TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule_ami$`microtrait_hmm-name`))
AMI_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(AMI_TRANSP_TOTAL.i$genome.size))
# Normalize by genome size
AMI_TRANSP_TOTAL_n = AMI_TRANSP_TOTAL.i %>% mutate(transp_total = rowSums(AMI_TRANSP_TOTAL.i[,4:ncol(AMI_TRANSP_TOTAL.i)])/AMI_TRANSP_TOTAL.i$genome.size)
# Aggregation
AMI_TRANSP_TOTAL.i_m = AMI_TRANSP_TOTAL_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
Aminoacids_T_Isolates = AMI_TRANSP_TOTAL.i_m %>% select(guild,genome.size_mean,transp_total_mean)

# Transport - Carbohydrate (Isolates)----
transp_rule_car    = transp_rule %>% filter(class==c("carbohydrate"))
transp_rule_car    = as.data.frame(rbind("id","guild","genome.size",transp_rule_car))
CAR_TRANSP_TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(transp_rule_car$`microtrait_hmm-name`))
CAR_TRANSP_TOTAL.i$genome.size = as.numeric(as.character(CAR_TRANSP_TOTAL.i$genome.size))
# Normalize by genome size
CAR_TRANSP_TOTAL_n   = CAR_TRANSP_TOTAL.i %>% mutate(transp_total = rowSums(CAR_TRANSP_TOTAL.i[,4:ncol(CAR_TRANSP_TOTAL.i)])/CAR_TRANSP_TOTAL.i$genome.size)
# Aggregation
CAR_TRANSP_TOTAL.i_m = CAR_TRANSP_TOTAL_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                   list(mean=mean), na.rm=TRUE))
CAR_TRANSP_TOTAL_Isolates     = CAR_TRANSP_TOTAL.i_m %>% select(guild,genome.size_mean,transp_total_mean)

# Osmolytes - genes (Isolates)----
OSMO_TOTAL.TOTAL.i = total_genes.guild.940_ISO %>% select(any_of(osmo_rule$`microtrait_hmm-name`))
OSMO_TOTAL.TOTAL.i$genome.size = as.numeric(as.character(OSMO_TOTAL.TOTAL.i$genome.size))
# Normalize by genome size
OSMO_TOTAL.TOTAL_n     = OSMO_TOTAL.TOTAL.i %>% mutate(OSMO_total = rowSums(OSMO_TOTAL.TOTAL.i[,4:ncol(OSMO_TOTAL.TOTAL.i)])/OSMO_TOTAL.TOTAL.i$genome.size)
# Aggregation
OSMO_TOTAL.TOTAL.i_m   = OSMO_TOTAL.TOTAL_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                     list(mean=mean), na.rm=TRUE))
Osmolyte_Isolates = OSMO_TOTAL.TOTAL.i_m %>% select(guild,genome.size_mean,OSMO_total_mean)

# Biofilm - genes (Isolates)----
BIO_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(biofilm_rule$`microtrait_hmm-name`))
BIO_TOTAL.MAG$genome.size = as.numeric(as.character(BIO_TOTAL.MAG$genome.size))
# Normalize by genome size
BIO_TOTAL.MAG_n           = BIO_TOTAL.MAG %>% mutate(BIO_total = rowSums(BIO_TOTAL.MAG[,4:ncol(BIO_TOTAL.MAG)])/BIO_TOTAL.MAG$genome.size)
# Aggregation
BIO_TOTAL.MAG_m   = BIO_TOTAL.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                             list(mean=mean), na.rm=TRUE))
Biofilm_Isolates  = BIO_TOTAL.MAG_m %>% select(guild,genome.size_mean,BIO_total_mean)

# High Temperature - genes (Isolates)----
TEMP_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(high.T_rule$`microtrait_hmm-name`))
TEMP_TOTAL.MAG$genome.size = as.numeric(as.character(TEMP_TOTAL.MAG$genome.size))
# Normalize by genome size
TEMP_TOTAL.MAG_n           = TEMP_TOTAL.MAG %>% mutate(TEMP_total = rowSums(TEMP_TOTAL.MAG[,4:ncol(TEMP_TOTAL.MAG)])/TEMP_TOTAL.MAG$genome.size)
# Aggregation
TEMP_TOTAL.MAG_m   = TEMP_TOTAL.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
TEMP_Isolates      = TEMP_TOTAL.MAG_m %>% select(guild,genome.size_mean,TEMP_total_mean)

# pH - genes (Isolates)----
pH_TOTAL.MAG = total_genes.guild.940_ISO %>% select(any_of(pH_rule$`microtrait_hmm-name`))
pH_TOTAL.MAG$genome.size = as.numeric(as.character(pH_TOTAL.MAG$genome.size))
# Normalize by genome size
pH_TOTAL.MAG_n           = pH_TOTAL.MAG %>% mutate(PH_total = rowSums(pH_TOTAL.MAG[,4:ncol(pH_TOTAL.MAG)])/pH_TOTAL.MAG$genome.size)
# Aggregation
pH_TOTAL.MAG_m   = pH_TOTAL.MAG_n %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                           list(mean=mean), na.rm=TRUE))
PH_Isolates      = pH_TOTAL.MAG_m %>% select(guild,genome.size_mean,PH_total_mean)

# Combine files----

Isolates_gen_trait_norm    = as.data.frame(cbind(Aminoacids_T_Isolates,PH_Isolates$PH_total_mean,
                                                 TEMP_Isolates$TEMP_total_mean,Biofilm_Isolates$BIO_total_mean,
                                                 Osmolyte_Isolates$OSMO_total_mean,CAR_TRANSP_TOTAL_Isolates$transp_total_mean,
                                                 Transporter_Isolates$transp_total_mean,Protein_Isolates$PR_total_mean,
                                                 GH_Isolates$GH_total_mean))
colnames(Isolates_gen_trait_norm) = c("guild","genome.size","amino.transport","pH","temp",
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
ISO_gen_trait.1        = merge(Isolates_gen_trait_norm,isolates.2.m,by = c("guild","genome.size"))
colnames(ISO_gen_trait.1) = c("Guild","Genome Size","Amino-transporter","pH-Tol",
                              "Temp-Tol", "Biofilm", "Osmolyte","GH-transporter",
                              "Total-transporter","Protein-enzyme","CAZy","Yield","MGT","OGT")
col_order = c("Guild", "Genome Size","Yield","OGT","MGT","Amino-transporter","GH-transporter",
              "Total-transporter","Protein-enzyme","CAZy","pH-Tol","Temp-Tol", 
              "Biofilm", "Osmolyte")
ISO_gen_trait.1   = ISO_gen_trait.1[, col_order]
ISO_gen_trait.1   = ISO_gen_trait.1 %>% filter(Yield <= 0.9)

write.csv(ISO_gen_trait.1, file = "Intermediate_Results/ISO_gen_trait.norm.csv")

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

#-------------------------------------------------------------------------------

# New Figures ----

# PCAS ----

# MAGs gene counts ----
MAG_gen_TOTAL_ab = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")

res.pca_AB      = prcomp(MAG_gen_TOTAL_ab[,3:14], scale = TRUE)
summary(res.pca_AB)
fviz_eig(res.pca_AB)

fviz_pca_var(res.pca_AB,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# MAGs gene counts + GC count ----
MAG_gen_TOTAL_GC = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")
MAG_gen_TOTAL_GC = as.data.frame(MAG_gen_TOTAL_GC[,-(15:17)])

res.pca_AB_CG    = prcomp(MAG_gen_TOTAL_GC[,3:15], scale = TRUE)
summary(res.pca_AB_CG)
fviz_eig(res.pca_AB_CG)

fviz_pca_var(res.pca_AB_CG,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# MAGs gene counts normalized ----
MAG_gen_TOTAL_no = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC_norm.csv",dec=".")

res.pca_norm      = prcomp(MAG_gen_TOTAL_no[,3:14], scale = TRUE)
summary(res.pca_norm)
fviz_eig(res.pca_norm)

fviz_pca_var(res.pca_norm,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# MAGs gene counts + GC count normalized ----
MAG_gen_TOTAL_no = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC_norm.csv",dec=".")

MAG_gen_TOTAL_no = as.data.frame(MAG_gen_TOTAL_no[,-(15:17)])

res.pca_no_CG    = prcomp(MAG_gen_TOTAL_no[,3:15], scale = TRUE)
summary(res.pca_no_CG)
fviz_eig(res.pca_no_CG)

fviz_pca_var(res.pca_no_CG,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Isolates gene counts ----
ISO_gen_trait.1  = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")

# PCA

res.pca_ISO      = prcomp(ISO_gen_trait.1[,2:ncol(ISO_gen_trait.1)], scale = TRUE)
summary(res.pca_ISO)
fviz_eig(res.pca_ISO)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Isolates gene counts normalized ----
ISO_gen_trait.1  = read.csv("Intermediate_Results/ISO_gen_trait.norm.csv",dec=".")

# PCA

res.pca_ISO_norm = prcomp(ISO_gen_trait.1[,3:ncol(ISO_gen_trait.1)], scale = TRUE)
summary(res.pca_ISO_norm)
fviz_eig(res.pca_ISO_norm)

fviz_pca_var(res.pca_ISO_norm,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# New correlation plots ----

# MAGs ----

MAG_gen_trait.2  = read.csv("Intermediate_Results/MAG_gen_trait_norm.csv",dec=".")
MAG_gen_trait.2  = MAG_gen_trait.2[, -1]
MAG_gen_trait.3  = read.csv("Intermediate_Results/MAG_gen_trait.5_norm.csv",dec=".")
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
write.csv(cor_1, file = "Output_Data/Fast.growing_correlation.MAG_norm.csv")

# Slow-growing bacteria

cor_2 = as.data.frame(cor(MAG_gen_trait.3[,3:13]))
write.csv(cor_2, file = "Output_Data/Slow.growing_correlation.MAG_norm.csv")

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

Figure_5_norm = ggpairs(MAG_gen_trait.2[,3:13], 
                        upper = list(continuous = my_fn),
                        lower = list(continuous = "smooth"))  
Figure_5_norm

# Figure S1 - Cross correlation MAGs ----

Figure_5.a_norm = ggpairs(MAG_gen_trait.3[,3:13], 
                          upper = list(continuous = my_fn),
                          lower = list(continuous = "smooth")) 
Figure_5.a_norm

# Isolates ----

ISO_gen_trait.1  = read.csv("Intermediate_Results/ISO_gen_trait.norm.csv",dec=".")
ISO_gen_trait.1  = ISO_gen_trait.1[, -1]

# Figure 2 - Cross correlation Isolates ----

ISO_gen_trait.1.f = ISO_gen_trait.1 %>% filter(MGT <= 5)

Figure_6_norm = ggpairs(ISO_gen_trait.1.f[,3:14], 
                        upper = list(continuous = my_fn),
                        lower = list(continuous = "smooth"))  
Figure_6_norm

# Figure S2 - Cross correlation Isolates ----

ISO_gen_trait.1.s = ISO_gen_trait.1 %>% filter(MGT > 5)

Figure_6.a_norm = ggpairs(ISO_gen_trait.1.s[,3:14], 
                          upper = list(continuous = my_fn),
                          lower = list(continuous = "smooth"))  
Figure_6.a_norm

# Correlation tables----

# Fast-growing bacteria

cor_1 = as.data.frame(cor(ISO_gen_trait.1.f[,3:14]))
write.csv(cor_1, file = "Output_Data/Fast.growing_correlation.isolates_norm.csv")

# Slow-growing bacteria

cor_2 = as.data.frame(cor(ISO_gen_trait.1.s[,3:14]))
write.csv(cor_2, file = "Output_Data/Slow.growing_correlation.isolates_norm.csv")

# Linear Regression Plots ----

# MAGs gene counts ----

MAG_gen_TOTAL_ab = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")
MAG_gen_TOTAL_ab = MAG_gen_TOTAL_ab %>% mutate(speed = case_when(mgr <= 5 ~ "fast",
                                                                 mgr  > 5 ~ "slow"))
# CAZy enzyme ----

# Fast
MAG_gen_TOTAL_ab_fast = MAG_gen_TOTAL_ab %>% filter(speed == "fast")

cazy_fast_basic = lm(CAZy ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(cazy_fast_basic)

cazy_fast_1     = lm(CAZy ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(cazy_fast_1)

cazy_fast_2     = lm(CAZy ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(cazy_fast_2)

cazy_fast_3     = lm(CAZy ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(cazy_fast_3)

cazy_power.law  = lm(log(CAZy) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(cazy_power.law)

AIC(cazy_fast_basic, cazy_fast_1, cazy_fast_2, cazy_fast_3, cazy_power.law)

# Slow
MAG_gen_TOTAL_ab_slow = MAG_gen_TOTAL_ab %>% filter(speed == "slow")

cazy_slow_basic = lm(CAZy ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(cazy_slow_basic)

cazy_slow_1     = lm(CAZy ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(cazy_slow_1)

cazy_slow_2     = lm(CAZy ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(cazy_slow_2)

cazy_slow_3     = lm(CAZy ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(cazy_slow_3)

cazy_power.law.slow  = lm(log(CAZy+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(cazy_power.law.slow)

AIC(cazy_slow_basic, cazy_slow_1, cazy_slow_2, cazy_slow_3, cazy_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$CAZy), 
     xlab="log(Genome size)", ylab="log(CAZy)", pch=18, col="black")
abline(cazy_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$CAZy), 
     xlab="log(Genome size)", ylab="log(CAZy)", pch=18, col="black")
abline(cazy_power.law.slow, col="gray", lty=6, lwd=2)

# Protein enzyme ----

# Fast

Protein_fast_basic = lm(Protein ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(Protein_fast_basic)

Protein_fast_1     = lm(Protein ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(Protein_fast_1)

Protein_fast_2     = lm(Protein ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(Protein_fast_2)

Protein_fast_3     = lm(Protein ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(Protein_fast_3)

Protein_power.law  = lm(log(Protein) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(Protein_power.law)

AIC(Protein_fast_basic, Protein_fast_1, Protein_fast_2, Protein_fast_3, Protein_power.law)

# Slow

Protein_slow_basic = lm(Protein ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(Protein_slow_basic)

Protein_slow_1     = lm(Protein ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(Protein_slow_1)

Protein_slow_2     = lm(Protein ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(Protein_slow_2)

Protein_slow_3     = lm(Protein ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(Protein_slow_3)

Protein_power.law.slow  = lm(log(Protein) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(Protein_power.law.slow)

AIC(Protein_slow_basic, Protein_slow_1, Protein_slow_2, Protein_slow_3, Protein_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$Protein), 
     xlab="log(Genome size)", ylab="log(Protein)", pch=18, col="black")
abline(Protein_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$Protein), 
     xlab="log(Genome size)", ylab="log(Protein)", pch=18, col="black")
abline(Protein_power.law.slow, col="gray", lty=6, lwd=2)

# Transport total Transporters ----

# Fast

tranport.total_fast_basic = lm(tranport.total ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(tranport.total_fast_basic)

tranport.total_fast_1     = lm(tranport.total ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(tranport.total_fast_1)

tranport.total_fast_2     = lm(tranport.total ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(tranport.total_fast_2)

tranport.total_fast_3     = lm(tranport.total ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(tranport.total_fast_3)

tranport.total_power.law  = lm(log(tranport.total) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(tranport.total_power.law)

AIC(tranport.total_fast_basic, tranport.total_fast_1, tranport.total_fast_2, tranport.total_fast_3, tranport.total_power.law)

# Slow

tranport.total_slow_basic = lm(tranport.total ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(tranport.total_slow_basic)

tranport.total_slow_1     = lm(tranport.total ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(tranport.total_slow_1)

tranport.total_slow_2     = lm(tranport.total ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(tranport.total_slow_2)

tranport.total_slow_3     = lm(tranport.total ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(tranport.total_slow_3)

tranport.total_power.law.slow  = lm(log(tranport.total+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(tranport.total_power.law.slow)

AIC(tranport.total_slow_basic, tranport.total_slow_1, tranport.total_slow_2, tranport.total_slow_3, tranport.total_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$tranport.total), 
     xlab="log(Genome size)", ylab="log(Total transporters)", pch=18, col="black")
abline(tranport.total_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$tranport.total), 
     xlab="log(Genome size)", ylab="log(Total transporters)", pch=18, col="black")
abline(tranport.total_power.law.slow, col="gray", lty=6, lwd=2)

# GH total Transporters ----

# Fast

GH.trasnport_fast_basic = lm(GH.trasnport ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(GH.trasnport_fast_basic)

GH.trasnport_fast_1     = lm(GH.trasnport ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(GH.trasnport_fast_1)

GH.trasnport_fast_2     = lm(GH.trasnport ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(GH.trasnport_fast_2)

GH.trasnport_fast_3     = lm(GH.trasnport ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(GH.trasnport_fast_3)

GH.trasnport_power.law  = lm(log(GH.trasnport+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(GH.trasnport_power.law)

AIC(GH.trasnport_fast_basic, GH.trasnport_fast_1, GH.trasnport_fast_2, GH.trasnport_fast_3, GH.trasnport_power.law)

# Slow

GH.trasnport_slow_basic = lm(GH.trasnport ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(GH.trasnport_slow_basic)

GH.trasnport_slow_1     = lm(GH.trasnport ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(GH.trasnport_slow_1)

GH.trasnport_slow_2     = lm(GH.trasnport ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(GH.trasnport_slow_2)

GH.trasnport_slow_3     = lm(GH.trasnport ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(GH.trasnport_slow_3)

GH.trasnport_power.law.slow  = lm(log(GH.trasnport+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(GH.trasnport_power.law.slow)

AIC(GH.trasnport_slow_basic, GH.trasnport_slow_1, GH.trasnport_slow_2, GH.trasnport_slow_3, GH.trasnport_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$GH.trasnport), 
     xlab="log(Genome size)", ylab="log(GH transporters)", pch=18, col="black")
abline(GH.trasnport_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$GH.trasnport), 
     xlab="log(Genome size)", ylab="log(GH transporters)", pch=18, col="black")
abline(GH.trasnport_power.law.slow, col="gray", lty=6, lwd=2)

# Amino total Transporters ----

# Fast

amino.transport_fast_basic = lm(amino.transport ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(amino.transport_fast_basic)

amino.transport_fast_1     = lm(amino.transport ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(amino.transport_fast_1)

amino.transport_fast_2     = lm(amino.transport ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(amino.transport_fast_2)

amino.transport_fast_3     = lm(amino.transport ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(amino.transport_fast_3)

amino.transport_power.law  = lm(log(amino.transport+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(amino.transport_power.law)

AIC(amino.transport_fast_basic, amino.transport_fast_1, amino.transport_fast_2, amino.transport_fast_3, amino.transport_power.law)

# Slow

amino.transport_slow_basic = lm(amino.transport ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(amino.transport_slow_basic)

amino.transport_slow_1     = lm(amino.transport ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(amino.transport_slow_1)

amino.transport_slow_2     = lm(amino.transport ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(amino.transport_slow_2)

amino.transport_slow_3     = lm(amino.transport ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(amino.transport_slow_3)

amino.transport_power.law.slow  = lm(log(amino.transport+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(amino.transport_power.law.slow)

AIC(amino.transport_slow_basic, amino.transport_slow_1, amino.transport_slow_2, amino.transport_slow_3, amino.transport_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$amino.transport), 
     xlab="log(Genome size)", ylab="log(Amino transporters)", pch=18, col="black")
abline(amino.transport_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$amino.transport), 
     xlab="log(Genome size)", ylab="log(Amino transporters)", pch=18, col="black")
abline(amino.transport_power.law.slow, col="gray", lty=6, lwd=2)

# Biofilm Transporters ----

# Fast

biofilm_fast_basic = lm(biofilm ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(biofilm_fast_basic)

biofilm_fast_1     = lm(biofilm ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(biofilm_fast_1)

biofilm_fast_2     = lm(biofilm ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(biofilm_fast_2)

biofilm_fast_3     = lm(biofilm ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(biofilm_fast_3)

biofilm_power.law  = lm(log(biofilm+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(biofilm_power.law)

AIC(biofilm_fast_basic, biofilm_fast_1, biofilm_fast_2, biofilm_fast_3, biofilm_power.law)

# Slow

biofilm_slow_basic = lm(biofilm ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(biofilm_slow_basic)

biofilm_slow_1     = lm(biofilm ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(biofilm_slow_1)

biofilm_slow_2     = lm(biofilm ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(biofilm_slow_2)

biofilm_slow_3     = lm(biofilm ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(biofilm_slow_3)

biofilm_power.law.slow  = lm(log(biofilm+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(biofilm_power.law.slow)

AIC(biofilm_slow_basic, biofilm_slow_1, biofilm_slow_2, biofilm_slow_3, biofilm_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$biofilm), 
     xlab="log(Genome size)", ylab="log(Biofilm)", pch=18, col="black")
abline(biofilm_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$biofilm), 
     xlab="log(Genome size)", ylab="log(Biofilm)", pch=18, col="black")
abline(biofilm_power.law.slow, col="gray", lty=6, lwd=2)

# Osmolyte Transporters ----

# Fast

osmolyte_fast_basic = lm(osmolyte ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(osmolyte_fast_basic)

osmolyte_fast_1     = lm(osmolyte ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(osmolyte_fast_1)

osmolyte_fast_2     = lm(osmolyte ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(osmolyte_fast_2)

osmolyte_fast_3     = lm(osmolyte ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(osmolyte_fast_3)

osmolyte_power.law  = lm(log(osmolyte+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(osmolyte_power.law)

AIC(osmolyte_fast_basic, osmolyte_fast_1, osmolyte_fast_2, osmolyte_fast_3, osmolyte_power.law)

# Slow

osmolyte_slow_basic = lm(osmolyte ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(osmolyte_slow_basic)

osmolyte_slow_1     = lm(osmolyte ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(osmolyte_slow_1)

osmolyte_slow_2     = lm(osmolyte ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(osmolyte_slow_2)

osmolyte_slow_3     = lm(osmolyte ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(osmolyte_slow_3)

osmolyte_power.law.slow  = lm(log(osmolyte+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(osmolyte_power.law.slow)

AIC(osmolyte_slow_basic, osmolyte_slow_1, osmolyte_slow_2, osmolyte_slow_3, osmolyte_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$osmolyte), 
     xlab="log(Genome size)", ylab="log(Osmolyte)", pch=18, col="black")
abline(osmolyte_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$osmolyte), 
     xlab="log(Genome size)", ylab="log(Osmolyte)", pch=18, col="black")
abline(osmolyte_power.law.slow, col="gray", lty=6, lwd=2)

# Heat Resistance Transporters ----

# Fast

temp_fast_basic = lm(temp ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(temp_fast_basic)

temp_fast_1     = lm(temp ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(temp_fast_1)

temp_fast_2     = lm(temp ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(temp_fast_2)

temp_fast_3     = lm(temp ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(temp_fast_3)

temp_power.law  = lm(log(temp) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(temp_power.law)

AIC(temp_fast_basic, temp_fast_1, temp_fast_2, temp_fast_3, temp_power.law)

# Slow

temp_slow_basic = lm(temp ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(temp_slow_basic)

temp_slow_1     = lm(temp ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(temp_slow_1)

temp_slow_2     = lm(temp ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(temp_slow_2)

temp_slow_3     = lm(temp ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(temp_slow_3)

temp_power.law.slow  = lm(log(temp) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(temp_power.law.slow)

AIC(temp_slow_basic, temp_slow_1, temp_slow_2, temp_slow_3, temp_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$temp), 
     xlab="log(Genome size)", ylab="log(Heat Resistance)", pch=18, col="black")
abline(temp_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$temp), 
     xlab="log(Genome size)", ylab="log(Heat Resistance)", pch=18, col="black")
abline(temp_power.law.slow, col="gray", lty=6, lwd=2)

# pH Resistance Transporters ----

# Fast

pH_fast_basic = lm(pH ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(pH_fast_basic)

pH_fast_1     = lm(pH ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(pH_fast_1)

pH_fast_2     = lm(pH ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(pH_fast_2)

pH_fast_3     = lm(pH ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(pH_fast_3)

pH_power.law  = lm(log(pH+1) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(pH_power.law)

AIC(pH_fast_basic, pH_fast_1, pH_fast_2, pH_fast_3, pH_power.law)

# Slow

pH_slow_basic = lm(pH ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(pH_slow_basic)

pH_slow_1     = lm(pH ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(pH_slow_1)

pH_slow_2     = lm(pH ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(pH_slow_2)

pH_slow_3     = lm(pH ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(pH_slow_3)

pH_power.law.slow  = lm(log(pH) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(pH_power.law.slow)

AIC(pH_slow_basic, pH_slow_1, pH_slow_2, pH_slow_3, pH_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$pH), 
     xlab="log(Genome size)", ylab="log(pH Resistance)", pch=18, col="black")
abline(pH_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$pH), 
     xlab="log(Genome size)", ylab="log(pH Resistance)", pch=18, col="black")
abline(pH_power.law.slow, col="gray", lty=6, lwd=2)

# Optimum Growth Temperature ----

# Fast

ogt_fast_basic = lm(ogt ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(ogt_fast_basic)

ogt_fast_1     = lm(ogt ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(ogt_fast_1)

ogt_fast_2     = lm(ogt ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(ogt_fast_2)

ogt_fast_3     = lm(ogt ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(ogt_fast_3)

ogt_power.law  = lm(log(ogt) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(ogt_power.law)

AIC(ogt_fast_basic, ogt_fast_1, ogt_fast_2, ogt_fast_3, ogt_power.law)

# Slow

ogt_slow_basic = lm(ogt ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(ogt_slow_basic)

ogt_slow_1     = lm(ogt ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(ogt_slow_1)

ogt_slow_2     = lm(ogt ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(ogt_slow_2)

ogt_slow_3     = lm(ogt ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(ogt_slow_3)

ogt_power.law.slow  = lm(log(ogt) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(ogt_power.law.slow)

AIC(ogt_slow_basic, ogt_slow_1, ogt_slow_2, ogt_slow_3, ogt_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$ogt), 
     xlab="log(Genome size)", ylab="log(Optimum Growth Temperature)", pch=18, col="black")
abline(ogt_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$ogt), 
     xlab="log(Genome size)", ylab="log(Optimum Growth Temperature)", pch=18, col="black")
abline(ogt_power.law.slow, col="gray", lty=6, lwd=2)

# Minimum Generation Time ----

# Fast

mgr_fast_basic = lm(mgr ~ genome.size, data = MAG_gen_TOTAL_ab_fast)
summary(mgr_fast_basic)

mgr_fast_1     = lm(mgr ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_fast)
summary(mgr_fast_1)

mgr_fast_2     = lm(mgr ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_fast)
summary(mgr_fast_2)

mgr_fast_3     = lm(mgr ~ GC_count_mean, data = MAG_gen_TOTAL_ab_fast)
summary(mgr_fast_3)

mgr_power.law  = lm(log(mgr) ~ log(genome.size), data = MAG_gen_TOTAL_ab_fast)
summary(mgr_power.law)

AIC(mgr_fast_basic, mgr_fast_1, mgr_fast_2, mgr_fast_3, mgr_power.law)

# Slow

mgr_slow_basic = lm(mgr ~ genome.size, data = MAG_gen_TOTAL_ab_slow)
summary(mgr_slow_basic)

mgr_slow_1     = lm(mgr ~ genome.size*completeness_mean, data = MAG_gen_TOTAL_ab_slow)
summary(mgr_slow_1)

mgr_slow_2     = lm(mgr ~ genome.size*contamination_mean, data = MAG_gen_TOTAL_ab_slow)
summary(mgr_slow_2)

mgr_slow_3     = lm(mgr ~ GC_count_mean, data = MAG_gen_TOTAL_ab_slow)
summary(mgr_slow_3)

mgr_power.law.slow  = lm(log(mgr) ~ log(genome.size), data = MAG_gen_TOTAL_ab_slow)
summary(mgr_power.law.slow)

AIC(mgr_slow_basic, mgr_slow_1, mgr_slow_2, mgr_slow_3, mgr_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(MAG_gen_TOTAL_ab_fast$genome.size), log(MAG_gen_TOTAL_ab_fast$mgr), 
     xlab="log(Genome size)", ylab="log(Minimum Generation Time)", pch=18, col="black")
abline(mgr_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(MAG_gen_TOTAL_ab_slow$genome.size), log(MAG_gen_TOTAL_ab_slow$mgr), 
     xlab="log(Genome size)", ylab="log(Minimum Generation Time)", pch=18, col="black")
abline(mgr_power.law.slow, col="gray", lty=6, lwd=2)

# Isolate gene counts ----

ISO_gen_TOTAL_ab = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_TOTAL_ab = ISO_gen_TOTAL_ab %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                                                 MGT  > 5 ~ "slow"))
# CAZy ----

# Fast

ISO_gen_TOTAL_ab_fast = ISO_gen_TOTAL_ab %>% filter(speed == "fast")

CAZy_fast_basic = lm(CAZy ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(CAZy_fast_basic)

CAZy_power.law  = lm(log(CAZy) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(CAZy_power.law)

AIC(CAZy_fast_basic, CAZy_power.law)

# Slow

ISO_gen_TOTAL_ab_slow = ISO_gen_TOTAL_ab %>% filter(speed == "slow")

CAZy_slow_basic = lm(CAZy ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(CAZy_slow_basic)

CAZy_power.law.slow  = lm(log(CAZy) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(CAZy_power.law.slow)

AIC(CAZy_slow_basic, CAZy_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$CAZy), 
     xlab="log(Genome size)", ylab="log(CAZy)", pch=18, col="black")
abline(CAZy_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$CAZy), 
     xlab="log(Genome size)", ylab="log(CAZy)", pch=18, col="black")
abline(CAZy_power.law.slow, col="gray", lty=6, lwd=2)

# Protein ----

# Fast

Protein.enzyme_fast_basic = lm(Protein.enzyme ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(Protein.enzyme_fast_basic)

Protein.enzyme_power.law  = lm(log(Protein.enzyme) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(Protein.enzyme_power.law)

AIC(Protein.enzyme_fast_basic, Protein.enzyme_power.law)

# Slow

Protein.enzyme_slow_basic = lm(Protein.enzyme ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(Protein.enzyme_slow_basic)

Protein.enzyme_power.law.slow  = lm(log(Protein.enzyme) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(Protein.enzyme_power.law.slow)

AIC(Protein.enzyme_slow_basic, Protein.enzyme_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$Protein.enzyme), 
     xlab="log(Genome size)", ylab="log(Protein)", pch=18, col="black")
abline(Protein.enzyme_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$Protein.enzyme), 
     xlab="log(Genome size)", ylab="log(Protein)", pch=18, col="black")
abline(Protein.enzyme_power.law.slow, col="gray", lty=6, lwd=2)

# Total transporter ----

# Fast

Total.transporter_fast_basic = lm(Total.transporter ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(Total.transporter_fast_basic)

Total.transporter_power.law  = lm(log(Total.transporter) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(Total.transporter_power.law)

AIC(Total.transporter_fast_basic, Total.transporter_power.law)

# Slow

Total.transporter_slow_basic = lm(Total.transporter ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(Total.transporter_slow_basic)

Total.transporter_power.law.slow  = lm(log(Total.transporter) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(Total.transporter_power.law.slow)

AIC(Total.transporter_slow_basic, Total.transporter_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$Total.transporter), 
     xlab="log(Genome size)", ylab="log(Total transporter)", pch=18, col="black")
abline(Total.transporter_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$Total.transporter), 
     xlab="log(Genome size)", ylab="log(Total transporter)", pch=18, col="black")
abline(Total.transporter_power.law.slow, col="gray", lty=6, lwd=2)

# Amino transporter ----

# Fast

Amino.transporter_fast_basic = lm(Amino.transporter ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(Amino.transporter_fast_basic)

Amino.transporter_power.law  = lm(log(Amino.transporter+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(Amino.transporter_power.law)

AIC(Amino.transporter_fast_basic, Amino.transporter_power.law)

# Slow

Amino.transporter_slow_basic = lm(Amino.transporter ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(Amino.transporter_slow_basic)

Amino.transporter_power.law.slow  = lm(log(Amino.transporter+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(Amino.transporter_power.law.slow)

AIC(Amino.transporter_slow_basic, Amino.transporter_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$Amino.transporter), 
     xlab="log(Genome size)", ylab="log(Amino transporter)", pch=18, col="black")
abline(Amino.transporter_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$Amino.transporter), 
     xlab="log(Genome size)", ylab="log(Amino transporter)", pch=18, col="black")
abline(Amino.transporter_power.law.slow, col="gray", lty=6, lwd=2)

# GH transporter ----

# Fast

GH.transporter_fast_basic = lm(GH.transporter ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(GH.transporter_fast_basic)

GH.transporter_power.law  = lm(log(GH.transporter+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(GH.transporter_power.law)

AIC(GH.transporter_fast_basic, GH.transporter_power.law)

# Slow

GH.transporter_slow_basic = lm(GH.transporter ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(GH.transporter_slow_basic)

GH.transporter_power.law.slow  = lm(log(GH.transporter+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(GH.transporter_power.law.slow)

AIC(GH.transporter_slow_basic, GH.transporter_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$GH.transporter), 
     xlab="log(Genome size)", ylab="log(GH transporter)", pch=18, col="black")
abline(GH.transporter_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$GH.transporter), 
     xlab="log(Genome size)", ylab="log(GH transporter)", pch=18, col="black")
abline(GH.transporter_power.law.slow, col="gray", lty=6, lwd=2)

# Osmolyte ----

# Fast

Osmolyte_fast_basic = lm(Osmolyte ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(Osmolyte_fast_basic)

Osmolyte_power.law  = lm(log(Osmolyte) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(Osmolyte_power.law)

AIC(Osmolyte_fast_basic, Osmolyte_power.law)

# Slow

Osmolyte_slow_basic = lm(Osmolyte ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(Osmolyte_slow_basic)

Osmolyte_power.law.slow  = lm(log(Osmolyte) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(Osmolyte_power.law.slow)

AIC(Osmolyte_slow_basic, Osmolyte_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$Osmolyte), 
     xlab="log(Genome size)", ylab="log(Osmolyte)", pch=18, col="black")
abline(Osmolyte_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$Osmolyte), 
     xlab="log(Genome size)", ylab="log(Osmolyte)", pch=18, col="black")
abline(Osmolyte_power.law.slow, col="gray", lty=6, lwd=2)

# Biofilm ----

# Fast

Biofilm_fast_basic = lm(Biofilm ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(Biofilm_fast_basic)

Biofilm_power.law  = lm(log(Biofilm+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(Biofilm_power.law)

AIC(Biofilm_fast_basic, Biofilm_power.law)

# Slow

Biofilm_slow_basic = lm(Biofilm ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(Biofilm_slow_basic)

Biofilm_power.law.slow  = lm(log(Biofilm) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(Biofilm_power.law.slow)

AIC(Biofilm_slow_basic, Biofilm_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$Biofilm), 
     xlab="log(Genome size)", ylab="log(Biofilm)", pch=18, col="black")
abline(Biofilm_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$Biofilm), 
     xlab="log(Genome size)", ylab="log(Biofilm)", pch=18, col="black")
abline(Biofilm_power.law.slow, col="gray", lty=6, lwd=2)

# Heat resistance ----

# Fast

Temp.Tol_fast_basic = lm(Temp.Tol ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(Temp.Tol_fast_basic)

Temp.Tol_power.law  = lm(log(Temp.Tol+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(Temp.Tol_power.law)

AIC(Temp.Tol_fast_basic, Temp.Tol_power.law)

# Slow

Temp.Tol_slow_basic = lm(Temp.Tol ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(Temp.Tol_slow_basic)

Temp.Tol_power.law.slow  = lm(log(Temp.Tol) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(Temp.Tol_power.law.slow)

AIC(Temp.Tol_slow_basic, Temp.Tol_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$Temp.Tol), 
     xlab="log(Genome size)", ylab="log(Heat resistance)", pch=18, col="black")
abline(Temp.Tol_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$Temp.Tol), 
     xlab="log(Genome size)", ylab="log(Heat resistance)", pch=18, col="black")
abline(Temp.Tol_power.law.slow, col="gray", lty=6, lwd=2)

# pH resistance ----

# Fast

pH.Tol_fast_basic = lm(pH.Tol ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(pH.Tol_fast_basic)

pH.Tol_power.law  = lm(log(pH.Tol) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(pH.Tol_power.law)

AIC(pH.Tol_fast_basic, pH.Tol_power.law)

# Slow

pH.Tol_slow_basic = lm(pH.Tol ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(pH.Tol_slow_basic)

pH.Tol_power.law.slow  = lm(log(pH.Tol) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(pH.Tol_power.law.slow)

AIC(pH.Tol_slow_basic, pH.Tol_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$pH.Tol), 
     xlab="log(Genome size)", ylab="log(pH resistance)", pch=18, col="black")
abline(pH.Tol_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$pH.Tol), 
     xlab="log(Genome size)", ylab="log(pH resistance)", pch=18, col="black")
abline(pH.Tol_power.law.slow, col="gray", lty=6, lwd=2)

# Minimum Generation Time ----

# Fast

MGT_fast_basic = lm(MGT ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(MGT_fast_basic)

MGT_power.law  = lm(log(MGT+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(MGT_power.law)

AIC(MGT_fast_basic, MGT_power.law)

# Slow

MGT_slow_basic = lm(MGT ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(MGT_slow_basic)

MGT_power.law.slow  = lm(log(MGT) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(MGT_power.law.slow)

AIC(MGT_slow_basic, MGT_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$MGT), 
     xlab="log(Genome size)", ylab="log(Minimum Generation Time)", pch=18, col="black")
abline(MGT_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$MGT), 
     xlab="log(Genome size)", ylab="log(Minimum Generation Time)", pch=18, col="black")
abline(MGT_power.law.slow, col="gray", lty=6, lwd=2)

# Optimum Growth Temperature ----

# Fast

OGT_fast_basic = lm(OGT ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(OGT_fast_basic)

OGT_power.law  = lm(log(OGT+1) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(OGT_power.law)

AIC(OGT_fast_basic, OGT_power.law)

# Slow

OGT_slow_basic = lm(OGT ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(OGT_slow_basic)

OGT_power.law.slow  = lm(log(OGT) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(OGT_power.law.slow)

AIC(OGT_slow_basic, OGT_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot(log(ISO_gen_TOTAL_ab_fast$Genome.Size), log(ISO_gen_TOTAL_ab_fast$OGT), 
     xlab="log(Genome size)", ylab="log(Optimum Growth Temperature)", pch=18, col="black")
abline(OGT_power.law, col="gray", lty=6, lwd=2)
# Slow Growing
plot(log(ISO_gen_TOTAL_ab_slow$Genome.Size), log(ISO_gen_TOTAL_ab_slow$OGT), 
     xlab="log(Genome size)", ylab="log(Optimum Growth Temperature)", pch=18, col="black")
abline(OGT_power.law.slow, col="gray", lty=6, lwd=2)

# Yield ----

# Fast

Yield_fast_basic = lm(Yield ~ Genome.Size, data = ISO_gen_TOTAL_ab_fast)
summary(Yield_fast_basic)

Yield_power.law  = lm(log(Yield) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_fast)
summary(Yield_power.law)

AIC(Yield_fast_basic, Yield_power.law)

# Slow

Yield_slow_basic = lm(Yield ~ Genome.Size, data = ISO_gen_TOTAL_ab_slow)
summary(Yield_slow_basic)

Yield_power.law.slow  = lm(log(Yield) ~ log(Genome.Size), data = ISO_gen_TOTAL_ab_slow)
summary(Yield_power.law.slow)

AIC(Yield_slow_basic, Yield_power.law.slow)

# Power Law Plot

# create a new plotting window and set the plotting area into a 2*1 array
par(mfrow = c(1, 2))
# Fast Growing
plot((ISO_gen_TOTAL_ab_fast$Genome.Size), (ISO_gen_TOTAL_ab_fast$Yield), 
     xlab="(Genome size)", ylab="(Yield)", pch=18, col="black")
abline(Yield_fast_basic, col="gray", lty=6, lwd=2)
# Slow Growing
plot((ISO_gen_TOTAL_ab_slow$Genome.Size), (ISO_gen_TOTAL_ab_slow$Yield), 
     xlab="(Genome size)", ylab="(Yield)", pch=18, col="black")
abline(Yield_slow_basic, col="gray", lty=6, lwd=2)

# Linear Regression Plots for LOMA ----

total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")
mag_stat                  = read.csv("Input_Data/LOMA_MAGs/mag_stats.csv",dec=".") 
mag_abun                  = read.csv("Input_Data/LOMA_MAGs/mag_adundance.csv",dec=".") 
mag_stat                  = mag_stat %>% full_join(mag_abun)
colnames(mag_stat)[7]     = "genome.size"

# Merge datasets
LOMA_MAGs                 = merge(x = mag_stat, y = total_genes.guild.940_MAG, 
                                  by = "genome.size")
# Remove duplicaes
LOMA_MAGs                 = LOMA_MAGs[!duplicated(LOMA_MAGs$genome.size), ]
LOMA_MAGs                 = LOMA_MAGs %>% select(id.y, guild, genome.size, completeness,
                                                 contamination, Average, 
                                                 Average.1, Average.2, Average.3)

# Erase gigantic MAG
LOMA_MAGs                 = subset(LOMA_MAGs, guild!= 247)

# Grassland
LOMA_MAGs_grass           = LOMA_MAGs %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size), completeness = mean(completeness),
            contamination = mean(contamination),
            abundance = sum(Average)*100/sum(LOMA_MAGs$Average)) %>% 
  mutate(condition = rep("Grassland Ambient" , length(unique(LOMA_MAGs$guild))))

# Grassland Drought
LOMA_MAGs_grass_Drought   = LOMA_MAGs %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size), completeness = mean(completeness),
            contamination = mean(contamination),
            abundance = sum(Average.1)*100/sum(LOMA_MAGs$Average.1)) %>% 
  mutate(condition = rep("Grassland Drought" , length(unique(LOMA_MAGs$guild))))

# Shrubland
LOMA_MAGs_shrub           = LOMA_MAGs %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size), completeness = mean(completeness),
            contamination = mean(contamination),
            abundance = sum(Average.2)*100/sum(LOMA_MAGs$Average.2)) %>% 
  mutate(condition = rep("Shrubland Ambient" , length(unique(LOMA_MAGs$guild))))

# Shrubland Drought
LOMA_MAGs_shrub_Drought   = LOMA_MAGs %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size), completeness = mean(completeness),
            contamination = mean(contamination),
            abundance = sum(Average.3)*100/sum(LOMA_MAGs$Average.3)) %>% 
  mutate(condition = rep("Shrubland Drought" , length(unique(LOMA_MAGs$guild))))

# Combined MAGs

LOMA_MAG_complete = as.data.frame(rbind(LOMA_MAGs_grass, LOMA_MAGs_grass_Drought,
                                        LOMA_MAGs_shrub, LOMA_MAGs_shrub_Drought))

write.csv(LOMA_MAG_complete, file = "Intermediate_Results/LOMA_MAG_complete.csv")

# Plot
LOMA_MAG_complete                 = read.csv("Intermediate_Results/LOMA_MAG_complete.csv",dec=".")

ggplot(LOMA_MAG_complete) +
  geom_point(aes(x = genome.size, y = abundance, size = abundance, color = condition))

ggplot(LOMA_MAG_complete, aes(x = condition, y = genome.size)) +
  geom_boxplot() + geom_jitter(position=position_jitter(0.1), 
                               aes(color = abundance), size = 2) + 
  scale_color_gradient(low="blue", high="red")

# Linear Regression Plots for FIRE project ----

total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")
hmm_fire    = read.csv("Input_Data/FIRE_MAGs/hmm_Fire.csv",dec=".")
mag_abun    = read.csv("Input_Data/FIRE_MAGs/mag_adundance_fire.csv",dec=".") 
mag_abun    = mag_abun[-c(440,546), ]
fire_meta   = read_excel("Input_Data/FIRE_MAGs/MAG_Dataset_BurnSeverity_ARNelson.xlsx")

# Dataframes from each treatment
Low_shallow   = hmm_fire %>% mutate(Rel.Abund = 100*mag_abun$Avg_low_shallow)
Low_shallow   = Low_shallow %>% mutate(treatment = rep("Low Shallow",each=nrow(Low_shallow)))
Low_deep      = hmm_fire %>% mutate(Rel.Abund = 100*mag_abun$Avg_low_deep)
Low_deep      = Low_deep %>% mutate(treatment = rep("Low Deep",each=nrow(Low_deep)))
High_shallow  = hmm_fire %>% mutate(Rel.Abund = 100*mag_abun$Avg_high_shallow)
High_shallow  = High_shallow %>% mutate(treatment = rep("High Shallow",each=nrow(High_shallow)))
High_deep     = hmm_fire %>% mutate(Rel.Abund = 100*mag_abun$Avg_high_deep)
High_deep     = High_deep %>% mutate(treatment = rep("High Deep",each=nrow(High_deep)))
data.fire     = as.data.frame(rbind(Low_shallow,Low_deep,High_shallow,High_deep))

# Select desired columns
FIRE_MAGs     = data.fire %>% select(id,Rel.Abund,treatment)

# Merge datasets
FIRE_MAGs                 = merge(x = FIRE_MAGs, y = total_genes.guild.940_MAG, 
                                  by = "id")
FIRE_MAGs                 = FIRE_MAGs %>% select(id, guild, genome.size, 
                                                 Rel.Abund, treatment)
# Low Shallow
FIRE_MAGs_LS              = FIRE_MAGs %>% filter(treatment == "Low Shallow")
FIRE_MAGs_LS              = FIRE_MAGs_LS %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size),
            abundance = sum(Rel.Abund)*100/sum(FIRE_MAGs_LS$Rel.Abund)) %>% 
  mutate(condition = rep("Low Shallow" , length(unique(FIRE_MAGs$guild))))

# Low Deep
FIRE_MAGs_LD              = FIRE_MAGs %>% filter(treatment == "Low Deep")
FIRE_MAGs_LD              = FIRE_MAGs_LD %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size),
            abundance = sum(Rel.Abund)*100/sum(FIRE_MAGs_LD$Rel.Abund)) %>% 
  mutate(condition = rep("Low Deep" , length(unique(FIRE_MAGs$guild))))

# High Shallow
FIRE_MAGs_SH              = FIRE_MAGs %>% filter(treatment == "High Shallow")
FIRE_MAGs_SH              = FIRE_MAGs_SH %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size),
            abundance = sum(Rel.Abund)*100/sum(FIRE_MAGs_SH$Rel.Abund)) %>% 
  mutate(condition = rep("High Shallow" , length(unique(FIRE_MAGs$guild))))

# High Deep
FIRE_MAGs_HD              = FIRE_MAGs %>% filter(treatment == "High Deep")
FIRE_MAGs_HD              = FIRE_MAGs_HD %>% group_by(guild) %>% 
  summarise(genome.size = mean(genome.size),
            abundance = sum(Rel.Abund)*100/sum(FIRE_MAGs_HD$Rel.Abund)) %>% 
  mutate(condition = rep("High Deep" , length(unique(FIRE_MAGs$guild))))

FIRE_MAGs     = as.data.frame(rbind(FIRE_MAGs_LS, FIRE_MAGs_LD, FIRE_MAGs_SH, FIRE_MAGs_HD))

write.csv(FIRE_MAGs, file = "Intermediate_Results/FIRE_MAGs_complete.csv")

# Plot 
ggplot(FIRE_MAGs) +
  geom_point(aes(x = genome.size, y = abundance, size = abundance, color = condition))

ggplot(FIRE_MAGs, aes(x = condition, y = genome.size)) +
  geom_boxplot() + geom_jitter(position=position_jitter(0.1), 
                               aes(color = abundance), size = 2) + 
  scale_color_gradient(low="blue", high="red")

ggplot(FIRE_MAGs, aes(x = condition, y = genome.size, color = abundance)) +
  geom_jitter(position=position_jitter(0.1), size = 2) + 
  scale_color_gradient(low="blue", high="red")

# All MAGs with biome data ----

metasoil_init = readRDS("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_emergent/MAG_metadata_soildb_biome_table.RDS")
total_genes.guild.940_MAG = read.csv("Intermediate_Results/total_genes.guild.940_MAG.csv",dec=".")

# Merge datasets
TOTAL_MAGs                = merge(x = total_genes.guild.940_MAG, y = metasoil_init, 
                                  by = "id")
TOTAL_MAGs_1              = TOTAL_MAGs %>% select(id, guild, genome.size, 
                                                  Biome)
TOTAL_MAGs_1              = TOTAL_MAGs_1 %>% group_by(guild,Biome) %>% 
  summarise(genome.size = mean(genome.size))
# Erase NAs
TOTAL_MAGs_2              = na.omit(TOTAL_MAGs_1)

# Plot

ggplot(TOTAL_MAGs_2, aes(x = Biome, y = genome.size/1e6)) +
  geom_boxplot() + geom_jitter(position=position_jitter(0.1), size = 2) + 
  scale_color_gradient(low="blue", high="red") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ----------------------------------------------------------------------------

# GGPLOTS - MAGs ----

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
                                             y = as.numeric(CAZy),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log(Genome Size)") + 
  ylab("Log(CAZy)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                   y = as.numeric(Protein),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log(Genome Size)") + 
  ylab("Log(Protein)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                                   y = as.numeric(tranport.total),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Transporters)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.5,5.5)
Figure_tranport.total

png("Output_Data/Figures/Figure_tranport.total.png",
    width=3500,height=3500*3/5,res=300)
print(Figure_tranport.total)
dev.off()

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
                                                             y = as.numeric(GH.trasnport),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (GH transporters)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                                   y = as.numeric(amino.transport),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4) + xlab("Log (Genome Size)") + 
  ylab("Log (Animo transporters)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                    y = as.numeric(osmolyte),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Osmolytes)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.7,3.5)
Figure_osmolyte

png("Output_Data/Figures/Figure_osmolyte.png",
    width=3500,height=3500*3/5,res=300)
print(Figure_osmolyte)
dev.off()

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
                                                  y = as.numeric(biofilm),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Biofilm)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.2,2.8)
Figure_biofilm

# Heat Resistance ----

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
                                                      y = as.numeric(temp_fast),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Heat Resistance)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.5,2.6)
Figure_temp_fast

# pH Resistance Transporters ----

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
                                             y = as.numeric(pH),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (pH Resistance)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(-2.7,2.85)
Figure_pH_fast

# Minimum Generation Time ----

mgr_data_fast = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_fast$genome.size), 
                                    log(MAG_gen_TOTAL_ab_fast$mgr), 
                                    rep("fast",nrow(MAG_gen_TOTAL_ab_fast))))
colnames(mgr_data_fast) = c("genome.size","mgr","speed")
mgr_data_slow = as.data.frame(cbind(log(MAG_gen_TOTAL_ab_slow$genome.size), 
                                    log(MAG_gen_TOTAL_ab_slow$mgr), 
                                    rep("slow",nrow(MAG_gen_TOTAL_ab_slow))))
colnames(mgr_data_slow) = c("genome.size","mgr","speed")

mgr_data     = as.data.frame(rbind(mgr_data_fast,
                                   mgr_data_slow))

Figure_mgr_fast  = ggplot(data = mgr_data, aes(x = as.numeric(genome.size), 
                                               y = as.numeric(mgr),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Minimum Generation Time)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
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
                                               y = as.numeric(ogt),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Optimum Growth Temperature)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(2.9,4.3)
Figure_ogt_fast

# GGPLOTS - Isolates ----

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
                                                                 y = as.numeric(CAZy),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log (Genome Size)") + 
  ylab("Log (CAZy)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                                                     y = as.numeric(Protein.enzyme),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 5) + xlab("Log (Genome Size)") + 
  ylab("Log (Protein)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                                     y = as.numeric(tranport.total),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 4, label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (Transporters)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(-0.5,5.5)
Figure_tranport.total_ISO

png("Output_Data/Figures/Figure_tranport.total_ISO.png",
    width=3500,height=3500*3/5,res=300)
print(Figure_tranport.total)
dev.off()

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
                                                                 y = as.numeric(GH.transporter),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "right") + xlab("Log (Genome Size)") + 
  ylab("Log (GH transporters)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                                       y = as.numeric(Amino.transporter),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4)) + xlab("Log (Genome Size)") + 
  ylab("Log (Amino transporters)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=24)) + 
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
                                                             y = as.numeric(osmolyte),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Osmolytes)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
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
                                                           y = as.numeric(Biofilm),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Biofilm)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(-3.2,2.8)
Figure_Biofilm_ISO

# Heat resistance ----

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
                                                             y = as.numeric(Temp.Tol),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Heat Resistance)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(0.5,2.6)
Figure_Temp.Tol_ISO

# pH resistance ----

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
                                                         y = as.numeric(pH.Tol),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (pH Resistance)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
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
                                                   y = as.numeric(Yield),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Yield)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5)
Figure_Yield_ISO

# Minimum Generation Time ----

mgr_data_fast_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_fast$Genome.Size), 
                                        log(ISO_gen_TOTAL_ab_fast$MGT), 
                                        rep("fast",nrow(ISO_gen_TOTAL_ab_fast))))
colnames(mgr_data_fast_ISO) = c("genome.size","mgr","speed")
mgr_data_slow_ISO = as.data.frame(cbind(log(ISO_gen_TOTAL_ab_slow$Genome.Size), 
                                        log(ISO_gen_TOTAL_ab_slow$MGT), 
                                        rep("slow",nrow(ISO_gen_TOTAL_ab_slow))))
colnames(mgr_data_slow_ISO) = c("genome.size","mgr","speed")

mgr_data_ISO     = as.data.frame(rbind(mgr_data_fast_ISO,
                                       mgr_data_slow_ISO))

Figure_mgr_ISO   = ggplot(data = mgr_data_ISO, aes(x = as.numeric(genome.size), 
                                                   y = as.numeric(mgr),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Minimum Generation Time)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
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
                                                   y = as.numeric(OGT),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "bottom", label.x = "left") + xlab("Log (Genome Size)") + 
  ylab("Log (Optimum Growth Temperature)") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(13,16.5) + ylim(2.9,4.3)
Figure_OGT_ISO

# S vs A tradeoffs (MAGs) ----

MAG_gen_TOTAL_ab = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")
MAG_gen_test     = MAG_gen_TOTAL_ab %>% mutate(A_trait = rowSums(MAG_gen_TOTAL_ab[,10:12]),
                                               S_trait = rowSums(MAG_gen_TOTAL_ab[,5:8]))
MAG_gen_test     = MAG_gen_test %>% mutate(A_S = A_trait/S_trait)
MAG_gen_test     = MAG_gen_test %>% mutate(speed = case_when(mgr <= 5 ~ "fast",
                                                             mgr  > 5 ~ "slow"))


Figure_genome_A_S_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((genome.size)/10^6), 
                                                           y = as.numeric((A_S)),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "right") + xlab("Mbp") + 
  ylab("A/S ratio") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(0.5,12.5) + ylim(0.8,8.1)

Figure_genome_A_S_MAGs

# S vs A tradeoffs (Isolates) ----

ISO_gen_TOTAL_ab = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_test     = ISO_gen_TOTAL_ab %>% mutate(A_trait = rowSums(ISO_gen_TOTAL_ab[,9:11]),
                                               S_trait = rowSums(ISO_gen_TOTAL_ab[,12:15]))
ISO_gen_test     = ISO_gen_test %>% mutate(A_S = A_trait/S_trait)
ISO_gen_test     = ISO_gen_test %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                                             MGT  > 5 ~ "slow"))


Figure_genome_A_S_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Genome.Size)/10^6), 
                                                          y = as.numeric((A_S)),color = speed)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "right") + xlab("Mbp") + 
  ylab("A/S ratio") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.position="none") + xlim(0.5,12.5) + ylim(0.8,8.1)

Figure_genome_A_S_ISO

# pH (MAGs) ----

Figure_Total_pH_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((tranport.total)), 
                                                         y = as.numeric((pH)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20) # + xlim(0.5,12.5) 
Figure_Total_pH_MAGs

Figure_amino_pH_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((amino.transport)), 
                                                         y = as.numeric((pH)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_pH_MAGs

Figure_amino_GH_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((GH.trasnport)), 
                                                         y = as.numeric((pH)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_GH_MAGs

Figure_protein_pH_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((Protein)), 
                                                           y = as.numeric((pH)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_pH_MAGs

Figure_protein_CAZy_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((CAZy)), 
                                                             y = as.numeric((pH)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_CAZy_MAGs

# pH (Isolates) ----

Figure_amino_pH_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Amino.transporter)), 
                                                        y = as.numeric((pH.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_pH_ISO

Figure_protein_pH_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Protein.enzyme)), 
                                                          y = as.numeric((pH.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_pH_ISO

Figure_amino_GH_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((GH.transporter)), 
                                                        y = as.numeric((pH.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_GH_ISO

Figure_Total_pH_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Total.transporter)), 
                                                        y = as.numeric((pH.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_Total_pH_ISO

Figure_CAZy_pH_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((CAZy)), 
                                                       y = as.numeric((pH.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("pH") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,20)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_CAZy_pH_ISO

# Temperature (MAGs) ----

Figure_Total_temp_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((tranport.total)), 
                                                           y = as.numeric((temp)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15) # + xlim(0.5,12.5) 
Figure_Total_temp_MAGs

Figure_amino_temp_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((amino.transport)), 
                                                           y = as.numeric((temp)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_temp_MAGs

Figure_GH_temp_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((GH.trasnport)), 
                                                        y = as.numeric((temp)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_GH_temp_MAGs

Figure_protein_temp_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((Protein)), 
                                                             y = as.numeric((temp)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_temp_MAGs

Figure_CAZy_temp_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((CAZy)), 
                                                          y = as.numeric((temp)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_CAZy_temp_MAGs

# Temperature (Isolates) ----

Figure_amino_temp_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Amino.transporter)), 
                                                          y = as.numeric((Temp.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_temp_ISO

Figure_protein_temp_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Protein.enzyme)), 
                                                            y = as.numeric((Temp.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_temp_ISO

Figure_GH_temp_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((GH.transporter)), 
                                                       y = as.numeric((Temp.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_GH_temp_ISO

Figure_Total_temp_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Total.transporter)), 
                                                          y = as.numeric((Temp.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_Total_temp_ISO

Figure_CAZy_temp_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((CAZy)), 
                                                         y = as.numeric((Temp.Tol)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("Temp") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,15)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_CAZy_temp_ISO

# Biofilm (MAGs) ----

Figure_Total_biofilm_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((tranport.total)), 
                                                              y = as.numeric((biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16) # + xlim(0.5,12.5) 
Figure_Total_biofilm_MAGs

Figure_amino_biofilm_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((amino.transport)), 
                                                              y = as.numeric((biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_biofilm_MAGs

Figure_GH_biofilm_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((GH.trasnport)), 
                                                           y = as.numeric((biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_GH_biofilm_MAGs

Figure_protein_biofilm_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((Protein)), 
                                                                y = as.numeric((biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_biofilm_MAGs

Figure_CAZy_biofilm_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((CAZy)), 
                                                             y = as.numeric((biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_CAZy_biofilm_MAGs

# Biofilm (Isolates) ----

Figure_amino_biofilm_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Amino.transporter)), 
                                                             y = as.numeric((Biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_biofilm_ISO

Figure_protein_biofilm_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Protein.enzyme)), 
                                                               y = as.numeric((Biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_biofilm_ISO

Figure_GH_biofilm_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((GH.transporter)), 
                                                          y = as.numeric((Biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_GH_biofilm_ISO

Figure_Total_biofilm_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Total.transporter)), 
                                                             y = as.numeric((Biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_Total_biofilm_ISO

Figure_CAZy_biofilm_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((CAZy)), 
                                                            y = as.numeric((Biofilm)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("Biofilm") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,16)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_CAZy_biofilm_ISO

# Osmolyte (MAGs) ----

Figure_Total_osmolyte_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((tranport.total)), 
                                                               y = as.numeric((osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30) # + xlim(0.5,12.5) 
Figure_Total_osmolyte_MAGs

Figure_amino_osmolyte_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((amino.transport)), 
                                                               y = as.numeric((osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_osmolyte_MAGs

Figure_GH_osmolyte_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((GH.trasnport)), 
                                                            y = as.numeric((osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_GH_osmolyte_MAGs

Figure_protein_osmolyte_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((Protein)), 
                                                                 y = as.numeric((osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_osmolyte_MAGs

Figure_CAZy_osmolyte_MAGs   = ggplot(data = MAG_gen_test, aes(x = as.numeric((CAZy)), 
                                                              y = as.numeric((osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_CAZy_osmolyte_MAGs

# Osmolyte (Isolates) ----

Figure_amino_osmolyte_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Amino.transporter)), 
                                                              y = as.numeric((Osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Amino transporters") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_amino_osmolyte_ISO

Figure_protein_osmolyte_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Protein.enzyme)), 
                                                                y = as.numeric((Osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Protein") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30) # + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_protein_osmolyte_ISO

Figure_GH_osmolyte_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((GH.transporter)), 
                                                           y = as.numeric((Osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("GH transporters") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_GH_osmolyte_ISO

Figure_Total_osmolyte_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((Total.transporter)), 
                                                              y = as.numeric((Osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("Total transporters") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_Total_osmolyte_ISO

Figure_CAZy_osmolyte_ISO   = ggplot(data = ISO_gen_test, aes(x = as.numeric((CAZy)), 
                                                             y = as.numeric((Osmolyte)),color = A_S)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"), size = 4), label.y = "top", label.x = "left") + 
  xlab("CAZy") + 
  ylab("Osmolyte") +
  geom_point(size = 3) + theme_classic() + theme(text = element_text(size=14)) + 
  theme() + scale_color_gradient(low="blue", high="red") + ylim(0,30)# + xlim(0.5,12.5) + ylim(0.8,8.1)
Figure_CAZy_osmolyte_ISO

# Additional models ----

library(domir)

# S vs A tradeoffs (MAGs) ----

MAG_gen_TOTAL_ab = read.csv("Intermediate_Results/MAG_gen_TOTAL_GC.csv",dec=".")
MAG_gen_test     = MAG_gen_TOTAL_ab %>% mutate(A_trait = rowSums(MAG_gen_TOTAL_ab[,10:12]),
                                               S_trait = rowSums(MAG_gen_TOTAL_ab[,5:8]))
MAG_gen_test     = MAG_gen_test %>% mutate(A_S = A_trait/S_trait)
MAG_gen_test     = MAG_gen_test %>% mutate(speed = case_when(mgr <= 5 ~ "fast",
                                                             mgr  > 5 ~ "slow"))

MAG_gen_test_log = as.data.frame(cbind(MAG_gen_test[,1:2],log(MAG_gen_test[3:14]),
                                       MAG_gen_test[,15:20],log(MAG_gen_test[21]),
                                       MAG_gen_test[,22]))
MAG_gen_test_log = subset(MAG_gen_test_log, amino.transport!="-Inf") 
MAG_gen_test_log = subset(MAG_gen_test_log, biofilm!="-Inf")
MAG_gen_test_log = subset(MAG_gen_test_log, GH.trasnport!="-Inf")
MAG_gen_test_log = subset(MAG_gen_test_log, CAZy!="-Inf")

# CAZy enzyme ----

cazy_fast_basic_total       = lm(CAZy ~ genome.size, data = MAG_gen_test)
summary(cazy_fast_basic_total)

cazy_fast_basic_total_A_S   = lm(CAZy ~ genome.size + A_S, data = MAG_gen_test)
summary(cazy_fast_basic_total_A_S)
domin(CAZy ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

cazy_fast_basic_total_A_S_I = lm(CAZy ~ genome.size * A_S, data = MAG_gen_test)
summary(cazy_fast_basic_total_A_S_I)
domin(CAZy ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

cazy_fast_basic_total_A_S_P = lm((CAZy) ~ (genome.size) + (A_S), 
                                 data = MAG_gen_test_log, na.action=na.exclude)
summary(cazy_fast_basic_total_A_S_P)
domin(CAZy ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Protein enzyme ----

Protein_fast_basic_total       = lm(Protein ~ genome.size, data = MAG_gen_test)
summary(Protein_fast_basic_total)

Protein_fast_basic_total_A_S   = lm(Protein ~ genome.size + A_S, data = MAG_gen_test)
summary(Protein_fast_basic_total_A_S)
domin(Protein ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

Protein_fast_basic_total_A_S_I = lm(Protein ~ genome.size * A_S, data = MAG_gen_test)
summary(Protein_fast_basic_total_A_S_I)
domin(Protein ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

Protein_fast_basic_total_A_S_P = lm((Protein) ~ (genome.size) + (A_S), 
                                    data = MAG_gen_test_log, na.action=na.exclude)
summary(Protein_fast_basic_total_A_S_P)
domin(Protein ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Transport total Transporters ----

tranport.total_fast_basic_total       = lm(tranport.total ~ genome.size, data = MAG_gen_test)
summary(tranport.total_fast_basic_total)

tranport.total_fast_basic_total_A_S   = lm(tranport.total ~ genome.size + A_S, data = MAG_gen_test)
summary(tranport.total_fast_basic_total_A_S)
domin(tranport.total ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

tranport.total_fast_basic_total_A_S_I = lm(tranport.total ~ genome.size * A_S, data = MAG_gen_test)
summary(tranport.total_fast_basic_total_A_S_I)
domin(tranport.total ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

tranport.total_fast_basic_total_A_S_P = lm((tranport.total) ~ (genome.size) + (A_S), 
                                           data = MAG_gen_test_log, na.action=na.exclude)
summary(tranport.total_fast_basic_total_A_S_P)
domin(tranport.total ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# GH total Transporters ----

GH.trasnport_fast_basic_total       = lm(GH.trasnport ~ genome.size, data = MAG_gen_test)
summary(GH.trasnport_fast_basic_total)

GH.trasnport_fast_basic_total_A_S   = lm(GH.trasnport ~ genome.size + A_S, data = MAG_gen_test)
summary(GH.trasnport_fast_basic_total_A_S)
domin(GH.trasnport ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

GH.trasnport_fast_basic_total_A_S_I = lm(GH.trasnport ~ genome.size * A_S, data = MAG_gen_test)
summary(GH.trasnport_fast_basic_total_A_S_I)
domin(GH.trasnport ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

GH.trasnport_fast_basic_total_A_S_P = lm((GH.trasnport) ~ (genome.size) + (A_S), 
                                         data = MAG_gen_test_log, na.action=na.exclude)
summary(GH.trasnport_fast_basic_total_A_S_P)
domin(GH.trasnport ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Amino total Transporters ----

amino.transport_fast_basic_total       = lm(amino.transport ~ genome.size, data = MAG_gen_test)
summary(amino.transport_fast_basic_total)

amino.transport_fast_basic_total_A_S   = lm(amino.transport ~ genome.size + A_S, data = MAG_gen_test)
summary(amino.transport_fast_basic_total_A_S)
domin(amino.transport ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

amino.transport_fast_basic_total_A_S_I = lm(amino.transport ~ genome.size * A_S, data = MAG_gen_test)
summary(amino.transport_fast_basic_total_A_S_I)
domin(amino.transport ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

amino.transport_fast_basic_total_A_S_P = lm((amino.transport) ~ (genome.size) + (A_S), 
                                            data = MAG_gen_test_log, na.action=na.exclude)
summary(amino.transport_fast_basic_total_A_S_P)
domin(amino.transport ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Osmolytes ----

osmolyte_fast_basic_total       = lm(osmolyte ~ genome.size, data = MAG_gen_test)
summary(osmolyte_fast_basic_total)

osmolyte_fast_basic_total_A_S   = lm(osmolyte ~ genome.size + A_S, data = MAG_gen_test)
summary(osmolyte_fast_basic_total_A_S)
domin(osmolyte ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

osmolyte_fast_basic_total_A_S_I = lm(osmolyte ~ genome.size * A_S, data = MAG_gen_test)
summary(osmolyte_fast_basic_total_A_S_I)
domin(osmolyte ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

osmolyte_fast_basic_total_A_S_P = lm((osmolyte) ~ (genome.size) + (A_S), 
                                     data = MAG_gen_test_log, na.action=na.exclude)
summary(osmolyte_fast_basic_total_A_S_P)
domin(osmolyte ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Biofilm ----

biofilm_fast_basic_total       = lm(biofilm ~ genome.size, data = MAG_gen_test)
summary(biofilm_fast_basic_total)

biofilm_fast_basic_total_A_S   = lm(biofilm ~ genome.size + A_S, data = MAG_gen_test)
summary(biofilm_fast_basic_total_A_S)
domin(biofilm ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

biofilm_fast_basic_total_A_S_I = lm(biofilm ~ genome.size * A_S, data = MAG_gen_test)
summary(biofilm_fast_basic_total_A_S_I)
domin(biofilm ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

biofilm_fast_basic_total_A_S_P = lm((biofilm) ~ (genome.size) + (A_S), 
                                    data = MAG_gen_test_log, na.action=na.exclude)
summary(biofilm_fast_basic_total_A_S_P)
domin(biofilm ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Heat Resistance ----

temp_fast_basic_total       = lm(temp ~ genome.size, data = MAG_gen_test)
summary(temp_fast_basic_total)

temp_fast_basic_total_A_S   = lm(temp ~ genome.size + A_S, data = MAG_gen_test)
summary(temp_fast_basic_total_A_S)
domin(temp ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

temp_fast_basic_total_A_S_I = lm(temp ~ genome.size * A_S, data = MAG_gen_test)
summary(temp_fast_basic_total_A_S_I)
domin(temp ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

temp_fast_basic_total_A_S_P = lm((temp) ~ (genome.size) + (A_S), 
                                 data = MAG_gen_test_log, na.action=na.exclude)
summary(temp_fast_basic_total_A_S_P)
domin(temp ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# pH Resistance Transporters ----

pH_fast_basic_total       = lm(pH ~ genome.size, data = MAG_gen_test)
summary(pH_fast_basic_total)

pH_fast_basic_total_A_S   = lm(pH ~ genome.size + A_S, data = MAG_gen_test)
summary(pH_fast_basic_total_A_S)
domin(pH ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

pH_fast_basic_total_A_S_I = lm(pH ~ genome.size * A_S, data = MAG_gen_test)
summary(pH_fast_basic_total_A_S_I)
domin(pH ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

pH_fast_basic_total_A_S_P = lm((pH) ~ (genome.size) + (A_S), 
                               data = MAG_gen_test_log, na.action=na.exclude)
summary(pH_fast_basic_total_A_S_P)
domin(pH ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Minimum Generation Time ----

mgr_fast_basic_total       = lm(mgr ~ genome.size, data = MAG_gen_test)
summary(mgr_fast_basic_total)

mgr_fast_basic_total_A_S   = lm(mgr ~ genome.size + A_S, data = MAG_gen_test)
summary(mgr_fast_basic_total_A_S)
domin(mgr ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

mgr_fast_basic_total_A_S_I = lm(mgr ~ genome.size * A_S, data = MAG_gen_test)
summary(mgr_fast_basic_total_A_S_I)
domin(mgr ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

mgr_fast_basic_total_A_S_P = lm((mgr) ~ (genome.size) + (A_S), 
                                data = MAG_gen_test_log, na.action=na.exclude)
summary(mgr_fast_basic_total_A_S_P)
domin(mgr ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# Optimum Growth Temperature ----

ogt_fast_basic_total       = lm(ogt ~ genome.size, data = MAG_gen_test)
summary(ogt_fast_basic_total)

ogt_fast_basic_total_A_S   = lm(ogt ~ genome.size + A_S, data = MAG_gen_test)
summary(ogt_fast_basic_total_A_S)
domin(ogt ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

ogt_fast_basic_total_A_S_I = lm(ogt ~ genome.size * A_S, data = MAG_gen_test)
summary(ogt_fast_basic_total_A_S_I)
domin(ogt ~ genome.size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test)

ogt_fast_basic_total_A_S_P = lm((ogt) ~ (genome.size) + (A_S), 
                                data = MAG_gen_test_log, na.action=na.exclude)
summary(ogt_fast_basic_total_A_S_P)
domin(ogt ~ genome.size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = MAG_gen_test_log)

# S vs A tradeoffs (Isolates) ----

ISO_gen_TOTAL_ab = read.csv("Intermediate_Results/ISO_gen_trait.1.csv",dec=".")
ISO_gen_test     = ISO_gen_TOTAL_ab %>% mutate(A_trait = rowSums(ISO_gen_TOTAL_ab[,9:11]),
                                               S_trait = rowSums(ISO_gen_TOTAL_ab[,12:15]))
ISO_gen_test     = ISO_gen_test %>% mutate(A_S = A_trait/S_trait)
ISO_gen_test     = ISO_gen_test %>% mutate(speed = case_when(MGT <= 5 ~ "fast",
                                                             MGT  > 5 ~ "slow"))

ISO_gen_test_log = as.data.frame(cbind(ISO_gen_test[,1:2],log(ISO_gen_test[3:15]),
                                       ISO_gen_test[,16:17],log(ISO_gen_test[18]),
                                       ISO_gen_test[,19]))
ISO_gen_test_log = subset(ISO_gen_test_log, Amino.transporter!="-Inf") 
ISO_gen_test_log = subset(ISO_gen_test_log, GH.transporter!="-Inf")
ISO_gen_test_log = subset(ISO_gen_test_log, Biofilm!="-Inf")

# CAZy enzyme ----

cazy_fast_basic_total       = lm(CAZy ~ Genome.Size, data = ISO_gen_test)
summary(cazy_fast_basic_total)

cazy_fast_basic_total_A_S   = lm(CAZy ~ Genome.Size + A_S, data = ISO_gen_test)
summary(cazy_fast_basic_total_A_S)
domin(CAZy ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

cazy_fast_basic_total_A_S_I = lm(CAZy ~ Genome.Size * A_S, data = ISO_gen_test)
summary(cazy_fast_basic_total_A_S_I)
domin(CAZy ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

cazy_fast_basic_total_A_S_P = lm((CAZy) ~ (Genome.Size) + (A_S), 
                                 data = ISO_gen_test_log, na.action=na.exclude)
summary(cazy_fast_basic_total_A_S_P)
domin(CAZy ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Protein enzyme ----

Protein_fast_basic_total       = lm(Protein.enzyme ~ Genome.Size, data = ISO_gen_test)
summary(Protein_fast_basic_total)

Protein_fast_basic_total_A_S   = lm(Protein.enzyme ~ Genome.Size + A_S, data = ISO_gen_test)
summary(Protein_fast_basic_total_A_S)
domin(Protein.enzyme ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

Protein_fast_basic_total_A_S_I = lm(Protein.enzyme ~ Genome.Size * A_S, data = ISO_gen_test)
summary(Protein_fast_basic_total_A_S_I)
domin(Protein.enzyme ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

Protein_fast_basic_total_A_S_P = lm((Protein.enzyme) ~ (Genome.Size) + (A_S), 
                                    data = ISO_gen_test_log, na.action=na.exclude)
summary(Protein_fast_basic_total_A_S_P)
domin(Protein.enzyme ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Transport total Transporters ----

tranport.total_fast_basic_total       = lm(Total.transporter ~ Genome.Size, data = ISO_gen_test)
summary(tranport.total_fast_basic_total)

tranport.total_fast_basic_total_A_S   = lm(Total.transporter ~ Genome.Size + A_S, data = ISO_gen_test)
summary(tranport.total_fast_basic_total_A_S)
domin(Total.transporter ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

tranport.total_fast_basic_total_A_S_I = lm(Total.transporter ~ Genome.Size * A_S, data = ISO_gen_test)
summary(tranport.total_fast_basic_total_A_S_I)
domin(Total.transporter ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

tranport.total_fast_basic_total_A_S_P = lm((Total.transporter) ~ (Genome.Size) + (A_S), 
                                           data = ISO_gen_test_log, na.action=na.exclude)
summary(tranport.total_fast_basic_total_A_S_P)
domin(Total.transporter ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# GH total Transporters ----

GH.trasnport_fast_basic_total       = lm(GH.transporter ~ Genome.Size, data = ISO_gen_test)
summary(GH.trasnport_fast_basic_total)

GH.trasnport_fast_basic_total_A_S   = lm(GH.transporter ~ Genome.Size + A_S, data = ISO_gen_test)
summary(GH.trasnport_fast_basic_total_A_S)
domin(GH.transporter ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

GH.trasnport_fast_basic_total_A_S_I = lm(GH.transporter ~ Genome.Size * A_S, data = ISO_gen_test)
summary(GH.trasnport_fast_basic_total_A_S_I)
domin(GH.transporter ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

GH.trasnport_fast_basic_total_A_S_P = lm((GH.transporter) ~ (Genome.Size) + (A_S), 
                                         data = ISO_gen_test_log, na.action=na.exclude)
summary(GH.trasnport_fast_basic_total_A_S_P)
domin(GH.transporter ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Amino total Transporters ----

amino.transport_fast_basic_total       = lm(Amino.transporter ~ Genome.Size, data = ISO_gen_test)
summary(amino.transport_fast_basic_total)

amino.transport_fast_basic_total_A_S   = lm(Amino.transporter ~ Genome.Size + A_S, data = ISO_gen_test)
summary(amino.transport_fast_basic_total_A_S)
domin(Amino.transporter ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

amino.transport_fast_basic_total_A_S_I = lm(Amino.transporter ~ Genome.Size * A_S, data = ISO_gen_test)
summary(amino.transport_fast_basic_total_A_S_I)
domin(Amino.transporter ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

amino.transport_fast_basic_total_A_S_P = lm((Amino.transporter) ~ (Genome.Size) + (A_S), 
                                            data = ISO_gen_test_log, na.action=na.exclude)
summary(amino.transport_fast_basic_total_A_S_P)
domin(Amino.transporter ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Osmolytes ----

osmolyte_fast_basic_total       = lm(Osmolyte ~ Genome.Size, data = ISO_gen_test)
summary(osmolyte_fast_basic_total)

osmolyte_fast_basic_total_A_S   = lm(Osmolyte ~ Genome.Size + A_S, data = ISO_gen_test)
summary(osmolyte_fast_basic_total_A_S)
domin(Osmolyte ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

osmolyte_fast_basic_total_A_S_I = lm(Osmolyte ~ Genome.Size * A_S, data = ISO_gen_test)
summary(osmolyte_fast_basic_total_A_S_I)
domin(Osmolyte ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

osmolyte_fast_basic_total_A_S_P = lm((Osmolyte) ~ (Genome.Size) + (A_S), 
                                     data = ISO_gen_test_log, na.action=na.exclude)
summary(osmolyte_fast_basic_total_A_S_P)
domin(Osmolyte ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Biofilm ----

biofilm_fast_basic_total       = lm(Biofilm ~ Genome.Size, data = ISO_gen_test)
summary(biofilm_fast_basic_total)

biofilm_fast_basic_total_A_S   = lm(Biofilm ~ Genome.Size + A_S, data = ISO_gen_test)
summary(biofilm_fast_basic_total_A_S)
domin(Biofilm ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

biofilm_fast_basic_total_A_S_I = lm(Biofilm ~ Genome.Size * A_S, data = ISO_gen_test)
summary(biofilm_fast_basic_total_A_S_I)
domin(Biofilm ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

biofilm_fast_basic_total_A_S_P = lm((Biofilm) ~ (Genome.Size) + (A_S), 
                                    data = ISO_gen_test_log, na.action=na.exclude)
summary(biofilm_fast_basic_total_A_S_P)
domin(Biofilm ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Heat Resistance ----

temp_fast_basic_total       = lm(Temp.Tol ~ Genome.Size, data = ISO_gen_test)
summary(temp_fast_basic_total)

temp_fast_basic_total_A_S   = lm(Temp.Tol ~ Genome.Size + A_S, data = ISO_gen_test)
summary(temp_fast_basic_total_A_S)
domin(Temp.Tol ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

temp_fast_basic_total_A_S_I = lm(Temp.Tol ~ Genome.Size * A_S, data = ISO_gen_test)
summary(temp_fast_basic_total_A_S_I)
domin(Temp.Tol ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

temp_fast_basic_total_A_S_P = lm((Temp.Tol) ~ (Genome.Size) + (A_S), 
                                 data = ISO_gen_test_log, na.action=na.exclude)
summary(temp_fast_basic_total_A_S_P)
domin(Temp.Tol ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# pH Resistance Transporters ----

pH_fast_basic_total       = lm(pH.Tol ~ Genome.Size, data = ISO_gen_test)
summary(pH_fast_basic_total)

pH_fast_basic_total_A_S   = lm(pH.Tol ~ Genome.Size + A_S, data = ISO_gen_test)
summary(pH_fast_basic_total_A_S)
domin(pH.Tol ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

pH_fast_basic_total_A_S_I = lm(pH.Tol ~ Genome.Size * A_S, data = ISO_gen_test)
summary(pH_fast_basic_total_A_S_I)
domin(pH.Tol ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

pH_fast_basic_total_A_S_P = lm((pH.Tol) ~ (Genome.Size) + (A_S), 
                               data = ISO_gen_test_log, na.action=na.exclude)
summary(pH_fast_basic_total_A_S_P)
domin(pH.Tol ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Minimum Generation Time ----

mgr_fast_basic_total       = lm(MGT ~ Genome.Size, data = ISO_gen_test)
summary(mgr_fast_basic_total)

mgr_fast_basic_total_A_S   = lm(MGT ~ Genome.Size + A_S, data = ISO_gen_test)
summary(mgr_fast_basic_total_A_S)
domin(MGT ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

mgr_fast_basic_total_A_S_I = lm(MGT ~ Genome.Size * A_S, data = ISO_gen_test)
summary(mgr_fast_basic_total_A_S_I)
domin(MGT ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

mgr_fast_basic_total_A_S_P = lm((MGT) ~ (Genome.Size) + (A_S), 
                                data = ISO_gen_test_log, na.action=na.exclude)
summary(mgr_fast_basic_total_A_S_P)
domin(MGT ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Optimum Growth Temperature ----

ogt_fast_basic_total       = lm(OGT ~ Genome.Size, data = ISO_gen_test)
summary(ogt_fast_basic_total)

ogt_fast_basic_total_A_S   = lm(OGT ~ Genome.Size + A_S, data = ISO_gen_test)
summary(ogt_fast_basic_total_A_S)
domin(OGT ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

ogt_fast_basic_total_A_S_I = lm(OGT ~ Genome.Size * A_S, data = ISO_gen_test)
summary(ogt_fast_basic_total_A_S_I)
domin(OGT ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

ogt_fast_basic_total_A_S_P = lm((OGT) ~ (Genome.Size) + (A_S), 
                                data = ISO_gen_test_log, na.action=na.exclude)
summary(ogt_fast_basic_total_A_S_P)
domin(OGT ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)

# Yield ----

Yield_fast_basic_total       = lm(Yield ~ Genome.Size, data = ISO_gen_test)
summary(Yield_fast_basic_total)

Yield_fast_basic_total_A_S   = lm(Yield ~ Genome.Size + A_S, data = ISO_gen_test)
summary(Yield_fast_basic_total_A_S)
domin(Yield ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

Yield_fast_basic_total_A_S_I = lm(Yield ~ Genome.Size * A_S, data = ISO_gen_test)
summary(Yield_fast_basic_total_A_S_I)
domin(Yield ~ Genome.Size * A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test)

Yield_fast_basic_total_A_S_P = lm((Yield) ~ (Genome.Size) + (A_S), 
                                  data = ISO_gen_test_log, na.action=na.exclude)
summary(Yield_fast_basic_total_A_S_P)
domin(Yield ~ Genome.Size + A_S, 
      lm, 
      list(summary, "r.squared"), 
      data = ISO_gen_test_log)