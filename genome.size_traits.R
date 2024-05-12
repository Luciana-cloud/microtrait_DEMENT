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

# CALLED DATA ----

total.guilds        = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/total.genes.guilds.selected.csv",dec=".")
loma_stat           = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_stats.txt") 
loma_stat           = as.data.frame(loma_stat[c(1,7)])
fire_stat           = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fire_metadata.txt")
fire_stat           = fire_stat[-c(440,546), ]
colnames(fire_stat) = c("id","size")
total_stat          = as.data.frame(rbind(fire_stat,loma_stat))
full_mat            = left_join(total.guilds, total_stat, by=c('id'))
a                   = as.data.frame(colnames(full_mat))
GH_rule             = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_GH.txt", header = TRUE)
prot_rule           = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_proteins.txt", header = TRUE)
osmo_rule           = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_osymolites.txt", header = TRUE)
biofilm_rule        = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_biofilm.txt", header = TRUE)
high.T_rule         = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_high_T.txt", header = TRUE)
low.T_rule          = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_low_T.txt", header = TRUE)
pH_rule             = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_pH_stress.txt", header = TRUE)
transp_rule         = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_transporters.txt", header = TRUE,sep = "\t")

# CAZy GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----

GH_TOTAL = full_mat %>% select(any_of(GH_rule$microtrait_hmm.name))
GH_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],GH_TOTAL))

GH_TOTAL_m = GH_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                               list(mean=mean), na.rm=TRUE))
GH_TOTAL_m = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,3:ncol(GH_TOTAL_m)]))

ggplot(GH_TOTAL_m, aes(size_mean,GH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") # It seems that group 62 is too big

ggplot(GH_TOTAL_m, aes(x = size_mean,y = GH_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

GH_TOTAL_m.1 = GH_TOTAL_m[-c(62), ]

ggplot(GH_TOTAL_m.1, aes(x = size_mean,y = GH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG")

ggplot(GH_TOTAL_m.1, aes(x = size_mean,y = GH_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# PROTEIN GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----

PROTEIN_TOTAL = full_mat %>% select(any_of(prot_rule$microtrait_rule.name))
PROTEIN_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],PROTEIN_TOTAL))

PROTEIN_TOTAL_m = PROTEIN_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                               list(mean=mean), na.rm=TRUE))
PROTEIN_TOTAL_m = PROTEIN_TOTAL_m %>% mutate(PROT_total = rowSums(PROTEIN_TOTAL_m[,3:ncol(PROTEIN_TOTAL_m)]))

ggplot(PROTEIN_TOTAL_m, aes(size_mean,PROT_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total Protein costs per FG") # It seems that group 62 is too big

ggplot(PROTEIN_TOTAL_m, aes(x = size_mean,y = PROT_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total Protein costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

PROTEIN_TOTAL_m.1 = PROTEIN_TOTAL_m[-c(62), ]

ggplot(PROTEIN_TOTAL_m.1, aes(x = size_mean,y = PROT_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total Protein costs per FG")

ggplot(PROTEIN_TOTAL_m.1, aes(x = size_mean,y = PROT_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total Protein costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# OSMOLYTE GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----

OSMO_TOTAL = full_mat %>% select(any_of(osmo_rule$microtrait_hmm.name))
OSMO_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],OSMO_TOTAL))

OSMO_TOTAL_m = OSMO_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
OSMO_TOTAL_m = OSMO_TOTAL_m %>% mutate(OSMO_total = rowSums(OSMO_TOTAL_m[,3:ncol(OSMO_TOTAL_m)]))

ggplot(OSMO_TOTAL_m, aes(size_mean,OSMO_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG") # It seems that group 62 is too big

ggplot(OSMO_TOTAL_m, aes(x = size_mean,y = OSMO_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

OSMO_TOTAL_m.1 = OSMO_TOTAL_m[-c(62), ]

ggplot(OSMO_TOTAL_m.1, aes(x = size_mean,y = OSMO_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG")

ggplot(OSMO_TOTAL_m.1, aes(x = size_mean,y = OSMO_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)





# ----------------------------------------------------------------------------- #

# Selecting GHs ###

GHs_guild       = as.data.frame(total.guilds[c(2,3,389,249,194,283,278,451,464,100,
                                             340,35,150,99,234,336,40,46,248,76,
                                             398,43,396)])
GHs_guild.1     = GHs_guild %>% mutate(sum   = rowSums(across(where(is.numeric)))/total_stat$size)
GHs_guild.2     = GHs_guild %>% mutate(sum.2 = rowSums(across(where(is.numeric))))

# Plot gene counts per functional groups - Normalized ###

ggplot(GHs_guild.1, aes(id,sum)) + geom_point() +
  xlab("MAG id") + ylab("Total GH costs")

ggplot(GHs_guild.1, aes(guild,sum)) + geom_point() +
  xlab("Functional group") + ylab("Total GH costs")

# Plot gene counts per functional groups - without Normalized ###

ggplot(GHs_guild.2, aes(id,sum.2)) + geom_point() +
  xlab("MAG id") + ylab("Total GH costs")

ggplot(GHs_guild.2, aes(guild,sum.2)) + geom_point() +
  xlab("Functional group") + ylab("Total GH costs")

# Guild size per functional group ###


GHs_full    = as.data.frame(full_mat[c(2,3,509,249,194,283,278,451,464,100,
                                               340,35,150,99,234,336,40,46,248,76,
                                               398,43,396)])
GHs_full.1  = GHs_full %>% mutate(sum   = rowSums(GHs_full[,4:23])/size)
GHs_full.1  = GHs_full.1 %>% mutate(sum.2 = rowSums(GHs_full[,4:23]))

ggplot(GHs_full.1, aes(size,sum,colour=factor(guild))) + geom_point() +
  xlab("Genome size") + ylab("Normalized Total GH costs")

ggplot(GHs_full.1, aes(size,sum.2,colour=factor(guild))) + geom_point() +
  xlab("Genome size") + ylab("Total GH costs")



# BIOFILM GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----

BIOFILM_TOTAL = full_mat %>% select(any_of(biofilm_rule$microtrait_hmm.name))
BIOFILM_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],BIOFILM_TOTAL))

BIOFILM_TOTAL_m = BIOFILM_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
BIOFILM_TOTAL_m = BIOFILM_TOTAL_m %>% mutate(BIOFILM_total = rowSums(BIOFILM_TOTAL_m[,3:ncol(BIOFILM_TOTAL_m)]))

ggplot(BIOFILM_TOTAL_m, aes(size_mean,BIOFILM_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total biofilm costs per FG") # It seems that group 62 is too big

ggplot(BIOFILM_TOTAL_m, aes(x = size_mean,y = BIOFILM_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total biofilm costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

BIOFILM_TOTAL_m.1 = BIOFILM_TOTAL_m[-c(62), ]

ggplot(BIOFILM_TOTAL_m.1, aes(x = size_mean,y = BIOFILM_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total biofilm costs per FG")

ggplot(BIOFILM_TOTAL_m.1, aes(x = size_mean,y = BIOFILM_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total biofilm costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# HEAT TOLERANCE GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----

HEAT.T_TOTAL = full_mat %>% select(any_of(high.T_rule$microtrait_hmm.name))
HEAT.T_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],HEAT.T_TOTAL))

HEAT.T_TOTAL_m = HEAT.T_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
HEAT.T_TOTAL_m = HEAT.T_TOTAL_m %>% mutate(HEAT.T_total = rowSums(HEAT.T_TOTAL_m[,3:ncol(HEAT.T_TOTAL_m)]))

ggplot(HEAT.T_TOTAL_m, aes(size_mean,HEAT.T_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total heat tolerance costs per FG") # It seems that group 62 is too big

ggplot(HEAT.T_TOTAL_m, aes(x = size_mean,y = HEAT.T_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total heat tolerance costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

HEAT.T_TOTAL_m.1 = HEAT.T_TOTAL_m[-c(62), ]

ggplot(HEAT.T_TOTAL_m.1, aes(x = size_mean,y = HEAT.T_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total heat tolerance costs per FG")

ggplot(HEAT.T_TOTAL_m.1, aes(x = size_mean,y = HEAT.T_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total heat tolerance costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)


# LOW TEMPERATURE TOLERANCE GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----

LOW.T_TOTAL = full_mat %>% select(any_of(low.T_rule$microtrait_hmm.name))
LOW.T_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],LOW.T_TOTAL))

LOW.T_TOTAL_m = LOW.T_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                       list(mean=mean), na.rm=TRUE))
LOW.T_TOTAL_m = LOW.T_TOTAL_m %>% mutate(LOW.T_total = rowSums(LOW.T_TOTAL_m[,3:ncol(LOW.T_TOTAL_m)]))

ggplot(LOW.T_TOTAL_m, aes(size_mean,LOW.T_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total low temperature tolerance costs per FG") # It seems that group 62 is too big

ggplot(LOW.T_TOTAL_m, aes(x = size_mean,y = LOW.T_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total low temperature tolerance costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

LOW.T_TOTAL_m.1 = LOW.T_TOTAL_m[-c(62), ]

ggplot(LOW.T_TOTAL_m.1, aes(x = size_mean,y = LOW.T_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total low temperature tolerance costs per FG")

ggplot(LOW.T_TOTAL_m.1, aes(x = size_mean,y = LOW.T_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total low temperature tolerance costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)


# PH TOLERANCE GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----

PH_TOTAL = full_mat %>% select(any_of(pH_rule$microtrait_hmm.name))
PH_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],PH_TOTAL))

PH_TOTAL_m = PH_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                     list(mean=mean), na.rm=TRUE))
PH_TOTAL_m = PH_TOTAL_m %>% mutate(PH_total = rowSums(PH_TOTAL_m[,3:ncol(PH_TOTAL_m)]))

ggplot(PH_TOTAL_m, aes(size_mean,PH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total pH tolerance costs per FG") # It seems that group 62 is too big

ggplot(PH_TOTAL_m, aes(x = size_mean,y = PH_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total pH tolerance costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

PH_TOTAL_m.1 = PH_TOTAL_m[-c(62), ]

ggplot(PH_TOTAL_m.1, aes(x = size_mean,y = PH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total pH tolerance costs per FG")

ggplot(PH_TOTAL_m.1, aes(x = size_mean,y = PH_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total pH tolerance costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# TRANSPORT GENE COST PER MEAN GENOME SIZE PER FUNCTIONAL GROUP ----
# TOTAL ----
length(unique(transp_rule$microtrait_hmm.name))

transp_rule  = transp_rule %>% filter(function.==c("transporter"))
TRANSP_TOTAL = full_mat %>% select(any_of(transp_rule$microtrait_hmm.name))
TRANSP_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],TRANSP_TOTAL))

TRANSP_TOTAL_m = TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                               list(mean=mean), na.rm=TRUE))
TRANSP_TOTAL_m = TRANSP_TOTAL_m %>% mutate(TRANSP_TOTAL = rowSums(TRANSP_TOTAL_m[,3:ncol(TRANSP_TOTAL_m)]))

ggplot(TRANSP_TOTAL_m, aes(size_mean,TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean transporter costs per FG") # It seems that group 62 is too big

ggplot(TRANSP_TOTAL_m, aes(x = size_mean,y = TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.x = 1e7,label.y = 85) + stat_regline_equation(label.x = 1e7,label.y = 80)

TRANSP_TOTAL_m.1 = TRANSP_TOTAL_m[-c(62), ]

ggplot(TRANSP_TOTAL_m.1, aes(x = size_mean,y = TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean transporter costs per FG")

ggplot(TRANSP_TOTAL_m.1, aes(x = size_mean,y = TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 75) + stat_regline_equation(label.y = 70)

# AMINOACIDS ----
transp_rule_ami  = transp_rule %>% filter(class==c("aminoacid","peptide"))
AMI_TRANSP_TOTAL = full_mat %>% select(any_of(transp_rule_ami$microtrait_hmm.name))
AMI_TRANSP_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],AMI_TRANSP_TOTAL))

AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                       list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL_m[,3:ncol(AMI_TRANSP_TOTAL_m)]))

ggplot(AMI_TRANSP_TOTAL_m, aes(size_mean,AMI_TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean aminoacid transporter costs per FG") # It seems that group 62 is too big

ggplot(AMI_TRANSP_TOTAL_m, aes(x = size_mean,y = AMI_TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean aminoacid transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

AMI_TRANSP_TOTAL_m.1 = AMI_TRANSP_TOTAL_m[-c(62), ]

ggplot(AMI_TRANSP_TOTAL_m.1, aes(x = size_mean,y = AMI_TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean aminoacid transporter costs per FG")

ggplot(AMI_TRANSP_TOTAL_m.1, aes(x = size_mean,y = AMI_TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean aminoacid transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)


# CARBOHYDRATE ----
transp_rule_carb = transp_rule %>% filter(class==c("carbohydrate"))
CARB_TRANSP_TOTAL = full_mat %>% select(any_of(transp_rule_carb$microtrait_hmm.name))
CARB_TRANSP_TOTAL = as_data_frame(cbind(full_mat[c("guild","id","size")],CARB_TRANSP_TOTAL))

CARB_TRANSP_TOTAL_m = CARB_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
CARB_TRANSP_TOTAL_m = CARB_TRANSP_TOTAL_m %>% mutate(CARB_TRANSP_TOTAL = rowSums(CARB_TRANSP_TOTAL_m[,3:ncol(CARB_TRANSP_TOTAL_m)]))

ggplot(CARB_TRANSP_TOTAL_m, aes(size_mean,CARB_TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean carbohydrate transporter costs per FG") # It seems that group 62 is too big

ggplot(CARB_TRANSP_TOTAL_m, aes(x = size_mean,y = CARB_TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean carbohydrate transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

CARB_TRANSP_TOTAL_m.1 = CARB_TRANSP_TOTAL_m[-c(62), ]

ggplot(CARB_TRANSP_TOTAL_m.1, aes(x = size_mean,y = CARB_TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean carbohydrate transporter costs per FG")

ggplot(CARB_TRANSP_TOTAL_m.1, aes(x = size_mean,y = CARB_TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean carbohydrate transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)


# LOMA + FIRE ANALYSIS ----

# CAZy GENE COST ----

GH_LOMA = full_mat[636:1168,] %>% select(any_of(GH_rule$microtrait_hmm.name))
GH_LOMA = as_data_frame(cbind(full_mat[636:1168,][c("guild","id","size")],GH_LOMA))

GH_LOMA_m = GH_LOMA %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                               list(mean=mean), na.rm=TRUE))
GH_LOMA_m = GH_LOMA_m %>% mutate(GH_total = rowSums(GH_LOMA_m[,3:ncol(GH_LOMA_m)]))

GH_LOMA_m = GH_LOMA_m[-c(44), ]

ggplot(GH_LOMA_m, aes(size_mean,GH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") # It seems that group 62 is too big

ggplot(GH_LOMA_m, aes(x = size_mean,y = GH_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

GH_FIRE = full_mat[1:635,] %>% select(any_of(GH_rule$microtrait_hmm.name))
GH_FIRE = as_data_frame(cbind(full_mat[1:635,][c("guild","id","size")],GH_FIRE))

GH_FIRE_m = GH_FIRE %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
GH_FIRE_m = GH_FIRE_m %>% mutate(GH_total = rowSums(GH_FIRE_m[,3:ncol(GH_FIRE_m)]))

ggplot(GH_FIRE_m, aes(size_mean,GH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") # It seems that group 62 is too big

ggplot(GH_FIRE_m, aes(x = size_mean,y = GH_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# PROTEIN GENE COST ----
PROTEIN_LOMA = full_mat[636:1168,] %>% select(any_of(prot_rule$microtrait_rule.name))
PROTEIN_LOMA = as_data_frame(cbind(full_mat[636:1168,][c("guild","id","size")],PROTEIN_LOMA))

PROTEIN_LOMA_m = PROTEIN_LOMA %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
PROTEIN_LOMA_m = PROTEIN_LOMA_m %>% mutate(PROTEIN_total = rowSums(PROTEIN_LOMA_m[,3:ncol(PROTEIN_LOMA_m)]))

PROTEIN_LOMA_m = PROTEIN_LOMA_m[-c(44), ]

ggplot(PROTEIN_LOMA_m, aes(size_mean,PROTEIN_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total Protein costs per FG") # It seems that group 62 is too big

ggplot(PROTEIN_LOMA_m, aes(x = size_mean,y = PROTEIN_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean Protein costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

PROTEIN_FIRE = full_mat[1:635,] %>% select(any_of(prot_rule$microtrait_rule.name))
PROTEIN_FIRE = as_data_frame(cbind(full_mat[1:635,][c("guild","id","size")],PROTEIN_FIRE))

PROTEIN_FIRE_m = PROTEIN_FIRE %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
PROTEIN_FIRE_m = PROTEIN_FIRE_m %>% mutate(PROTEIN_total = rowSums(PROTEIN_FIRE_m[,3:ncol(PROTEIN_FIRE_m)]))

ggplot(PROTEIN_FIRE_m, aes(size_mean,PROTEIN_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total Protein costs per FG") # It seems that group 62 is too big

ggplot(PROTEIN_FIRE_m, aes(x = size_mean,y = PROTEIN_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total Protein costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# OSMOLYTE GENE COST ----
OSMOLYTE_LOMA = full_mat[636:1168,] %>% select(any_of(osmo_rule$microtrait_hmm.name))
OSMOLYTE_LOMA = as_data_frame(cbind(full_mat[636:1168,][c("guild","id","size")],OSMOLYTE_LOMA))

OSMOLYTE_LOMA_m = OSMOLYTE_LOMA %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                       list(mean=mean), na.rm=TRUE))
OSMOLYTE_LOMA_m = OSMOLYTE_LOMA_m %>% mutate(OSMOLYTE_total = rowSums(OSMOLYTE_LOMA_m[,3:ncol(OSMOLYTE_LOMA_m)]))

OSMOLYTE_LOMA_m = OSMOLYTE_LOMA_m[-c(44), ]

ggplot(OSMOLYTE_LOMA_m, aes(size_mean,OSMOLYTE_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG") # It seems that group 62 is too big

ggplot(OSMOLYTE_LOMA_m, aes(x = size_mean,y = OSMOLYTE_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean osmolyte costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

OSMOLYTE_FIRE = full_mat[1:635,] %>% select(any_of(osmo_rule$microtrait_hmm.name))
OSMOLYTE_FIRE = as_data_frame(cbind(full_mat[1:635,][c("guild","id","size")],OSMOLYTE_FIRE))

OSMOLYTE_FIRE_m = OSMOLYTE_FIRE %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                       list(mean=mean), na.rm=TRUE))
OSMOLYTE_FIRE_m = OSMOLYTE_FIRE_m %>% mutate(OSMOLYTE_total = rowSums(OSMOLYTE_FIRE_m[,3:ncol(OSMOLYTE_FIRE_m)]))

ggplot(OSMOLYTE_FIRE_m, aes(size_mean,OSMOLYTE_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG") # It seems that group 62 is too big

ggplot(OSMOLYTE_FIRE_m, aes(x = size_mean,y = OSMOLYTE_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# TRANSPORTER GENE COST ----
transp_rule  = transp_rule %>% filter(function.==c("transporter"))
TRANSPORTER_LOMA = full_mat[636:1168,] %>% select(any_of(transp_rule$microtrait_hmm.name))
TRANSPORTER_LOMA = as_data_frame(cbind(full_mat[636:1168,][c("guild","id","size")],TRANSPORTER_LOMA))

TRANSPORTER_LOMA_m = TRANSPORTER_LOMA %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSPORTER_LOMA_m = TRANSPORTER_LOMA_m %>% mutate(TRANSPORTER_total = rowSums(TRANSPORTER_LOMA_m[,3:ncol(TRANSPORTER_LOMA_m)]))

TRANSPORTER_LOMA_m = TRANSPORTER_LOMA_m[-c(44), ]

ggplot(TRANSPORTER_LOMA_m, aes(size_mean,TRANSPORTER_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total transporter costs per FG") # It seems that group 62 is too big

ggplot(TRANSPORTER_LOMA_m, aes(x = size_mean,y = TRANSPORTER_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean transporter costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 85) + stat_regline_equation(label.y = 80)

TRANSPORTER_FIRE = full_mat[1:635,] %>% select(any_of(transp_rule$microtrait_hmm.name))
TRANSPORTER_FIRE = as_data_frame(cbind(full_mat[1:635,][c("guild","id","size")],TRANSPORTER_FIRE))

TRANSPORTER_FIRE_m = TRANSPORTER_FIRE %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                         list(mean=mean), na.rm=TRUE))
TRANSPORTER_FIRE_m = TRANSPORTER_FIRE_m %>% mutate(TRANSPORTER_total = rowSums(TRANSPORTER_FIRE_m[,3:ncol(TRANSPORTER_FIRE_m)]))

ggplot(TRANSPORTER_FIRE_m, aes(size_mean,TRANSPORTER_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total transporter costs per FG") # It seems that group 62 is too big

ggplot(TRANSPORTER_FIRE_m, aes(x = size_mean,y = TRANSPORTER_total)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total transporter costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 85) + stat_regline_equation(label.y = 80)


# GC content relationship ----

fire_metadata = read_excel("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/MAG_Dataset_BurnSeverity_ARNelson.xlsx")
fire_metadata = fire_metadata[-c(440,546), ]
fire_stat     =as.data.frame(cbind(fire_stat,fire_metadata$GC))
colnames(fire_stat) = c("id","size","CG")

GH_FIRE = full_mat[1:635,] %>% select(any_of(GH_rule$microtrait_hmm.name))
GH_FIRE = as_data_frame(cbind(full_mat[1:635,][c("guild","id","size")],GH_FIRE))
GH_FIRE = left_join(GH_FIRE, fire_stat, by=c('id',"size"))
GH_FIRE = GH_FIRE %>% select(CG, everything())

GH_FIRE_m = GH_FIRE %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
GH_FIRE_m = GH_FIRE_m %>% select(`fire_metadata$GC`, everything())
GH_FIRE_m = GH_FIRE_m %>% mutate(GH_total = rowSums(GH_FIRE_m[,4:ncol(GH_FIRE_m)]))

ggplot(GH_FIRE_m, aes(CG_mean,GH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean GC content") + ylab("Mean total GH costs per FG")

ggplot(GH_FIRE_m, aes(x = CG_mean,y = GH_total)) + geom_point() +
  xlab("Mean GC content") + ylab("Mean total GH costs per FG") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

ggplot(GH_FIRE_m, aes(x = CG_mean,y = (size_mean))) + geom_point() +
  xlab("Mean GC content") + ylab("Mean genome size") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2)) +
  geom_point() + stat_cor(label.y = 10e6) + stat_regline_equation(label.y = 0)

ggplot(GH_FIRE, aes(x = CG,y = size)) + geom_point() +
  xlab("GC content") + ylab("Genome size") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2)) +
  geom_point() + stat_cor(label.y = 10e6) + stat_regline_equation(label.y = 30)

# Isolates ----

# Calling data
isolates      = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/hmm_isolates.csv",dec=".")
df_1          = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils_unpublished.tsv",sep="\t")
df_1          = subset(df_1, select = -c(23))
df_2          = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils.tsv",sep="\t")
df_2          = subset(df_2, select = -c(14,16,18))
iso.metadata  = as.data.frame(rbind(df_1,df_2))
iso.metadata  = filter(iso.metadata, High.Quality %in% c("Yes"))

# Clean data
test          = as.data.frame(isolates[3:2300])
sums          = as.data.frame(rowSums(test, na.rm = FALSE, dims = 1))
temp          = cbind(sums,isolates$id,1:nrow(isolates))
isolates      = isolates[isolates$id %in% iso.metadata$taxon_oid, ]

# Other data
hmm_loma      = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/hmm_Loma_full.csv",dec=".")
gene_loma     = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/litter_mags_metadata.txt",dec=".")
mag_stat      = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_stats.txt") 
mag_abun.loma = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_adundance.txt") 
mag_stat      = mag_stat %>% full_join(mag_abun.loma)
hmm_fire      = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/hmm_Fire_full.csv",dec=".")
mag_abun.fire = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/mag_adundance_fire.txt") 
mag_abun.fire = mag_abun.fire[-c(440,546), ]

# Selecting genes based on Loma and Fire examples
mat.gene.loma = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/loma.genes.selected.csv",dec=".")
mat.gene.fire = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fire.genes.selected.csv",dec=".")
T.gene.select = as.data.frame(unique(c(colnames(mat.gene.loma[3:293]),
                                       colnames(mat.gene.fire[3:310]))))
colnames(T.gene.select)      = c("gene")
temp.fire                    = as.data.frame(cbind(hmm_fire$id,hmm_fire %>% select(T.gene.select$gene)))
temp.loma                    = as.data.frame(cbind(hmm_loma$id,hmm_loma %>% select(T.gene.select$gene)))
colnames(temp.fire)[1]       = c("id")
colnames(temp.loma)[1]       = c("id")
temp.isolates                = as.data.frame(cbind(isolates$id,isolates %>% select(T.gene.select$gene)))
colnames(temp.isolates)[1]   = c("id")
total.mags                   = as.data.frame(rbind(temp.fire,temp.loma,temp.isolates))

# write.csv(total.mags, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/total.mags.csv")
total.mags                   = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/total.mags.csv",dec=".")

# Grouping (I am using the preexisting functional groups from Loma and the MAGs)
set.seed(1)
distance.total  = vegdist(total.mags[3:507], method = "jaccard", binary = TRUE)
cluster.total   = hclust(distance.total, method="ward.D2")

# similarity within guilds
my_vec       = c()
nguilds.1    = seq(2, nrow(total.mags), 2)

for(i in nguilds.1) {
  v                      = cutree(cluster.total,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distance.total ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1
  my_vec   <- c(my_vec, my_out)   
}
mat_r2_1  = as.data.frame(cbind(nguilds.1,my_vec))

# Similarity among guilds ####
v                      = cutree(cluster.total,k=148)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
test.clus.1            = adonis.pair(distance.total, genome2guild[,"guild"], nper = 1, 
                                     corr.method = "fdr")
adonis_2               = vegan::adonis2(distance.total ~ guild, data = genome2guild, perm = 1)
mat.gene.total         = as.data.frame(cbind(genome2guild,total.mags))
# write.csv(mat.gene.total, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/mat.gene.total.csv")
mat.gene.total         = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/mat.gene.total.csv",dec=".")

# Select isolates ####
mat.gene.isolates      = mat.gene.total[1169:6469,]
iso.metadata.1         = iso.metadata[iso.metadata$taxon_oid %in% isolates$id, ]
names(iso.metadata.1)[names(iso.metadata.1) == 'taxon_oid'] <- 'id'
iso.metadata.gsize     = iso.metadata.1 %>% select(id,Genome.Size....assembled)
names(iso.metadata.gsize)[names(iso.metadata.gsize) == 'Genome.Size....assembled'] <- 'Genome.Size'
mat.gene.isolates.1    = merge(mat.gene.isolates, iso.metadata.gsize, by="id")

# CAZy gene cost per functional group of isolates ----

GH_rule     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_GH.txt", header = TRUE)
GH_TOTAL    = mat.gene.isolates.1 %>% select(any_of(GH_rule$microtrait_hmm.name))
GH_TOTAL.1  = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                          GH_TOTAL))
GH_TOTAL_m  = GH_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                               list(mean=mean), na.rm=TRUE))
GH_TOTAL_m  = GH_TOTAL_m %>% mutate(GH_total = rowSums(GH_TOTAL_m[,3:ncol(GH_TOTAL_m)]))

ggplot(GH_TOTAL_m, aes(Genome.Size_mean,GH_total,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total GH costs per FG") 

ggplot(GH_TOTAL_m, aes(x = Genome.Size_mean,y = GH_total)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean total GH costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# Protein gene cost per functional group of isolates ----

prot_rule   = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_proteins.txt", header = TRUE)
PT_TOTAL    = mat.gene.isolates.1 %>% select(any_of(prot_rule$microtrait_rule.name))
PT_TOTAL.1  = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                  PT_TOTAL))
PT_TOTAL_m  = PT_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                  list(mean=mean), na.rm=TRUE))
PT_TOTAL_m  = PT_TOTAL_m %>% mutate(PT_TOTAL = rowSums(PT_TOTAL_m[,3:ncol(PT_TOTAL_m)]))

ggplot(PT_TOTAL_m, aes(Genome.Size_mean,PT_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total protein costs per FG") 

ggplot(PT_TOTAL_m, aes(x = Genome.Size_mean,y = PT_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean total protein costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# Osmolyte gene cost per functional group of isolates ----

osmo_rule   = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_osymolites.txt", header = TRUE)
OS_TOTAL    = mat.gene.isolates.1 %>% select(any_of(osmo_rule$microtrait_hmm.name))
OS_TOTAL.1  = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                  OS_TOTAL))
OS_TOTAL_m  = OS_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                  list(mean=mean), na.rm=TRUE))
OS_TOTAL_m  = OS_TOTAL_m %>% mutate(OS_TOTAL = rowSums(OS_TOTAL_m[,3:ncol(OS_TOTAL_m)]))

ggplot(OS_TOTAL_m, aes(Genome.Size_mean,OS_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total osmolyte costs per FG") 

ggplot(OS_TOTAL_m, aes(x = Genome.Size_mean,y = OS_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean total osmolyte costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# Biofilm gene cost per functional group of isolates ----

biofilm_rule = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_biofilm.txt", header = TRUE)
BI_TOTAL     = mat.gene.isolates.1 %>% select(any_of(biofilm_rule$microtrait_hmm.name))
BI_TOTAL.1   = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                   BI_TOTAL))
BI_TOTAL_m   = BI_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
BI_TOTAL_m   = BI_TOTAL_m %>% mutate(BI_TOTAL = rowSums(BI_TOTAL_m[,3:ncol(BI_TOTAL_m)]))

ggplot(BI_TOTAL_m, aes(Genome.Size_mean,BI_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total biofilm costs per FG") 

ggplot(BI_TOTAL_m, aes(x = Genome.Size_mean,y = BI_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean total biofilm costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 10) + stat_regline_equation(label.y = 9)

# High tolerance gene cost per functional group of isolates ----

high.T_rule  = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_high_T.txt", header = TRUE)
HT_TOTAL     = mat.gene.isolates.1 %>% select(any_of(high.T_rule$microtrait_hmm.name))
HT_TOTAL.1   = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                   HT_TOTAL))
HT_TOTAL_m   = HT_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
HT_TOTAL_m   = HT_TOTAL_m %>% mutate(HT_TOTAL = rowSums(HT_TOTAL_m[,3:ncol(HT_TOTAL_m)]))

ggplot(HT_TOTAL_m, aes(Genome.Size_mean,HT_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total high temperature gene costs per FG") 

ggplot(HT_TOTAL_m, aes(x = Genome.Size_mean,y = HT_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean total high temperature gene costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 10) + stat_regline_equation(label.y = 9)

# Low tolerance gene cost per functional group of isolates ----

low.T_rule   = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_low_T.txt", header = TRUE)
LT_TOTAL     = mat.gene.isolates.1 %>% select(any_of(low.T_rule$microtrait_hmm.name))
LT_TOTAL.1   = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                   LT_TOTAL))
LT_TOTAL_m   = LT_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
LT_TOTAL_m   = LT_TOTAL_m %>% mutate(LT_TOTAL = rowSums(LT_TOTAL_m[,3:ncol(LT_TOTAL_m)]))

ggplot(LT_TOTAL_m, aes(Genome.Size_mean,LT_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total low temperature gene costs per FG") 

ggplot(LT_TOTAL_m, aes(x = Genome.Size_mean,y = LT_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean total low temperature gene costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 6) + stat_regline_equation(label.y = 5)

# pH gene cost per functional group of isolates ----

pH_rule      = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_pH_stress.txt", header = TRUE)
pH_TOTAL     = mat.gene.isolates.1 %>% select(any_of(pH_rule$microtrait_hmm.name))
pH_TOTAL.1   = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                   pH_TOTAL))
pH_TOTAL_m   = pH_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
pH_TOTAL_m   = pH_TOTAL_m %>% mutate(pH_TOTAL = rowSums(pH_TOTAL_m[,3:ncol(pH_TOTAL_m)]))

ggplot(pH_TOTAL_m, aes(Genome.Size_mean,pH_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total pH gene costs per FG") 

ggplot(pH_TOTAL_m, aes(x = Genome.Size_mean,y = pH_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean total pH gene costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 10) + stat_regline_equation(label.y = 9)

# Transporters gene cost per functional group of isolates ----

transp_rule  = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/microtrait_transporters.txt", header = TRUE,sep = "\t")
transp_rule  = transp_rule %>% filter(function.==c("transporter"))
Ts_TOTAL     = mat.gene.isolates.1 %>% select(any_of(transp_rule$microtrait_hmm.name))
Ts_TOTAL.1   = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],
                                   Ts_TOTAL))
Ts_TOTAL_m   = Ts_TOTAL.1 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
Ts_TOTAL_m   = Ts_TOTAL_m %>% mutate(Ts_TOTAL = rowSums(Ts_TOTAL_m[,3:ncol(Ts_TOTAL_m)]))

ggplot(Ts_TOTAL_m, aes(Genome.Size_mean,Ts_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean total transporters gene costs per FG") 

ggplot(Ts_TOTAL_m, aes(x = Genome.Size_mean,y = Ts_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG (isolates)") + 
  ylab("Mean transporters gene costs per FG (isolates)") + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 100) + stat_regline_equation(label.y = 90)

# AMINOACIDS ----

transp_rule_ami  = transp_rule %>% filter(class==c("aminoacid","peptide"))
AMI_TRANSP_TOTAL = mat.gene.isolates.1 %>% select(any_of(transp_rule_ami$microtrait_hmm.name))
AMI_TRANSP_TOTAL = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],AMI_TRANSP_TOTAL))

AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                               list(mean=mean), na.rm=TRUE))
AMI_TRANSP_TOTAL_m = AMI_TRANSP_TOTAL_m %>% mutate(AMI_TRANSP_TOTAL = rowSums(AMI_TRANSP_TOTAL_m[,3:ncol(AMI_TRANSP_TOTAL_m)]))

ggplot(AMI_TRANSP_TOTAL_m, aes(Genome.Size_mean,AMI_TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean aminoacid transporter costs per FG") # It seems that group 62 is too big

ggplot(AMI_TRANSP_TOTAL_m, aes(x = Genome.Size_mean,y = AMI_TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean aminoacid transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)

# CARBOHYDRATE ----
transp_rule_carb = transp_rule %>% filter(class==c("carbohydrate"))
CARB_TRANSP_TOTAL = mat.gene.isolates.1 %>% select(any_of(transp_rule_carb$microtrait_hmm.name))
CARB_TRANSP_TOTAL = as_data_frame(cbind(mat.gene.isolates.1[c("guild","id","Genome.Size")],CARB_TRANSP_TOTAL))

CARB_TRANSP_TOTAL_m = CARB_TRANSP_TOTAL %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                                 list(mean=mean), na.rm=TRUE))
CARB_TRANSP_TOTAL_m = CARB_TRANSP_TOTAL_m %>% mutate(CARB_TRANSP_TOTAL = rowSums(CARB_TRANSP_TOTAL_m[,3:ncol(CARB_TRANSP_TOTAL_m)]))

ggplot(CARB_TRANSP_TOTAL_m, aes(Genome.Size_mean,CARB_TRANSP_TOTAL,colour=factor(guild))) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean carbohydrate transporter costs per FG") # It seems that group 62 is too big

ggplot(CARB_TRANSP_TOTAL_m, aes(x = Genome.Size_mean,y = CARB_TRANSP_TOTAL)) + geom_point() +
  xlab("Mean genome size per FG") + ylab("Mean carbohydrate transporter costs per FG") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point() + stat_cor(label.y = 35) + stat_regline_equation(label.y = 30)