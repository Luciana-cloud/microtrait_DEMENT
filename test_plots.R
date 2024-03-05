library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

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
GH_rule             = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/microtrait_GH.txt", header = TRUE)
prot_rule           = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/microtrait_proteins.txt", header = TRUE)
osmo_rule           = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/microtrait_osymolites.txt", header = TRUE)

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


