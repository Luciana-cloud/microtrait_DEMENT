# Packages ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stats)
library(corrr)
# devtools::install_github("rsquaredacademy/olsrr")
library(olsrr)

# Calling data ####

mat_ori    = read.csv("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_trait_matrixatgranularity.csv",dec=".")
gen_size   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/litter_mags_metadata.txt",dec=".") 

mag_stat   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/mag_stats.txt") 
mag_abun   = read.delim("C:/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/mag_adundance.txt") 

mag_stat   = mag_stat %>% full_join(mag_abun)

# Erase weird MAGs, which were found in preliminary analysis: 116, 175

mat_ori    = mat_ori[-c(116,175),] 
gen_size   = gen_size[-c(116,175),] 
mag_stat   = mag_stat[-c(116,175),]

mat_trait    = mat_ori %>% select(2:191)/gen_size$length
mat_trait    = as.data.frame(cbind(mat_ori$id,mat_trait))
mat_trait    = mat_trait %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))

# Objective 1: Traits responsible for changes in abundance ####

# MAG abundance Analysis  ####

# Grass-ambient vs grass-drought  ####

# Changes in MAGs abundance ####

mag_stat = mag_stat %>% mutate(ratio1 = log10(mag_stat$Average.1/mag_stat$Average))
mag_stat = mag_stat %>% mutate(pos = ratio1 >= 0)
mag_stat = mag_stat %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p<-ggplot(data=mag_stat, aes(x=id_2, y=ratio1, fill = pos)) +
   geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
   xlab("MAG ID")
p

# Changes in MAGs weighted abundance ####

mag_stat = mag_stat %>% mutate(ratio1.1 = log10(mag_stat$Average.1/mag_stat$Average)*
                                 abs(mag_stat$Average.1-mag_stat$Average))
mag_stat = mag_stat %>% mutate(pos.1 = ratio1.1 >= 0)
mag_stat = mag_stat %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p<-ggplot(data=mag_stat, aes(x=id_2, y=ratio1.1, fill = pos.1)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p

# Correlation results with the highest granularity ####

mat_trait   = mat_trait %>% mutate(rat_am_dro_gra = mag_stat$ratio1)
mat_trait   = mat_trait %>% mutate(rat_am_dro_gra = ifelse(is.na(rat_am_dro_gra),0,rat_am_dro_gra),
                                   rat_am_dro_gra = ifelse(rat_am_dro_gra==Inf,0,rat_am_dro_gra),
                                   rat_am_dro_gra = ifelse(rat_am_dro_gra ==-Inf, 0, rat_am_dro_gra))
mat_trait   = mat_trait %>% mutate(rat_am_dro_gra.1 = mag_stat$ratio1.1)
mat_trait   = mat_trait %>% mutate(rat_am_dro_gra.1 = ifelse(is.na(rat_am_dro_gra.1),0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1==Inf,0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1 ==-Inf, 0, rat_am_dro_gra.1))

mat_trait_1 = as.data.frame(mat_trait[,-192])
mat = (as.matrix(mat_trait_1[2:193]))

corre   = as.data.frame(cor(mat[,1:190], mat[,191]))
corre.1 = as.data.frame(cor(mat[,1:190], mat[,192]))

# Better prediction attempt  ####

test.1 = as.data.frame(cbind(mat_trait_1[2:191]))

# Erase traits based on correlation matrix

seq.1   = seq(1,190)
temp.1  = as.data.frame(cbind((colSums(test.1)),seq.1))
erase.1 = temp.1 %>% filter(V1==0)
erase.1$row_names <- row.names(erase.1)

test.2 = subset(test.1, select =-((erase.1$row_names)))
test.2 = test.1[-(erase.1$seq.1)]
temp.2 = as.data.frame(colSums(test.2))

corre.1.1 = as.data.frame(cor(as.matrix(test.2), mat[,192]))

m.all = lm(mat[,192]~.,data=test.2[1:50])
temp.1.50 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,192]~.,data=test.2[51:100])
temp.51.100 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,192]~.,data=test.2[101:115])
temp.101.115 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,192]~.,data=test.2[117:127])
temp.117.127 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,192]~.,data=test.2[128:140])
temp.128.140 = ols_step_forward_p(m.all, details = FALSE)

erase.2 = as.data.frame(c(temp.1.50[["metrics"]][["variable"]],
                          temp.51.100[["metrics"]][["variable"]],
                          temp.101.115[["metrics"]][["variable"]],
                          temp.117.127[["metrics"]][["variable"]],
                          temp.128.140[["metrics"]][["variable"]]))
colnames(erase.2) = c("trait")

test.3 = test.2 %>% select((erase.2$trait))
corre.1.2 = as.data.frame(cor(as.matrix(test.3), mat[,192]))
test.3 = test.3[-59]
m.all.1 = lm(mat[,192]~.,data=test.3[1:55])
test.3  = test.3[-55]
m.all.1 = lm(mat[,192]~.,data=test.3)
test.3  = test.3[-58]
m.all.1 = lm(mat[,192]~.,data=test.3)
temp.58 = ols_step_forward_p(m.all.1, details = FALSE)

erase.3 = as.data.frame(temp.58[["metrics"]][["variable"]])
colnames(erase.3) = c("trait")

test.4 = test.3 %>% select((erase.3$trait))
write.csv(test.4, file = "test.4.csv")
corre.1.3 = as.data.frame(cor(as.matrix(test.4), mat[,192]))
m.all.2 = lm(mat[,192]~.,data=test.4)
temp.37 = ols_step_forward_p(m.all.2, details = FALSE)

mat.1 = (as.data.frame(cbind(mat_trait_1$rat_am_dro_gra.1,test.4)))

# Leaps test ####

x = model.matrix(mat.1$`mat_trait_1$rat_am_dro_gra.1`~.-1,data=mat.1[1:30])
y = mat.1$`mat_trait_1$rat_am_dro_gra.1`
bestmods = leaps(x,y,nbest=1)

#a = as.data.frame(colnames(test.3))

#best.mod = ols_step_best_subset(model=m.all.2)
#str(best.mod)
#plot(best.mod)

#k.1  <- ols_step_all_possible(m.all)
#k.1
#pdf("k.1.pdf")
#plot(k.1)
#dev.off() 

# Correlation plots  ####

x <- mat.1 %>% 
  correlate() %>% 
  focus(`mat[, 192]`)

x %>% 
  mutate(rowname = factor(term, levels = term[order(`mat[, 192]`)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = `mat[, 192]`)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_gra") +
  xlab("")

# Correlation results with the lowest granularity ####

resource_acquisition1 = rowSums(mat_trait_1 %>% select(2:84), na.rm=FALSE)/83
resource_use1         = rowSums(mat_trait_1 %>% select(85:161), na.rm=FALSE)/77
stress_tolerance1     = rowSums(mat_trait_1 %>% select(162:190), na.rm=FALSE)/29

mat_trait_1_low       = cbind(resource_acquisition1,resource_use1,stress_tolerance1,
                              mat_trait_1 %>% select(191:192))

mat_low = (as.matrix(mat_trait_1_low))

x_low <- mat_low %>% 
  correlate() %>% 
  focus(rat_am_dro_gra)

x_low %>% 
  mutate(rowname = factor(term, levels = term[order(rat_am_dro_gra)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_am_dro_gra)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_gra") +
  xlab("")

# Correlation results with the middle granularity ####

RA_substrate_uptake       = rowSums(mat_trait_1 %>% select(2:42), na.rm=FALSE)/41
RA_substrate_degradation  = rowSums(mat_trait_1 %>% select(43:61), na.rm=FALSE)/19
RA_substrate_assimilation = rowSums(mat_trait_1 %>% select(62:84), na.rm=FALSE)/23
RU_chemotrophy            = rowSums(mat_trait_1 %>% select(85:131), na.rm=FALSE)/47
RU_phototrophy            = rowSums(mat_trait_1 %>% select(132:161), na.rm=FALSE)/30
ST_general                = rowSums(mat_trait_1 %>% select(162:165), na.rm=FALSE)/4
ST_specific               = rowSums(mat_trait_1 %>% select(166:190), na.rm=FALSE)/25

mat_trait_1_mid        = cbind(RA_substrate_uptake,RA_substrate_degradation,
                               RA_substrate_assimilation,RU_chemotrophy,RU_phototrophy,
                               ST_general,ST_specific,mat_trait_1 %>% select(191:192))

mat_mid = (as.matrix(mat_trait_1_mid))

x_mid <- mat_mid %>% 
  correlate() %>% 
  focus(rat_am_dro_gra)

x_mid %>% 
  mutate(rowname = factor(term, levels = term[order(rat_am_dro_gra)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_am_dro_gra)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_gra") +
  xlab("")

# Correlation results with the middle 1 granularity ####

RA_uptake_aromatic       = rowSums(mat_trait_1 %>% select(2), na.rm=FALSE)/1
RA_uptake_biopolymer     = rowSums(mat_trait_1 %>% select(3), na.rm=FALSE)/1
RA_uptake_carbohydrate   = rowSums(mat_trait_1 %>% select(4:14), na.rm=FALSE)/11
RA_uptake_carboxylate    = rowSums(mat_trait_1 %>% select(15:17), na.rm=FALSE)/3
RA_uptake_aminoacids     = rowSums(mat_trait_1 %>% select(18), na.rm=FALSE)/1
RA_uptake_ions           = rowSums(mat_trait_1 %>% select(19:20), na.rm=FALSE)/2
RA_uptake_lipid          = rowSums(mat_trait_1 %>% select(21:23), na.rm=FALSE)/3
RA_uptake_N_compound     = rowSums(mat_trait_1 %>% select(24:30), na.rm=FALSE)/7
RA_uptake_nucleid_acid   = rowSums(mat_trait_1 %>% select(31:34), na.rm=FALSE)/4
RA_uptake_organoP        = rowSums(mat_trait_1 %>% select(35), na.rm=FALSE)/1
RA_uptake_osmolyte       = rowSums(mat_trait_1 %>% select(36), na.rm=FALSE)/1
RA_uptake_other          = rowSums(mat_trait_1 %>% select(37), na.rm=FALSE)/1
RA_uptake_peptide        = rowSums(mat_trait_1 %>% select(38), na.rm=FALSE)/1
RA_uptake_S_compound     = rowSums(mat_trait_1 %>% select(39), na.rm=FALSE)/1
RA_uptake_siderophore    = rowSums(mat_trait_1 %>% select(40), na.rm=FALSE)/1
RA_uptake_vitamin        = rowSums(mat_trait_1 %>% select(41:42), na.rm=FALSE)/2
RA_degradation_complex   = rowSums(mat_trait_1 %>% select(43:48), na.rm=FALSE)/6
RA_degradation_simple    = rowSums(mat_trait_1 %>% select(49:61), na.rm=FALSE)/13
RA_assimilation_C_comp   = rowSums(mat_trait_1 %>% select(62:76), na.rm=FALSE)/15
RA_assimilation_N_comp   = rowSums(mat_trait_1 %>% select(77:81), na.rm=FALSE)/5
RA_assimilation_S_comp   = rowSums(mat_trait_1 %>% select(82), na.rm=FALSE)/1
RA_assimilation_P_comp   = rowSums(mat_trait_1 %>% select(83:84), na.rm=FALSE)/2
RU_chemotrophy            = rowSums(mat_trait_1 %>% select(85:131), na.rm=FALSE)/47
RU_phototrophy            = rowSums(mat_trait_1 %>% select(132:161), na.rm=FALSE)/30
ST_general                = rowSums(mat_trait_1 %>% select(162:165), na.rm=FALSE)/4
ST_specific_high_T        = rowSums(mat_trait_1 %>% select(166:168), na.rm=FALSE)/3
ST_specific_low_T         = rowSums(mat_trait_1 %>% select(169:173), na.rm=FALSE)/5
ST_specific_desiccation   = rowSums(mat_trait_1 %>% select(174:177), na.rm=FALSE)/4
ST_specific_pH_stress     = rowSums(mat_trait_1 %>% select(178:184), na.rm=FALSE)/7
ST_specific_oxidative     = rowSums(mat_trait_1 %>% select(185:186), na.rm=FALSE)/2
ST_specific_O_limit       = rowSums(mat_trait_1 %>% select(187:188), na.rm=FALSE)/2
ST_specific_envelope      = rowSums(mat_trait_1 %>% select(189:190), na.rm=FALSE)/2

mat_trait_1_mid_1         = cbind(RA_uptake_aromatic,RA_uptake_biopolymer,
                               RA_uptake_carbohydrate,RA_uptake_carboxylate,
                               RA_uptake_aminoacids,RA_uptake_ions,RA_uptake_lipid,
                               RA_uptake_N_compound,RA_uptake_nucleid_acid,
                               RA_uptake_organoP,RA_uptake_osmolyte,RA_uptake_other,
                               RA_uptake_peptide,RA_uptake_S_compound,
                               RA_uptake_siderophore,RA_uptake_vitamin,
                               RA_degradation_complex,RA_degradation_simple,
                               RA_assimilation_C_comp,RA_assimilation_N_comp,
                               RA_assimilation_S_comp,RA_assimilation_P_comp,
                               RU_chemotrophy,RU_phototrophy,ST_general,
                               ST_specific_high_T,ST_specific_low_T,
                               ST_specific_desiccation,ST_specific_pH_stress,
                               ST_specific_oxidative,ST_specific_O_limit,
                               ST_specific_envelope,mat_trait_1 %>% select(191:192))

mat_mid_1 = (as.matrix(mat_trait_1_mid_1))

x_mid_1 <- mat_mid_1 %>% 
  correlate() %>% 
  focus(rat_am_dro_gra)

x_mid_1 %>% 
  mutate(rowname = factor(term, levels = term[order(rat_am_dro_gra)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_am_dro_gra)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_gra") +
  xlab("")

# Shrub-ambient vs shrub-drought  ####

# Changes in MAGs abundance

mag_stat.1 = mag_stat %>% mutate(ratio1 = log10(mag_stat$Average.3/mag_stat$Average.2))
mag_stat.1 = mag_stat.1 %>% mutate(pos = ratio1 >= 0)
mag_stat.1 = mag_stat.1 %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p<-ggplot(data=mag_stat.1, aes(x=id_2, y=ratio1, fill = pos)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p
# Correlation plots

mat_traits   = mat_trait %>% mutate(rat_am_dro_shr = mag_stat.1$ratio1)
mat_traits   = mat_traits %>% mutate(rat_am_dro_shr = ifelse(is.na(rat_am_dro_shr),0,rat_am_dro_shr),
                                   rat_am_dro_shr = ifelse(rat_am_dro_shr==Inf,0,rat_am_dro_shr),
                                   rat_am_dro_shr = ifelse(rat_am_dro_shr ==-Inf, 0, rat_am_dro_shr))
mat_trait_1s = as.data.frame(mat_traits[,-192])
mat_trait_1s = as.data.frame(mat_trait_1s[,-192])
mats = (as.matrix(mat_trait_1s[2:192]))

corre1 = as.data.frame(cor(mats[,1:190], mats[,191]))

# Correlation plots  ####

x1 <- mats %>% 
  correlate() %>% 
  focus(rat_am_dro_shr)

x1 %>% 
  mutate(rowname = factor(term, levels = term[order(rat_am_dro_shr)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_am_dro_shr)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_shr") +
  xlab("")

# Correlation results with the lowest granularity ####

resource_acquisition1 = rowSums(mat_trait_1s %>% select(2:84), na.rm=FALSE)/83
resource_use1         = rowSums(mat_trait_1s %>% select(85:161), na.rm=FALSE)/77
stress_tolerance1     = rowSums(mat_trait_1s %>% select(162:190), na.rm=FALSE)/29

mat_trait_1s_low      = cbind(resource_acquisition1,resource_use1,stress_tolerance1,
                              mat_trait_1s %>% select(191:192))

mat_low_s = (as.matrix(mat_trait_1s_low))

x_low_s <- mat_low_s %>% 
  correlate() %>% 
  focus(rat_am_dro_shr)

x_low_s %>% 
  mutate(rowname = factor(term, levels = term[order(rat_am_dro_shr)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_am_dro_shr)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_shr") +
  xlab("")

# Correlation results with the middle granularity ####

RA_substrate_uptake       = rowSums(mat_trait_1s %>% select(2:42), na.rm=FALSE)/41
RA_substrate_degradation  = rowSums(mat_trait_1s %>% select(43:61), na.rm=FALSE)/19
RA_substrate_assimilation = rowSums(mat_trait_1s %>% select(62:84), na.rm=FALSE)/23
RU_chemotrophy            = rowSums(mat_trait_1s %>% select(85:131), na.rm=FALSE)/47
RU_phototrophy            = rowSums(mat_trait_1s %>% select(132:161), na.rm=FALSE)/30
ST_general                = rowSums(mat_trait_1s %>% select(162:165), na.rm=FALSE)/4
ST_specific               = rowSums(mat_trait_1s %>% select(166:190), na.rm=FALSE)/25

mat_trait_1s_mid          = cbind(RA_substrate_uptake,RA_substrate_degradation,
                               RA_substrate_assimilation,RU_chemotrophy,RU_phototrophy,
                               ST_general,ST_specific,mat_trait_1s %>% select(191:192))

mat_mid_1s = (as.matrix(mat_trait_1s_mid))

x_mid_1s <- mat_mid_1s %>% 
  correlate() %>% 
  focus(rat_am_dro_shr)

x_mid_1s %>% 
  mutate(rowname = factor(term, levels = term[order(rat_am_dro_shr)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_am_dro_shr)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_shr") +
  xlab("")

# Correlation results with the middle 1 granularity ####

RA_uptake_aromatic       = rowSums(mat_trait_1s %>% select(2), na.rm=FALSE)/1
RA_uptake_biopolymer     = rowSums(mat_trait_1s %>% select(3), na.rm=FALSE)/1
RA_uptake_carbohydrate   = rowSums(mat_trait_1s %>% select(4:14), na.rm=FALSE)/11
RA_uptake_carboxylate    = rowSums(mat_trait_1s %>% select(15:17), na.rm=FALSE)/3
RA_uptake_aminoacids     = rowSums(mat_trait_1s %>% select(18), na.rm=FALSE)/1
RA_uptake_ions           = rowSums(mat_trait_1s %>% select(19:20), na.rm=FALSE)/2
RA_uptake_lipid          = rowSums(mat_trait_1s %>% select(21:23), na.rm=FALSE)/3
RA_uptake_N_compound     = rowSums(mat_trait_1s %>% select(24:30), na.rm=FALSE)/7
RA_uptake_nucleid_acid   = rowSums(mat_trait_1s %>% select(31:34), na.rm=FALSE)/4
RA_uptake_organoP        = rowSums(mat_trait_1s %>% select(35), na.rm=FALSE)/1
RA_uptake_osmolyte       = rowSums(mat_trait_1s %>% select(36), na.rm=FALSE)/1
RA_uptake_other          = rowSums(mat_trait_1s %>% select(37), na.rm=FALSE)/1
RA_uptake_peptide        = rowSums(mat_trait_1s %>% select(38), na.rm=FALSE)/1
RA_uptake_S_compound     = rowSums(mat_trait_1s %>% select(39), na.rm=FALSE)/1
RA_uptake_siderophore    = rowSums(mat_trait_1s %>% select(40), na.rm=FALSE)/1
RA_uptake_vitamin        = rowSums(mat_trait_1s %>% select(41:42), na.rm=FALSE)/2
RA_degradation_complex   = rowSums(mat_trait_1s %>% select(43:48), na.rm=FALSE)/6
RA_degradation_simple    = rowSums(mat_trait_1s %>% select(49:61), na.rm=FALSE)/13
RA_assimilation_C_comp   = rowSums(mat_trait_1s %>% select(62:76), na.rm=FALSE)/15
RA_assimilation_N_comp   = rowSums(mat_trait_1s %>% select(77:81), na.rm=FALSE)/5
RA_assimilation_S_comp   = rowSums(mat_trait_1s %>% select(82), na.rm=FALSE)/1
RA_assimilation_P_comp   = rowSums(mat_trait_1s %>% select(83:84), na.rm=FALSE)/2
RU_chemotrophy            = rowSums(mat_trait_1s %>% select(85:131), na.rm=FALSE)/47
RU_phototrophy            = rowSums(mat_trait_1s %>% select(132:161), na.rm=FALSE)/30
ST_general                = rowSums(mat_trait_1s %>% select(162:165), na.rm=FALSE)/4
ST_specific_high_T        = rowSums(mat_trait_1s %>% select(166:168), na.rm=FALSE)/3
ST_specific_low_T         = rowSums(mat_trait_1s %>% select(169:173), na.rm=FALSE)/5
ST_specific_desiccation   = rowSums(mat_trait_1s %>% select(174:177), na.rm=FALSE)/4
ST_specific_pH_stress     = rowSums(mat_trait_1s %>% select(178:184), na.rm=FALSE)/7
ST_specific_oxidative     = rowSums(mat_trait_1s %>% select(185:186), na.rm=FALSE)/2
ST_specific_O_limit       = rowSums(mat_trait_1s %>% select(187:188), na.rm=FALSE)/2
ST_specific_envelope      = rowSums(mat_trait_1s %>% select(189:190), na.rm=FALSE)/2

mat_trait_s_mid_1         = cbind(RA_uptake_aromatic,RA_uptake_biopolymer,
                                  RA_uptake_carbohydrate,RA_uptake_carboxylate,
                                  RA_uptake_aminoacids,RA_uptake_ions,RA_uptake_lipid,
                                  RA_uptake_N_compound,RA_uptake_nucleid_acid,
                                  RA_uptake_organoP,RA_uptake_osmolyte,RA_uptake_other,
                                  RA_uptake_peptide,RA_uptake_S_compound,
                                  RA_uptake_siderophore,RA_uptake_vitamin,
                                  RA_degradation_complex,RA_degradation_simple,
                                  RA_assimilation_C_comp,RA_assimilation_N_comp,
                                  RA_assimilation_S_comp,RA_assimilation_P_comp,
                                  RU_chemotrophy,RU_phototrophy,ST_general,
                                  ST_specific_high_T,ST_specific_low_T,
                                  ST_specific_desiccation,ST_specific_pH_stress,
                                  ST_specific_oxidative,ST_specific_O_limit,
                                  ST_specific_envelope,mat_trait_1s %>% select(191:192))

mat_mid_s = (as.matrix(mat_trait_s_mid_1))

x_mid_s <- mat_trait_s_mid_1 %>% 
  correlate() %>% 
  focus(rat_am_dro_shr)

x_mid_s %>% 
  mutate(rowname = factor(term, levels = term[order(rat_am_dro_shr)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_am_dro_shr)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_am_dro_shr") +
  xlab("")

# Grass-ambient vs Shrubland-ambient  ####

# Changes in MAGs abundance

mag_stat.2 = mag_stat %>% mutate(ratio1 = log10(mag_stat$Average.2/mag_stat$Average))
mag_stat.2 = mag_stat.2 %>% mutate(pos = ratio1 >= 0)
mag_stat.2 = mag_stat.2 %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p<-ggplot(data=mag_stat.2, aes(x=id_2, y=ratio1, fill = pos)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p
# Correlation plots

mat_traitsg   = mat_trait %>% mutate(rat_shr_gras = mag_stat.2$ratio1)
mat_traitsg   = mat_traitsg %>% mutate(rat_shr_gras = ifelse(is.na(rat_shr_gras),0,rat_shr_gras),
                                     rat_shr_gras = ifelse(rat_shr_gras==Inf,0,rat_shr_gras),
                                     rat_shr_gras = ifelse(rat_shr_gras ==-Inf, 0, rat_shr_gras))
mat_trait_1gs = as.data.frame(mat_traitsg[,-192])
mat_trait_1gs = as.data.frame(mat_trait_1gs[,-192])
matgs = (as.matrix(mat_trait_1gs[2:192]))

corre2 = as.data.frame(cor(matgs[,1:190], matgs[,191]))


# Correlation plots  ####

x2 <- matgs %>% 
  correlate() %>% 
  focus(rat_shr_gras)

x2 %>% 
  mutate(rowname = factor(term, levels = term[order(rat_shr_gras)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_shr_gras)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_shr_gras") +
  xlab("")

# Correlation results with the lowest granularity ####

resource_acquisition1 = rowSums(mat_trait_1gs %>% select(2:84), na.rm=FALSE)/83
resource_use1         = rowSums(mat_trait_1gs %>% select(85:161), na.rm=FALSE)/77
stress_tolerance1     = rowSums(mat_trait_1gs %>% select(162:190), na.rm=FALSE)/29

mat_trait_2s_low      = cbind(resource_acquisition1,resource_use1,stress_tolerance1,
                              mat_trait_1gs %>% select(191:192))

mat_low_sg = (as.matrix(mat_trait_2s_low))

x_low_sg <- mat_low_sg %>% 
  correlate() %>% 
  focus(rat_shr_gras)

x_low_sg %>% 
  mutate(rowname = factor(term, levels = term[order(rat_shr_gras)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_shr_gras)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_shr_gras") +
  xlab("")

# Correlation results with the middle granularity ####

RA_substrate_uptake       = rowSums(mat_trait_1gs %>% select(2:42), na.rm=FALSE)/41
RA_substrate_degradation  = rowSums(mat_trait_1gs %>% select(43:61), na.rm=FALSE)/19
RA_substrate_assimilation = rowSums(mat_trait_1gs %>% select(62:84), na.rm=FALSE)/23
RU_chemotrophy            = rowSums(mat_trait_1gs %>% select(85:131), na.rm=FALSE)/47
RU_phototrophy            = rowSums(mat_trait_1gs %>% select(132:161), na.rm=FALSE)/30
ST_general                = rowSums(mat_trait_1gs %>% select(162:165), na.rm=FALSE)/4
ST_specific               = rowSums(mat_trait_1gs %>% select(166:190), na.rm=FALSE)/25

mat_trait_sg_mid          = cbind(RA_substrate_uptake,RA_substrate_degradation,
                                  RA_substrate_assimilation,RU_chemotrophy,RU_phototrophy,
                                  ST_general,ST_specific,mat_trait_1gs %>% select(191:192))

mat_mid_sg = (as.matrix(mat_trait_sg_mid))

x_mid_sg <- mat_mid_sg %>% 
  correlate() %>% 
  focus(rat_shr_gras)

x_mid_sg %>% 
  mutate(rowname = factor(term, levels = term[order(rat_shr_gras)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_shr_gras)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_shr_gras") +
  xlab("")

# Correlation results with the middle 1 granularity ####

RA_uptake_aromatic       = rowSums(mat_trait_1gs %>% select(2), na.rm=FALSE)/1
RA_uptake_biopolymer     = rowSums(mat_trait_1gs %>% select(3), na.rm=FALSE)/1
RA_uptake_carbohydrate   = rowSums(mat_trait_1gs %>% select(4:14), na.rm=FALSE)/11
RA_uptake_carboxylate    = rowSums(mat_trait_1gs %>% select(15:17), na.rm=FALSE)/3
RA_uptake_aminoacids     = rowSums(mat_trait_1gs %>% select(18), na.rm=FALSE)/1
RA_uptake_ions           = rowSums(mat_trait_1gs %>% select(19:20), na.rm=FALSE)/2
RA_uptake_lipid          = rowSums(mat_trait_1gs %>% select(21:23), na.rm=FALSE)/3
RA_uptake_N_compound     = rowSums(mat_trait_1gs %>% select(24:30), na.rm=FALSE)/7
RA_uptake_nucleid_acid   = rowSums(mat_trait_1gs %>% select(31:34), na.rm=FALSE)/4
RA_uptake_organoP        = rowSums(mat_trait_1gs %>% select(35), na.rm=FALSE)/1
RA_uptake_osmolyte       = rowSums(mat_trait_1gs %>% select(36), na.rm=FALSE)/1
RA_uptake_other          = rowSums(mat_trait_1gs %>% select(37), na.rm=FALSE)/1
RA_uptake_peptide        = rowSums(mat_trait_1gs %>% select(38), na.rm=FALSE)/1
RA_uptake_S_compound     = rowSums(mat_trait_1gs %>% select(39), na.rm=FALSE)/1
RA_uptake_siderophore    = rowSums(mat_trait_1gs %>% select(40), na.rm=FALSE)/1
RA_uptake_vitamin        = rowSums(mat_trait_1gs %>% select(41:42), na.rm=FALSE)/2
RA_degradation_complex   = rowSums(mat_trait_1gs %>% select(43:48), na.rm=FALSE)/6
RA_degradation_simple    = rowSums(mat_trait_1gs %>% select(49:61), na.rm=FALSE)/13
RA_assimilation_C_comp   = rowSums(mat_trait_1gs %>% select(62:76), na.rm=FALSE)/15
RA_assimilation_N_comp   = rowSums(mat_trait_1gs %>% select(77:81), na.rm=FALSE)/5
RA_assimilation_S_comp   = rowSums(mat_trait_1gs %>% select(82), na.rm=FALSE)/1
RA_assimilation_P_comp   = rowSums(mat_trait_1gs %>% select(83:84), na.rm=FALSE)/2
RU_chemotrophy            = rowSums(mat_trait_1gs %>% select(85:131), na.rm=FALSE)/47
RU_phototrophy            = rowSums(mat_trait_1gs %>% select(132:161), na.rm=FALSE)/30
ST_general                = rowSums(mat_trait_1gs %>% select(162:165), na.rm=FALSE)/4
ST_specific_high_T        = rowSums(mat_trait_1gs %>% select(166:168), na.rm=FALSE)/3
ST_specific_low_T         = rowSums(mat_trait_1gs %>% select(169:173), na.rm=FALSE)/5
ST_specific_desiccation   = rowSums(mat_trait_1gs %>% select(174:177), na.rm=FALSE)/4
ST_specific_pH_stress     = rowSums(mat_trait_1gs %>% select(178:184), na.rm=FALSE)/7
ST_specific_oxidative     = rowSums(mat_trait_1gs %>% select(185:186), na.rm=FALSE)/2
ST_specific_O_limit       = rowSums(mat_trait_1gs %>% select(187:188), na.rm=FALSE)/2
ST_specific_envelope      = rowSums(mat_trait_1gs %>% select(189:190), na.rm=FALSE)/2

mat_trait_sg_mid_1         = cbind(RA_uptake_aromatic,RA_uptake_biopolymer,
                                  RA_uptake_carbohydrate,RA_uptake_carboxylate,
                                  RA_uptake_aminoacids,RA_uptake_ions,RA_uptake_lipid,
                                  RA_uptake_N_compound,RA_uptake_nucleid_acid,
                                  RA_uptake_organoP,RA_uptake_osmolyte,RA_uptake_other,
                                  RA_uptake_peptide,RA_uptake_S_compound,
                                  RA_uptake_siderophore,RA_uptake_vitamin,
                                  RA_degradation_complex,RA_degradation_simple,
                                  RA_assimilation_C_comp,RA_assimilation_N_comp,
                                  RA_assimilation_S_comp,RA_assimilation_P_comp,
                                  RU_chemotrophy,RU_phototrophy,ST_general,
                                  ST_specific_high_T,ST_specific_low_T,
                                  ST_specific_desiccation,ST_specific_pH_stress,
                                  ST_specific_oxidative,ST_specific_O_limit,
                                  ST_specific_envelope,mat_trait_1gs %>% select(191:192))

mat_mid_sg = (as.matrix(mat_trait_sg_mid_1))

x_mid_sg <- mat_trait_sg_mid_1 %>% 
  correlate() %>% 
  focus(rat_shr_gras)

x_mid_sg %>% 
  mutate(rowname = factor(term, levels = term[order(rat_shr_gras)])) %>%  # Order by correlation strength
  ggplot(aes(x = term, y = rat_shr_gras)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with rat_shr_gras") +
  xlab("")

# Objective 2: Grouping based on principal traits  ####






# Objective 3: Traits-tradeoff analysis  ####





###############################################################################
###############################################################################
###############################################################################




###############################################################################
# Normalized data
###############################################################################

mat_trait    = mat_ori %>% select(2:191)/gen_size$length
mat_trait    = as.data.frame(cbind(mat_ori$id,mat_trait))
trait_sd     = mat_trait %>% select(1,43:61) # substrate degradation
trait_tp     = mat_trait %>% select(1:42) # transporters or substrate uptake
trait_st     = mat_trait %>% select(1,162:190) # stress related traits

###############################################################################
# Exploratory plots - gene costs (normalized genes) #
###############################################################################

# complex compounds
boxplot(trait_sd%>% select(2:7,8),ylab="gene cost",
        names=c("cellulose","chitin", "heteromannan", "linkage-glucan", "xylan",
                "xyloglucan","protein"))
# protein + complex compounds
boxplot(trait_sd%>% select(2,6,8),ylab="gene cost",
        names=c("cellulose","xylan","protein"))
# DEMENT-related transporters
boxplot(trait_tp%>% select(8,4,7,18,24:27,31:35,38),ylab="gene cost",
        names=c("monosacc","amino sugar","carb_phos",
                "f_aminoac","amine","ammonium","nitrate","nitrite",
                "nucleobase","nucleoside","nucleotide","ribonucle",
                "organP","peptide"))
# DEMENT-related stress-related traits
boxplot(trait_st%>% select(13:16,28),ylab="gene cost",
        names=c("sol_transp","sol_synt","EPS biosyn (S)","osmo_sensing",
                "misfolded proteins"))
# Heat-related traits
boxplot(trait_st%>% select(5:7),ylab="gene cost",
        names=c("heat shock proteins","ATP proteases","transcription factor"))

# Resource - Use
boxplot(mat_trait%>% select(93:99,125),ylab="gene cost",
        names=c("oxi.pentose","nonoxi.pentose","glyoxylate.cycle",
                "TCA.cycle","ETC.complex.I","ETC.complex.II","ETC.complex.III",
                "nitrite.oxidation"))

###############################################################################
# Targeted traits for DEMENT modeling
###############################################################################

trait_tar = mat_trait %>% select(4,6:8,18,23:25,31:35,38,43:49,174:177,93:99,125) # 
genomes   = as.list(mat_trait[1])
test      = trait_tar %>% summarize_if(is.numeric, sum, na.rm=TRUE)

###############################################################################
# Forming functional guilds 
# (https://www.youtube.com/watch?v=GPOUGpF-Sno)
# https://github.com/ukaraoz/microtrait
###############################################################################

library(vegan)
set.seed(1)

# Gene cost
###############################################################################
distanceg  = vegdist(trait_tar, method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar), 2)
plot(clusterg)

# Similiarity within guilds
###############################################################################
my_vec = c()

for(i in nguilds) {
  v                      = cutree(clusterg,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distanceg ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1        
  my_vec <- c(my_vec, my_out)    
}

# Plotting
mat_r2_c  = as.data.frame(cbind(nguilds,my_vec))
ggplot(data=mat_r2_c,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 20)) + 
  geom_hline(yintercept=0.8349532, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 100, linetype="dashed", color = "red")

###############################################################################
# Finding the optimal number of guilds
# https://www.r-bloggers.com/2019/01/10-tips-for-choosing-the-optimal-number-of-clusters/
# https://r-tastic.co.uk/post/optimal-number-of-clusters/
###############################################################################

# Manual Test 
# Silhouette width : average dissimilarity between and object and all other 
# objects within a cluster to which it belongs
# https://www.youtube.com/watch?v=ezn6bFz-Phk&t=2s (19:00)
###############################################################################
clusterg_2 = clusterg
Si = numeric(nrow(trait_tar))
for (k in 2:(nrow(trait_tar) - 1))
{
  sil = silhouette(cutree(clusterg_2,k=k),distanceg)
  Si[k] = summary(sil)$avg.width
}

k.best = which.max(Si)
plot(1:nrow(trait_tar),Si,type="h",main="Silhouette-optimal number of clusters",
     xlab = "Number of clusters",ylab = "Average silhouette width")
axis(1,k.best,paste("optimum",k.best,sep="\n"),col="red",font=2,col.axis="red")
points(k.best,max(Si),pch=16,col="red",cex=1.5)

dend = as.dendrogram(clusterg_2)
heatmap(as.matrix(distanceg),Rowv = dend, symm = T, margin = c(3,3))

or = tabasco(trait_tar, scale="log",dend)
or = tabasco(trait_tar,dend)

# CALINSKY CRITERION
###############################################################################
cal_fit2 <- cascadeKM(trait_tar, 1, 10, iter = 1000)
plot(cal_fit2, sortg = TRUE, grpmts.plot = TRUE)
calinski.best2 <- as.numeric(which.max(cal_fit2$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best2, "\n")

# Elbow method
###############################################################################
library(factoextra)
# install.packages("NbClust")
library(NbClust)

fviz_nbclust(trait_tar, kmeans, method = "wss", k.max = 10) + 
  theme_minimal() + ggtitle("the Elbow Method")

# Gap Statistic
###############################################################################
library(cluster)
gap_stat <- clusGap(trait_tar, kmeans, 100, B = 100, verbose = interactive())
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

# silhouette
###############################################################################
fviz_nbclust(trait_tar, kmeans, method = "silhouette", k.max = 100) + 
  theme_minimal() + ggtitle("The Silhouette Plot")

###############################################################################
# Working with 3 guilds
###############################################################################

v                      = cutree(clusterg,k=4) # Clusters
genome2guild_5         = data.frame(guild = factor(v))
rownames(genome2guild_5) = names(v)
mat_trait_5     = as.data.frame(cbind(genome2guild_5,mat_ori$id,trait_tar))  
temp = mat_trait_5 %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:34]))
a   = as.data.frame(colnames(mat))
colnames(mat) = c("amino.sug.trans","glycosi.trans",
                  "FOS.fayt.trans","monosac.trans","aminoac.trans",
                  "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                  "nucleos.trans","nucleot.trans","ribonucle.trans",
                  "organophos.trans","peptide.trans","cellulose","chitin",
                  "heteroman","glucan.mix","xylan","xyloglucan","protein",
                  "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                  "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                  "ETC.II","ETC.III","nitrite")

# Ordination
###############################################################################
set.seed(2)
colnames(trait_tar) = c("amino.sug.trans","glycosi.trans",
                        "FOS.fayt.trans","monosac.trans","aminoac.trans",
                        "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                        "nucleos.trans","nucleot.trans","ribonucle.trans",
                        "organophos.trans","peptide.trans","cellulose","chitin",
                        "heteroman","glucan.mix","xylan","xyloglucan","protein",
                        "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                        "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                        "ETC.II","ETC.III","nitrite")
temp_1 = metaMDS(trait_tar,autotransform = F,trymax = 2000)
ordiplot(temp_1)
ordiplot(temp_1,type = "t")
library(ggvegan)

# One Panel
###############################################################################
fort = fortify(temp_1)

ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,size = 5),
                      alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(temp_1$stress,2)))

# Substrate
###############################################################################
substrate = as.data.frame(substring(mat_ori$id,1,1)) 
colnames(substrate) = c("subst")

ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,
                                    shape = substrate$subst,size = 5),alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(nmds$stress,2)))

# Other codes
###############################################################################
other = as.data.frame(substring(mat_ori$id,20)) 
colnames(other) = c("other")
other = other %>% mutate(other = replace(other,other == ".orig","orig"))
other = other %>% mutate(other = replace(other,other == ".permissive","permissive"))
other = other %>% mutate(other = replace(other,other == ".strict","strict"))

ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,
                                    shape = other$other,size = 5),alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(nmds$stress,2)))

# Spider plots for YAS framework traits
###############################################################################
b   = as.data.frame(colnames(temp)) # extract columns to extract the data for each category
r_acqui = rowSums(temp %>% select(2:22), na.rm=FALSE)
s_tol   = rowSums(temp %>% select(23:26), na.rm=FALSE)
r_use   = rowSums(temp %>% select(27:34), na.rm=FALSE)

library(fmsb)

temp1  = as.data.frame(cbind(temp$guild,r_acqui,s_tol,r_use))
total  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2  = as.data.frame(rbind(rep(max(total),3),
                             rep(min(r_use),3),temp1[2:4]))
colnames(temp2) = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
radarchart(temp2) # Spider plots using 100 guilds and total costs

temp1a  = as.data.frame(cbind(temp$guild,r_acqui/21,s_tol/4,r_use/8))
totala  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2a  = as.data.frame(rbind(rep(max(temp1a$V2),3),
                              rep(min(temp1a$V4),3),temp1a[2:4]))
colnames(temp2a) = c("Resource Acquisition","Stress Tolerance",
                     "Resource Use")
radarchart(temp2a) # Spider plots using 100 guilds and costs per trait

# Binary plots scaled by trait number
###############################################################################
r_acqui1 = rowSums(temp %>% select(2:22), na.rm=FALSE)/21
s_tol1   = rowSums(temp %>% select(23:26), na.rm=FALSE)/4
r_use1   = rowSums(temp %>% select(27:34), na.rm=FALSE)/8

par(mfrow=c(1,3))
plot(r_acqui1,s_tol1,xlab = "Resource Acquisition", ylab = "Stress Tolerance",
     col = "red",pch = 15,cex.lab = 1.5,cex = 3)
plot(s_tol1,r_use1,xlab = "Stress Tolerance", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5,cex = 3)
plot(r_acqui1,r_use1,xlab = "Resource Acquisition", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5,cex = 3)
par(mfrow=c(1,1))


############

###############################################################################
# Working with 100 guilds
###############################################################################

v                      = cutree(clusterg,k=100) # Clusters
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
mat_trait_f     = as.data.frame(cbind(genome2guild,mat_ori$id,trait_tar))  
temp = mat_trait_f %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat = (as.matrix(temp[2:34]))
a   = as.data.frame(colnames(mat))
colnames(mat) = c("amino.sug.trans","glycosi.trans",
                  "FOS.fayt.trans","monosac.trans","aminoac.trans",
                  "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                  "nucleos.trans","nucleot.trans","ribonucle.trans",
                  "organophos.trans","peptide.trans","cellulose","chitin",
                  "heteroman","glucan.mix","xylan","xyloglucan","protein",
                  "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                  "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                  "ETC.II","ETC.III","nitrite")

# Heatmap of trait distribution
###############################################################################
heatmap(mat,Colv = NA) # including main cluster
heatmap(mat) # including both clusters

###############################################################################
# Ordination
###############################################################################
nmds = metaMDS(distanceg,trymax = 2000)
nmds
stressplot(nmds)
nmds$points

NMDS1 <- nmds$points[,1] 
NMDS2 <- nmds$points[,2]
df_nor2.plot<-cbind(mat_trait_f, NMDS1, NMDS2)

ggplot(data=df_nor2.plot,aes(x=NMDS1,y=NMDS2,color =guild))+geom_point(size = 4) + 
  annotate("text", x=-0.4, y=-0.6, label=paste('Stress =',round(nmds$stress,2)))

###############################################################################
# Binary plots with total values
###############################################################################
b   = as.data.frame(colnames(temp)) # extract columns to extract the data for each cathegory
r_acqui = rowSums(temp %>% select(2:22), na.rm=FALSE)
s_tol   = rowSums(temp %>% select(23:26), na.rm=FALSE)
r_use   = rowSums(temp %>% select(27:34), na.rm=FALSE)

par(mfrow=c(1,3))
plot(r_acqui,s_tol,xlab = "Resource Acquisition", ylab = "Stress Tolerance",
     col = "red",pch = 15,cex.lab = 1.5)
plot(s_tol,r_use,xlab = "Stress Tolerance", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
plot(r_acqui,r_use,xlab = "Resource Acquisition", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

###############################################################################
# Binary plots scaled by trait number
###############################################################################
r_acqui1 = rowSums(temp %>% select(2:22), na.rm=FALSE)/21
s_tol1   = rowSums(temp %>% select(23:26), na.rm=FALSE)/4
r_use1   = rowSums(temp %>% select(27:34), na.rm=FALSE)/8

par(mfrow=c(1,3))
plot(r_acqui1,s_tol1,xlab = "Resource Acquisition", ylab = "Stress Tolerance",
     col = "red",pch = 15,cex.lab = 1.5)
plot(s_tol1,r_use1,xlab = "Stress Tolerance", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
plot(r_acqui1,r_use1,xlab = "Resource Acquisition", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

###############################################################################
# Spider plots for YAS framework traits
###############################################################################

library(fmsb)

temp1  = as.data.frame(cbind(temp$guild,r_acqui,s_tol,r_use))
total  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2  = as.data.frame(rbind(rep(max(total),3),
                             rep(min(r_use),3),temp1[2:4]))
colnames(temp2) = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
radarchart(temp2) # Spider plots using 100 guilds and total costs

temp1a  = as.data.frame(cbind(temp$guild,r_acqui/21,s_tol/4,r_use/8))
totala  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2a  = as.data.frame(rbind(rep(max(r_acqui/21),3),
                             rep(min(r_use),3),temp1a[2:4]))
colnames(temp2a) = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
radarchart(temp2a) # Spider plots using 100 guilds and costs per trait

###############################################################################
# Spider plots for within YAS framework traits
###############################################################################

# Resource acquisition
###############################################################################
resource_acquisition1 = temp %>% select(2:22)
totala  = rowSums(resource_acquisition1, na.rm=FALSE)
colnames(resource_acquisition1) <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                                     "13","14","15","16","17","18","19","20","21")
RA_total   = as.data.frame(rbind(rep(max(totala),3),
                                 rep(min(resource_acquisition1),3),resource_acquisition1))
radarchart(RA_total)

par(mfrow=c(1,2))
plot(temp$Resource.Acquisition.Substrate.degradation.simple.compound.degradation.protein.degradation,
     temp$Resource.Acquisition.Substrate.uptake.free.amino.acids.transport,xlab = "protein degradation", 
     ylab = "free aminoacids transporters",
     col = "red",pch = 15,cex.lab = 1.5)
plot(temp$Resource.Acquisition.Substrate.degradation.simple.compound.degradation.protein.degradation,
     temp$Resource.Acquisition.Substrate.degradation.complex.carbohydrate.depolymerization.xylan.and.heteroxylan.breakdown,xlab = "protein degradation", 
     ylab = "xylan-heteroxylan degradation",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

# Stress Tolerance
###############################################################################
stress_tolerance1 = temp %>% select(23:26)
totala  = rowSums(stress_tolerance1, na.rm=FALSE)
colnames(stress_tolerance1) <- c("solute.transport",
                                 "solute.synthesis","EPS.biosynthesis(S)",
                                 "osmotic.sensors")
st_total   = as.data.frame(rbind(rep(max(totala),3),
                                 rep(min(stress_tolerance1),3),stress_tolerance1))
radarchart(st_total)

par(mfrow=c(1,2))
plot(temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.accumulation.of.compatible.solutes.synthesis.of.compatible.solutes,
     temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.accumulation.of.compatible.solutes.transport.of.compatible.solutes,xlab = "solute synthesis", 
     ylab = "solute transport",
     col = "red",pch = 15,cex.lab = 1.5)
plot(temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.accumulation.of.compatible.solutes.synthesis.of.compatible.solutes,
     temp$Stress.Tolerance.Specific.desiccation.osmotic.salt.stress.enhancement.of.microenvironment.EPS.biosynthesis.export,xlab = "solute synthesis", 
     ylab = "EPS biosynthesis",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

# Resource use
###############################################################################
resource_use1 = temp %>% select(27:34)
totala  = rowSums(resource_use1, na.rm=FALSE)
colnames(resource_use1) <- c("1","2","3","4","5","6","7","8")
st_total   = as.data.frame(rbind(rep(max(totala),3),
                                 rep(min(resource_use1),3),resource_use1))
radarchart(st_total)

par(mfrow=c(1,2))
plot(temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.electron.transport.chain..ETC.complex.III,
     temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.electron.transport.chain..ETC.complex.II,xlab = "ETC complex III", 
     ylab = "ETC complex II",
     col = "red",pch = 15,cex.lab = 1.5)
plot(temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.electron.transport.chain..ETC.complex.III,
     temp$Resource.Use.Chemotrophy.chemoorganoheterotrophy.aerobic.respiration.oxidative.pentose.phosphate.pathway,xlab = "ETC complex III", 
     ylab = "oxidative pentose phosphate pathway ",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

###############################################################################
# Guild 65 (invest the most in resource acquisition)
###############################################################################

guild_65 = mat_trait_f %>% filter(guild == 65)

a = as.data.frame(colnames(guild_65))
resource_acquisition1 = rowSums(guild_65 %>% select(3:23), na.rm=FALSE)
stress_tolerance1     = rowSums(guild_65 %>% select(24:27), na.rm=FALSE)
resource_use1         = rowSums(guild_65 %>% select(28:35), na.rm=FALSE)

temp1a  = as.data.frame(cbind(guild_65$guild,resource_acquisition1,
                              resource_use1,
                              stress_tolerance1))
temp2a  = as.data.frame(rbind(rep(max(resource_acquisition1),3),
                              rep(min(resource_use1),3),temp1a[2:4]))
radarchart(temp1a)

# Stress tolerance
###############################################################################
guild_6st = guild_65 %>% select(24:27)
colnames(guild_6st) <- c("solute.transport",
                         "solute.synthesis","EPS.biosynthesis(S)",
                         "osmotic.sensors")
temp6st   = as.data.frame(rbind(rep(max(guild_6st),3),
                                rep(min(guild_6st),3),guild_6st))
radarchart(temp6st)

# Resource Acquisition
###############################################################################
guild_6ra = guild_65 %>% select(3:23)
colnames(guild_6ra) <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                         "13","14","15","16","16","18","19","20","21")
temp6ra   = as.data.frame(rbind(rep(max(guild_6ra),3),
                                rep(min(guild_6ra),3),guild_6ra))
radarchart(temp6ra)

# Resource use
###############################################################################
guild_6ru = guild_65 %>% select(28:35)
colnames(guild_6ru) <- c("1","2","3","4","5","6","7","8")
temp6ru   = as.data.frame(rbind(rep(max(guild_6ru),3),
                                 rep(min(guild_6ru),3),guild_6ru))
radarchart(temp6ru)

###############################################################################
# Correlation plots
###############################################################################

library(ggcorrplot)

# Correlation plots - all samples
###############################################################################

colnames(trait_tar) = c("amino.sug.trans","glycosi.trans",
                        "FOS.fayt.trans","monosac.trans","aminoac.trans",
                        "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                        "nucleos.trans","nucleot.trans","ribonucle.trans",
                        "organophos.trans","peptide.trans","cellulose","chitin",
                        "heteroman","glucan.mix","xylan","xyloglucan","protein",
                        "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                        "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                        "ETC.II","ETC.III","nitrite")
corr <- round(cor(trait_tar), 1)
p.mat <- cor_pmat(trait_tar)
ggcorrplot(corr)
# using hierarchical clustering
ggcorrplot(corr, hc.order = TRUE, outline.col = "white")
# Leave blank on no significant coefficient
ggcorrplot(corr, p.mat = p.mat, hc.order = TRUE,
           type = "lower", insig = "blank")

# Correlation plots - 100 guilds
###############################################################################

temp_100 = temp %>% select(2:34)
colnames(temp_100) = c("amino.sug.trans","glycosi.trans",
                        "FOS.fayt.trans","monosac.trans","aminoac.trans",
                        "lipid.trans","amide.trans","NH4.trans","nucleob.trans",
                        "nucleos.trans","nucleot.trans","ribonucle.trans",
                        "organophos.trans","peptide.trans","cellulose","chitin",
                        "heteroman","glucan.mix","xylan","xyloglucan","protein",
                        "sol.trans","sol.synt","EPS.synt","osmotic.sensor",
                        "oxi.pentose","nonoxi.pentose","glyoxylate","TCA","ETC.I",
                        "ETC.II","ETC.III","nitrite")
corr100 <- round(cor(temp_100), 1)
p.mat100 <- cor_pmat(temp_100)
ggcorrplot(corr100)
ggcorrplot(corr100, hc.order = TRUE, outline.col = "white")
# Leave blank on no significant coefficient
ggcorrplot(corr100, p.mat = p.mat100, hc.order = TRUE,
           type = "lower", insig = "blank")

###############################################################################
# CAZymes-related traits
###############################################################################

trait_tarCAZ = mat_trait %>% select(6,8,43:48) # 
genomes      = as.list(mat_trait[1])
test         = trait_tar %>% summarize_if(is.numeric, sum, na.rm=TRUE)
names(trait_tarCAZ)
colnames(trait_tarCAZ) = c("glycoside.trans","monosac.trans","cellulose","chitin",
                           "heteromannan","mixed.glucan","xylan","xyloglucan")
boxplot(trait_tarCAZ,ylab="gene cost")

###############################################################################
# Forming functional guilds 
###############################################################################

distanceg  = vegdist(trait_tarCAZ, method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar), 2)
plot(clusterg)

# Similiarity within guilds
###############################################################################
my_vec = c()

for(i in nguilds) {
  v                      = cutree(clusterg,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distanceg ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1        
  my_vec <- c(my_vec, my_out)    
}

# Plotting
mat_r2_c  = as.data.frame(cbind(nguilds,my_vec))
ggplot(data=mat_r2_c,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 20)) + 
  geom_hline(yintercept=0.9200740, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 100, linetype="dashed", color = "red")

barplot(c(29,67,4),ylab="Number of guilds",names.arg=c("Grass","Both","Shrub"))