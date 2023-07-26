# Final Analysis

# General goal: fitness traits tradeoffs and microbial life history strategies 
# using genome-scale data

# Packages ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stats)
library(corrr)
library(olsrr)

# Calling data and preprocessing ####

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

# Aim 1 ####
# Traits responsible for adaptations under drought in grassland in Loma ####

# Ratio of abundances drought grassland ####
mag_stat.1 = mag_stat %>% mutate(ratio1.1 = log10((mag_stat$Average.1/mag_stat$Average)*
                                                  abs(mag_stat$Average.1-mag_stat$Average)))
mag_stat.1 = mag_stat.1 %>% mutate(pos.1 = ratio1.1 >= 0)
mag_stat.1 = mag_stat.1 %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p.1<-ggplot(data=mag_stat.1, aes(x=id_2, y=ratio1.1, fill = pos.1)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p.1

# Preliminary filtering ####

mat_trait   = mat_trait %>% mutate(rat_am_dro_gra.1 = mag_stat.1$ratio1.1)
mat_trait   = mat_trait %>% mutate(rat_am_dro_gra.1 = ifelse(is.na(rat_am_dro_gra.1),0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1==Inf,0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1 ==-Inf, 0, rat_am_dro_gra.1))

mat_trait_1 = as.data.frame(mat_trait[,-192])
mat = (as.matrix(mat_trait_1[2:192]))

corre.1 = as.data.frame(cor(mat[,1:190], mat[,191]))

# Erase traits based on correlation matrix

test.1    = as.data.frame(cbind(mat_trait[2:191]))
seq.1     = seq(1,190)
temp.1    = as.data.frame(cbind((colSums(test.1)),seq.1))
erase.1   = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
test.1    = test.1[-(erase.1$seq.1)]
temp.1    = as.data.frame(colSums(test.1))
corre.1.1 = as.data.frame(cor(as.matrix(test.1), mat[,191]))
#write.csv(test.1, file = "test.1.csv")
test.1    = read.csv("test.1.csv") 

# Best predictors statistics ####

m.all        = lm(mat[,191]~.,data=test.1[2:50])
temp.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat[,191]~.,data=test.1[51:100])
temp.51.100  = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat[,191]~.,data=test.1[101:116])
temp.101.116 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat[,191]~.,data=test.1[117:128])
temp.117.128 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat[,191]~.,data=test.1[129:138])
temp.129.138 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat[,191]~.,data=test.1[139:141])
temp.139.141 = ols_step_forward_p(m.all, details = FALSE)
erase.1a     = as.data.frame(c(temp.1.50[["metrics"]][["variable"]],
                               temp.51.100[["metrics"]][["variable"]],
                               temp.101.116[["metrics"]][["variable"]],
                               temp.117.128[["metrics"]][["variable"]],
                               temp.129.138[["metrics"]][["variable"]],
                               temp.139.141[["metrics"]][["variable"]]))
colnames(erase.1a) = c("trait")
test.1             = test.1 %>% select((erase.1a$trait))
corre.1.2          = as.data.frame(cor(as.matrix(test.1),mat[,191]))
a                  = as.data.frame(colnames(test.1))
test.1             = test.1[c(-58,-60,-65)]
m.all              = lm(mat[,191]~.,data=test.1)
best.predictors    = ols_step_forward_p(m.all, details = FALSE)
erase.1a           = as.data.frame(best.predictors[["metrics"]][["variable"]])
colnames(erase.1a) = c("trait")
test.1             = test.1 %>% select((erase.1a$trait))
m.all              = lm(mat[,191]~.,data=test.1)
#write.csv(test.1, file = "grassland.drought.best.predictors.csv")
test.1             = read.csv("grassland.drought.best.predictors.csv") 

matgs = as.matrix(cbind(mat[,191],test.1[2:43]))
x2 <- matgs %>% 
  correlate() %>% 
  focus(`mat[, 191]`)
x2 %>% 
  mutate(rowname = factor(term, levels = term[order(`mat[, 191]`)])) 

gran_1 = c("S.Tolerance","R.Acquisition","R.Acquisition","R.Acquisition",
           "S.Tolerance","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Use","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Use","S.Tolerance","R.Acquisition","R.Acquisition","R.Acquisition",
           "S.Tolerance","R.Acquisition","S.Tolerance","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Use","R.Acquisition","R.Acquisition","S.Tolerance",
           "R.Use","S.Tolerance","S.Tolerance","R.Acquisition","S.Tolerance",
           "S.Tolerance","R.Acquisition","R.Acquisition")

gran_2 = c("Desiccation","S.Uptake","S.Uptake","S.Assimilation","T.Tolerance",
           "S.Uptake","S.Uptake","S.Degradation","S.Assimilation","S.Uptake",
           "Chemotrophy","S.Assimilation","S.Assimilation","S.Degradation",
           "S.Assimilation","S.Uptake","Chemotrophy","pH.Stress","S.Uptake",
           "S.Uptake","S.Uptake","T.Tolerance","S.Degradation","Phage.Resistance",
           "S.Degradation","S.Degradation","S.Degradation","S.Uptake","S.Uptake",
           "S.Assimilation","Chemotrophy","S.Uptake","S.Uptake","O.Stress",
           "Phototrophy","T.Tolerance","T.Tolerance","S.Degradation",
           "Desiccation","O.Stress","S.Assimilation","S.Degradation")

x2 = as.data.frame(cbind(gran_1,gran_2,x2))
colnames(x2) = c("trait.level.1","trait.level.2","term","pearson")
g <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.1)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g + scale_fill_manual(values = c("#fc8d62","#33a02c","#2166ac"))
g1 <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.2)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g1 + scale_fill_manual(values=c("#a6d96a","#c6dbef",
                                "#6baed6","#41b6c4","#0570b0","#006d2c",
                                "#a50026","#f46d43","#d73027","#08306b"))

# Traits responsible for adaptations under drought in shrubland in Loma ####

# Ratio of abundances drought shrubland ####
mag_stat.2 = mag_stat %>% mutate(ratio1.1 = log10((mag_stat$Average.3/mag_stat$Average.2)*
                                                  abs(mag_stat$Average.3-mag_stat$Average.2)))
mag_stat.2 = mag_stat.2 %>% mutate(pos.1 = ratio1.1 >= 0)
mag_stat.2 = mag_stat.2 %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p.2<-ggplot(data=mag_stat.2, aes(x=id_2, y=ratio1.1, fill = pos.1)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p.2

# Preliminary filtering ####

mat_trait.2   = mat_trait   %>% mutate(rat_am_dro_gra.1 = mag_stat.2$ratio1.1)
mat_trait.2   = mat_trait.2 %>% mutate(rat_am_dro_gra.1 = ifelse(is.na(rat_am_dro_gra.1),0,rat_am_dro_gra.1),
                                       rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1==Inf,0,rat_am_dro_gra.1),
                                       rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1 ==-Inf, 0, rat_am_dro_gra.1))
mat_trait.2   = as.data.frame(mat_trait.2[,-192])
mat.2         = (as.matrix(mat_trait.2[2:192]))
corre.2       = as.data.frame(cor(mat.2[,1:190], mat.2[,191]))

# Erase traits based on correlation matrix

test.2    = as.data.frame(cbind(mat_trait.2[2:191]))
seq.2     = seq(1,190)
temp.2    = as.data.frame(cbind((colSums(test.2)),seq.2))
erase.2   = temp.2 %>% filter(V1==0)
erase.2$row_names = row.names(erase.2)
test.2    = test.2[-(erase.2$seq.2)]
temp.2    = as.data.frame(colSums(test.2))
corre.2.1 = as.data.frame(cor(as.matrix(test.2), mat.2[,191]))
#write.csv(test.2, file = "test.2.csv")
test.2    = read.csv("test.2.csv") 

# Best predictors statistics ####

m.all        = lm(mat.2[,191]~.,data=test.2[2:50])
temp.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.2[,191]~.,data=test.2[51:94])
temp.51.94   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.2[,191]~.,data=test.2[95:109])
temp.95.109  = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.2[,191]~.,data=test.2[110:112])
temp.110.112 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.2[,191]~.,data=test.2[113:126])
temp.113.126 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.2[,191]~.,data=test.2[127:141])
temp.127.141 = ols_step_forward_p(m.all, details = FALSE)
erase.2a     = as.data.frame(c(temp.1.50[["metrics"]][["variable"]],
                               temp.51.94[["metrics"]][["variable"]],
                               temp.95.109[["metrics"]][["variable"]],
                               temp.110.112[["metrics"]][["variable"]],
                               temp.113.126[["metrics"]][["variable"]],
                               temp.127.141[["metrics"]][["variable"]]))
colnames(erase.2a) = c("trait")
test.2             = test.2 %>% select((erase.2a$trait))
corre.2.2          = as.data.frame(cor(as.matrix(test.2),mat.2[,191]))
a                  = as.data.frame(colnames(test.2))
test.2             = test.2[c(-44,-46,-47,-48,-60)]
m.all              = lm(mat.2[,191]~.,data=test.2)
best.predictors    = ols_step_forward_p(m.all, details = FALSE)
erase.2a           = as.data.frame(best.predictors[["metrics"]][["variable"]])
colnames(erase.2a) = c("trait")
test.2             = test.2 %>% select((erase.2a$trait))
m.all              = lm(mat.2[,191]~.,data=test.2)
#write.csv(test.2, file = "shrubland.drought.best.predictors.csv")
test.2    = read.csv("shrubland.drought.best.predictors.csv") 

matgs = as.matrix(cbind(mat.2[,191],test.2[2:42]))
x2 <- matgs %>% 
  correlate() %>% 
  focus(`mat.2[, 191]`)
x2 %>% 
  mutate(rowname = factor(term, levels = term[order(`mat.2[, 191]`)])) 

gran_1 = c("R.Acquisition","R.Acquisition","S.Tolerance","R.Acquisition",
           "R.Acquisition","S.Tolerance","R.Acquisition","R.Acquisition",
           "R.Acquisition","S.Tolerance","R.Acquisition","R.Use","R.Acquisition",
           "R.Use","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Use",
           "R.Acquisition","R.Use","S.Tolerance","S.Tolerance","R.Use",
           "R.Acquisition","Mean.Generation.Time","R.Acquisition","R.Use",
           "S.Tolerance","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","S.Tolerance","R.Acquisition")

gran_2 = c("S.Uptake","S.Degradation","T.Tolerance","S.Degradation","S.Uptake",
           "T.Tolerance","S.Assimilation","S.Assimilation","S.Uptake","pH.Stress",
           "S.Uptake","Chemotrophy","S.Degradation","Phototrophy","S.Assimilation",
           "S.Uptake","S.Uptake","S.Uptake","S.Uptake","S.Degradation","Chemotrophy",
           "S.Assimilation","Chemotrophy","O.Stress","pH.Stress","Chemotrophy","S.Uptake",
           "Mean.Generation.Time","S.Uptake","Chemotrophy","Desiccation","S.Degradation",
           "S.Degradation","S.Uptake","S.Uptake","S.Uptake","S.Degradation","S.Degradation",
           "S.Uptake","pH.Stress","S.Assimilation")

x2 = as.data.frame(cbind(gran_1,gran_2,x2))
colnames(x2) = c("trait.level.1","trait.level.2","term","pearson")
g <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.1)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g + scale_fill_manual(values = c("#999999","#fc8d62","#33a02c","#2166ac"))
g1 <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.2)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g1 + scale_fill_manual(values=c("#a6d96a","#c6dbef",
                                "#999999","#6baed6","#41b6c4","#006d2c",
                                "#a50026","#f46d43","#d73027","#08306b"))

# Traits responsible for adaptations under in shrubland in Loma ####

# Ratio of abundances shrubland ####
mag_stat.3 = mag_stat %>% mutate(ratio1.1 = log10((mag_stat$Average.2/mag_stat$Average)*
                                                  abs(mag_stat$Average.2-mag_stat$Average)))
mag_stat.3 = mag_stat.3 %>% mutate(pos.1 = ratio1.1 >= 0)
mag_stat.3 = mag_stat.3 %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))
p.3<-ggplot(data=mag_stat.3, aes(x=id_2, y=ratio1.1, fill = pos.1)) +
  geom_bar(stat="identity") + theme(legend.position="none") + ylab("Log-transformed Ratio") + 
  xlab("MAG ID")
p.3

# Preliminary filtering ####

mat_trait.3   = mat_trait   %>% mutate(rat_am_dro_gra.1 = mag_stat.3$ratio1.1)
mat_trait.3   = mat_trait.3 %>% mutate(rat_am_dro_gra.1 = ifelse(is.na(rat_am_dro_gra.1),0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1==Inf,0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1 ==-Inf, 0, rat_am_dro_gra.1))

mat_trait.3   = as.data.frame(mat_trait.3[,-192])
mat.3         = (as.matrix(mat_trait.3[2:192]))

corre.3       = as.data.frame(cor(mat.3[,1:190], mat.3[,191]))

# Erase traits based on correlation matrix

test.3    = as.data.frame(cbind(mat_trait.3[2:191]))
seq.3     = seq(1,190)
temp.3    = as.data.frame(cbind((colSums(test.3)),seq.3))
erase.3   = temp.3 %>% filter(V1==0)
erase.3$row_names = row.names(erase.3)
test.3    = test.3[-(erase.3$seq.3)]
temp.3    = as.data.frame(colSums(test.3))
corre.3.1 = as.data.frame(cor(as.matrix(test.3), mat.3[,191]))
#write.csv(test.3, file = "test.3.csv")
test.3    = read.csv("test.3.csv") 

# Best predictors statistics (Stepwise forward regression) ####

m.all        = lm(mat.3[,191]~.,data=test.3[2:50])
temp.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.3[,191]~.,data=test.3[51:94])
temp.51.94   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.3[,191]~.,data=test.3[95:116])
temp.95.116  = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat.3[,191]~.,data=test.3[117:141])
temp.117.141 = ols_step_forward_p(m.all, details = FALSE)
erase.3a     = as.data.frame(c(temp.1.50[["metrics"]][["variable"]],
                               temp.51.94[["metrics"]][["variable"]],
                               temp.95.116[["metrics"]][["variable"]],
                               temp.117.141[["metrics"]][["variable"]]))
colnames(erase.3a) = c("trait")
test.3             = test.3 %>% select((erase.3a$trait))
corre.3.2          = as.data.frame(cor(as.matrix(test.3),mat.3[,191]))
a                  = as.data.frame(colnames(test.3))
test.3             = test.3[c(-36,-45,-50)]
m.all              = lm(mat.3[,191]~.,data=test.3)
best.predictors    = ols_step_forward_p(m.all, details = FALSE)
erase.3a           = as.data.frame(best.predictors[["metrics"]][["variable"]])
colnames(erase.3a) = c("trait")
test.3             = test.3 %>% select((erase.3a$trait))
m.all              = lm(mat.3[,191]~.,data=test.3)
#write.csv(test.3, file = "shrubland.best.predictors.csv")
test.3    = read.csv("shrubland.best.predictors.csv") 

matgs = as.matrix(cbind(mat.3[,191],test.3[2:33]))
x2 <- matgs %>% 
  correlate() %>% 
  focus(`mat.3[, 191]`)
x2 %>% 
  mutate(rowname = factor(term, levels = term[order(`mat.3[, 191]`)])) 

gran_1 = c("S.Tolerance","S.Tolerance","R.Acquisition","S.Tolerance",
           "Mean.Generation.Time","R.Acquisition","R.Acquisition","R.Use",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Use",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Use","R.Acquisition","R.Acquisition",
           "R.Acquisition","S.Tolerance","R.Acquisition","S.Tolerance",
           "S.Tolerance","R.Acquisition","R.Acquisition","R.Use","R.Use",
           "S.Tolerance","R.Acquisition","R.Use")

gran_2 = c("T.Tolerance","T.Tolerance","S.Degradation","EPS.biosynthesis",
           "Mean.Generation.Time","S.Uptake","S.Uptake","Chemotrophy",
           "S.Assimilation","S.Uptake","S.Assimilation","Phototrophy",
           "S.Uptake","S.Uptake","S.Assimilation","S.Uptake","S.Degradation",
           "Chemotrophy","S.Assimilation","S.Uptake","S.Uptake","pH.Stress",
           "S.Uptake","T.Tolerance","O.Stress","S.Degradation",
           "S.Assimilation","Phototrophy","Chemotrophy","T.Tolerance",
           "S.Degradation","Phototrophy")

x2 = as.data.frame(cbind(gran_1,gran_2,x2))
colnames(x2) = c("trait.level.1","trait.level.2","term","pearson")
g <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.1)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g + scale_fill_manual(values = c("#999999","#fc8d62","#33a02c","#2166ac"))
g1 <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.2)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g1 + scale_fill_manual(values=c("#a6d96a","#c6dbef",
                                "#999999","#6baed6","#41b6c4","#006d2c",
                                "#a50026","#f46d43","#d73027","#08306b"))

# Aim 2 ####

# Functional groups under different conditions

ratio     = as.data.frame((mag_stat$Average.1/mag_stat$Average))
colnames(ratio) = "ratio"
test.1    = as.data.frame(cbind(test.1,ratio))
test.1    = as.data.frame(test.1[c(-163,-164,-263,-299,-387,-396,-399,-425,-479,
                                   -485,-487,-491,-499,-506,-511),])

