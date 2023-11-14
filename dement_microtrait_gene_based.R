# FUNCIONAL GROUPS BASED ON GENE INFORMATION ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stats)
library(corrr)
library(olsrr)
library(vegan)
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
# install.packages("remotes")
# remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(fmsb)

# LOMA RIDGE ####

# Calling data and preprocessing ####
hmm_loma    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/hmm_Loma.csv",dec=".")
gene_loma   = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/litter_mags_metadata.txt",dec=".")
mag_stat    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_stats.txt") 
mag_abun    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_adundance.txt") 
mag_stat    = mag_stat %>% full_join(mag_abun)

# Gene selection ####

# Grassland ####
mat_trait_g   = hmm_loma %>% mutate(grass_abund = mag_stat$Average)
mat_g         = (as.matrix(mat_trait_g[3:1299]))
# Best predictors statistics - step 1 ####
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[3:75])
temp.1.75      = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[76:140])
temp.76.140    = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[141:145])
temp.141.145   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[146:148])
temp.146.148   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[149:300])
temp.149.300   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[301:500])
temp.301.500   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[501:650])
temp.501.650   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[651:800])
temp.651.800   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[801:920])
temp.801.920   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[921:1030])
temp.921.1030  = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[1031:1200])
temp.1031.1200 = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_g[,1297]~.,data=mat_trait_g[1201:1296])
temp.1201.1296 = ols_step_forward_p(m.all, details = FALSE)
erase.1a     = as.data.frame(c(temp.1.75[["predictors"]],
                               temp.76.140[["predictors"]],
                               temp.141.145[["predictors"]],
                               temp.146.148[["predictors"]],
                               temp.149.300[["predictors"]],
                               temp.301.500[["predictors"]],
                               temp.501.650[["predictors"]],
                               temp.651.800[["predictors"]],
                               temp.801.920[["predictors"]],
                               temp.921.1030[["predictors"]],
                               temp.1031.1200[["predictors"]],
                               temp.1201.1296[["predictors"]]))
colnames(erase.1a) = c("trait")
test.g             = test.1 %>% select((erase.1a$trait))
# Best predictors statistics - step 2 ####
m.all              = lm(mat_g[,1297]~.,data=test.g[1:100])
temp.1.100.2       = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=test.g[101:200])
temp.101.200.2     = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=test.g[201:300])
temp.201.300.2     = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=test.g[301:400])
temp.301.400.2     = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=test.g[401:440])
temp.401.440.2     = ols_step_forward_p(m.all, details = FALSE)
erase.2a     = as.data.frame(c(temp.1.100.2[["predictors"]],
                               temp.101.200.2[["predictors"]],
                               temp.201.300.2[["predictors"]],
                               temp.301.400.2[["predictors"]],
                               temp.401.440.2[["predictors"]]))
colnames(erase.2a) = c("trait")
test.g             = test.g %>% select((erase.2a$trait))
# Best predictors statistics - step 3 ####
m.all              = lm(mat_g[,1297]~.,data=test.g[1:70])
temp.1.70.3        = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=test.g[71:140])
temp.71.140.3      = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=test.g[141:220])
temp.141.220.3     = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=test.g[221:300])
temp.221.300.3     = ols_step_forward_p(m.all, details = FALSE)
erase.3a     = as.data.frame(c(temp.1.70.3[["predictors"]],
                               temp.71.140.3[["predictors"]],
                               temp.141.220.3[["predictors"]],
                               temp.221.300.3[["predictors"]]))
colnames(erase.3a) = c("trait")
test.g             = test.g %>% select((erase.3a$trait))

# Best predictors statistics - step 4 ####
mat.g    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/temp_grass_ambient.csv",dec=".")

m.all              = lm(mat_g[,1297]~.,data=mat.g[2:100])
temp.1.100.4       = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,1297]~.,data=mat.g[101:236])
temp.101.236.4     = ols_step_forward_p(m.all, details = FALSE)
erase.4a           = as.data.frame(c(temp.1.100.4[["predictors"]],
                                     temp.101.236.4[["predictors"]]))
colnames(erase.4a) = c("trait")
mat.g              = mat.g %>% select((erase.4a$trait))

# Best predictors statistics - step 5 ####
m.all              = lm(mat_g[,1297]~.,data=mat.g)
temp.5             = ols_step_forward_p(m.all, details = FALSE)
erase.5a           = as.data.frame((temp.5[["predictors"]]))
colnames(erase.5a) = c("trait")
mat.g              = mat.g %>% select((erase.5a$trait))

write.csv(mat.g, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/grassland.ambient.best.predictors.csv")

# Shrubland ####
mat_trait_s   = hmm_loma %>% mutate(shrub_abund = mag_stat$Average.2)
mat_s         = (as.matrix(mat_trait_s[3:1299]))
# Best predictors statistics - step 1 ####
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[3:150])
temp.s.3.150     = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[151:300])
temp.s.151.300   = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[301:500])
temp.s.301.500   = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[501:700])
temp.s.501.700   = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[701:850])
temp.s.701.850   = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[851:950])
temp.s.851.950   = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[951:990])
temp.s.951.990   = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[991:1100])
temp.s.991.1100  = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[1101:1200])
temp.s.1101.1200 = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_s[,1297]~.,data=mat_trait_s[1201:1296])
temp.s.1201.1296 = ols_step_forward_p(m.all, details = FALSE)
erase.1s     = as.data.frame(c(temp.s.3.150[["predictors"]],
                               temp.s.151.300[["predictors"]],
                               temp.s.301.500[["predictors"]],
                               temp.s.501.700[["predictors"]],
                               temp.s.701.850[["predictors"]],
                               temp.s.851.950[["predictors"]],
                               temp.s.951.990[["predictors"]],
                               temp.s.991.1100[["predictors"]],
                               temp.s.1101.1200[["predictors"]],
                               temp.s.1201.1296[["predictors"]]))
colnames(erase.1s) = c("trait")
mat_s              = as.data.frame(mat_s) %>% select((erase.1s$trait))

# Best predictors statistics - step 2 ####
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s[1:100])
temp.s1.1.100    = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s[101:200])
temp.s1.101.200  = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s[201:300])
temp.s1.201.300  = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s[301:406])
temp.s1.301.406  = ols_step_forward_p(m.all, details = FALSE)
erase.2s         = as.data.frame(c(temp.s1.1.100[["predictors"]],
                                   temp.s1.101.200[["predictors"]],
                                   temp.s1.201.300[["predictors"]],
                                   temp.s1.301.406[["predictors"]]))
colnames(erase.2s) = c("trait")
mat_s              = as.data.frame(mat_s) %>% select((erase.2s$trait))

# Best predictors statistics - step 3 ####
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s[1:100])
temp.s2.1.100    = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s[101:200])
temp.s2.101.200  = ols_step_forward_p(m.all, details = FALSE)
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s[201:267])
temp.s2.201.267  = ols_step_forward_p(m.all, details = FALSE)
erase.3s         = as.data.frame(c(temp.s2.1.100[["predictors"]],
                                   temp.s2.101.200[["predictors"]],
                                   temp.s2.201.267[["predictors"]]))
colnames(erase.3s) = c("trait")
mat_s              = as.data.frame(mat_s) %>% select((erase.3s$trait))

# Best predictors statistics - step 4 ####
m.all            = lm(mat_trait_s[,1299]~.,data=mat_s)
temp.s4          = ols_step_forward_p(m.all, details = FALSE)
erase.4s         = as.data.frame(c(temp.s4[["predictors"]]))
colnames(erase.4s) = c("trait")
mat_s              = as.data.frame(mat_s) %>% select((erase.4s$trait))

write.csv(mat_s, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/shrubland.ambient.best.predictors.csv")

# Grassland Drought ####

mat_trait_gd   = hmm_loma %>% mutate(grass.abund.dro = mag_stat$Average.1)
mat.g.dro      = (as.matrix(mat_trait_gd[3:1299]))

# Best predictors statistics - step 1 ####
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[3:75])
temp.gd.3.75      = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[76:150])
temp.gd.76.150    = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[151:300])
temp.gd.151.300   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[301:450])
temp.gd.301.450   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[451:600])
temp.gd.451.600   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[601:700])
temp.gd.601.700   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[701:730])
temp.gd.701.730   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[731:740])
temp.gd.731.740   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[741:800])
temp.gd.741.800   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[801:900])
temp.gd.801.900   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[901:1000])
temp.gd.901.1000  = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat.g.dro[,1297]~.,data=mat_trait_gd[1001:1296])
temp.gd.1001.1296 = ols_step_forward_p(m.all, details = FALSE)

erase.1gd           = as.data.frame(c(temp.gd.3.75[["predictors"]],
                                      temp.gd.76.150[["predictors"]],
                                      temp.gd.151.300[["predictors"]],
                                      temp.gd.301.450[["predictors"]],
                                      temp.gd.451.600[["predictors"]],
                                      temp.gd.601.700[["predictors"]],
                                      temp.gd.701.730[["predictors"]],
                                      temp.gd.731.740[["predictors"]],
                                      temp.gd.741.800[["predictors"]],
                                      temp.gd.801.900[["predictors"]],
                                      temp.gd.901.1000[["predictors"]],
                                      temp.gd.1001.1296[["predictors"]]))
colnames(erase.1gd) = c("trait")
mat_gd              = mat_trait_gd %>% select((erase.1gd$trait))

# Best predictors statistics - step 2 #### 

mat_gd    = read.csv("temp_1.csv",dec=".")

m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[2:100])
temp.gd.1.2.100     = ols_step_forward_p(m.all, details = FALSE)
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[101:200])
temp.gd.1.101.200   = ols_step_forward_p(m.all, details = FALSE)
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[201:300])
temp.gd.1.201.300   = ols_step_forward_p(m.all, details = FALSE)
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[301:410])
temp.gd.1.301.410   = ols_step_forward_p(m.all, details = FALSE)

erase.2gd           = as.data.frame(c(temp.gd.1.2.100[["predictors"]],temp.gd.1.101.200[["predictors"]],
                                      temp.gd.1.201.300[["predictors"]],temp.gd.1.301.410[["predictors"]]))
colnames(erase.2gd) = c("trait")
mat_gd              = mat_gd %>% select((erase.2gd$trait))

# Best predictors statistics - step 3 #### 
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[1:60])
temp.gd.2.1.60      = ols_step_forward_p(m.all, details = FALSE)
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[61:200])
temp.gd.2.61.200    = ols_step_forward_p(m.all, details = FALSE)
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[201:268])
temp.gd.2.101.268   = ols_step_forward_p(m.all, details = FALSE)
erase.3gd           = as.data.frame(c(temp.gd.2.1.60[["predictors"]],temp.gd.2.61.200[["predictors"]],
                                      temp.gd.2.101.268[["predictors"]]))
colnames(erase.3gd) = c("trait")
mat_gd              = mat_gd %>% select((erase.3gd$trait))

# Best predictors statistics - step 4 #### 
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[1:100])
temp.gd.4.1.100     = ols_step_forward_p(m.all, details = FALSE)
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[101:199])
temp.gd.4.101.199   = ols_step_forward_p(m.all, details = FALSE)
erase.4gd           = as.data.frame(c(temp.gd.4.1.100[["predictors"]],temp.gd.4.101.199[["predictors"]]))
colnames(erase.4gd) = c("trait")
mat_gd              = mat_gd %>% select((erase.4gd$trait))

# Best predictors statistics - step 5 #### 
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[1:100])
temp.gd.5.1.100     = ols_step_forward_p(m.all, details = FALSE)
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd[101:152])
temp.gd.5.101.152   = ols_step_forward_p(m.all, details = FALSE)
erase.5gd           = as.data.frame(c(temp.gd.5.1.100[["predictors"]],temp.gd.5.101.152[["predictors"]]))
colnames(erase.5gd) = c("trait")
mat_gd              = mat_gd %>% select((erase.5gd$trait))

# Best predictors statistics - step 6 #### 
m.all               = lm(mat.g.dro[,1297]~.,data=mat_gd)
temp.gd.5           = ols_step_forward_p(m.all, details = FALSE)
erase.6gd           = as.data.frame(c(temp.gd.5[["predictors"]]))
colnames(erase.6gd) = c("trait")
mat_gd              = mat_gd %>% select((erase.6gd$trait))

write.csv(mat_gd, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/grassland.drought.best.predictors.csv")

# Shrubland Drought ####

mat_trait_sd   = hmm_loma %>% mutate(shrub.abund.dro = mag_stat$Average.3)
mat.s.dro      = (as.matrix(mat_trait_gd[3:1299]))
# Best predictors statistics - step 1 ####
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[3:80])
temp.sd.3.80       = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[81:150])
temp.sd.81.150     = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[151:300])
temp.sd.151.300    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[301:450])
temp.sd.301.450    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[451:600])
temp.sd.451.600    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[601:700])
temp.sd.601.700    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[701:730])
temp.sd.701.730    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[731:740])
temp.sd.731.740    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[741:800])
temp.sd.741.800    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[801:900])
temp.sd.801.900    = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[901:1000])
temp.sd.901.1000   = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[1001:1100])
temp.sd.1001.1100  = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[1101:1200])
temp.sd.1101.1200  = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat.s.dro[,1297]~.,data=mat_trait_sd[1201:1296])
temp.sd.1201.1296  = ols_step_forward_p(m.all, details = FALSE)
erase.1sd          = as.data.frame(c(temp.sd.3.80[["predictors"]],
                                     temp.sd.81.150[["predictors"]],
                                     temp.sd.151.300[["predictors"]],
                                     temp.sd.301.450[["predictors"]],
                                     temp.sd.451.600[["predictors"]],
                                     temp.sd.601.700[["predictors"]],
                                     temp.sd.701.730[["predictors"]],
                                     temp.sd.731.740[["predictors"]],
                                     temp.sd.741.800[["predictors"]],
                                     temp.sd.801.900[["predictors"]],
                                     temp.sd.901.1000[["predictors"]],
                                     temp.sd.1001.1100[["predictors"]],
                                     temp.sd.1101.1200[["predictors"]],
                                     temp.sd.1201.1296[["predictors"]]))
colnames(erase.1sd) = c("trait")
mat_sd              = mat_trait_sd %>% select((erase.1sd$trait))

# Best predictors statistics - step 2 ####
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[1:100])
temp.sd.1.1.100      = ols_step_forward_p(m.all, details = FALSE)
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[101:200])
temp.sd.1.101.200    = ols_step_forward_p(m.all, details = FALSE)
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[201:300])
temp.sd.1.201.300    = ols_step_forward_p(m.all, details = FALSE)
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[301:415])
temp.sd.1.301.415    = ols_step_forward_p(m.all, details = FALSE)
erase.2sd            = as.data.frame(c(temp.sd.1.1.100[["predictors"]],
                                       temp.sd.1.101.200[["predictors"]],
                                       temp.sd.1.201.300[["predictors"]],
                                       temp.sd.1.301.415[["predictors"]]))
colnames(erase.2sd) = c("trait")
mat_sd              = mat_trait_sd %>% select((erase.2sd$trait))

# Best predictors statistics - step 3 ####
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[1:50])
temp.sd.2.1.50       = ols_step_forward_p(m.all, details = FALSE)
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[51:100])
temp.sd.2.51.100     = ols_step_forward_p(m.all, details = FALSE)
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[101:200])
temp.sd.2.101.200    = ols_step_forward_p(m.all, details = FALSE)
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[201:242])
temp.sd.2.201.242    = ols_step_forward_p(m.all, details = FALSE)
erase.3sd            = as.data.frame(c(temp.sd.2.1.50[["predictors"]],
                                       temp.sd.2.51.100[["predictors"]],
                                       temp.sd.2.101.200[["predictors"]],
                                       temp.sd.2.201.242[["predictors"]]))
colnames(erase.3sd) = c("trait")
mat_sd              = mat_sd %>% select((erase.3sd$trait))

# Best predictors statistics - step 4 ####
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[1:50])
temp.sd.3.1.50       = ols_step_forward_p(m.all, details = FALSE)
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd[51:199])
temp.sd.3.51.199     = ols_step_forward_p(m.all, details = FALSE)
erase.4sd            = as.data.frame(c(temp.sd.3.1.50[["predictors"]],
                                       temp.sd.3.51.199[["predictors"]]))
colnames(erase.4sd) = c("trait")
mat_sd              = mat_sd %>% select((erase.4sd$trait))                                       

# Best predictors statistics - step 5 ####
m.all                = lm(mat.s.dro[,1297]~.,data=mat_sd)
temp.sd.5            = ols_step_forward_p(m.all, details = FALSE)
erase.5sd            = as.data.frame(c(temp.sd.5[["predictors"]]))
colnames(erase.5sd) = c("trait")
mat_sd              = mat_sd %>% select((erase.5sd$trait)) 

write.csv(mat_sd, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/shrubland.drought.best.predictors.csv")

# FIRE DATA ####

# Calling data and preprocessing ####
hmm_fire    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/hmm_Fire.csv",dec=".")
mag_abun    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/mag_adundance_fire.txt") 
mag_abun    = mag_abun[-c(440,546), ]

# Low shallow ####
mat_trait_ls   = hmm_fire %>% mutate(L.shallow_abund = mag_abun$Avg_low_shallow)
mat_ls         = (as.matrix(mat_trait_ls[3:1215]))

# Best predictors statistics - step 1 ####
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[3:100])
temp.1.100     = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[101:300])
temp.101.300   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[301:400])
temp.301.400   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[401:500])
temp.401.500   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[501:700])
temp.501.700   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[701:750])
temp.701.750   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[751:850])
temp.751.850   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[851:950])
temp.851.950   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[951:1100])
temp.951.1100  = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat_trait_ls[1101:1212])
temp.1101.1212 = ols_step_forward_p(m.all, details = FALSE)
erase.1ls           = as.data.frame(c(temp.1.100[["predictors"]],
                                      temp.101.300[["predictors"]],
                                      temp.301.400[["predictors"]],
                                      temp.401.500[["predictors"]],
                                      temp.501.700[["predictors"]],
                                      temp.701.750[["predictors"]],
                                      temp.751.850[["predictors"]],
                                      temp.851.950[["predictors"]],
                                      temp.951.1100[["predictors"]],
                                      temp.1101.1212[["predictors"]]))
colnames(erase.1ls) = c("trait")
mat.ls              = mat_trait_ls %>% select((erase.1ls$trait))

# Best predictors statistics - step 2 ####
m.all          = lm(mat_ls[,1213]~.,data=mat.ls[1:100])
temp.1.1.100   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat.ls[101:200])
temp.1.101.200 = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat.ls[201:244])
temp.1.201.244 = ols_step_forward_p(m.all, details = FALSE)
erase.2ls      = as.data.frame(c(temp.1.1.100[["predictors"]],
                                 temp.1.101.200[["predictors"]],
                                 temp.1.201.244[["predictors"]]))
colnames(erase.2ls) = c("trait")
mat.ls              = mat_trait_ls %>% select((erase.2ls$trait))

# Best predictors statistics - step 3 ####
m.all          = lm(mat_ls[,1213]~.,data=mat.ls[1:100])
temp.2.1.100   = ols_step_forward_p(m.all, details = FALSE)
m.all          = lm(mat_ls[,1213]~.,data=mat.ls[101:159])
temp.2.101.159 = ols_step_forward_p(m.all, details = FALSE)
erase.3ls      = as.data.frame(c(temp.2.1.100[["predictors"]],
                                 temp.2.101.159[["predictors"]]))
colnames(erase.3ls) = c("trait")
mat.ls              = mat_trait_ls %>% select((erase.3ls$trait))

# Best predictors statistics - step 4 ####
m.all          = lm(mat_ls[,1213]~.,data=mat.ls)
temp.3         = ols_step_forward_p(m.all, details = FALSE)
erase.4ls      = as.data.frame(c(temp.3[["predictors"]]))
colnames(erase.4ls) = c("trait")
mat.ls              = mat_trait_ls %>% select((erase.4ls$trait))

write.csv(mat.ls, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/low.shallow.best.predictors.csv")

# Low deep ####
mat_trait_ld   = hmm_fire %>% mutate(L.deep_abund = mag_abun$Avg_low_deep)
mat_ld         = (as.matrix(mat_trait_ld[3:1215]))

# Best predictors statistics - step 1 ####
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[3:200])
temp.ld.1.200     = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[201:400])
temp.ld.201.400   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[401:430])
temp.ld.401.430   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[431:500])
temp.ld.431.500   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[501:700])
temp.ld.501.700   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[701:900])
temp.ld.701.900   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[901:1100])
temp.ld.901.1100  = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_ld[,1213]~.,data=mat_trait_ld[1101:1212])
temp.ld.1101.1212 = ols_step_forward_p(m.all, details = FALSE)
erase.1ld           = as.data.frame(c(temp.ld.1.200[["predictors"]],
                                      temp.ld.201.400[["predictors"]],
                                      temp.ld.401.430[["predictors"]],
                                      temp.ld.431.500[["predictors"]],
                                      temp.ld.501.700[["predictors"]],
                                      temp.ld.701.900[["predictors"]],
                                      temp.ld.901.1100[["predictors"]],
                                      temp.ld.1101.1212[["predictors"]]))
colnames(erase.1ld) = c("trait")
mat_ld              = mat_trait_ld %>% select((erase.1ld$trait))  

# Best predictors statistics - step 2 ####
m.all             = lm(mat_trait_ld[,1215]~.,data=mat_ld[1:100])
temp.ld.1.1.100   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_trait_ld[,1215]~.,data=mat_ld[101:220])
temp.ld.1.101.220 = ols_step_forward_p(m.all, details = FALSE)
erase.2ld         = as.data.frame(c(temp.ld.1.1.100[["predictors"]],
                                    temp.ld.1.101.220[["predictors"]]))
colnames(erase.2ld) = c("trait")
mat.ld              = mat_ld %>% select((erase.2ld$trait))

# Best predictors statistics - step 3 ####
m.all             = lm(mat_trait_ld[,1215]~.,data=mat.ld)
temp.ld.2         = ols_step_forward_p(m.all, details = FALSE)
erase.3ld         = as.data.frame(c(temp.ld.2[["predictors"]]))
colnames(erase.3ld) = c("trait")
mat.ld              = mat_ld %>% select((erase.3ld$trait))

write.csv(mat.ld, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/low.deep.best.predictors.csv")

# High shallow ####
mat_trait_hs   = hmm_fire %>% mutate(H.shallow_abund = mag_abun$Avg_high_shallow)
mat_hs         = (as.matrix(mat_trait_hs[3:1215]))

# Best predictors statistics - step 1 ####
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[3:200])
temp.hs.1.200     = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[201:300])
temp.hs.201.300   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[301:400])
temp.hs.301.400   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[401:500])
temp.hs.401.500   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[501:505])
temp.hs.501.505   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[506:509])
temp.hs.506.509   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[510])
temp.hs.510       = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[511:600])
temp.hs.511.600   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[601:700])
temp.hs.601.700   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[701:750])
temp.hs.701.750   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[751:800])
temp.hs.751.800   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[801:900])
temp.hs.801.900   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[901:1000])
temp.hs.901.1000  = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[1001:1100])
temp.hs.1001.1100 = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat_trait_hs[1101:1212])
temp.hs.1101.1212 = ols_step_forward_p(m.all, details = FALSE)
erase.1hs         = as.data.frame(c(temp.hs.1.200[["predictors"]],
                                    temp.hs.201.300[["predictors"]],
                                    temp.hs.301.400[["predictors"]],
                                    temp.hs.401.500[["predictors"]],
                                    temp.hs.501.505[["predictors"]],
                                    temp.hs.506.509[["predictors"]],
                                    temp.hs.511.600[["predictors"]],
                                    temp.hs.601.700[["predictors"]],
                                    temp.hs.701.750[["predictors"]],
                                    temp.hs.751.800[["predictors"]],
                                    temp.hs.801.900[["predictors"]],
                                    temp.hs.901.1000[["predictors"]],
                                    temp.hs.1001.1100[["predictors"]],
                                    temp.hs.1101.1212[["predictors"]]))
colnames(erase.1hs) = c("trait")
mat.ld              = mat_trait_hs %>% select((erase.1hs$trait))

# Best predictors statistics - step 2 ####
m.all             = lm(mat_hs[,1213]~.,data=mat.ld[1:100])
temp.hs.1.1.100   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.ld[101:110])
temp.hs.1.101.110 = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.ld[111:115])
temp.hs.1.111.115 = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.ld[116:118])
temp.hs.1.116.118 = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.ld[119:268])
temp.hs.1.119.268 = ols_step_forward_p(m.all, details = FALSE)
erase.2hs         = as.data.frame(c(temp.hs.1.1.100[["predictors"]],
                                    temp.hs.1.101.110[["predictors"]],
                                    temp.hs.1.111.115[["predictors"]],
                                    temp.hs.1.116.118[["predictors"]],
                                    temp.hs.1.119.268[["predictors"]]))
colnames(erase.2hs) = c("trait")
mat.hs              = mat_trait_hs %>% select((erase.2hs$trait))

# Best predictors statistics - step 3 ####
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[1:50])
temp.hs.2.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[51:80])
temp.hs.2.51.80   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[81:84])
temp.hs.2.81.84   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[85:155])
temp.hs.2.85.155  = ols_step_forward_p(m.all, details = FALSE)
erase.3hs         = as.data.frame(c(temp.hs.2.1.50[["predictors"]],
                                    temp.hs.2.51.80[["predictors"]],
                                    temp.hs.2.81.84[["predictors"]],
                                    temp.hs.2.85.155[["predictors"]]))
colnames(erase.3hs) = c("trait")
mat.hs              = mat_trait_hs %>% select((erase.3hs$trait))

# Best predictors statistics - step 4 ####
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[1:50])
temp.hs.3.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[51:70])
temp.hs.3.51.70   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[71:75])
temp.hs.3.71.75   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[76:145])
temp.hs.3.76.145  = ols_step_forward_p(m.all, details = FALSE)
erase.4hs         = as.data.frame(c(temp.hs.3.1.50[["predictors"]],
                                    temp.hs.3.51.70[["predictors"]],
                                    temp.hs.3.71.75[["predictors"]],
                                    temp.hs.3.76.145[["predictors"]]))
colnames(erase.4hs) = c("trait")
mat.hs              = mat.hs %>% select((erase.4hs$trait))

write.csv(mat.hs, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.shallow.best.predictors.csv")

# High deep ####
mat_trait_dh   = hmm_fire %>% mutate(H.deep_abund = mag_abun$Avg_high_deep)
mat_hd         = (as.matrix(mat_trait_dh[3:1215]))
# Best predictors statistics - step 1 ####



write.csv(mat.hs, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.deep.best.predictors.csv")



# NEON DATA - DECOMPOSERS ####

# Calling data and preprocessing ####
hmm_neon    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_neon/hmm_neon.csv",dec=".")






