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
library(factoextra)
library(igraph)

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

# Selected genes ####

shrubland.drought = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/shrubland.drought.best.predictors.csv")
shrubland.ambient = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/shrubland.ambient.best.predictors.csv")
grassland.drought = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/grassland.drought.best.predictors.csv")
grassland.ambient = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/grassland.ambient.best.predictors.csv")
gene.select       = as.data.frame(unique(c(colnames(shrubland.drought[2:109]),
                                    colnames(shrubland.ambient[2:95]),
                                    colnames(grassland.drought[2:100]),
                                    colnames(grassland.ambient[2:127]))))
colnames(gene.select)        = c("gene")
mat.gene.loma                = as.data.frame(cbind(hmm_loma$id,
                                            hmm_loma %>% select(gene.select$gene)))
colnames(mat.gene.loma)[1]   = "loma.id"

write.csv(mat.gene.loma, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/loma.genes.selected.csv")

# Functional groups ####

# Hierarchical Clustering 
set.seed(1)
distance.loma  = vegdist(mat.gene.loma[2:292], method = "jaccard", binary = TRUE)
cluster.loma   = hclust(distance.loma, method="ward.D2")

# similarity within guilds ####
my_vec       = c()
nguilds.1    = seq(2, nrow(mat.gene.loma), 2)

for(i in nguilds.1) {
  v                      = cutree(cluster.loma,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distance.loma ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1
  my_vec   <- c(my_vec, my_out)   
}
mat_r2_1  = as.data.frame(cbind(nguilds.1,my_vec))

# Similarity among guilds ####
v                      = cutree(cluster.loma,k=33)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.loma,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
test.clus.1            = adonis.pair(distance.loma, genome2guild[,"guild"], nper = 999, 
                                     corr.method = "fdr")
adonis_2               = vegan::adonis2(distance.loma ~ guild, data = genome2guild, perm = 999)

# Plotting clustering ####
fviz_dend(cluster.loma, k = 33,                 # Cut in four groups
          cex = 0.25,
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
)

# Plotting similarity ####
ggplot(data=mat_r2_1,aes(x=nguilds.1,y=my_vec)) + geom_line() +
  xlab("Number of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 16)) + 
  geom_hline(yintercept =0.5148247, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 33, linetype="dashed", color = "red") +
  theme_classic()

# Ordination of the 43 functional groups ####
mat.gene.loma          = as.data.frame(cbind(genome2guild,mat.gene.loma))
write.csv(mat.gene.loma, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/loma.genes.guilds.selected.csv")
set.seed(1992)
ordination_1           = metaMDS(mat.gene.loma[3:293],autotransform = F,trymax = 500)
fort.1                 = fortify(ordination_1)
write.csv(fort.1, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/loma.fort.1.csv")

# Ordination 1 ####
fort.1                 = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/loma.fort.1.csv")
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                           mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.loma$guild),size = 2),
                           alpha=0.5) + 
  geom_segment(data=subset(fort.1,score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort.1,score=="species"), # "species"
            mapping=aes(label=label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2)))

# Ordination 2 ####
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.loma$guild),size = 2),
                      alpha=0.5) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2))) + 
  stat_ellipse(data = subset(fort.1,score=="sites"), 
               aes(x = NMDS1, y = NMDS2, color = as.factor(mat.gene.loma$guild)))

# Abundance of functional groups under conditions ####
fg_abundance                = as.data.frame(cbind(mat.gene.loma$guild,mag_stat))
colnames(fg_abundance)[1]   = "guild"

grassland.a    = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Average)/sum(fg_abundance$Average)) %>% 
  mutate(condition = rep("grassland.ambient" , nrow(grassland.a)))
grassland.d    = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Average.1)/sum(fg_abundance$Average.1)) %>% 
  mutate(condition = rep("grassland.drought" , nrow(grassland.a)))
shrubland.a    = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Average.2)/sum(fg_abundance$Average.2)) %>% 
  mutate(condition = rep("shrubland.ambient" , nrow(grassland.a)))
shrubland.d    = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Average.3)/sum(fg_abundance$Average.3)) %>% 
  mutate(condition = rep("shrubland.drought" , nrow(grassland.a)))

fg_ab.fig      =  as.data.frame(rbind(grassland.a,grassland.d,shrubland.a,
                                      shrubland.d))
colnames(fg_ab.fig) = c("guild","abundance","condition")
write.csv(fg_ab.fig, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/fg_ab.loma.csv")

ggplot(fg_ab.fig, aes(fill=as.factor(guild), y=abundance, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  theme(text = element_text(size=25)) + 
  labs(y="Abundance",x = element_blank()) + 
  theme(legend.title = element_blank())

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

# Best predictors statistics - step 5 ####
mat.hs            = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.shallow.best.predictors.csv",dec=".")
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[2:50])
temp.hs.5         = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[51:76])
temp.hs.5.1       = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[77:144])
temp.hs.5.2       = ols_step_forward_p(m.all, details = FALSE)
erase.5hs         = as.data.frame(c(temp.hs.5[["predictors"]],
                                    temp.hs.5.1[["predictors"]],
                                    temp.hs.5.2[["predictors"]]))
colnames(erase.5hs) = c("trait")
mat.hs              = mat.hs %>% select((erase.5hs$trait))

# Best predictors statistics - step 6 ####
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[1:62])
temp.hs.6         = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hs[,1213]~.,data=mat.hs[63:128])
temp.hs.6.1       = ols_step_forward_p(m.all, details = FALSE)
erase.6hs         = as.data.frame(c(temp.hs.6[["predictors"]],
                                    temp.hs.6.1[["predictors"]]))
colnames(erase.6hs) = c("trait")
mat.hs              = mat.hs %>% select((erase.6hs$trait))

write.csv(mat.hs, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.shallow.best.predictors.csv")

# High deep ####
mat_trait_dh   = hmm_fire %>% mutate(H.deep_abund = mag_abun$Avg_high_deep)
mat_hd         = (as.matrix(mat_trait_dh[3:1215]))

# Best predictors statistics - step 1 ####
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[3:50])
temp.hd.1.50      = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[51:100])
temp.hd.51.100    = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[101:110])
temp.hd.101.110   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[111:115])
temp.hd.111.115   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[116:150])
temp.hd.116.150   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[151:160])
temp.hd.151.160   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[161:170])
temp.hd.161.170   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[171:174])
temp.hd.171.174   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[175:176])
temp.hd.175.176   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[177:180])
temp.hd.177.180   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[181:200])
temp.hd.181.200   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[201:300])
temp.hd.201.300   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[301:400])
temp.hd.301.400   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[401:410])
temp.hd.401.410   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[411:414])
temp.hd.411.414   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[415:500])
temp.hd.415.500   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[501:520])
temp.hd.501.520   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[521:650])
temp.hd.521.650   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[651:800])
temp.hd.651.800   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[801:900])
temp.hd.801.900   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[901:950])
temp.hd.901.950   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[951:1050])
temp.hd.951.1050  = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat_trait_dh[1051:1212])
temp.hd.1051.1212 = ols_step_forward_p(m.all, details = FALSE)
erase.1hd         = as.data.frame(c(temp.hd.1.50[["predictors"]],
                                    temp.hd.51.100[["predictors"]],
                                    temp.hd.101.110[["predictors"]],
                                    temp.hd.111.115[["predictors"]],
                                    temp.hd.116.150[["predictors"]],
                                    temp.hd.151.160[["predictors"]],
                                    temp.hd.161.170[["predictors"]],
                                    temp.hd.171.174[["predictors"]],
                                    temp.hd.175.176[["predictors"]],
                                    temp.hd.177.180[["predictors"]],
                                    temp.hd.181.200[["predictors"]],
                                    temp.hd.201.300[["predictors"]],
                                    temp.hd.301.400[["predictors"]],
                                    temp.hd.401.410[["predictors"]],
                                    temp.hd.411.414[["predictors"]],
                                    temp.hd.415.500[["predictors"]],
                                    temp.hd.501.520[["predictors"]],
                                    temp.hd.521.650[["predictors"]],
                                    temp.hd.651.800[["predictors"]],
                                    temp.hd.801.900[["predictors"]],
                                    temp.hd.901.950[["predictors"]],
                                    temp.hd.951.1050[["predictors"]],
                                    temp.hd.1051.1212[["predictors"]]))
colnames(erase.1hd) = c("trait")
mat.hd              = mat_trait_dh %>% select((erase.1hd$trait))

# Best predictors statistics - step 2 ####
mat.hd            = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.deep.best.predictors.csv",dec=".")
m.all             = lm(mat_hd[,1213]~.,data=mat.hd[2:100])
temp.hd.1.1.100   = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat.hd[101:200])
temp.hd.1.101.200 = ols_step_forward_p(m.all, details = FALSE)
m.all             = lm(mat_hd[,1213]~.,data=mat.hd[201:274])
temp.hd.1.201.274 = ols_step_forward_p(m.all, details = FALSE)
erase.2hd         = as.data.frame(c(temp.hd.1.1.100[["predictors"]],
                                    temp.hd.1.101.200[["predictors"]],
                                    temp.hd.1.201.274[["predictors"]]))
colnames(erase.2hd) = c("trait")
mat.hd              = mat.hd %>% select((erase.2hd$trait))

# Best predictors statistics - step 3 ####
m.all             = lm(mat_hd[,1213]~.,data=mat.hd[1:163])
temp.hd.2         = ols_step_forward_p(m.all, details = FALSE)
erase.3hd         = as.data.frame(c(temp.hd.2[["predictors"]]))
colnames(erase.3hd) = c("trait")
mat.hd              = mat.hd %>% select((erase.3hd$trait))

write.csv(mat.hd, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.deep.best.predictors.csv")

# Selected genes ####

low.shallow  = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/low.shallow.best.predictors.csv")
low.deep     = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/low.deep.best.predictors.csv")
high.shallow = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.deep.best.predictors.csv")
high.deep    = read.csv(file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/high.shallow.best.predictors.csv")
gene.select  = as.data.frame(unique(c(colnames(low.shallow[2:109]),
                                           colnames(low.deep[2:84]),
                                           colnames(high.shallow[2:76]),
                                           colnames(high.deep[2:121]))))
colnames(gene.select)        = c("gene")
mat.gene.fire                = as.data.frame(cbind(hmm_fire$id,
                                                   hmm_fire %>% select(gene.select$gene)))
colnames(mat.gene.fire)[1]   = "loma.id"

write.csv(mat.gene.fire, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fire.genes.selected.csv")

# Functional groups ####

# Hierarchical Clustering 
set.seed(1)
distance.fire  = vegdist(mat.gene.fire[2:309], method = "jaccard", binary = TRUE)
cluster.fire   = hclust(distance.fire, method="ward.D2")

# similarity within guilds ####
my_vec       = c()
nguilds.1    = seq(2, nrow(mat.gene.fire), 2)

for(i in nguilds.1) {
  v                      = cutree(cluster.fire,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distance.fire ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1
  my_vec   <- c(my_vec, my_out)   
}
mat_r2_1  = as.data.frame(cbind(nguilds.1,my_vec))

# Similarity among guilds ####
v                      = cutree(cluster.fire,k=42)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.fire,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
test.clus.1            = adonis.pair(distance.fire, genome2guild[,"guild"], nper = 999, 
                                     corr.method = "fdr")
adonis_2               = vegan::adonis2(distance.fire ~ guild, data = genome2guild, perm = 999)

# Plotting clustering ####
fviz_dend(cluster.fire, k = 33,                 # Cut in four groups
          cex = 0.25,
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
)

# Plotting similarity ####
ggplot(data=mat_r2_1,aes(x=nguilds.1,y=my_vec)) + geom_line() +
  xlab("Number of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 16)) + 
  geom_hline(yintercept =0.5560235, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 42, linetype="dashed", color = "red") +
  theme_classic()

# Ordination of the 43 functional groups ####
mat.gene.fire          = as.data.frame(cbind(genome2guild,mat.gene.fire))
write.csv(mat.gene.fire, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fire.genes.guilds.selected.csv")
set.seed(1)
ordination_1           = metaMDS(mat.gene.fire[3:310],autotransform = T,trymax = 500)
fort.1                 = fortify(ordination_1)
write.csv(fort.1, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fire.fort.1.csv")

# Ordination 1 ####
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.fire$guild),size = 2),
                      alpha=0.5) + 
  geom_segment(data=subset(fort.1,score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort.1,score=="species"), # "species"
            mapping=aes(label=label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2)))

# Ordination 2 ####
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.fire$guild),size = 2),
                      alpha=0.5) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2))) + 
  stat_ellipse(data = subset(fort.1,score=="sites"), 
               aes(x = NMDS1, y = NMDS2, color = as.factor(mat.gene.fire$guild)))

# Abundance of functional groups under conditions ####
fg_abundance                = as.data.frame(cbind(mat.gene.fire$guild,mag_abun))
colnames(fg_abundance)[1]   = "guild"

low_shallow    = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Avg_low_shallow)/sum(fg_abundance$Avg_low_shallow)) %>% 
  mutate(condition = rep("low_shallow" , nrow(low_shallow)))
low_deep       = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Avg_low_deep)/sum(fg_abundance$Avg_low_deep)) %>% 
  mutate(condition = rep("low_deep" , nrow(low_deep)))
high_deep      = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Avg_high_deep)/sum(fg_abundance$Avg_high_deep)) %>% 
  mutate(condition = rep("high_deep" , nrow(high_deep)))
high_shallow   = fg_abundance %>% group_by(guild) %>% 
  summarise(abundance = sum(Avg_high_shallow)/sum(fg_abundance$Avg_high_shallow)) %>% 
  mutate(condition = rep("high_shallow" , nrow(high_shallow)))

fg_ab.fig      =  as.data.frame(rbind(low_shallow,low_deep,high_deep,
                                      high_shallow))
colnames(fg_ab.fig) = c("guild","abundance","condition")
write.csv(fg_ab.fig, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fg_ab.fire.csv")

ggplot(fg_ab.fig, aes(fill=as.factor(guild), y=abundance, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  theme(text = element_text(size=25)) + 
  labs(y="Abundance",x = element_blank()) + 
  theme(legend.title = element_blank())

# LOMA-FIRE ####

# Calling data and preprocessing ####
hmm_loma         = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/hmm_Loma_full.csv",dec=".")
gene_loma        = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/litter_mags_metadata.txt",dec=".")
mag_stat         = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_stats.txt") 
mag_abun.loma    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_adundance.txt") 
mag_stat         = mag_stat %>% full_join(mag_abun.loma)
hmm_fire         = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/hmm_Fire_full.csv",dec=".")
mag_abun.fire    = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/mag_adundance_fire.txt") 
mag_abun.fire    = mag_abun.fire[-c(440,546), ]

# Selected genes ####
mat.gene.loma = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/loma.genes.selected.csv",dec=".")
mat.gene.fire = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fire.genes.selected.csv",dec=".")
T.gene.select = as.data.frame(unique(c(colnames(mat.gene.loma[3:293]),
                                       colnames(mat.gene.fire[3:310]))))
colnames(T.gene.select)      = c("gene")
temp.fire                    = as.data.frame(cbind(hmm_fire$id,hmm_fire %>% select(T.gene.select$gene)))
temp.loma                    = as.data.frame(cbind(hmm_loma$id,hmm_loma %>% select(T.gene.select$gene)))
colnames(temp.fire)[1]       = c("id")
colnames(temp.loma)[1]       = c("id")
total.mags                   = as.data.frame(rbind(temp.fire,temp.loma))

# Hierarchical Clustering 
set.seed(1)
distance.total  = vegdist(total.mags[2:506], method = "jaccard", binary = TRUE)
cluster.total   = hclust(distance.total, method="ward.D2")

# similarity within guilds ####
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
v                      = cutree(cluster.total,k=63)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.total,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
test.clus.1            = adonis.pair(distance.total, genome2guild[,"guild"], nper = 999, 
                                     corr.method = "fdr")
adonis_2               = vegan::adonis2(distance.total ~ guild, data = genome2guild, perm = 999)

# Plotting clustering ####
fviz_dend(cluster.total, k = 63,                 # Cut in four groups
          cex = 0.25,
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
)

# Plotting similarity ####
ggplot(data=mat_r2_1,aes(x=nguilds.1,y=my_vec)) + geom_line() +
  xlab("Number of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 16)) + 
  geom_hline(yintercept =0.5455844, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 63, linetype="dashed", color = "red") +
  theme_classic()

# Ordination of the 63 functional groups ####
mat.gene.total         = as.data.frame(cbind(genome2guild,total.mags))
write.csv(mat.gene.total, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/total.genes.guilds.selected.csv")
set.seed(1)
ordination_1           = metaMDS(mat.gene.total[3:507],autotransform = T,trymax = 5000)
fort.1                 = fortify(ordination_1)
write.csv(fort.1, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/total.fort.1.csv")

# Ordination 1 ####
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.total$guild),size = 2),
                      alpha=0.5) + 
  geom_segment(data=subset(fort.1,score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort.1,score=="species"), # "species"
            mapping=aes(label=label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2)))

# Ordination 2 ####
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.total$guild),size = 2),
                      alpha=0.5) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2))) + 
  stat_ellipse(data = subset(fort.1,score=="sites"), 
               aes(x = NMDS1, y = NMDS2, color = as.factor(mat.gene.total$guild)))

# Ordination 3 ####

mat.gene.total         = mat.gene.total %>% mutate(condition = c(rep("fire" , nrow(mag_abun.fire)),
                                                                     rep("loma" , nrow(mag_abun.loma))))
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.total$guild),shape=as.factor(mat.gene.total$condition),size = 2),
                      alpha=0.5) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2)))

# Abundance of functional groups under conditions ####

# Microbial network ####
myedgeslist <- data.frame(to = mat.gene.total$guild,
                          from = mat.gene.total$condition)
mygraph <- myedgeslist %>% igraph::graph_from_data_frame(directed = T) 
mygraph %>% igraph::plot.igraph()


# NEON DATA - DECOMPOSERS ####

# Calling data and preprocessing ####
hmm_neon       = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_neon/hmm_neon.csv",dec=".")
col.names      = as.data.frame(colnames(hmm_neon))

# Selecting decomposers genes ####
hmm_neon_decom = as.data.frame(cbind(hmm_neon$id,hmm_neon[,1793:1833]))
colnames(hmm_neon_decom)[1]   = "id"
temp.1         = as.data.frame(cbind((colSums(hmm_neon_decom[,2:42]))))
erase.1        = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
hmm_neon_decom = hmm_neon_decom[ , !(names(hmm_neon_decom) %in% erase.1$row_names)]
a              = as.data.frame(rowSums(hmm_neon_decom[2:40]))
hmm_neon_decom = hmm_neon_decom[-c(76,308,378,473,531),]

# Hierarchical Clustering 
set.seed(1)
distance.neon  = vegdist(hmm_neon_decom[2:40], method = "jaccard", binary = TRUE)
cluster.neon   = hclust(distance.neon, method="ward.D2")

# similarity within guilds ####
my_vec       = c()
nguilds.1    = seq(2, nrow(hmm_neon_decom), 2)

for(i in nguilds.1) {
  v                      = cutree(cluster.neon,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distance.neon ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1
  my_vec   <- c(my_vec, my_out)   
}
mat_r2_1  = as.data.frame(cbind(nguilds.1,my_vec))

# Similarity among guilds ####
v                      = cutree(cluster.neon,k=48)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distance.neon,genome2guild$guild,
                                                         perm = 999,p.adjust.m='BH')
test.clus.1            = adonis.pair(distance.neon, genome2guild[,"guild"], nper = 999, 
                                     corr.method = "fdr")
adonis_2               = vegan::adonis2(distance.neon ~ guild, data = genome2guild, perm = 999)

# Plotting clustering ####
fviz_dend(cluster.neon, k = 48,                 # Cut in four groups
          cex = 0.25,
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
)

# Plotting similarity ####
ggplot(data=mat_r2_1,aes(x=nguilds.1,y=my_vec)) + geom_line() +
  xlab("Number of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 16)) + 
  geom_hline(yintercept =0.6104834, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 48, linetype="dashed", color = "red") +
  theme_classic()

# Ordination of the 48 functional groups ####
mat.gene.neon         = as.data.frame(cbind(genome2guild,hmm_neon_decom))
write.csv(mat.gene.neon, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_neon/total.genes.guilds.selected.csv")
set.seed(1)
ordination_1           = metaMDS(mat.gene.neon[3:41],autotransform = T,trymax = 100)
fort.1                 = fortify(ordination_1)

# Ordination 1 ####
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.neon$guild),size = 2),
                      alpha=0.5) + 
  geom_segment(data=subset(fort.1,score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort.1,score=="species"), # "species"
            mapping=aes(label=label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2)))

# Ordination 2 ####
ggplot() + geom_point(data=subset(fort.1,score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat.gene.neon$guild),size = 2),
                      alpha=0.5) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2))) + 
  stat_ellipse(data = subset(fort.1,score=="sites"), 
               aes(x = NMDS1, y = NMDS2, color = as.factor(mat.gene.neon$guild)))

