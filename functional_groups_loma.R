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
# which traits are responsible for adaptations under drought in grassland in Loma

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

mat_trait   = mat_trait %>% mutate(rat_am_dro_gra.1 = mag_stat$ratio1.1)
mat_trait   = mat_trait %>% mutate(rat_am_dro_gra.1 = ifelse(is.na(rat_am_dro_gra.1),0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1==Inf,0,rat_am_dro_gra.1),
                                   rat_am_dro_gra.1 = ifelse(rat_am_dro_gra.1 ==-Inf, 0, rat_am_dro_gra.1))

mat_trait_1 = as.data.frame(mat_trait[,-192])
mat = (as.matrix(mat_trait_1[2:192]))

corre.1 = as.data.frame(cor(mat[,1:190], mat[,191]))

# Erase traits based on correlation matrix ####

test.1    = as.data.frame(cbind(mat_trait_1[2:191]))
seq.1     = seq(1,190)
temp.1    = as.data.frame(cbind((colSums(test.1)),seq.1))
erase.1   = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
test.2    = test.1[-(erase.1$seq.1)]
temp.2    = as.data.frame(colSums(test.2))
corre.1.1 = as.data.frame(cor(as.matrix(test.2), mat[,191]))

# Better prediction attempt  ####

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
corre.1.2 = as.data.frame(cor(as.matrix(test.3), mat[,191]))
test.3 = test.3[-59]
m.all.1 = lm(mat[,192]~.,data=test.3[1:55])
test.3  = test.3[-55]
m.all.1 = lm(mat[,192]~.,data=test.3)
test.3  = test.3[-58]
m.all.1 = lm(mat[,192]~.,data=test.3)
temp.58 = ols_step_forward_p(m.all.1, details = FALSE)

erase.3 = as.data.frame(temp.58[["metrics"]][["variable"]])
colnames(erase.3) = c("trait")

test.4  = test.3 %>% select((erase.3$trait))
# write.csv(test.4, file = "test.4.csv")
test.4  = read.csv("test.4.csv") 

corre.1.3 = as.data.frame(cor(as.matrix(test.4[2:38]), mat[,191]))
# 37 traits out of 192 are the best predictors for changes in MAGs abundance 
# under drought

# Aim 2 ####
# Do we see any functional groups in Loma under drought and which traits drive
# these functional groups

trait_tar = test.4[2:38]*mag_stat$Average.1
genomes   = as.list(mat_trait[1])
test      = as.data.frame(rowSums(trait_tar))
test$row_names = row.names(test)
erase.n   = test %>% filter(`rowSums(trait_tar)`==0)
trait_tar$row_names = row.names(trait_tar)
trait_tar = trait_tar[!(trait_tar$row_names %in% erase.n$row_names),]

# Forming functional guilds 

library(vegan)
set.seed(1)

# Gene cost
distanceg  = vegdist(trait_tar[1:37], method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar), 2)
plot(clusterg)

# Similiarity within guilds
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
  geom_hline(yintercept=0.6503766, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 6, linetype="dashed", color = "red")

# Working with 6 guilds





