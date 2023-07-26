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

mag_stat = mag_stat %>% mutate(ratio1.1 = log10((mag_stat$Average.1/mag_stat$Average)*
                                 abs(mag_stat$Average.1-mag_stat$Average)))
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

m.all = lm(mat[,191]~.,data=test.2[1:50])
temp.1.50 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,191]~.,data=test.2[51:100])
temp.51.100 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,191]~.,data=test.2[101:115])
temp.101.115 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,191]~.,data=test.2[117:127])
temp.117.127 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,191]~.,data=test.2[128:137])
temp.128.137 = ols_step_forward_p(m.all, details = FALSE)
m.all = lm(mat[,191]~.,data=test.2[138:140])
temp.138.140 = ols_step_forward_p(m.all, details = FALSE)

erase.2 = as.data.frame(c(temp.1.50[["metrics"]][["variable"]],
                          temp.51.100[["metrics"]][["variable"]],
                          temp.101.115[["metrics"]][["variable"]],
                          temp.117.127[["metrics"]][["variable"]],
                          temp.128.137[["metrics"]][["variable"]],
                          temp.138.140[["metrics"]][["variable"]]))
colnames(erase.2) = c("trait")

test.3 = test.2 %>% select((erase.2$trait))
corre.1.2 = as.data.frame(cor(as.matrix(test.3), mat[,191]))
test.3 = test.3[-56]
m.all.1 = lm(mat[,191]~.,data=test.3)
temp.72 = ols_step_forward_p(m.all.1, details = FALSE)
erase.3 = as.data.frame(temp.72[["metrics"]][["variable"]])
colnames(erase.3) = c("trait")

test.4  = test.3 %>% select((erase.3$trait))
#write.csv(test.4, file = "test.4.csv")
test.4  = read.csv("test.4.csv") 

corre.1.3 = as.data.frame(cor(as.matrix(test.4[2:43]), mat[,191]))
# 42 traits out of 192 are the best predictors for changes in MAGs abundance 
# under drought
# write.csv(corre.1.3, file = "corre.1.3.csv")
corre.1.3  = read.csv("corre.1.3.csv") 

matgs = as.matrix(cbind(mat[,191],test.4[2:43]))
x2 <- matgs %>% 
  correlate() %>% 
  focus(`mat[, 191]`)
x2 %>% 
  mutate(rowname = factor(term, levels = term[order(`mat[, 191]`)])) 

gran_1 = c("S.Tolerance","R.Acquisition","R.Acquisition","R.Acquisition",
           "S.Tolerance","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","S.Tolerance",
           "R.Use","R.Acquisition","R.Use","S.Tolerance","S.Tolerance",
           "R.Acquisition","R.Use","R.Acquisition","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Use",
           "R.Acquisition","S.Tolerance","R.Use","R.Acquisition",
           "R.Acquisition","R.Acquisition","R.Acquisition","R.Acquisition",
           "R.Use","S.Tolerance","R.Acquisition","S.Tolerance","R.Acquisition")

gran_2 = c("Desiccation","S.uptake","S.uptake","S.assimilation",
           "Temp.low","S.uptake","S.uptake","S.degradation",
           "S.assimilation","S.assimilation","S.uptake","S.assimilation",
           "S.uptake","S.assimilation","S.degradation","Protection",
           "Phototrophy","S.uptake","Chemotrophy","Temp.low","pH.stress",
           "S.uptake","Chemotrophy","S.degradation","S.degradation",
           "S.assimilation","S.degradation","S.assimilation","Chemotrophy",
           "S.uptake","O.limitation","Chemotrophy","S.uptake","S.uptake",
           "S.degradation","S.uptake","S.degradation","Phototrophy",
           "Phage.resistance","S.uptake","Temp.low","S.degradation")

x2 = as.data.frame(cbind(gran_1,gran_2,x2))
colnames(x2) = c("trait.level.1","trait.level.2","term","pearson")
g <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.1)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g
g1 <- ggplot(x2, aes(x = term, y = pearson,fill = trait.level.2)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank())
g1 + scale_fill_manual(values=c("#a6d96a","#c6dbef",
                                "#6baed6","#41b6c4","#0570b0","#006d2c",
                                "#253494","#f46d43","#d73027","#a50026",
                                "#08306b"))
# Aim 2 ####
# Do we see any functional groups in Loma under drought and which traits drive
# these functional groups: average MAGs abundance
#library(ggpubr)
test.4    = read.csv("test.4.csv") 
ratio     = as.data.frame((mag_stat$Average.1/mag_stat$Average))
colnames(ratio) = "ratio"
test.4    = as.data.frame(cbind(test.4,ratio))
test.4    = as.data.frame(test.4[c(-163,-164,-263,-299,-387,-396,-399,-425,-479,
                  -485,-487,-491,-499,-506,-511),])
#ggqqplot(log10(test.4[,44]+1))
# Z-transformation ####
trait_tar = scale(test.4[,2:44])
genomes   = as.list(mat_trait[1])
test      = as.data.frame(rowSums(trait_tar))
test$row_names = row.names(test)

# Forming functional guilds 

library(vegan)
set.seed(1)

# Gene cost
distanceg  = vegdist(trait_tar, method = "euclidean", binary = FALSE)
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
# write.csv(mat_r2_c, file = "mat_r2_c.csv")
# mat_r2_c  = read.csv("mat_r2_c.csv") 
ggplot(data=mat_r2_c,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 10)) + 
  geom_hline(yintercept=0.32582104, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 10, linetype="dashed", color = "red")

# Working with 6 guilds

v                        = cutree(clusterg,k=10) # Clusters
genome2guild_5           = data.frame(guild = factor(v))
rownames(genome2guild_5) = names(v)
mat_ori$row_names        = row.names(mat_ori)
mat_ori                  = as.data.frame(mat_ori[c(-163,-164,-263,-299,-387,-396,-399,-425,-479,
                                                  -485,-487,-491,-499,-506,-511),])
mat_ori.1                = as.data.frame(cbind(mat_ori,test.4$ratio))
mat_trait_5              = as.data.frame(cbind(genome2guild_5,mat_ori.1$id,trait_tar)) 

temp.6 = mat_trait_5 %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat.6  = (as.matrix(temp.6[2:44]))
a      = as.data.frame(colnames(mat.6))
library(ggvegan)
set.seed(2)
temp_1 = metaMDS(trait_tar)

fort = fortify(temp_1)
# write.csv(fort, file = "fort.csv")
# fort  = read.csv("fort.csv") 
ggplot() + geom_point(data=subset(fort,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,size = 4),
                      alpha=0.5) + 
  geom_segment(data=subset(fort,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score==""), # "species"
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(temp_1$stress,2)))

# Aim 3 ####
# Do we see changes in functional groups in Loma under drought over time

# Time 1: grassRedDry1 ####

trait_tar_0 = test.4[2:38]*mag_stat$grassRedDry1
genomes   = as.list(mat_trait[1])
test      = as.data.frame(rowSums(trait_tar_0))
test$row_names = row.names(test)
erase.n   = test %>% filter(`rowSums(trait_tar_0)`==0)
trait_tar_0$row_names = row.names(trait_tar_0)
trait_tar_0 = trait_tar_0[!(trait_tar_0$row_names %in% erase.n$row_names),]

# Forming functional guilds 

library(vegan)
set.seed(1)

# Gene cost
distanceg  = vegdist(trait_tar_0[1:37], method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar_0), 2)
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
mat_r2_c.0  = as.data.frame(cbind(nguilds,my_vec))
# write.csv(mat_r2_c.0, file = "mat_r2_c.0.csv")
mat_r2_c.0  = read.csv("mat_r2_c.0.csv")

ggplot(data=mat_r2_c.0,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 20)) + 
  geom_hline(yintercept=0.6678939, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 6, linetype="dashed", color = "red")

# Working with 6 guilds

v                        = cutree(clusterg,k=6) # Clusters
genome2guild_5           = data.frame(guild = factor(v))
rownames(genome2guild_5) = names(v)
mat_ori$row_names        = row.names(mat_ori)
mat_ori.1                = mat_ori[!(mat_ori$row_names %in% erase.n$row_names),]
mat_trait_5              = as.data.frame(cbind(genome2guild_5,mat_ori.1$id,trait_tar_0)) 

temp.6 = mat_trait_5 %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat.6  = (as.matrix(temp.6[2:38]))
a      = as.data.frame(colnames(mat.6))
library(ggvegan)
set.seed(2)
temp_1 = metaMDS(trait_tar_0[1:37],autotransform = T,
                 trymax = 1000,noshare=0.1)

fort.0 = fortify(temp_1)
# write.csv(fort.0, file = "fort.0.csv")
fort.0  = read.csv("fort.0.csv")

ggplot() + geom_point(data=subset(fort.0,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =mat_trait_5$guild,size = 4),
                      alpha=0.5) + 
  geom_segment(data=subset(fort.0,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort,Score=="species"), # "species"
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(temp_1$stress,2)))

# Time 2: grassRedWet1 ####

trait_tar_1 = test.4[2:38]*mag_stat$grassRedWet1
genomes   = as.list(mat_trait[1])
test      = as.data.frame(rowSums(trait_tar_1))
test$row_names = row.names(test)
erase.n   = test %>% filter(`rowSums(trait_tar_1)`==0)
trait_tar_1$row_names = row.names(trait_tar_1)
trait_tar_1 = trait_tar_1[!(trait_tar_1$row_names %in% erase.n$row_names),]

# Forming functional guilds 

library(vegan)
set.seed(1)

# Gene cost
distanceg  = vegdist(trait_tar_1[1:37], method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar_1), 2)
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
  geom_hline(yintercept=0.6687614, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 6, linetype="dashed", color = "red")

# Working with 6 guilds

v                        = cutree(clusterg,k=6) # Clusters
genome2guild_5           = data.frame(guild = factor(v))
rownames(genome2guild_5) = names(v)
mat_ori$row_names        = row.names(mat_ori)
mat_ori.1                = mat_ori[!(mat_ori$row_names %in% erase.n$row_names),]
mat_trait_5              = as.data.frame(cbind(genome2guild_5,mat_ori.1$id,trait_tar_1)) 

temp.6 = mat_trait_5 %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat.6  = (as.matrix(temp.6[2:38]))
a      = as.data.frame(colnames(mat.6))
library(ggvegan)
set.seed(2)
temp_1 = metaMDS(trait_tar_1[1:37],autotransform = T,
                 trymax = 1000,noshare=0.1)

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
  geom_text(data=subset(fort,Score==""), # "species"
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(temp_1$stress,2)))

# Time 3: grassRedDry2 ####

trait_tar_2 = test.4[2:38]*mag_stat$grassRedDry2
genomes   = as.list(mat_trait[1])
test      = as.data.frame(rowSums(trait_tar_2))
test$row_names = row.names(test)
erase.n   = test %>% filter(`rowSums(trait_tar_2)`==0)
trait_tar_2$row_names = row.names(trait_tar_2)
trait_tar_2 = trait_tar_2[!(trait_tar_2$row_names %in% erase.n$row_names),]

# Forming functional guilds 

library(vegan)
set.seed(1)

# Gene cost
distanceg  = vegdist(trait_tar_2[1:37], method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar_2), 2)
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
mat_r2_c.4  = as.data.frame(cbind(nguilds,my_vec))
# write.csv(mat_r2_c.4, file = "mat_r2_c.4.csv")
mat_r2_c.4  = read.csv("mat_r2_c.4.csv")
ggplot(data=mat_r2_c.4,aes(x=nguilds,y=my_vec)) + geom_line() +
  xlab("# of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 20)) + 
  geom_hline(yintercept=0.6695554, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 6, linetype="dashed", color = "red")

# Working with 6 guilds

v                        = cutree(clusterg,k=6) # Clusters
genome2guild_5           = data.frame(guild = factor(v))
rownames(genome2guild_5) = names(v)
mat_ori$row_names        = row.names(mat_ori)
mat_ori.1                = mat_ori[!(mat_ori$row_names %in% erase.n$row_names),]
mat_trait_5              = as.data.frame(cbind(genome2guild_5,mat_ori.1$id,trait_tar_2)) 

temp.6 = mat_trait_5 %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat.6  = (as.matrix(temp.6[2:38]))
a      = as.data.frame(colnames(mat.6))
library(ggvegan)
set.seed(2)
temp_1 = metaMDS(trait_tar_2[1:37],autotransform = T,
                 trymax = 1000,noshare=0.1)

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
  geom_text(data=subset(fort,Score==""), # "species"
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(temp_1$stress,2)))

# Time 4: grassRedWet2

trait_tar_3 = test.4[2:38]*mag_stat$grassRedWet2
genomes   = as.list(mat_trait[1])
test      = as.data.frame(rowSums(trait_tar_3))
test$row_names = row.names(test)
erase.n   = test %>% filter(`rowSums(trait_tar_3)`==0)
trait_tar_3$row_names = row.names(trait_tar_3)
trait_tar_3 = trait_tar_3[!(trait_tar_3$row_names %in% erase.n$row_names),]

# Forming functional guilds 

library(vegan)
set.seed(1)

# Gene cost
distanceg  = vegdist(trait_tar_3[1:37], method = "bray", binary = FALSE)
clusterg   = hclust(distanceg, method="ward")
nguilds    = seq(2, nrow(trait_tar_3), 2)
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
  geom_hline(yintercept=0.6389547, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 6, linetype="dashed", color = "red")

# Working with 6 guilds

v                        = cutree(clusterg,k=6) # Clusters
genome2guild_5           = data.frame(guild = factor(v))
rownames(genome2guild_5) = names(v)
mat_ori$row_names        = row.names(mat_ori)
mat_ori.1                = mat_ori[!(mat_ori$row_names %in% erase.n$row_names),]
mat_trait_5              = as.data.frame(cbind(genome2guild_5,mat_ori.1$id,trait_tar_3)) 

temp.6 = mat_trait_5 %>% group_by(guild) %>% summarize_if(is.numeric, mean, na.rm=TRUE)
mat.6  = (as.matrix(temp.6[2:38]))
a      = as.data.frame(colnames(mat.6))
library(ggvegan)
set.seed(2)
temp_1 = metaMDS(trait_tar_3[1:37],autotransform = T,
                 trymax = 100,noshare=0.1)

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
  geom_text(data=subset(fort,Score==""), # "species"
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(temp_1$stress,2)))

# Aim 4 ####
# Do we see changes in trait tradeoffs in Loma under drought over time

# Trait tradeoffs binary ####

b   = as.data.frame(colnames(temp.6)) # extract columns to extract the data for each cathegory
r_acqui = rowSums(temp.6 %>% select(11,27,13,10,15,29,
                                    5,25,38,26,43,16,28,9,
                                    36,7,3,35,37,4,14,19,
                                    41,8,43,34,31,12), na.rm=FALSE)/28
s_tol   = rowSums(temp.6 %>% select(17,2,40,32,22,6,42,21), na.rm=FALSE)/8
r_use   = rowSums(temp.6 %>% select(30,24,33,20,18,39), na.rm=FALSE)/6

par(mfrow=c(1,3))
plot(r_acqui,s_tol,xlab = "Resource Acquisition", ylab = "Stress Tolerance",
     col = "red",pch = 15,cex.lab = 1.5)
plot(s_tol,r_use,xlab = "Stress Tolerance", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
plot(r_acqui,r_use,xlab = "Resource Acquisition", ylab = "Resource Use",
     col = "red",pch = 15,cex.lab = 1.5)
par(mfrow=c(1,1))

# Trait tradeoffs spider ####

library(fmsb)

temp1  = as.data.frame(cbind(temp.6$guild,r_acqui,s_tol,r_use))
total  = rowSums(temp1 %>% select(2:4), na.rm=FALSE)
temp2  = as.data.frame(rbind(rep(max(total),3),
                             rep(min(r_use),3),temp1[2:4]))
colnames(temp2) = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
radarchart(temp2,plwd=3) # Spider plots using 100 guilds and total costs
