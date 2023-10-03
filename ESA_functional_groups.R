# Final Analysis

setwd("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT")

# General goal: fitness traits tradeoffs and microbial life history strategies 
# using genome-scale data

# Packages ####

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

# Calling data and preprocessing ####

mat_ori    = read.csv("litter_mags_trait_matrixatgranularity.csv",dec=".")
gen_size   = read.delim("litter_mags_metadata.txt",dec=".") 

mag_stat   = read.delim("mag_stats.txt") 
mag_abun   = read.delim("mag_adundance.txt") 

mag_stat   = mag_stat %>% full_join(mag_abun)

# Erase weird MAGs, which were found in preliminary analysis: 116, 175

mat_ori    = mat_ori[-c(116,175),] 
gen_size   = gen_size[-c(116,175),] 
mag_stat   = mag_stat[-c(116,175),]

mat_trait  = as.data.frame(cbind((mat_ori %>% select(2:49)/gen_size$length),
                                  mat_ori %>% select(50:161),
                                  mat_ori %>% select(162:190)/gen_size$length,
                                  mat_ori %>% select(191)))
mat_trait    = as.data.frame(cbind(mat_ori$id,mat_trait))
mat_trait    = mat_trait %>% mutate(id_2 = seq(from = 1, to = 531, by = 1))

# Aim 1: Trait selection ####

# Grassland ####

mat_trait_g   = mat_trait %>% mutate(grass_abund = mag_stat$Average)
mat_trait_g   = as.data.frame(mat_trait_g[,-192])
mat_g         = (as.matrix(mat_trait_g[2:192]))
# corre.g       = as.data.frame(cor(mat_g[,1:190], mat_g[,191]))
# Erase traits based on correlation matrix

test.1        = as.data.frame(cbind(mat_trait_g[2:191]))
seq.1         = seq(1,190)
temp.1        = as.data.frame(cbind((colSums(test.1)),seq.1))
erase.1       = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
test.1        = test.1[-(erase.1$seq.1)]
temp.1        = as.data.frame(colSums(test.1))
# corre.1.1     = as.data.frame(cor(as.matrix(test.1), mat_g[,191]))
# Best predictors statistics

m.all        = lm(mat_g[,191]~.,data=test.1[1:50])
temp.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[51:95])
temp.51.95   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[96:99])
temp.96.99   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[100:102])
temp.100.102 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[103:109])
temp.103.109 = ols_step_forward_p(m.all, details = FALSE)

# Erasing 110 and 111 because the trait values are the same as trait 109

m.all        = lm(mat_g[,191]~.,data=test.1[112:115])
temp.112.115 = ols_step_forward_p(m.all, details = FALSE)

# Erasing 116 because the trait value is the same as trait 115

m.all        = lm(mat_g[,191]~.,data=test.1[117:125])
temp.117.125 = ols_step_forward_p(m.all, details = FALSE)

# Erasing 126 because the trait value is the same as trait 113

m.all        = lm(mat_g[,191]~.,data=test.1[127:134])
temp.127.134 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[135:140])
temp.135.140 = ols_step_forward_p(m.all, details = FALSE)

erase.1a     = as.data.frame(c(temp.1.50[["predictors"]],
                               temp.51.95[["predictors"]],
                               temp.96.99[["predictors"]],
                               temp.100.102[["predictors"]],
                               temp.103.109[["predictors"]],
                               temp.112.115[["predictors"]],
                               temp.117.125[["predictors"]],
                               temp.127.134[["predictors"]],
                               temp.135.140[["predictors"]]))
colnames(erase.1a) = c("trait")
test.g             = test.1 %>% select((erase.1a$trait))
# corre.1.2          = as.data.frame(cor(as.matrix(test.g),mat_g[,191]))
a                  = as.data.frame(colnames(test.g))
test.g             = test.g[c(-48,-57,-60)]
m.all              = lm(mat_g[,191]~.,data=test.g)
best.predictors    = ols_step_forward_p(m.all, details = FALSE)
erase.1a           = as.data.frame(best.predictors[["predictors"]])
colnames(erase.1a) = c("trait")
test.g             = test.g %>% select((erase.1a$trait))
m.all              = lm(mat_g[,191]~.,data=test.g)
write.csv(erase.1a, file = "grassland.best.predictors.csv")

# Shrubland ####

mat_trait_s   = mat_trait %>% mutate(shrub_abund = mag_stat$Average.2)
mat_trait_s   = as.data.frame(mat_trait_s[,-192])
mat_s         = (as.matrix(mat_trait_s[2:192]))
# corre.g       = as.data.frame(cor(mat_g[,1:190], mat_g[,191]))
# Erase traits based on correlation matrix

test.1        = as.data.frame(cbind(mat_trait_s[2:191]))
seq.1         = seq(1,190)
temp.1        = as.data.frame(cbind((colSums(test.1)),seq.1))
erase.1       = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
test.1        = test.1[-(erase.1$seq.1)]
temp.1        = as.data.frame(colSums(test.1))
# corre.1.1     = as.data.frame(cor(as.matrix(test.1), mat_s[,191]))
# Best predictors statistics

m.all        = lm(mat_s[,191]~.,data=test.1[1:50])
temp.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_s[,191]~.,data=test.1[51:95])
temp.51.95   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_s[,191]~.,data=test.1[96:99])
temp.96.99  = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_s[,191]~.,data=test.1[100:113])
temp.100.113 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_s[,191]~.,data=test.1[114:115])
temp.114.115 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_s[,191]~.,data=test.1[116:140])
temp.116.140 = ols_step_forward_p(m.all, details = FALSE)

erase.1a     = as.data.frame(c(temp.1.50[["predictors"]],
                               temp.51.95[["predictors"]],
                               temp.96.99[["predictors"]],
                               temp.100.113[["predictors"]],
                               temp.114.115[["predictors"]],
                               temp.116.140[["predictors"]]))
colnames(erase.1a) = c("trait")
test.s             = test.1 %>% select((erase.1a$trait))
# corre.1.2          = as.data.frame(cor(as.matrix(test.g),mat_g[,191]))
a                  = as.data.frame(colnames(test.s))
test.s             = test.s[c(-60,-61)]
m.all              = lm(mat_s[,191]~.,data=test.s)
best.predictors    = ols_step_forward_p(m.all, details = FALSE)
erase.1a           = as.data.frame(best.predictors[["predictors"]])
colnames(erase.1a) = c("trait")
test.s             = test.s %>% select((erase.1a$trait))
m.all              = lm(mat_s[,191]~.,data=test.s)
write.csv(erase.1a, file = "shrubland.best.predictors.csv")

# Grassland Drought ####

mat_trait_gd   = mat_trait %>% mutate(grass_abund = mag_stat$Average.1)
mat_trait_gd   = as.data.frame(mat_trait_gd[,-192])
mat_gd         = (as.matrix(mat_trait_gd[2:192]))
# corre.gd       = as.data.frame(cor(mat_gd[,1:190], mat_gd[,191]))
# Erase traits based on correlation matrix

test.1        = as.data.frame(cbind(mat_trait_gd[2:191]))
seq.1         = seq(1,190)
temp.1        = as.data.frame(cbind((colSums(test.1)),seq.1))
erase.1       = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
test.1        = test.1[-(erase.1$seq.1)]
temp.1        = as.data.frame(colSums(test.1))
# corre.1.1     = as.data.frame(cor(as.matrix(test.1), mat_gd[,191]))
# Best predictors statistics

m.all        = lm(mat_gd[,191]~.,data=test.1[1:50])
temp.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_gd[,191]~.,data=test.1[51:93])
temp.51.93   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_gd[,191]~.,data=test.1[94:99])
temp.94.99   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_gd[,191]~.,data=test.1[100:115])
temp.100.115 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_gd[,191]~.,data=test.1[116:140])
temp.116.140 = ols_step_forward_p(m.all, details = FALSE)
erase.1a     = as.data.frame(c(temp.1.50[["predictors"]],
                               temp.51.93[["predictors"]],
                               temp.96.99[["predictors"]],
                               temp.94.99[["predictors"]],
                               temp.100.115[["predictors"]],
                               temp.116.140[["predictors"]]))
colnames(erase.1a) = c("trait")
test.gd            = test.1 %>% select((erase.1a$trait))
# corre.1.2          = as.data.frame(cor(as.matrix(test.g),mat_gd[,191]))
a                  = as.data.frame(colnames(test.gd))
m.all              = lm(mat_gd[,191]~.,data=test.gd[1:46])
temp.1.46          = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_gd[,191]~.,data=test.gd[47:65])
temp.47.65         = ols_step_forward_p(m.all, details = FALSE)
erase.1a           = as.data.frame(c(temp.1.46[["predictors"]],
                                     temp.47.65[["predictors"]]))
colnames(erase.1a) = c("trait")
test.gd            = test.gd %>% select((erase.1a$trait))
a                  = as.data.frame(colnames(test.gd))
m.all              = lm(mat_gd[,191]~.,data=test.gd)
best.predictors    = ols_step_forward_p(m.all, details = FALSE)
erase.1a           = as.data.frame(best.predictors[["predictors"]])
colnames(erase.1a) = c("trait")
test.gd            = test.gd %>% select((erase.1a$trait))
m.all              = lm(mat_gd[,191]~.,data=test.gd)
write.csv(erase.1a, file = "grassland.drought.best.predictors.csv")

# Shrubland Drought ####

mat_trait_g   = mat_trait %>% mutate(grass_abund = mag_stat$Average.3)
mat_trait_g   = as.data.frame(mat_trait_g[,-192])
mat_g         = (as.matrix(mat_trait_g[2:192]))
# corre.g       = as.data.frame(cor(mat_g[,1:190], mat_g[,191]))
# Erase traits based on correlation matrix

test.1        = as.data.frame(cbind(mat_trait_g[2:191]))
seq.1         = seq(1,190)
temp.1        = as.data.frame(cbind((colSums(test.1)),seq.1))
erase.1       = temp.1 %>% filter(V1==0)
erase.1$row_names = row.names(erase.1)
test.1        = test.1[-(erase.1$seq.1)]
temp.1        = as.data.frame(colSums(test.1))
# corre.1.1     = as.data.frame(cor(as.matrix(test.1), mat_g[,191]))
# Best predictors statistics

m.all        = lm(mat_g[,191]~.,data=test.1[1:50])
temp.1.50    = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[51:93])
temp.51.93   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[94:99])
temp.94.99   = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[100:102])
temp.100.102 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[103:115])
temp.103.115 = ols_step_forward_p(m.all, details = FALSE)
m.all        = lm(mat_g[,191]~.,data=test.1[116:140])
temp.116.140 = ols_step_forward_p(m.all, details = FALSE)
erase.1a     = as.data.frame(c(temp.1.50[["predictors"]],
                               temp.51.93[["predictors"]],
                               temp.94.99[["predictors"]],
                               temp.100.102[["predictors"]],
                               temp.103.115[["predictors"]],
                               temp.116.140[["predictors"]]))
colnames(erase.1a) = c("trait")
test.g             = test.1 %>% select((erase.1a$trait))
# corre.1.2          = as.data.frame(cor(as.matrix(test.g),mat_g[,191]))
a                  = as.data.frame(colnames(test.g))

m.all              = lm(mat_g[,191]~.,data=test.g[1:45])
temp.1.45          = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,191]~.,data=test.g[46:50])
temp.46.50         = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,191]~.,data=test.g[51:60])
temp.51.60         = ols_step_forward_p(m.all, details = FALSE)
m.all              = lm(mat_g[,191]~.,data=test.g[61:63])
temp.61.63         = ols_step_forward_p(m.all, details = FALSE)
erase.1a     = as.data.frame(c(temp.1.45[["predictors"]],
                               temp.46.50[["predictors"]],
                               temp.51.60[["predictors"]],
                               temp.61.63[["predictors"]]))
colnames(erase.1a) = c("trait")
test.g             = test.g %>% select((erase.1a$trait))
a                  = as.data.frame(colnames(test.g))
test.g1            = test.g[c(-38,-45,-50)]
a                  = as.data.frame(colnames(test.g1))
test.g1            = test.g1[c(-22)] # erasing trait 22
m.all              = lm(mat_g[,191]~.,data=test.g1)
best.predictors    = ols_step_forward_p(m.all, details = FALSE)
erase.1a           = as.data.frame(best.predictors[["predictors"]])
colnames(erase.1a) = c("trait")
test.g1            = test.g1 %>% select((erase.1a$trait))
m.all              = lm(mat_g[,191]~.,data=test.g1)
write.csv(erase.1a, file = "shrubland.drought.best.predictors.csv")

# Selected traits ####

shrubland.drought = read.csv(file = "shrubland.drought.best.predictors.csv")
shrubland.ambient = read.csv(file = "shrubland.best.predictors.csv")
grassland.drought = read.csv(file = "grassland.drought.best.predictors.csv")
grassland.ambient = read.csv(file = "grassland.best.predictors.csv")
trait.selec       = as.data.frame(unique(c(shrubland.drought$trait,shrubland.ambient$trait,
                          grassland.drought$trait,grassland.ambient$trait)))
colnames(trait.selec) = c("trait")
mat_trait.fin         = as.data.frame(cbind(mat_trait$`mat_ori$id`,
                                            mat_trait %>% select((trait.selec$trait))))
# write.csv(mat_trait.fin, file = "selected.trait.values.csv")

# Aim 2 : Functional groups ####

# Data scaling ####
mat_trait.fin = read.csv(file = "selected.trait.values.csv")
trait_tar     = scale(mat_trait.fin[,3:90],center = FALSE) # scaling with center=False to avoid breaking the sparsity structure of the data

# Hierarchical Clustering ####
set.seed(1)
distanceg.1  = vegdist(trait_tar, method = "euclidean", binary = FALSE)
clusterg.1   = hclust(distanceg.1, method="ward")
nguilds.1    = seq(2, nrow(trait_tar), 2)
plot(clusterg.1)

library(factoextra)
fviz_dend(clusterg.1, cex = 0.5)

fviz_dend(clusterg.1, k = 15,                 # Cut in four groups
          cex = 0.25,                 # label size
          k_colors = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                       "#3F4921", "#7FDCC0", "#CBD588", "#5F7FC7",
                       "#673770", "#D3D93E", "#38333E", "#508578", 
                       "#D7C1B1", "#689030", "#AD6F3B"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
)

# Similarity within guilds ####
my_vec = c()

for(i in nguilds.1) {
  v                      = cutree(clusterg.1,k=i)
  genome2guild           = data.frame(guild = factor(v))
  rownames(genome2guild) = names(v)
  adonis_2               = vegan::adonis2(distanceg.1 ~ guild, data = genome2guild, perm = 1)
  temp                   = adonis_2$R2
  temp_1                 = temp[1]
  my_out                 = temp_1
  my_vec   <- c(my_vec, my_out)   
}
mat_r2_1  = as.data.frame(cbind(nguilds.1,my_vec))

# Similarity among guilds ####
v                      = cutree(clusterg.1,k=15)
genome2guild           = data.frame(guild = factor(v))
rownames(genome2guild) = names(v)
adonis_1               = pairwiseAdonis::pairwise.adonis(distanceg.1,genome2guild$guild,perm = 999)
test.clus.1            = adonis.pair(distanceg.1, genome2guild[,"guild"], nper = 1000, 
                                     corr.method = "fdr")
adonis_2               = vegan::adonis2(distanceg.1 ~ guild, data = genome2guild, perm = 1)

# Plotting similarity ####
pdf("Test.pdf", height = 3, width = 5)
ggplot(data=mat_r2_1,aes(x=nguilds.1,y=my_vec)) + geom_line() +
  xlab("Number of guilds") + ylab("Similarity within guilds") +
  theme(text = element_text(size = 16)) + 
  geom_hline(yintercept=0.4100782, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 15, linetype="dashed", color = "red") +
  theme_classic()
dev.off()

# Ordination of the 15 functional groups ####
mat_trait_1a           = as.data.frame(cbind(genome2guild,mat_ori$id,trait_tar))
write.csv(mat_trait_1a, file = "mat_trait_1a.csv")
mat_trait_1a           = read.csv(file = "mat_trait_1a.csv")

set.seed(16)
ordination_1           = metaMDS(trait_tar,autotransform = T,trymax = 500)
fort.1                 = fortify(ordination_1)

p1 = ggplot() + geom_point(data=subset(fort.1,Score=="sites"),
                      mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat_trait_1a$guild),size = 2),
                      alpha=0.5) + 
  geom_segment(data=subset(fort.1,Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="darkgray",
               linewidth=0.8) + 
  geom_text(data=subset(fort.1,Score==""), # "species"
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1)) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2)))

ggplot() + geom_point(data=subset(fort.1,Score=="sites"),
                           mapping = aes(x=NMDS1,y=NMDS2,color =as.factor(mat_trait_1a$guild),size = 2),
                           alpha=0.5) + 
  geom_abline(intercept=0,slope=0,linetype="dashed",linewidth=0.8,colour="gray") + 
  geom_vline(aes(xintercept=0),linetype="dashed",linewidth=0.8,colour="gray") + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) + 
  annotate("text", x=-1, y=-1, label=paste('Stress =',round(ordination_1$stress,2))) + 
  stat_ellipse(data = subset(fort.1,Score=="sites"), 
            aes(x = NMDS1, y = NMDS2, color = as.factor(mat_trait_1a$guild))) + 
  scale_colour_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                               "#3F4921", "#7FDCC0", "#CBD588", "#5F7FC7",
                               "#673770", "#D3D93E", "#38333E", "#508578", 
                               "#D7C1B1", "#689030", "#AD6F3B"))

# Aim 3 : Life history strategies ####

final_trait = as.data.frame(cbind(genome2guild$guild,mat_trait$`mat_ori$id`,
                                  mat_trait.fin[,3:90]))
a.1         = as.data.frame(colnames(final_trait))
write.csv(final_trait, file = "final_trait.csv")
# Totals ####
r_acqui.t   = c(rowSums(final_trait %>% select(27,45,3,50,72,73,39,53,78,15,81,
                                               17,67,14,77,42,44,32,13,25,88,51,
                                               89,85,18,29,79,16,35,47,34,59,57,
                                               11,90,40,24,65,43,58,7,75,20,12,
                                               10,8,41,38,49,74,28,68,70,19,71),
                        na.rm=FALSE)/55)
s_tol.t     = c(rowSums(final_trait %>% select(4,5,64,6,36,80,37,46,33,82,66,69,
                                               54,31,22,52),na.rm=FALSE)/16)
r_use.t     = c(rowSums(final_trait %>% select(83,62,63,86,87,23,84,9,76,21,55,
                                               56,60,26,30,61,48), na.rm=FALSE)/17)
temp1.t     = as.data.frame(cbind(r_acqui.t,s_tol.t,r_use.t))

# FGn (1:15) ####
group_n   = final_trait %>% filter(`genome2guild$guild` == 15)
r_acqui   = c(rowSums(group_n %>% select(27,45,3,50,72,73,39,53,78,15,81,
                                               17,67,14,77,42,44,32,13,25,88,51,
                                               89,85,18,29,79,16,35,47,34,59,57,
                                               11,90,40,24,65,43,58,7,75,20,12,
                                               10,8,41,38,49,74,28,68,70,19,71),
                        na.rm=FALSE)/55)
s_tol     = c(rowSums(group_n %>% select(4,5,64,6,36,80,37,46,33,82,66,69,
                                               54,31,22,52),na.rm=FALSE)/16)
r_use     = c(rowSums(group_n %>% select(83,62,63,86,87,23,84,9,76,21,55,
                                               56,60,26,30,61,48), na.rm=FALSE)/17)
temp1     = as.data.frame(cbind(r_acqui,s_tol,r_use))

temp1.1   = as.data.frame(cbind(mean(r_acqui),mean(s_tol),mean(r_use)))
temp2     = as.data.frame(rbind(rep(max(temp1.t),3),
                                    rep(min(temp1.t),3),temp1))
temp2.1   = as.data.frame(rbind(rep(max(temp1.t),3),
                                    rep(min(temp1.t),3),temp1.1))
colnames(temp2)   = c("Resource Acquisition","Stress Tolerance",
                    "Resource Use")
colnames(temp2.1) = c("Resource Acquisition","Stress Tolerance",
                      "Resource Use")
radarchart(temp2,plwd=1) 
radarchart(temp2.1,plwd=2) 

# Aim 3 : Abundance of functional groups under conditions ####

fg_abundance = as.data.frame(cbind(genome2guild$guild,mat_trait$`mat_ori$id`,
                                  mag_stat[,7:27]))
# write.csv(fg_abundance, file = "fg_abundance.csv")
fg_abundance = read.csv(file = "fg_abundance.csv")

grassland.a    = fg_abundance %>% group_by(`genome2guild$guild`) %>% 
  summarise(abundance = sum(Average)/sum(fg_abundance$Average)) %>% 
  mutate(condition = rep("grassland.ambient" , nrow(grassland.a)))
grassland.d    = fg_abundance %>% group_by(`genome2guild$guild`) %>% 
  summarise(abundance = sum(Average.1)/sum(fg_abundance$Average.1)) %>% 
  mutate(condition = rep("grassland.drought" , nrow(grassland.a)))
shrubland.a    = fg_abundance %>% group_by(`genome2guild$guild`) %>% 
  summarise(abundance = sum(Average.2)/sum(fg_abundance$Average.2)) %>% 
  mutate(condition = rep("shrubland.ambient" , nrow(grassland.a)))
shrubland.d    = fg_abundance %>% group_by(`genome2guild$guild`) %>% 
  summarise(abundance = sum(Average.3)/sum(fg_abundance$Average.3)) %>% 
  mutate(condition = rep("shrubland.drought" , nrow(grassland.a)))

fg_ab.fig      =  as.data.frame(rbind(grassland.a,grassland.d,shrubland.a,
                                      shrubland.d))
colnames(fg_ab.fig) = c("guild","abundance","condition")
# write.csv(fg_ab.fig, file = "fg_ab.fig.csv")
fg_ab.fig = read.csv(file = "fg_ab.fig.csv")

ggplot(fg_ab.fig, aes(fill=as.factor(guild), y=abundance, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                             "#3F4921", "#7FDCC0", "#CBD588", "#5F7FC7",
                             "#673770", "#D3D93E", "#38333E", "#508578", 
                             "#D7C1B1", "#689030", "#AD6F3B")) + 
  theme(text = element_text(size=25)) + 
  labs(y="Abundance",x = element_blank()) + 
  theme(legend.title = element_blank())
