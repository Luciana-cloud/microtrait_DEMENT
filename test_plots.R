# Test Figures ###

total.guilds    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/total.genes.guilds.selected.csv",dec=".")
a               = as.data.frame(colnames(total.guilds))
loma_stat       = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/mag_stats.txt") 
loma_stat       = as.data.frame(loma_stat[c(1,7)])
fire_stat       = read.delim("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/fire_metadata.txt")
fire_stat       = fire_stat[-c(440,546), ]
colnames(fire_stat) = c("id","size")
total_stat      = as.data.frame(rbind(fire_stat,loma_stat))

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

full_mat    = left_join(total.guilds, total_stat, by=c('id'))

ggplot(full_mat, aes(id,size)) + geom_point() +
  xlab("MAG id") + ylab("Genome size")

ggplot(full_mat, aes(guild,size)) + geom_point() +
  xlab("Functional group") + ylab("Genome size")

full_mat.2    = left_join(total_stat, GHs_guild.2, by=c('id'))

ggplot(full_mat.2, aes(size,sum.2)) + geom_point() +
  xlab("Genome size") + ylab("Total GH costs")


