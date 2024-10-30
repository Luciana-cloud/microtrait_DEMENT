
total_data = as.data.frame(rbind(total_genes.guild.940_ISO,total_genes.guild.940_MAG))
# select columns

total_data     = total_data %>% select(c("id","genome.size","guild"))
hmm_trait      = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/IMG_MAG_results/trait_granularity_3_MAG_IMG.csv",dec=".")
hmm_trait      = hmm_trait %>% select(-c(X))
names(hmm_trait)[names(hmm_trait) == 'growthrate_d'] = 'mgt'     
burnt_trait    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAGs_burnt/trait_matrixatgranularity3_Fire.csv",dec=".")
burnt_trait    = burnt_trait %>% select(-c(X))
names(burnt_trait)[names(burnt_trait) == 'growthrate_d'] = 'mgt' 
loma_trait     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/MAG_Loma/litter_mags_trait_matrixatgranularity.csv",dec=".")
loma_trait     = loma_trait %>% add_column(ogt = NA)
loma_trait     = loma_trait %>% relocate(mgt, .after = ogt)
isolate_trait  = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/isolates_rds/trait_granularity_3_isolates.csv",dec=".")
isolate_trait  = isolate_trait %>% add_column(ogt = NA)
isolate_trait  = isolate_trait %>% add_column(mgt = NA)
isolate_trait  = isolate_trait %>% select(-c(X))

total_trait    = as.data.frame(rbind(hmm_trait,loma_trait,burnt_trait,isolate_trait))

# Merge datasets----

total.granularity.3 = total_data %>% left_join(total_trait, by='id')
write.csv(total.granularity.3, file = "C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/total.granularity.3_datasets.csv")

sheet_write(total.granularity.3,
            ss = "https://docs.google.com/spreadsheets/d/1t_DpBn0Ra256ZgaZUtYlgRGT9JBAdR1lGi8nbIgRO4o/edit?gid=0#gid=0",
            sheet = "total.granularity.3")

# CUE

isolates     = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Data/dement_isolates_CUE.csv",dec=".")
isolatest    = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/total.granularity.3_datasets.csv",dec=".")
isolates.1   = subset(isolates, isolates$CUE !='NaN')
total_genes.guild.940_ISO = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MAG_database/Intermediate_results/total_genes.guild.940_ISO.csv",dec=".")

# Merge

isolates.2   = isolates.1 %>% left_join(total_genes.guild.940_ISO, by='id')
isolates.2.m = isolates.2 %>% group_by(guild) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
Figure_test_1.CUE = ggplot(data = isolates.2.m, aes(x = genome.size_mean, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CUE") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,1) + labs(title = "Functional Groups")
Figure_test_1.CUE

# Metadata

meta_iso   = readr::read_tsv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils_unpublished.tsv")
meta_iso   = as.data.frame(cbind(meta_iso$`IMG Genome ID`,meta_iso$Phylum,
                                 meta_iso$Class,meta_iso$Order,meta_iso$Family,
                                 meta_iso$Genus))
meta_iso.1 = readr::read_tsv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/isolates_GOLD/soils.tsv")
meta_iso.1 = as.data.frame(cbind(meta_iso.1$`IMG Genome ID`,meta_iso.1$Phylum,
                                 meta_iso.1$Class,meta_iso.1$Order,meta_iso.1$Family,
                                 meta_iso.1$Genus))

# final meta
final_meta   = as.data.frame(rbind(meta_iso,meta_iso.1))
colnames(final_meta) = c("id","Phylum","Class","Order","Family","Genus")

merged_df    = merge(isolates.1,final_meta,by = "id")

# Phylum

phylum = merged_df %>% group_by(Phylum) %>% summarise(across(where(is.numeric), 
                                                                   list(mean=mean), na.rm=TRUE))
Figure_test_1.phylum = ggplot(data = phylum, aes(x = genome_length_mean, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CUE") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,1) + labs(title = "Phylum")
Figure_test_1.phylum

# Class

Class = merged_df %>% group_by(Class) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
Figure_test_1.Class = ggplot(data = Class, aes(x = genome_length_mean, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CUE") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,1) + labs(title = "Class")
Figure_test_1.Class

# Order

Order = merged_df %>% group_by(Order) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
Figure_test_1.Order = ggplot(data = Order, aes(x = genome_length_mean, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CUE") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,1) + labs(title = "Order")
Figure_test_1.Order

# Family

Family = merged_df %>% group_by(Family) %>% summarise(across(where(is.numeric), 
                                                           list(mean=mean), na.rm=TRUE))
Figure_test_1.Family = ggplot(data = Family, aes(x = genome_length_mean, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CUE") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,1) + labs(title = "Family")
Figure_test_1.Family

# Genus

Genus = merged_df %>% group_by(Genus) %>% summarise(across(where(is.numeric), 
                                                             list(mean=mean), na.rm=TRUE))
Figure_test_1.Genus = ggplot(data = Genus, aes(x = genome_length_mean, y = CUE_mean)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) + xlab("Genome Size") + 
  ylab("CUE") +
  geom_point() + theme_classic() + theme(text = element_text(size=14)) + 
  xlim(0,1.25e7) + ylim(0,1) + labs(title = "Genus")
Figure_test_1.Genus
