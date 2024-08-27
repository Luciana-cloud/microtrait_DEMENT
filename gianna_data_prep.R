
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


