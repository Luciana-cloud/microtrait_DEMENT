# DEBmicroTrait

# Isolate Model ----

# Required input files:
#   1. isolates_microtrait.csv   — microTrait output for isolates (wide format,
#                                  columns separated by ":" in names)
#   2. isolates_metadata.csv     — isolate metadata file with genome_length,
#                                  mingentime, optimumT, rrncopies, CUE, Phylum,
#                                  Gene Count, etc.
# Output:
#   isolates2traits_new.csv      — formatted input for DEBmicroTrait

# 1. LOAD FILES

isolates_microtrait = read_csv("Input_Data/SOIL_ISOLATES/dement_isolates_CUE.csv")
df_1            = read.delim("Input_Data/SOIL_ISOLATES/metadata_2.tsv",sep="\t")
df_1            = subset(df_1, select = -c(23))
df_2            = read.delim("Input_Data/SOIL_ISOLATES/metadata_1.tsv",sep="\t")
df_2            = subset(df_2, select = -c(14,16,18))
isolates_metadata  = as.data.frame(rbind(df_1,df_2))

# Join the two files
df = isolates_microtrait[isolates_microtrait$id %in% isolates_metadata$taxon_oid, ]
colnames(isolates_metadata)[1] = "id"
df = left_join(df, isolates_metadata, by = "id")

# 2. GRAM STAIN LOOKUP FROM PHYLUM 

# Based on established microbiology — Gram stain is highly conserved at phylum level
# Adjust or expand this list if you have phyla not listed here

gram_lookup = tibble(
  Phylum = c(
    # Gram-positive
    "Firmicutes", "Bacillota", "Actinobacteria", "Actinomycetota",
    "Tenericutes", "Chloroflexi",
    # Gram-negative
    "Proteobacteria", "Pseudomonadota", "Bacteroidetes", "Bacteroidota",
    "Acidobacteria", "Acidobacteriota", "Cyanobacteria", "Cyanobacteriota",
    "Verrucomicrobia", "Verrucomicrobiota", "Planctomycetes", "Planctomycetota",
    "Chlorobi", "Chlorobiota", "Spirochaetes", "Spirochaetota",
    "Deinococcus-Thermus", "Deinococcota", "Fusobacteria", "Fusobacteriota",
    "Nitrospirae", "Nitrospirota", "Gemmatimonadetes", "Gemmatimonadota",
    "Myxococcota", "Desulfobacterota", "Campylobacterota",
    # ... existing entries ...
    # Add these:
    "Bdellovibrionota", "Thermodesulfobacteriota", "Aquificota",
    "Candidatus Saccharibacteria", "Chloroflexota", "Armatimonadota",
    "Candidatus Microgenomates", "Balneolota", "Thermomicrobiota"
  ),
  gram_stain = c(
    # Gram-positive
    "(+)", "(+)", "(+)", "(+)", "(+)", "(-)",
    # Gram-negative
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)",
    # ... existing entries ...
    # Add these (all negative):
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)"
  )
)

df = df %>% left_join(gram_lookup, by = "Phylum")

# Check how many are unmatched — these will need manual assignment
n_missing_gram = sum(is.na(df$gram_stain))
cat("Isolates with unmatched Gram stain (check Phylum names):", n_missing_gram, "\n")
if (n_missing_gram > 0) {
  cat("Unmatched phyla:\n")
  print(unique(df$Phylum[is.na(df$gram_stain)]))
}

# 3. MAP Z-TRAIT CATEGORIES 

# Helper: safely select and sum columns matching a pattern
sum_cols = function(data, pattern) {
  cols = grep(pattern, colnames(data), value = TRUE)
  if (length(cols) == 0) {
    warning(paste("No columns found matching pattern:", pattern))
    return(rep(0, nrow(data)))
  }
  rowSums(data[, cols, drop = FALSE], na.rm = TRUE)
}

# z_sugars: all carbohydrate transport subcategories
df$z_sugars = sum_cols(df, "carbohydrate transport")

# z_organic_acids: all carboxylate transport subcategories
df$z_organic_acids = sum_cols(df, "carboxylate transport")

# z_amino_acids: free amino acid transport
df$z_amino_acids = sum_cols(df, "free amino acids transport")

# z_fatty_acids: lipid transport subcategories
df$z_fatty_acids = sum_cols(df, "lipid transport")

# z_nucleotides: nucleic acid component transport subcategories
df$z_nucleotides = sum_cols(df, "nucleic acid component transport")

# z_hydrolases: complex carbohydrate depolymerization + protein degradation
df$z_hydrolases = sum_cols(df, "complex carbohydrate depolymerization") +
  sum_cols(df, "protein degradation")

# z_auxins: not present in microTrait output — set to 0
df$z_auxins = 0

# 4. BUILD FINAL OUTPUT TABLE 

# Adjust column name references below if your metadata file uses different names

isolates2traits = df %>%
  transmute(
    Isolate     = id,
    Genome_size = genome_length,           # bp
    Gene_count  = `Gene.Count....assembled`,
    GC          = NA,                      # not available — fill if you have it
    rRNA_genes  = rrncopies,
    tRNA_genes  = NA,                      # not available — fill if you have it
    Min_gen_time = mingentime,             # hours
    OGT         = optimumT,               # degrees C
    gram_stain  = gram_stain,
    z_sugars    = z_sugars,
    z_organic_acids = z_organic_acids,
    z_fatty_acids   = z_fatty_acids,
    z_auxins        = z_auxins,
    z_nucleotides   = z_nucleotides,
    z_amino_acids   = z_amino_acids,
    z_hydrolases    = z_hydrolases,
    CUE_observed    = CUE                 # keep for validation later
  )

# 6. SAVE OUTPUT 

write.csv(isolates2traits, file = "Input_data/DEB_microtrait/isolates2traits_new.csv", row.names = FALSE)

# Save DEBmicroTrait-ready version (without CUE_observed)
isolates2traits %>% select(-CUE_observed) %>% 
  write.csv("Input_data/DEB_microtrait/isolates2traits_DEBinput.csv", row.names = FALSE)

# 7. VALIDATION

sim = read_csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/DEBmicroTrait.jl-main/DEBmicroTrait.jl-main/files/my_isolates_BGE_glucose.csv")

# Join simulated BGE with observed CUE

validation = sim %>%
  inner_join(isolates2traits %>% select(Isolate, CUE_observed), 
             by = c("isolate_id" = "Isolate")) %>%
  filter(!is.na(BGE) & !is.na(CUE_observed))

# Linear model
model = lm(BGE ~ CUE_observed, data = validation)
r2    = summary(model)$r.squared
cat("R²:", r2, "\n")

# Scatter plot
ggplot(validation, aes(x = CUE_observed, y = BGE)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 0.3, y = 0.95, 
           label = paste0("R² = ", round(r2, 3))) +
  xlab("Observed CUE") +
  ylab("Simulated BGE (glucose)") +
  theme_classic()

# MAG Model ----

# 1. LOAD FILES

microtrait = read_csv("Input_Data/IMG_JGI_MAGs/total.granularity.3_datasets.csv")
metadata   = read_csv("Input_Data/IMG_JGI_MAGs/IMG_bindata_withmeta_norestricted.csv")

# 2. FILTER TO JGI MAGs ONLY

microtrait_jgi = microtrait %>%
  filter(!grepl("^RYN|^shrub|^grass", id, ignore.case = TRUE))

# 3. JOIN WITH METADATA

df = microtrait_jgi %>% left_join(metadata, by = c("id" = "Bin.ID"))

cat("MAGs after join:", nrow(df), "\n")
cat("MAGs with missing metadata:", sum(is.na(df$Phylum)), "\n")

# Remove MAGs with missing critical metadata
df = df %>%
  filter(!is.na(Phylum) & !is.na(X16s.rRNA) & !is.na(mgt) & !is.na(ogt))

cat("MAGs with complete metadata:", nrow(df), "\n")

# 4. GRAM STAIN LOOKUP FROM PHYLUM

gram_lookup = tibble(
  Phylum = c(
    # Gram-positive
    "Firmicutes", "Bacillota", "Actinobacteria", "Actinomycetota",
    "Tenericutes", "Chloroflexi",
    # Gram-negative
    "Proteobacteria", "Pseudomonadota", "Bacteroidetes", "Bacteroidota",
    "Acidobacteria", "Acidobacteriota", "Cyanobacteria", "Cyanobacteriota",
    "Verrucomicrobia", "Verrucomicrobiota", "Planctomycetes", "Planctomycetota",
    "Chlorobi", "Chlorobiota", "Spirochaetes", "Spirochaetota",
    "Deinococcus-Thermus", "Deinococcota", "Fusobacteria", "Fusobacteriota",
    "Nitrospirae", "Nitrospirota", "Gemmatimonadetes", "Gemmatimonadota",
    "Myxococcota", "Desulfobacterota", "Campylobacterota",
    # Additional phyla common in soil MAGs
    "Bdellovibrionota", "Thermodesulfobacteriota", "Aquificota",
    "Candidatus Saccharibacteria", "Chloroflexota", "Armatimonadota",
    "Candidatus Microgenomates", "Balneolota", "Thermomicrobiota",
    "Patescibacteria", "Dependentiae", "Eremiobacterota",
    "Methylomirabilota", "Modulibacteria", "Omnitrophota",
    "Calditrichota", "Caldisericota", "Chrysiogenota",
    "Deferribacterota", "Elusimicrobiota", "Fibrobacterota",
    "Halobacterota", "Hydrogenedentota", "Kiritimatiellota",
    "Latescibacterota", "Marinisomatota", "Sumerlaeota",
    "Synergistota", "Thermotogota", "WPS-2",
    "Actinobacteriota", 
    "Firmicutes_A", "Firmicutes_B", "Firmicutes_C",
    "Firmicutes_D", "Firmicutes_E", "Firmicutes_F",
    "Firmicutes_G", "Firmicutes_H"
  ),
  gram_stain = c(
    # Gram-positive
    "(+)", "(+)", "(+)", "(+)", "(+)", "(-)",
    # Gram-negative
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)",
    # Additional — all Gram-negative
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(-)", "(-)", "(-)", "(-)", "(-)", "(-)",
    "(+)", "(+)", "(+)", "(+)", 
    "(+)", "(+)", "(+)", "(+)", "(+)"
  )
)

df = df %>%
  left_join(gram_lookup, by = "Phylum")

# Remove MAGs with unmatched gram stain instead of defaulting
df = df %>% filter(!is.na(gram_stain))
cat("MAGs after removing unmatched gram stain:", nrow(df), "\n")

# Check unmatched phyla
n_missing_gram = sum(is.na(df$gram_stain))
cat("\nMAGs with unmatched Gram stain:", n_missing_gram, "\n")

# 5. MAP Z-TRAIT CATEGORIES

sum_cols = function(data, pattern) {
  cols = grep(pattern, colnames(data), value = TRUE)
  if (length(cols) == 0) {
    warning(paste("No columns found matching pattern:", pattern))
    return(rep(0, nrow(data)))
  }
  rowSums(data[, cols, drop = FALSE], na.rm = TRUE)
}

df$z_sugars       = sum_cols(df, "carbohydrate.transport")
df$z_organic_acids = sum_cols(df, "carboxylate.transport")
df$z_amino_acids  = sum_cols(df, "free.amino.acids.transport")
df$z_fatty_acids  = sum_cols(df, "lipid.transport")
df$z_nucleotides  = sum_cols(df, "nucleic.acid.component.transport")
df$z_hydrolases   = sum_cols(df, "complex.carbohydrate.depolymerization") +
  sum_cols(df, "protein.degradation")
df$z_auxins       = 0  # not captured in microTrait

# 6. BUILD FINAL OUTPUT TABLE

MAGs2traits = df %>%
  transmute(
    Isolate      = id,
    Genome_size  = genome.size,
    Gene_count   = `Gene.Count.....assembled`,
    GC           = `GC.....assembled`,
    rRNA_genes   = X16s.rRNA,
    tRNA_genes   = tRNA.Genes,
    Min_gen_time = mgt,
    OGT          = ogt,
    gram_stain   = gram_stain,
    z_sugars     = z_sugars,
    z_organic_acids = z_organic_acids,
    z_fatty_acids   = z_fatty_acids,
    z_auxins        = z_auxins,
    z_nucleotides   = z_nucleotides,
    z_amino_acids   = z_amino_acids,
    z_hydrolases    = z_hydrolases,
    # Keep quality metrics for reference
    Completeness    = Bin.Completeness,
    Contamination   = Bin.Contamination,
    Phylum          = Phylum
  )


# rRNA copy number imputed as 1 for MAGs where X16s.rRNA = 0
# due to known underrepresentation of rRNA genes in MAG assemblies
# Reference to be added

MAGs2traits = MAGs2traits %>%
  mutate(rRNA_genes = ifelse(rRNA_genes == 0, 1, rRNA_genes))

# 8. SAVE OUTPUT

# Full version with metadata for reference
write.csv(MAGs2traits, "Input_data/DEB_microtrait/MAGs2traits_full.csv", row.names = FALSE)

# DEBmicroTrait input version
MAGs2traits %>%
  select(-Completeness, -Contamination, -Phylum) %>%
  write.csv("Input_data/DEB_microtrait/MAGs2traits_DEBinput.csv", row.names = FALSE)

# 7. VALIDATION

sim_1000  = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/DEBmicroTrait.jl-main/DEBmicroTrait.jl-main/files/my_MAGs_BGE_glucose.csv")
sim_10000 = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/DEBmicroTrait.jl-main/DEBmicroTrait.jl-main/files/my_MAGs_BGE_glucose_10000.csv")

cat("Valid BGE in tspan=1000:", sum(!is.na(sim_1000$BGE)), "\n")
cat("Valid BGE in tspan=10000:", sum(!is.na(sim_10000$BGE)), "\n")

# For MAGs that were NaN in tspan=1000 but valid in tspan=10000
rescued = sim_10000 %>%
  filter(!is.na(BGE)) %>%
  filter(is.na(sim_1000$BGE[match(MAG_id, sim_1000$MAG_id)]))

cat("MAGs rescued by longer simulation:", nrow(rescued), "\n")
cat("Median BGE of rescued MAGs:", median(rescued$BGE), "\n")
cat("Median mingt of rescued MAGs:", median(rescued$mingt), "\n")

# Linear regression

# Join guild from microTrait
sim_valid = sim_10000 %>% 
  filter(!is.na(BGE) & BGE > 0) %>%
  left_join(microtrait_jgi %>% select(id, guild), 
            by = c("MAG_id" = "id"))

# Average by guild
sim_guild = sim_valid %>%
  group_by(guild) %>%
  summarise(
    BGE_mean = mean(BGE, na.rm = TRUE),
    genomesize_mean = mean(genomesize, na.rm = TRUE),
    n = n()
  ) %>%
  filter(!is.na(guild))

cat("Number of guilds with valid BGE:", nrow(sim_guild), "\n")

# Log-log linear model
model_MAG_guild = lm(log(BGE_mean) ~ log(genomesize_mean), data = sim_guild)
summary(model_MAG_guild)

# Bigger dataset

sim_50000 = read.csv("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/DEBmicroTrait.jl-main/DEBmicroTrait.jl-main/files/my_MAGs_BGE_glucose_50000.csv")

sim_valid_50000 = sim_50000 %>%
  filter(!is.na(BGE) & BGE > 0) %>%
  left_join(microtrait_jgi %>% select(id, guild),
            by = c("MAG_id" = "id"))

sim_guild_50000 = sim_valid_50000 %>%
  group_by(guild) %>%
  summarise(BGE_mean = mean(BGE, na.rm = TRUE),
            genomesize_mean = mean(genomesize, na.rm = TRUE),
            n = n()) %>%
  filter(!is.na(guild))

cat("Guilds with valid BGE (tspan=10000):", 1267, "\n")
cat("Guilds with valid BGE (tspan=50000):", nrow(sim_guild_50000), "\n")

# Regression
model_50000 = lm(log(BGE_mean) ~ log(genomesize_mean), 
                  data = sim_guild_50000)
summary(model_50000)

