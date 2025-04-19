# Supplementary Material 

# Calling packages----
# library(googledrive)
# library(googlesheets4)
library(readxl)

# Set directory----

setwd("C:/Users/lucia/OneDrive - Wageningen University & Research/UCI_projects/project_1_Microtrait_DEMENT/Manuscript/Dataset")

# Calling main data----
dat     = read_excel("Input_Data/microTrait_Rules/Table 1.xlsx", sheet = "ST2.microtrait_hmms")

# Cazymes----
GH_rule = read_excel("Input_Data/microTrait_Rules/microtrait_GH.xlsx")
GH_rule_combined = GH_rule %>% left_join(dat, by='microtrait_hmm-name')
GH_rule_combined = select(GH_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))

# Protein----
PR_rule = read_excel("Input_Data/microTrait_Rules/microtrait_proteins.xlsx")
colnames(PR_rule) = "microtrait_hmm-name"
PR_rule_combined = PR_rule %>% left_join(dat, by='microtrait_hmm-name')
PR_rule_combined = select(PR_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))

# Transporters----
transp_rule = read_excel("Input_Data/microTrait_Rules/microtrait_transporters.xlsx")
colnames(transp_rule)[1] = "microtrait_hmm-name"
transp_rule  = transp_rule %>% filter(function_gene == c("transporter"))
transp_rule_combined = transp_rule %>% left_join(dat, by='microtrait_hmm-name')
transp_rule_combined = select(transp_rule_combined, c('microtrait_hmm-name',"class",'microtrait_hmm-description'))

# Osmolyte----
osmo_rule   = read_excel("Input_Data/microTrait_Rules/microtrait_osmolytes.xlsx")
colnames(osmo_rule) = "microtrait_hmm-name"
osmo_rule_combined = osmo_rule %>% left_join(dat, by='microtrait_hmm-name')
osmo_rule_combined = select(osmo_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))

# Biofilm----
biofilm_rule = read_excel("Input_Data/microTrait_Rules/microtrait_biofilm.xlsx")
colnames(biofilm_rule) = "microtrait_hmm-name"
biofilm_rule_combined = biofilm_rule %>% left_join(dat, by='microtrait_hmm-name')
biofilm_rule_combined = select(biofilm_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))

# High Temperature----
high.T_rule  = read_excel("Input_Data/microTrait_Rules/microtrait_high_T.xlsx")
colnames(high.T_rule) = "microtrait_hmm-name"
high.T_rule_combined = high.T_rule %>% left_join(dat, by='microtrait_hmm-name')
high.T_rule_combined = select(high.T_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))

# pH----
pH_rule      = read_excel("Input_Data/microTrait_Rules/microtrait_pH_stress.xlsx")
colnames(pH_rule) = "microtrait_hmm-name"
pH_rule_combined = pH_rule %>% left_join(dat, by='microtrait_hmm-name')
pH_rule_combined = select(pH_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))


