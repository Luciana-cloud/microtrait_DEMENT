# Supplementary Material 

# Calling packages----
library(googledrive)
library(googlesheets4)
library(readxl)

# Calling main data----
dat     = read_excel("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/others/Table 1.xlsx", sheet = "ST2.microtrait_hmms")

# Cazymes----
GH_rule = read_sheet("https://docs.google.com/spreadsheets/d/1U6qWJJHossiK3kIwV2XYNKgsXkVrUc42QHiigMhX8Wg/edit?gid=1910000517#gid=1910000517")
GH_rule_combined = GH_rule %>% left_join(dat, by='microtrait_hmm-name')
GH_rule_combined = select(GH_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))
sheet_write(GH_rule_combined,
            ss = "https://docs.google.com/spreadsheets/d/1MtOeY6CLdFbHVeYMAvoRgEd2e3DpWp_og0aQaXBctI0/edit?gid=0#gid=0",
            sheet = "GH_rule_combined")

# Protein----
PR_rule = read_sheet("https://docs.google.com/spreadsheets/d/1y8kqtT9mBdf-34wDUO9LN11kN0Ugz7vk9G2ZlZSc31Y/edit?gid=1394473486#gid=1394473486")
colnames(PR_rule) = "microtrait_hmm-name"
PR_rule_combined = PR_rule %>% left_join(dat, by='microtrait_hmm-name')
PR_rule_combined = select(PR_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))
sheet_write(PR_rule_combined,
            ss = "https://docs.google.com/spreadsheets/d/1MtOeY6CLdFbHVeYMAvoRgEd2e3DpWp_og0aQaXBctI0/edit?gid=0#gid=0",
            sheet = "PR_rule_combined")

# Transporters----
transp_rule = read_sheet("https://docs.google.com/spreadsheets/d/1NJGOHKHM8IpKEZNAs69R_XlftkZ_GrhzOhQvTZfPz_w/edit?gid=1971062019#gid=1971062019")
colnames(transp_rule)[1] = "microtrait_hmm-name"
transp_rule  = transp_rule %>% filter(function_gene == c("transporter"))
transp_rule_combined = transp_rule %>% left_join(dat, by='microtrait_hmm-name')
transp_rule_combined = select(transp_rule_combined, c('microtrait_hmm-name',"class",'microtrait_hmm-description'))
sheet_write(transp_rule_combined,
            ss = "https://docs.google.com/spreadsheets/d/1MtOeY6CLdFbHVeYMAvoRgEd2e3DpWp_og0aQaXBctI0/edit?gid=0#gid=0",
            sheet = "transp_rule_combined")

# Osmolyte----
osmo_rule   = read_sheet("https://docs.google.com/spreadsheets/d/1WQ27I2Hd9jtCOv3z3cZnSB_A3xTt9bP5tZlOjiU4wVg/edit?gid=379499472#gid=379499472")
colnames(osmo_rule) = "microtrait_hmm-name"
osmo_rule_combined = osmo_rule %>% left_join(dat, by='microtrait_hmm-name')
osmo_rule_combined = select(osmo_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))
sheet_write(osmo_rule_combined,
            ss = "https://docs.google.com/spreadsheets/d/1MtOeY6CLdFbHVeYMAvoRgEd2e3DpWp_og0aQaXBctI0/edit?gid=0#gid=0",
            sheet = "osmo_rule_combined")

# Biofilm----
biofilm_rule = read_sheet("https://docs.google.com/spreadsheets/d/1-FR1s9-txuPZWg8uJ21wCn8S0hYP3JDAZpF3Qmx93QU/edit?gid=1459756284#gid=1459756284")
colnames(biofilm_rule) = "microtrait_hmm-name"
biofilm_rule_combined = biofilm_rule %>% left_join(dat, by='microtrait_hmm-name')
biofilm_rule_combined = select(biofilm_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))
sheet_write(biofilm_rule_combined,
            ss = "https://docs.google.com/spreadsheets/d/1MtOeY6CLdFbHVeYMAvoRgEd2e3DpWp_og0aQaXBctI0/edit?gid=0#gid=0",
            sheet = "biofilm_rule_combined")

# High Temperature----
high.T_rule  = read_sheet("https://docs.google.com/spreadsheets/d/1-PeZ-F2RnXVFMa0Hkmg8Q-zlp38bCf4lxC-gpI9jEu4/edit?gid=123098806#gid=123098806")
colnames(high.T_rule) = "microtrait_hmm-name"
high.T_rule_combined = high.T_rule %>% left_join(dat, by='microtrait_hmm-name')
high.T_rule_combined = select(high.T_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))
sheet_write(high.T_rule_combined,
            ss = "https://docs.google.com/spreadsheets/d/1MtOeY6CLdFbHVeYMAvoRgEd2e3DpWp_og0aQaXBctI0/edit?gid=0#gid=0",
            sheet = "high.T_rule_combined")

# pH----
pH_rule      = read_sheet("https://docs.google.com/spreadsheets/d/1kANfYGvbb8tDiYEkJ_Whb9kocfhv8oDWpbNITMdZFQg/edit?gid=1661588961#gid=1661588961")
colnames(pH_rule) = "microtrait_hmm-name"
pH_rule_combined = pH_rule %>% left_join(dat, by='microtrait_hmm-name')
pH_rule_combined = select(pH_rule_combined, c('microtrait_hmm-name','microtrait_hmm-description'))
sheet_write(pH_rule_combined,
            ss = "https://docs.google.com/spreadsheets/d/1MtOeY6CLdFbHVeYMAvoRgEd2e3DpWp_og0aQaXBctI0/edit?gid=0#gid=0",
            sheet = "pH_rule_combined")

