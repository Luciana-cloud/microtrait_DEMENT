# ============================================================
# Find the anti-mode (local minimum) between the two BGE populations,
# then compare genome.size~CUE regression with vs without exclusion.
# ============================================================

# MAG ----

library(dplyr)

# --- Find the local minimum in the density, as an objective threshold ---
d = density(sim_10000$BGE[sim_10000$BGE > 0.5 & sim_10000$BGE < 1], na.rm = TRUE, n = 512)
# look for the local min between the two peaks (roughly x in 0.78-0.90 based on the histogram)
window = d$x > 0.75 & d$x < 0.92
local_min_x = d$x[window][which.min(d$y[window])]
cat("Data-driven anti-mode threshold (local minimum in density):", round(local_min_x, 3), "\n\n")

cat("How many rows would be excluded at this threshold?\n")
cat("  BGE >", round(local_min_x,3), ":", sum(sim_10000$BGE > local_min_x, na.rm=TRUE), "\n\n")

# --- Genome size check at this threshold ---
sim_10000 %>%
  mutate(near_one = BGE > local_min_x) %>%
  filter(!is.na(near_one)) %>%
  group_by(near_one) %>%
  summarise(mean_genome_size = mean(genomesize, na.rm=TRUE), n = n())

# --- Sensitivity comparison: guild-level regression WITH vs WITHOUT exclusion ---
guild_lookup = total_genes.guild.940_MAG %>% select(id, guild)

cue_with = sim_10000 %>%
  left_join(guild_lookup, by = c("MAG_id" = "id")) %>%
  filter(!is.na(guild))

# WITH the near-1 values (current approach)
guild_cue_with = cue_with %>% group_by(guild) %>%
  summarise(CUE_mean = mean(BGE, na.rm=TRUE), genome.size = mean(genomesize, na.rm=TRUE))
lm_with = lm(CUE_mean ~ log(genome.size), data = guild_cue_with)

# WITHOUT the near-1 values (excluded at the data-driven threshold)
guild_cue_without = cue_with %>% filter(BGE <= local_min_x) %>% group_by(guild) %>%
  summarise(CUE_mean = mean(BGE, na.rm=TRUE), genome.size = mean(genomesize, na.rm=TRUE))
lm_without = lm(CUE_mean ~ log(genome.size), data = guild_cue_without)

cat("=== WITH near-1 values, n guilds =", nrow(guild_cue_with), "===\n")
print(summary(lm_with)$coefficients)
cat("Adj R2:", summary(lm_with)$adj.r.squared, "\n\n")

cat("=== WITHOUT near-1 values (BGE >", round(local_min_x,3), "excluded), n guilds =", nrow(guild_cue_without), "===\n")
print(summary(lm_without)$coefficients)
cat("Adj R2:", summary(lm_without)$adj.r.squared, "\n")

# Artifact diagnosis -----

sim_10000 = sim_10000 %>% mutate(near_one = BGE > 0.848)

cat("Comparing other simulation outputs: near-1 vs rest\n\n")

sim_10000 %>%
  filter(!is.na(near_one)) %>%
  group_by(near_one) %>%
  summarise(
    n = n(),
    mean_rgrowth = mean(rgrowth, na.rm=TRUE),
    mean_BP      = mean(BP, na.rm=TRUE),
    mean_BR      = mean(BR, na.rm=TRUE),
    mean_mingt   = mean(mingt, na.rm=TRUE),
    mean_rrn     = mean(rrn, na.rm=TRUE),
    pct_BR_negative = mean(BR < 0, na.rm=TRUE) * 100,
    pct_BR_near_zero = mean(abs(BR) < 1e-6, na.rm=TRUE) * 100
  )

# Distribution of BR (respiration) specifically for near-1 rows - if respiration
# is suspiciously near zero or negative, that's a strong signal of non-convergence
cat("\nBR (respiration) summary for near-1 rows only:\n")
print(summary(sim_10000$BR[sim_10000$near_one]))

cat("\nBR (respiration) summary for rest:\n")
print(summary(sim_10000$BR[!sim_10000$near_one]))

# ISOLATES ----

my_isolates_BGE_glucose = read_csv("Output_Data/my_isolates_BGE_glucose.csv")  # adjust path if different

hist(my_isolates_BGE_glucose$BGE, breaks = 50)

# check respiration near zero, same signature as the MAG artifact
summary(my_isolates_BGE_glucose$BR[my_isolates_BGE_glucose$BGE > 0.9])
summary(my_isolates_BGE_glucose$BR[my_isolates_BGE_glucose$BGE <= 0.9])

sum(my_isolates_BGE_glucose$BR < 0, na.rm = TRUE)

# Look for the local minimum just before the big spike (visually around 0.75-0.82)
d = density(my_isolates_BGE_glucose$BGE[my_isolates_BGE_glucose$BGE > 0.5 & 
                                          my_isolates_BGE_glucose$BGE < 0.9], 
            na.rm = TRUE, n = 512)

window = d$x > 0.70 & d$x < 0.82
local_min_x_iso = d$x[window][which.min(d$y[window])]
cat("Isolate-specific anti-mode threshold:", round(local_min_x_iso, 3), "\n\n")

cat("How many isolates excluded at this threshold?\n")
cat("  BGE >", round(local_min_x_iso,3), ":", 
    sum(my_isolates_BGE_glucose$BGE > local_min_x_iso, na.rm=TRUE), 
    "out of", nrow(my_isolates_BGE_glucose), "\n\n")

# Genome size skew check, same as we did for MAGs
my_isolates_BGE_glucose %>%
  mutate(near_ceiling = BGE > local_min_x_iso) %>%
  filter(!is.na(near_ceiling)) %>%
  group_by(near_ceiling) %>%
  summarise(mean_genome_size = mean(genomesize, na.rm=TRUE), n = n())

# Check BR across finer BGE bands - where does it actually collapse?
my_isolates_BGE_glucose %>%
  mutate(BGE_band = case_when(
    BGE <= 0.70 ~ "<=0.70",
    BGE <= 0.75 ~ "0.70-0.75",
    BGE <= 0.80 ~ "0.75-0.80",
    BGE <= 0.85 ~ "0.80-0.85",
    BGE <= 0.90 ~ "0.85-0.90",
    TRUE ~ ">0.90"
  )) %>%
  group_by(BGE_band) %>%
  summarise(n = n(),
            mean_BR = mean(BR, na.rm=TRUE),
            median_BR = median(BR, na.rm=TRUE),
            pct_BR_zero = mean(BR == 0, na.rm=TRUE)*100)

# Also check: how many isolates are in the big spike specifically (0.79-0.82, per the histogram)?
cat("\nIsolates in the 0.79-0.82 spike specifically:", 
    sum(my_isolates_BGE_glucose$BGE > 0.79 & my_isolates_BGE_glucose$BGE <= 0.82, na.rm=TRUE), "\n")
cat("Isolates in 0.70-0.79 (before the spike):", 
    sum(my_isolates_BGE_glucose$BGE > 0.70 & my_isolates_BGE_glucose$BGE <= 0.79, na.rm=TRUE), "\n")

