# Linear Regression Tests

# MAG-GH----

linear_model <- lm((GH_total) ~ genome.size_mean, data=GH_TOTAL_m.1) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# MAG-Protein----

linear_model <- lm((PR_total) ~ genome.size_mean, data=PR_TOTAL_m.1) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# MAG-Transporters----

linear_model <- lm((transp_total) ~ genome.size_mean, data=TRANSP_TOTAL_m.1) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# MAG-Osmolytes----

linear_model <- lm((OSMO_total) ~ genome.size_mean, data=OSMO_TOTAL.MAG_m.1) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# ISOLATES-GH----

linear_model <- lm((GH_total) ~ genome.size_mean, data=GH_TOTAL_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# ISOLATES-Protein----

linear_model <- lm((PR_total) ~ genome.size_mean, data=PH_TOTAL_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# ISOLATES-Transporters----

linear_model <- lm((transp_total) ~ genome.size_mean, data=TRANSP_TOTAL_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# ISOLATES-Osmolytes----

linear_model <- lm((OSMO_total) ~ genome.size_mean, data=OSMO_TOTAL_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

