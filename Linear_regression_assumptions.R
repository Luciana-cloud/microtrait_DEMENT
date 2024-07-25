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

# MAG-Biofilm----

linear_model <- lm((BIO_total) ~ genome.size_mean, data=BIO_TOTAL.MAG_m.1) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# MAG-High Temp----

linear_model <- lm((TEMP_total) ~ genome.size_mean, data=TEMP.MAG_m.1) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# MAG-Aminoacid-Transporters----

linear_model <- lm(log(AMI_TRANSP_TOTAL+1) ~ genome.size_mean, data=AMI_TRANSP_TOTAL_m.1) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# MAG-Carbohydrate-Transporters----

linear_model <- lm(log(CAR_TRANSP_TOTAL+1) ~ genome.size_mean, data=CAR_TRANSP_TOTAL_m.1) 
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

# ISOLATES-Aminoacid-Transporters----

linear_model <- lm((AMI_TRANSP_TOTAL) ~ genome.size_mean, data=AMI_TRANSP_TOTAL.i_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# ISOLATES-Carbohydrate-Transporters----

linear_model <- lm((CAR_TRANSP_TOTAL) ~ genome.size_mean, data=CAR_TRANSP_TOTAL.i_m) 
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

# ISOLATES-Biofilm----

linear_model <- lm((BIO_total) ~ genome.size_mean, data=BIO_TOTAL.MAG_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# ISOLATES-Heat----

linear_model <- lm((TEMP_total) ~ genome.size_mean, data=TEMP_TOTAL.MAG_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))

# ISOLATES-pH----

linear_model <- lm((PH_total) ~ genome.size_mean, data=pH_TOTAL.MAG_m) 
summary(linear_model)            

# QQ-Plot
par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(linear_model, which = 3)
# 2. Normality
plot(linear_model, which = 2)
par(mfrow = c(1, 1))


