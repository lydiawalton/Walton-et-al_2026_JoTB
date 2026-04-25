#######################################################i
# Modelling Pisaster body temperatures in the field
#
# Author: Lydia Walton
#Date updated: 14-Jan-2026
#######################################################i

#dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Packages for plotting
library(ggplot2)
library(patchwork)
library(kableExtra)
library(EnvStats)
library(ggpubr)
library(EnvStats)
#Packages for statistics and modelling
library(glmmTMB)
library(lme4)
library(DHARMa)
library(sjPlot)
library(ggpmisc)

# Setup ----

#------------------------------------------Load data files

PO.bodytemps.Bamfield <- read.csv("MS Data/Pisaster body temperatures_Barkley Sound_2024.csv")

PO.bodytemps.Sidney <- read.csv("MS Data/Pisaster body temperatures_Sidney_2024.csv")

environ.parameters <- read.csv("MS Data/Site environmental parameters_2024.csv")

#-----------------------------------------Create theme for plotting

#Lydia's theme
LW_theme <- theme_classic() +
  theme(axis.title.x = element_text(vjust = -1, size = 14),
        axis.title.y = element_text(vjust = 2, size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = "black"))

## Data manipulation ----

######################i
## BAMFIELD DATA
######################i

#Make a subset of the data
## Exclude individuals with no central disc temperature and/or no size
PO.bodytemp.Bamfield.sub <- PO.bodytemps.Bamfield %>% 
  drop_na(disc.temp) %>% 
  drop_na(arm.length.cm)


#Change name for transect number
PO.bodytemp.Bamfield.sub <- PO.bodytemp.Bamfield.sub %>% 
  mutate(transect.number = recode(transect.number, "1" = "Low tide line", "2" = "High tide line"))

#Rename df
Bamfield.master.df <- PO.bodytemp.Bamfield.sub


######################i
## SIDNEY DATA
######################i

#Make a subset of the data
## Exclude individuals with no central disc temperature and/or no size
PO.bodytemp.Sidney.sub <- PO.bodytemps.Sidney %>% 
  drop_na(disc.temp) %>% 
  drop_na(arm.length.cm)


#Change name of transect number
PO.bodytemp.Sidney.sub <- PO.bodytemp.Sidney.sub %>% 
  mutate(transect.number = recode(transect.number, "1" = "Low tide line", "2" = "High tide line",
                                  "1A" = "Low tide line", "1B" = "Low tide line",
                                  "2B" = "High tide line"))

#Rename df
Sidney.master.df <- PO.bodytemp.Sidney.sub

#Change name for animals that were measured but were not on the survey transect
Sidney.outside.sub <- Sidney.master.df %>% 
  filter(transect.number != "Outside transect")

#Make df with both localities for COMPARISONS
Allsites.master <- rbind(Bamfield.master.df, Sidney.master.df)

#Add environmental parameters to master df
environ.parameters <- environ.parameters %>% 
  drop_na(sst)

Allsites.master <- merge(Allsites.master, environ.parameters, by = c("locality", "site", "transect.number"))

#Use rescaled body size (scales and centres the data for modelling)
Allsites.master$arm.length_z <- arm::rescale(Allsites.master$arm.length.cm)

#Use rescaled substrate temp (scales and centres the data for modedlling)
Allsites.master$substrate_z <- arm::rescale(Allsites.master$substrate.ave)


#######################i
## Visualization  ----
#######################i

# Make a split plot with each panel showing the relationship between body temp and each predictor
## Use colors to represent different sites (does the relationship between predictor and response change (different slope) across sites?)

### arm length
arm.length.plot <- ggplot(Allsites.master, aes(x = arm.length.cm, y = disc.temp, color = site)) +
  geom_smooth(method = lm) +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Arm length (cm)", y = "Body Temperature")

### aggregated
aggregated.plot <- ggplot(Allsites.master, aes(x = aggregated, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Aggregation", y = "Body Temperature")

### substrate temp
substrate.temp.plot <- ggplot(Allsites.master, aes(x = substrate.ave, y = disc.temp, color = site)) +
  geom_smooth(method = lm) +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Substrate temperature", y = "Body Temperature")

### sst
sst.plot <- ggplot(Allsites.master, aes(x = sst, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "SST", y = "Body Temperature")

### air
air.plot <- ggplot(Allsites.master, aes(x = air, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Air", y = "Body Temperature")

### humidity
humidity.plot <- ggplot(Allsites.master, aes(x = humidity, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Humidity", y = "Body Temperature")

### wind (daily average)
wind.ave.plot <- ggplot(Allsites.master, aes(x = wind.ave.km, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Average daily wind", y = "Body Temperature")

### wind (daily max)
wind.max.plot <- ggplot(Allsites.master, aes(x = wind.max.km, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Maximum daily wind", y = "Body Temperature")

### transect
transect.plot <- ggplot(Allsites.master, aes(x = transect.number, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Transect", y = "Body Temperature")

### microhabitat
microhabitat.plot <- ggplot(Allsites.master, aes(x = microhabitat, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Microhabitat", y = "Body Temperature")

### site
site.plot <- ggplot(Allsites.master, aes(x = site, y = disc.temp, color = site)) +
  geom_boxplot() +
  scale_color_manual(values = c("#77B6EA", "#1C7293", "#8FD694", "#355070",
                                "#295135", "#065A82", "#7DBA84", "#0E402D")) +
  LW_theme +
  theme(legend.position = "none") +
  labs(x = "Site", y = "Body Temperature")

### Split plot ----

split.plot <- arm.length.plot + aggregated.plot + substrate.temp.plot + sst.plot + air.plot +
  humidity.plot + wind.ave.plot + wind.max.plot + transect.plot + microhabitat.plot + site.plot +
  plot_layout(ncol = 3)

split.plot

####################i
## Modelling ----
####################i

### Dredge  ----

library(MuMIn)

# Predictor variables

## Continuous: Arm length, Substrate temperature
## Categorical: Aggregation, Transect, Microhabitat
## One measurement/day: SST, Air, Humidity, Avg wind speed, Max wind speed

## Random effect: Site

# Fit the full model with all predictors
global.model <- glmmTMB(disc.temp ~ arm.length_z + substrate_z +
                          aggregated + transect.number + microhabitat +
                          sst + air + humidity + wind.ave.km + wind.max.km
                          + (1 | site),
                        family = gaussian, data = Allsites.master,
                        na.action = na.fail)

# Perform model selection using dredge (NOTE: This step takes a while in R)
dredged.model <- dredge(global.model)

# View the model selection table (NOTE: This step also takes a while to run)
#print(model_selection)

# Get the best model
best.model <- get.models(dredged.model, 1)[[1]]

# Print the best model
summary(best.model)

tab_model(best.model)

# If there are several good models (close in AIC), you can perform model averaging using model.avg

#Model averaging for a subset of models (NOTE: This step does not make sense for interactions)
avg.model <- model.avg(dredged.model, subset = delta < 2)  #Models with ΔAIC < 2
summary(avg.model)

#### Top dredge models ----

# Select the top models whose cumulative AIC weights sum to 95%

# Extract the AIC weights
dredged.model$CumWeight <- cumsum(dredged.model$weight)

# Subset models that are within the top 95% cumulative AIC weights
top.95.percent <- subset(dredged.model, CumWeight <= 0.95)

# View the top 95% models
print(top.95.percent)

# Change this into a df so we can look at the model selection more easily
top.95.percent.df <- as.data.frame(top.95.percent)

#write.csv(top.95.percent.df, file = "top.95.percent.df.csv", row.names = FALSE)

### Q: How many models are in the top 95%? Which predictors were most commonly included? 

# Look at the models that are within 2 AIC units of the best model
top.80.percent <- subset(dredged.model, CumWeight <= 0.80)

# View the top 80% models
print(top.80.percent)

top.80.percent.df <- as.data.frame(top.80.percent)


#### Collinearity diagnostics ----

# Collinearity refers to the degree of correlation between predictor variables.
# If predictor variables are highly collinear, it can lead to unreliable estimates and poor model performance. 

# You can use diagnostic tools (such as Variance Inflation Factor (VIF)) to check for collinearity among predictors

library(car)

# Refit the model without random effects
fixed.effects.model <- lm(disc.temp ~ arm.length_z + aggregated +  
                            substrate_z + sst + air + humidity +
                            wind.ave.km + wind.max.km + 
                            transect.number + microhabitat,
                          data = Allsites.master)

# Calculate VIF for the fixed-effects model
vif(fixed.effects.model)

## Air and sst are potentially collinear - not an issue based on above
## Wind avg and max are also potentially collinear - keeping only average wind anyways (based on best model)

################################i
# Summary of best dredged model
################################i

# 1) The random effects from site explain some variation in disc.temp, while locality contributes negligible random variation.

# 2) The fixed effects show that many predictors significantly affect disc.temp, 
#    Positive relationship = solitary, air temp, humidity, substrate temp
#    Negative relationship = microhabitat (PE, Sh, TP, UR), SST
#    Not significant = Wind ave

# 3) Some fixed effects were not included in the "best" model
#    dropped terms = arm.length, transect.number, and wind max


### GLMM ----
#### Table S4 ----

## Note: bt stands for body temp

bt.glmm <- glmmTMB(disc.temp ~ aggregated + air + humidity + sst + substrate_z 
                   + wind.ave.km + microhabitat +
                     (1 | site),
                   data = Allsites.master, family = gaussian)

summary(bt.glmm)

#visualize the model output
tab_model(bt.glmm,
          show.est = TRUE)

#Simulate the residuals using dharma
bt.dharma <- simulateResiduals(fittedModel = bt.glmm, plot = T)

#Look at dispersion of residuals
testDispersion(bt.glmm, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa") #High p-value = no evidence of dispersion issues

#Look at the effects of predictors on body temp
plot(predict(bt.glmm), Allsites.master$disc.temp, 
     xlab = "Predicted", ylab = "Observed")
abline(0, 1, col = "red")

#Plot all effects on a split plot
library(effects)
plot(allEffects(bt.glmm)) #Effect line represents how the predictor affects the response


### GLMM interaction ----
bt.glmm.interaction <- glmmTMB(disc.temp ~ aggregated + air*humidity + sst + substrate_z 
                   + wind.ave.km + microhabitat +
                     (1 | site),
                   data = Allsites.master, family = gaussian)

summary(bt.glmm.interaction)

#visualize the model output
tab_model(bt.glmm.interaction,
          show.est = TRUE)

#Simulate the residuals using dharma
bt.dharma.interaction <- simulateResiduals(fittedModel = bt.glmm.interaction, plot = T)
## Residuals seem to deviate more with the interaction than without

#Look at dispersion of residuals
testDispersion(bt.glmm.interaction, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa") #High p-value = no evidence of dispersion issues

#Look at the effects of predictors on body temp
plot(predict(bt.glmm.interaction), Allsites.master$disc.temp, 
     xlab = "Predicted", ylab = "Observed")
abline(0, 1, col = "red")

#Compare the AIC of the glmm with and without the interaction
AIC(bt.glmm, bt.glmm.interaction) #Not very different; go with the less complex model (without the interaction)


### GAM ----

library(gamm4)

#Change site to a factor
Allsites.master$site <- as.factor(Allsites.master$site)

#Make the gam 
bt.gam <- gam(disc.temp ~ aggregated + microhabitat + substrate_z +
                s(air, k = 4) + s(humidity, k = 4) + s(sst, k = 4) + s(wind.ave.km, k = 6) +
                s(site, bs = "re"),
              data = Allsites.master)

summary(bt.gam)

tab_model(bt.gam)
## edf (estimated df) shows linearity - if edf is close to 1 the relationship is nearly linear
##     if edf is higher than 1, this indicates more flexibility; modelling a nonlinear relationship
## air and sst are linear; humidity is almost linear (probably a product of sampling effort); wind is not linear

#plot the smooth terms
plot(bt.gam, pages = 1, shade = TRUE)

#check the residuals
gam.check(bt.gam)
## Histogram looks normally distributed
## Residuals aren't completely random - more at higher body temps
## Q-Q plot is decent fit until temps above 20 are reached

#test concurvity (analogous to multicollinearity in linear models)
concurvity(bt.gam, full = TRUE)
## Values close to 1 indicate high concurvity, suggesting predictors are highly correlated
## likely means the k needs to be adjusted...

#Visualize the outputs
#https://mfasiolo.github.io/mgcViz/reference/plot.gamViz.html

library(mgcViz)

bt.viz <- getViz(bt.gam)

#Look at the smooth terms
print(plot(bt.viz), ask = FALSE)

#As sst increases, body temperature did not increases - body temp was not tracking sst, body temperature was relatively lower than sst
#Air and humidity likely driving more of the body temp at low tide 

#Animals keeping body temp within a certain window regardless of sst
#Body temp on average was varying less than sst, as sst was going up, body temp was remained lower than sst


### GAM interaction ----

#Try an interaction in the model
bt.gam.interaction <- gam(disc.temp ~ aggregated + microhabitat + substrate_z +
                            s(air, humidity, k = 4) + s(sst, k = 4) + s(wind.ave.km, k = 6) +
                            s(site, bs = "re"),
                          data = Allsites.master, method = "REML")

summary(bt.gam.interaction)

#Look at the model output
tab_model(bt.gam.interaction)

#Visualize the smoothed terms
bt.smooth.plot <- getViz(bt.gam.interaction)
print(plot(bt.smooth.plot), ask = FALSE)

#Plot only the first smooth, looking at the interaction between humidity and air
plot(bt.smooth.plot, select = 1) + l_dens(type = "cond") + l_fitLine() + l_ciLine()






