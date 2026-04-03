#######################################################################################################i
# Title: Describing Pisaster body temperatures in the rocky shore intertidal of Vancouver Island
#
# Author: Lydia Walton
# Last updated: 04-Apr-2026
######################################################################################################i

#graphics.off() #Clear open plots
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
library(see)
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

#Drop anything with NA in SST column - gets rid of extra rows in csv
environ.parameters <- environ.parameters %>% 
  drop_na(sst)

Allsites.master <- merge(Allsites.master, environ.parameters, by = c("locality", "site", "transect.number"))

#Use rescaled body size (scales and centres the data for modelling)
Allsites.master$arm.length_z <- arm::rescale(Allsites.master$arm.length.cm)


##########################i
# Body temperature  ----
##########################i

#Min, max and mean body temperatures at each location
bt.summarytable.location <- Allsites.master %>% 
  group_by(locality) %>% 
  summarise(min.bt = min(disc.temp), max.bt = max(disc.temp),
            mean.bt = mean(disc.temp), SD.bt = sd(disc.temp)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

bt.summarytable.location

#Min, max and mean body temperatures at each site and transect
bt.summarytable.site <- Allsites.master %>% 
  group_by(locality, site, date, transect.number) %>% 
  summarise(min.bt = min(disc.temp), max.bt = max(disc.temp),
            mean.bt = mean(disc.temp), SD.bt = sd(disc.temp)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

bt.summarytable.site

#Looking for the extremes (cooler than the lowest SST and hotter than the highest air temp)
##What is the % of animals that were cooler than 14.5C (lowest SST recorded at any site)
coolest.bt <- Allsites.master %>% 
  filter(disc.temp < 14.5) #86 individuals cooler than 14.5

86/740*100 #11.6%

##What is the % of animals hotter than 26.6C (highest air temp recorded at any site)
hottest.bt <- Allsites.master %>% 
  filter(disc.temp > 26.6) #zero individuals hotter than the highest air temp

#What is the % of animals equal to or hotter than 23C (threshold for negative feeding effects; Pincebourde et al. 2008)
over23.bt <- Allsites.master %>% 
  filter(disc.temp > 22.9) #6 animals

6/740*100 #0.81% equal to or over 23C


##########################i
# SST and Air temp ----
##########################i

#Look at min, mean and max air across locality
air.summarytable.locality <- Allsites.master %>% 
  group_by(locality, site, date) %>% 
  summarise(min = min(air), max = max(air),
            mean = mean(air), SD = sd(air)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

air.summarytable.locality

#Look at min, mean and max sst across locality
sst.summarytable.locality <- Allsites.master %>% 
  group_by(locality, site, date) %>% 
  summarise(min = min(sst), max = max(sst),
            mean = mean(sst), SD = sd(sst)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

sst.summarytable.locality

#Look at min, mean and max humidity across locality
humidity.summarytable.locality <- Allsites.master %>% 
  group_by(locality, site, date) %>% 
  summarise(min = min(humidity), max = max(humidity),
            mean = mean(humidity), SD = sd(humidity)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

humidity.summarytable.locality


#############################i
# Substrate temperature ----
#############################i

#Look at min, mean and max body size across sites
substrate.summarytable.site <- Allsites.master %>% 
  group_by(locality, site, transect.number) %>% 
  summarise(min.temp = min(substrate.ave), max.temp = max(substrate.ave),
            mean.temp = mean(substrate.ave), SD.temp = sd(substrate.ave)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

substrate.summarytable.site


####################i
# Body sizes ----
####################i

#Look at min, mean and max body size across sites
size.summarytable.site <- Allsites.master %>% 
  group_by(locality, site, transect.number) %>% 
  summarise(min.size = min(arm.length.cm), max.size = max(arm.length.cm),
            mean.size = mean(arm.length.cm), SD.size = sd(arm.length.cm)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

size.summarytable.site

#Look at min, mean and max body size across locality
size.summarytable.locality <- Allsites.master %>% 
  group_by(locality) %>% 
  summarise(min.size = min(arm.length.cm), max.size = max(arm.length.cm),
            mean.size = mean(arm.length.cm), SD.size = sd(arm.length.cm)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

size.summarytable.locality

#Histogram of body sizes
hist(Allsites.master$arm.length.cm)

#Bin body sizes for plotting
## <=6 cm is juvenile; >6 cm is adult (based on Gooding and Harley 2015)
Bodysize.plotting <- Allsites.master %>% 
  mutate(life.stage = case_when(
    arm.length.cm <= 6 ~ "Juvenile", 
    arm.length.cm > 6 ~ "Adult",
    TRUE ~ "OTHER"
  ))
  
## Plot the distribution across location
Bodysize.plot.location <- Bodysize.plotting %>% 
  ggplot(aes(x = locality, fill = life.stage)) +
  geom_bar() +
  scale_fill_manual(values = c("#5E4B56", "#BBC8CA")) +
  scale_y_continuous(breaks = seq(from = 0, to = 600, by = 50)) +
  LW_theme

Bodysize.plot.location

## Across microhabitats
Bodysize.plot.microhabitat <- Bodysize.plotting %>% 
  ggplot(aes(x = microhabitat, fill = life.stage)) +
  geom_bar() +
  scale_fill_manual(values = c("#5E4B56", "#BBC8CA")) +
  scale_y_continuous(breaks = seq(from = 0, to = 600, by = 50)) +
  LW_theme

Bodysize.plot.microhabitat

#Table showing number of juveniles vs adults across locality
lifestage.table <- Bodysize.plotting %>% 
  group_by(locality, life.stage) %>% 
  summarise(n = length(life.stage)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

lifestage.table

#Table showing the number of juveniles vs adults across microhabitats and locality
lifestage.MH.table <- Bodysize.plotting %>% 
  group_by(locality, microhabitat, life.stage) %>% 
  summarise(n = length(life.stage)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

lifestage.MH.table

##################i
# Figure ----
##################i

## Statistics ----

### Arm length ----
## Table S1

bt.size.glmm <- glmmTMB(disc.temp ~ arm.length_z*locality + (1|site),
                        data = Allsites.master, family = gaussian)

summary(bt.size.glmm)

#Check the model fit
bt.size.dharma <- simulateResiduals(fittedModel = bt.size.glmm, plot = F)

plot(bt.size.dharma) #Some slight issues with fit, no outliers, dispersion good

tab_model(bt.size.glmm,
          show.est = TRUE) #visualizing model outputs

### Aggregation ----
## Table S2

bt.aggregation.glmm <- glmmTMB(disc.temp ~ aggregated*locality
                               + (1|site),
                               data = Allsites.master, family = gaussian)

summary(bt.aggregation.glmm)


bt.aggregation.dharma <- simulateResiduals(fittedModel = bt.aggregation.glmm, plot = F)

plot(bt.aggregation.dharma) #Within-group deviations detected, no outliers, dispersion good

tab_model(bt.aggregation.glmm,
          show.est = TRUE) #visualizing model outputs


## Plots ----

### Arm length ----
bt.size <- Allsites.master %>% 
  ggplot(aes(x = arm.length.cm, y = disc.temp, color = locality)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = gaussian)) +
  stat_poly_eq(use_label(c("adj.R2", "n")),
               label.x = 0.95) +
  scale_color_manual(name = "Location",
                     values = c("grey75","grey25")) +
  labs(y = "Surface body temperature (°C)", x = "Arm length (cm)") +
  scale_y_continuous(breaks = seq(from = 10, to = 26, by = 2)) +
  scale_x_continuous(breaks = seq(from = 0, to = 18, by = 2)) +
  LW_theme +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 6, linetype = "dotted", color = "blue", linewidth = 1)

bt.size

### Aggregation ----

#Create sample sizes
n_values <- Allsites.master %>%
  group_by(locality, aggregated) %>%
  summarize(n = n())

#Plot
bt.aggregation <- Allsites.master %>% 
  ggplot(aes(x = locality, y = disc.temp, fill = aggregated)) +
  geom_boxplot() +
  scale_fill_manual(labels = c("Aggregated", "Solitary"),
                    values = c("grey90","grey30")) +
  labs(x = "Location", y = "Surface body temperature (°C)",
       fill = "Aggregation behaviour") +
  LW_theme +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(10, 27), breaks = seq(from = 10, to = 27, by = 2)) +
  #Add sample size
  geom_text(data = n_values, aes(x = locality, y = 11, label = paste("n =", n)), 
            position = position_dodge(width = 0.75), vjust = -0.5) +
  #Add significance
  geom_text(x = 1, y = 27, label = "ns") +
  geom_text(x = 2, y = 27, label = "**")


bt.aggregation


#### MS plot (thermal inertia) ----
thermal.inertia <- bt.size + bt.aggregation +
  plot_annotation(tag_levels = "A")

thermal.inertia

#Save figure
#png("MS Figures/Fig.Thermal_Inertia.png", width = 9, height = 7, units = "in", res = 600)
#thermal.inertia
#dev.off()


##################i
# Figure ----
##################i

## Statistics ----

### Table S3 ----

#New df for modelling removing the ordered factors
Allsites.modelling <- Allsites.master

Allsites.modelling$transect.number <- as.character(Allsites.modelling$transect.number)

#glmm
bt.transect.glmm <- glmmTMB(disc.temp ~ transect.number*site
                            + (1|locality),
                            data = Allsites.modelling, family = gaussian)

summary(bt.transect.glmm)

#Check the model fit
bt.transect.dharma <- simulateResiduals(fittedModel = bt.transect.glmm, plot = F)

plot(bt.transect.dharma) #Some slight issues with fit, no outliers, dispersion good

tab_model(bt.transect.glmm,
          show.est = TRUE) #visualizing model outputs


## Plots ----

#Add a new column with the temperature difference between body and sst and air
Allsites.master <- Allsites.master %>% 
  mutate(delta.temp.SST = disc.temp - sst)

Allsites.master <- Allsites.master %>% 
  mutate(delta.temp.air = disc.temp - air)

Allsites.master <- Allsites.master %>% 
  mutate(delta.temp.sub = disc.temp - substrate.ave) 

### Percent diff ----

#Which animals were hotter than the ambient air?
bt.delta.air.hot <- Allsites.master %>% 
  filter(delta.temp.air > 0) #10 animals (80% found in exposed habitats)

#Which animals were cooler than ambient air?
bt.delta.air.cold <- Allsites.master %>% 
  filter(delta.temp.air < 0) #727 animals

#Min, max, mean and SD temp difference
bt.delta.air <- Allsites.master

bt.delta.air.table <- bt.delta.air %>% 
  group_by(locality, site) %>% 
  summarise(min = min(delta.temp.air), max = max(delta.temp.air),
            mean = mean(delta.temp.air), SD = sd(delta.temp.air)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

bt.delta.air.table

#Which animals were colder than the ambient seawater?
bt.delta.sst.cold <- Allsites.master %>% 
  filter(delta.temp.SST < 0) #341 animals (52% of total animals)

#Min, max and mean temp difference
bt.delta.SST.cold.table <- bt.delta.sst.cold %>% 
  group_by(locality, site) %>% 
  summarise(min = min(delta.temp.SST), max = max(delta.temp.SST),
            mean = mean(delta.temp.SST), SD = sd(delta.temp.SST)) %>% 
  kbl() %>% 
  kable_classic(full_width = F, html = "Times New Roman")

bt.delta.SST.cold.table


#What is the proportion of animals cooler than SST across microhabitats
bt.delta.sst.habitat <- bt.delta.sst.cold %>% 
  filter(microhabitat == "under rock") #231 animals 

231/341*100 #67.7% cooler than SST Pisaster found underneath rocks


### Bamfield SST ----

Bamfield.plotting <- Allsites.master %>% 
  filter(locality == "Bamfield")

Bamfield.plotting$transect.number <- factor(Bamfield.plotting$transect.number,
                                            ordered = TRUE,
                                            levels = c("Low tide line", "High tide line"))

#Calculate n
n_values.transect <- Bamfield.plotting %>%
  group_by(locality, site, transect.number) %>%
  summarize(n = n())

#PLOT
bt.SST.bam <- Bamfield.plotting %>% 
  ggplot(aes(x = site, y = delta.temp.SST, fill = transect.number)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(name = "Transect",
                    values = c("#A0DDE6", "#4357AD")) +
  facet_wrap(~locality) +
  LW_theme +
  theme(legend.position = "bottom",
        axis.title.x = element_blank()) +
  #Add zero line (no difference between sst and body temperature)
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", linewidth=1) +
  labs(x = "Site", y = expression(Delta ~ "temperature (body - SST)")) +
  scale_y_continuous(breaks = seq(from = -5, to = 10, by = 1)) +
  #Add sample size
  geom_text(data = n_values.transect, aes(x = site, y = -5, label = paste("n =", n)), 
            position = position_dodge(width = 1), vjust = -0.5,
            show.legend = FALSE)

bt.SST.bam

### Sidney SST ----

Sidney.plotting <- Allsites.master %>% 
  filter(locality == "Sidney")

Sidney.plotting$transect.number <- factor(Sidney.plotting$transect.number,
                                          ordered = TRUE,
                                          levels = c("Low tide line", "High tide line"))

##Calculate n
n_values.transect.sid <- Sidney.plotting %>%
  group_by(locality, site, transect.number) %>%
  summarize(n = n())

#PLOT
bt.SST.sid <- Sidney.plotting %>% 
  ggplot(aes(x = site, y = delta.temp.SST, fill = transect.number)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(name = "Transect",
                    values = c("#A0DDE6", "#4357AD")) +
  facet_wrap(~locality) +
  LW_theme +
  theme(legend.position = "bottom") +
  #Add zero line (no difference in temperature)
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", linewidth=1) +
  labs(x = "Site", y = expression(Delta ~ "temperature (body - SST)")) +
  scale_y_continuous(breaks = seq(from = -7, to = 6, by = 1)) +
  #Add sample size
  geom_text(data = n_values.transect.sid, aes(x = site, y = -6.7, label = paste("n =", n)), 
            position = position_dodge(width = 1), vjust = -0.5,
            show.legend = FALSE) 

bt.SST.sid

### Bamfield Air ----

#PLOT
bt.Air.bam <- Bamfield.plotting %>% 
  ggplot(aes(x = site, y = delta.temp.air, fill = transect.number)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(name = "Transect",
                    values = c("#9E2A2B", "#540B0E")) +
  facet_wrap(~locality) +
  LW_theme +
  theme(legend.position = "bottom",
        axis.title.x = element_blank()) +
  #Add zero line for no difference in temp
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", linewidth=1) +
  labs(x = "Site", y = expression(Delta ~ "temperature (body - air)")) +
  scale_y_continuous(breaks = seq(from = -16, to = 3, by = 2)) +
  #Add sample size
  geom_text(data = n_values.transect, aes(x = site, y = -16.7, label = paste("n =", n)), 
            position = position_dodge(width = 1), vjust = -0.5,
            show.legend = FALSE) 

bt.Air.bam


### Sidney Air ----

#PLOT
bt.Air.sid <- Sidney.plotting %>% 
  ggplot(aes(x = site, y = delta.temp.air, fill = transect.number)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(name = "Transect",
                    values = c("#9E2A2B", "#540B0E")) +
  facet_wrap(~locality) +
  LW_theme +
  theme(legend.position = "bottom") +
  #Add zero line for no difference in temperature
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", linewidth=1) +
  labs(x = "Site", y = expression(Delta ~ "temperature (body - air)")) +
  scale_y_continuous(breaks = seq(from = -18, to = 2, by = 2)) +
  #Add sample size
  geom_text(data = n_values.transect.sid, aes(x = site, y = -16.7, label = paste("n =", n)), 
            position = position_dodge(width = 1), vjust = -0.5,
            show.legend = FALSE) 

bt.Air.sid

#### MS plot ----

deltaT.env <- bt.SST.bam + bt.Air.bam + bt.SST.sid + bt.Air.sid +
  plot_layout(nrow = 2, ncol = 2,
              guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

deltaT.env

#Save figure
#png("MS Figures/Fig.deltaT_environment.png", width = 12, height = 10, units = "in", res = 600)
#deltaT.env
#dev.off()


# Figure ----

## Plots ----

### Bamfield ----

##Calculate n
n_values.bam <- Bamfield.master.df %>%
  group_by(microhabitat) %>%
  summarize(n = n())

#PLOT
bt.substrate.bam <- Bamfield.plotting %>% 
  ggplot() +
  geom_violinhalf(aes(x = microhabitat, y = delta.temp.sub, fill = microhabitat)) +
  coord_flip() +
  scale_fill_manual(name = "Microhabitat",
                    values = c("#89043D", "#682D63", "#0A1045", "#84C0C6", "#C17817")) +
  #facet_wrap(~transect.number) +
  LW_theme +
  theme(legend.position = "top",
        legend.box = "vertical") +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", linewidth=1) +
  guides(color = "none") +
  labs(x = "Microhabitat", y = expression(Delta ~ "temperature (body - substrate)")) +
  scale_y_continuous(breaks = seq(from = -8, to = 4, by = 1)) +
  #Add sample size
  geom_text(data = n_values.bam, aes(x = microhabitat, y = -7.4, label = paste("n =", n)), 
            position = position_dodge(width = 0.75), vjust = -0.5) 


bt.substrate.bam

### Sidney ----

##Calculate n
n_values.sid <- Sidney.master.df %>%
  group_by(microhabitat) %>%
  summarize(n = n())

#PLOT
bt.substrate.sid <- Sidney.plotting %>% 
  ggplot() +
  geom_violinhalf(aes(x = microhabitat, y = delta.temp.sub, fill = microhabitat)) +
  coord_flip() +
  scale_fill_manual(name = "Microhabitat",
                    values = c("#89043D", "#682D63", "#0A1045", "#C17817")) +
  #facet_wrap(~transect.number) +
  LW_theme +
  theme(legend.position = "none",
        legend.box = "vertical") +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", linewidth=1) +
  guides(color = "none") +
  labs(x = "Microhabitat", y = expression(Delta ~ "temperature (body - substrate)")) +
  scale_y_continuous(breaks = seq(from = -7, to = 2, by = 1)) +
  #Add sample size
  geom_text(data = n_values.sid, aes(x = microhabitat, y = -7.4, label = paste("n =", n)), 
            position = position_dodge(width = 0.75), vjust = -0.5) 

bt.substrate.sid

#### MS plot ----

bt.sub.microhabit <- bt.substrate.bam + bt.substrate.sid +
  plot_layout(ncol = 1, nrow = 2) +
  plot_annotation(tag_levels = "A") 

bt.sub.microhabit

#Save figure
#png("MS Figures/Fig.DeltaT_Microhabitat.png", width = 8, height = 9, units = "in", res = 600)
#bt.sub.microhabit
#dev.off()

