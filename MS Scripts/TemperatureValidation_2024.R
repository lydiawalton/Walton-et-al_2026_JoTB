##############################################################################i
# Title: Validation of thermal imaging camera using a digital thermometer
#
# Author: Lydia Walton
# Last updated: 04-Apr-2026
##############################################################################i

graphics.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Packages for plotting
library(ggplot2)

# Setup ----

#------------------------------------------Load data files

temp.validation.df <- read.csv("MS Data/Temperature-validation_2024.csv")

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


#Look at the data
hist(temp.validation.df$IR.temp)

hist(temp.validation.df$thermocouple.temp)

#pivot to long format
temp.validation.long <- temp.validation.df %>%
  pivot_longer(
    cols = c(IR.temp, thermocouple.temp),
    names_to = "Measurement_Type",
    values_to = "Temperature")

temp.validation.plot <- temp.validation.long %>% 
  ggplot(aes(x = Measurement_Type, y = Temperature, fill = Measurement_Type)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(name = "Measurement Type",
                    values = c("#A0DDE6", "#4357AD")) +
  LW_theme

temp.validation.plot

#Paired t-test
t.test(temp.validation.df$IR.temp, temp.validation.df$thermocouple.temp, paired = TRUE)

#Test for normality
hist(temp.validation.df$temp.diff)
qqnorm(temp.validation.df$temp.diff)
qqline(temp.validation.df$temp.diff)

shapiro.test(temp.validation.df$temp.diff)

#Regression
temp.lm <- lm(IR.temp ~ thermocouple.temp, data = temp.validation.df)

summary(temp.lm)

plot(temp.validation.df$thermocouple.temp, temp.validation.df$IR.temp,
     xlab = "Thermocouple Temperature",
     ylab = "IR Temperature")

abline(temp.lm, col = "blue", lwd = 2)

