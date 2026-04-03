################################################################i
# Title: Site map for Chapter 3
#
# Author: Lydia Walton
# Last updated: 14-Jan-2026
################################################################i

# Setup ----
rm(list = ls(all=T)) #Clear environment
#Set working directory
library(here)
#Data manipulation
library(dplyr)
library(tidyr)
#Data visualization
library(ggplot2)
library(kableExtra)
library(patchwork)
#Making the site map
library(maps) 
library(mapdata)
library(tidyverse)
library(PBSmapping)
library(cowplot)
library(mapproj)
library(sf)
library(ggspatial)

#Create theme for plotting

#Font change to Times New Roman as "A"
#windowsFonts(A = windowsFont("Times New Roman"))

#Lydia's theme
LW_theme <- theme_classic() +
  theme(text = element_text(family = "A"),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = "black"))


###################i
# Site map ----
###################i

## Van Isle ----

#Use this data for the Van Isle map 
VanIsle.sf <- read_sf("MS Data/bc-coast.shp")

#Figure out the coordinates 
ggplot() +
  geom_sf(data = VanIsle.sf, fill = "grey75") +
  coord_sf(xlim = c(-129, -122), ylim = c(47.5, 51.5)) +
  LW_theme

# crop map to show the west coast of Canada
VanIsle.Inset <- ggplot() +
  geom_sf(data = VanIsle.sf, fill = "grey75", color = "white") +
  coord_sf(xlim = c(-128.5, -122), ylim = c(48, 51.5)) +
  LW_theme +
  labs(x = "Longitude", y = "Latitude") +
  #theme(axis.title.x = element_blank(),
  #      axis.title.y = element_blank(),
  #      axis.text = element_blank(),
  #     axis.ticks = element_blank(),
  #     plot.margin=grid::unit(c(0,0,0,0), "mm"),
  #      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #Add text for Canada and Vancouver Island
  annotate("text", x = -123, y = 51, label = "Canada", 
           size = 5, fontface = "bold") +
  annotate("text", x = -126.3, y = 50.1, label = "Vancouver \nIsland", 
           size = 4, fontface = "bold") +
  #Add text for Bamfiels and Sidney panels
  annotate("text", x = -125.25, y = 48.7, label = "A", 
           size = 3, fontface = "bold") +
  annotate("text", x = -123.45, y = 48.45, label = "B", 
           size = 3, fontface = "bold") +
  #Add box where Bamfield and Sidney sites will be
  annotate("rect",xmin = -124.9, xmax = -125.55, ymin = 48.75, ymax = 49.1, 
           col = "red", fill = NA, size = 0.8) +
  annotate("rect",xmin = -123.3, xmax = -123.6, ymin = 48.5, ymax = 48.75, 
           col = "red", fill = NA, size = 0.8)

VanIsle.Inset

#png("MS Figures/Fig_VanIsle.png", width = 10, height = 9, units = "in", res = 600)
#VanIsle.Inset
#dev.off()



## Bamfield ----

#Read in sites
POsites <- read.csv("MS Data/Summer2024_SiteMap.csv")

#Filter for Bamfield (sites to be included in this map)
POsites.bamfield <- POsites %>% 
  filter(Locality == "Bamfield")


#Hakai basemap
hakai.basemap <- read_sf("MS Data/COAST_TEST2.shp")
st_crs(hakai.basemap)


barkley.sound.surveys <- st_crop(hakai.basemap,
                                 c(xmin = -125.1, xmax = -125.18, ymin = 48.81, ymax = 48.853))


#Map
SiteMap.Bamfield <- ggplot() + 
  LW_theme +
  geom_sf(data = barkley.sound.surveys, color= "white", fill='grey75') +
  annotation_scale(location = "bl", width_hint = 0.3) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  #Add site locations
  geom_point(data=POsites.bamfield, 
             aes(x=Lon, y=Lat, 
                 color = Site.description, shape = Site.description), 
             size=4) +
  scale_colour_manual(values = c("#FA9F42", "#3185FC")) +
  scale_shape_manual(values = c(17,16)) +
  labs(shape = "Site description", color = "Site description") +
  theme(legend.position = "top",
        #legend.justification = "top",
        axis.title = element_text(size = 12)) +
  annotation_north_arrow(location = "tl", 
                         pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_nautical(
                           fill = c("grey40", "white"),
                           line_col = "grey20"))
#Add site labels 
SiteMap.Bamfield <- SiteMap.Bamfield +
  geom_label(aes(label= "BMSC", x=-125.1352, y= 48.837), size=4) +
  geom_label(aes(label= "EB", x=-125.147, y= 48.8352), size=4) +
  geom_label(aes(label= "GN", x=-125.1109, y= 48.836), size=4) +
  geom_label(aes(label= "BB", x=-125.153518, y= 48.829), size=4) +
  geom_label(aes(label= "BI", x=-125.136202, y= 48.8335), size=4) +
  geom_text(aes(label= "Bamfield", x=-125.125, y= 48.814), size=5, fontface = "bold") +
  geom_text(aes(label= "Trevor \nChannel", x=-125.154, y= 48.841), size=5,fontface = "bold")


SiteMap.Bamfield


#png("MS Figures/Fig_Bamfield-Sites.png", width = 10, height = 9, units = "in", res = 600)
#SiteMap.Bamfield
#dev.off()



## Sidney ----

#Filter for Sidney (sites to be included in this map)
POsites.sidney <- POsites %>% 
  filter(Locality == "Sidney")


sidney.surveys <- st_crop(hakai.basemap,
                          c(xmin = -123.2, xmax = -123.7, ymin = 48.5, ymax = 48.75))

## Map
SiteMap.Sidney <- ggplot() + 
  LW_theme +
  geom_sf(data = sidney.surveys, color= "white", fill='grey75') +
  annotation_scale(location = "bl", width_hint = 0.4) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  #Add site locations
  geom_point(data=POsites.sidney, 
             aes(x=Lon, y=Lat, 
                 color = Site.description, shape = Site.description), 
             size=4) +
  scale_colour_manual(values = c("#FA9F42", "#3185FC")) +
  scale_shape_manual(values = c(17,16)) +
  labs(shape = "Site description", color = "Site description") +
  theme(legend.position = "top",
        #legend.justification = "top",
        axis.title = element_text(size = 12)) +
  annotation_north_arrow(location = "tl", 
                         pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_nautical(
                           fill = c("grey40", "white"),
                           line_col = "grey20"))
#Add site labels 
SiteMap.Sidney <- SiteMap.Sidney +
  geom_label(aes(label= "YYJ", x=-123.430729, y= 48.65), size=4) +
  geom_label(aes(label= "CB", x=-123.47369, y= 48.676), size=4) +
  geom_label(aes(label= "GB", x=-123.393472, y= 48.657), size=4) +
  geom_label(aes(label= "MP", x=-123.485151, y= 48.70), size=4) +
  geom_label(aes(label= "SR", x=-123.446249, y= 48.71), size=4) +
  geom_text(aes(label= "North \nSaanich", x=-123.425, y= 48.68), size=4, fontface = "bold") +
  geom_text(aes(label= "Saanich \nInlet", x=-123.502569, y= 48.627076), size=4.5,fontface = "bold")


SiteMap.Sidney


#png("MS Figures/Fig_Sidney-Sites.png", width = 10, height = 9, units = "in", res = 600)
#SiteMap.Sidney
#dev.off()








