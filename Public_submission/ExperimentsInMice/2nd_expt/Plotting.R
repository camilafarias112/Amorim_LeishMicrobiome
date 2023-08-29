# In this script I modifed the GraphPad file sent to me by Tej Singh - I just out in more tidy-text file and here I plot the results
# EXPT 2 - POST-REVIEWERS

# Libraries ----
library(ggpubr)
library(ggbreak)
library(ggforce)
library(ggthemes)
library(patchwork)
library(tidyverse)

# Import colors ----
load("../../../../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/Dark24")

# Import files and plotting----
# Flow cytometry results and CFUs:
Flow_CFUs_Parasite <- read_delim("Flow_CFUs_Parasite.txt", 
                                 delim = "\t", escape_double = FALSE, 
                                 col_types = cols(group = col_factor(levels = c("PBS","Sa_IgG", "Sa_aIL1B", "Sa_aIL1R"))), 
                                 trim_ws = TRUE)

Flow_CFUs_Parasite %>%
  pivot_longer(!sampleID & !group, names_to = "measurement", values_to = "value") %>%
  filter(measurement %in% c("CD90.2p","CD45p",
                            "CD4","CD8",
                            "Neutrophils")) %>%
  filter(group != "Sa_aIL1A") %>%
  mutate(value = value/1000) %>%
ggplot(., aes(x=group, y=value,
              color=group)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=3) +
  theme_classic() +
  scale_color_manual(values = c("gray",
                                Dark24[1],
                                #Dark24[6],
                                Dark24[7],
                                Dark24[16])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 17), legend.title = element_blank(),
        axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("N cells 10^3")) +
  facet_wrap(. ~fct_inorder(measurement), scales = "free")

# Supplemental image:
Flow_CFUs_Parasite %>%
  pivot_longer(!sampleID & !group, names_to = "measurement", values_to = "value") %>%
  filter(measurement %in% c("CD90.2n","abT")) %>%
  filter(group != "Sa_aIL1A") %>%
  mutate(value = value/1000) %>%
  ggplot(., aes(x=group, y=value,
                color=group)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=3) +
  theme_classic() +
  scale_color_manual(values = c("gray",
                                Dark24[1],
                                #Dark24[6],
                                Dark24[7],
                                Dark24[16])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="right", legend.text = element_text(size = 14), legend.title = element_blank(),
        axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("N cells 10^3")) +
  facet_wrap(. ~ factor(measurement, levels = c("CD90.2n","abT","CD4")), scales = "free")

# CFUs:
cfus <- Flow_CFUs_Parasite %>%
  pivot_longer(!sampleID & !group, names_to = "measurement", values_to = "value") %>%
  filter(measurement %in% c("CFUs")) %>%
  mutate(value = value+1) %>%
  ggplot(., aes(x=group, y=value,
                color=group)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=3) +
  theme_classic() +
  scale_color_manual(values = c("gray",
                                Dark24[1],
                                #Dark24[6],
                                Dark24[7],
                                Dark24[16])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 17), legend.title = element_blank(),
        axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("CFUs/ear (S. aureus)")) +
  scale_y_break(c(10,400000), scales = 10) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + annotation_logticks(sides = "left", outside = TRUE) +
  coord_cartesian(clip = "off")

# Parasite counts:
parasite_counts <- Flow_CFUs_Parasite %>%
  pivot_longer(!sampleID & !group, names_to = "measurement", values_to = "value") %>%
  filter(measurement %in% c("Parasite_load")) %>%
  ggplot(., aes(x=group, y=value,
                color=group)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=3) +
  theme_classic() +
  scale_color_manual(values = c("gray",
                                Dark24[1],
                                #Dark24[6],
                                Dark24[7],
                                Dark24[16])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="none", legend.text = element_text(size = 17), legend.title = element_blank(),
        axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  xlab("") + ylab(paste("Parasite load")) +
  ylim(2,5)

parasite_counts / cfus
  
# Skin thickness and Pathology score:
Thickness_Pathology <- read_delim("Thickness_Pathology.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_types = cols(group = col_factor(levels = c("PBS","Sa_IgG", "Sa_aIL1B", "Sa_aIL1R"))), 
                                  trim_ws = TRUE)

# Sa infection:
Sa_infection <- Thickness_Pathology %>%
  filter(measurement %in% "Skin_tickness") %>%
  filter(group %in% c("PBS","Sa_IgG")) %>%
  pivot_longer(!sampleID & !group & !measurement, names_to = "weeks", values_to = "value") %>%
  mutate(weeks = case_when(
    weeks %in% "week0" ~ 0,
    weeks %in% "week1" ~ 1,
    weeks %in% "week2" ~ 2,
    weeks %in% "week3" ~ 3,
    weeks %in% "week4" ~ 4,
    weeks %in% "week5" ~ 5,
    weeks %in% "week6" ~ 6)) %>%
  ggplot(., aes(x=weeks, y=value,
                color=group)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="right", strip.background = element_blank(), strip.text = element_text(size=17)) +
  stat_summary(aes(group=group, color=group), fun = mean, geom = "line", size=1) +
  stat_summary(aes(group=group, color=group), fun = mean, geom = "point", size=2.5) +
  scale_color_manual(values = c("gray",
                                Dark24[1],
                                #Dark24[6],
                                Dark24[7],
                                Dark24[16])) + 
  #stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.2) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", alpha=0.6) +
  xlab("Weeks")

# Treatment:
treatments <- Thickness_Pathology %>%
  filter(measurement %in% "Skin_tickness") %>%
  filter(group %in% c("Sa","Sa_IgG","Sa_aIL1B","Sa_aIL1R")) %>%
  pivot_longer(!sampleID & !group & !measurement, names_to = "weeks", values_to = "value") %>%
  mutate(weeks = case_when(
    weeks %in% "week0" ~ 0,
    weeks %in% "week1" ~ 1,
    weeks %in% "week2" ~ 2,
    weeks %in% "week3" ~ 3,
    weeks %in% "week4" ~ 4,
    weeks %in% "week5" ~ 5,
    weeks %in% "week6" ~ 6)) %>%
  ggplot(., aes(x=weeks, y=value,
                color=group)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17), 
        axis.text.y = element_text(size = 17), axis.title = element_text(size = 17),
        legend.position="right", strip.background = element_blank(), strip.text = element_text(size=17)) +
  stat_summary(aes(group=group, color=group), fun = mean, geom = "line", size=1) +
  stat_summary(aes(group=group, color=group), fun = mean, geom = "point", size=2.5) +
  scale_color_manual(values = c(Dark24[1],
                                #Dark24[6],
                                Dark24[7],
                                Dark24[16])) + 
  #stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.2) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", alpha=0.6) +
  xlab("Weeks")

Sa_infection / treatments
