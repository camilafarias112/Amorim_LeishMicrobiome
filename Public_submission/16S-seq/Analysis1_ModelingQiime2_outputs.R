# In this script, I bring in Qiime2 outputs and explore the microbiome profiles in CL lesions.

# Libraries ----
library(ggpubr)
library(ggthemes)
library(patchwork)
library(ggrepel)
library(gplots)
#library(qiime2R)
#library(DESeq2)
#library(phyloseq)
library(ggforce)
library(ggdist)
library(scales)

library(tidyverse)
theme_set(theme_classic())

# Color palette ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/Dark24")

# Load study design ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/targets_all")

# Load metadata and modeling ----
metadata <- read_q2metadata("Qiime2/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples_wTypesCol.tsv")
metadata$Time_point <- factor(metadata$Time_point, levels = c("day0", "day30","day60","day90",
                                                              "day120", "day150", "day180", "day210", "day240", "day270"))
colnames(metadata)[6] <- colnames(targets_all)[1]
metadata$LTCP_patient_ID <- as.character(metadata$LTCP_patient_ID)

metadata <- metadata %>%
  left_join(targets_all %>%
              select(LTCP_patient_ID, enviromental_samples, # do I need to include more?
                     sex, age, illness_duration_days,
                     local_of_lesion, size_lesion_mm2,
                     lymphadenopathy, DTH_mm2, treatment_outcome,
                     healing_time_days, treatment_other_drug, sample_RNAseq,
                     group_RNAseq), "LTCP_patient_ID")

# Shannon alpha diversity ----
shannon <- read_qza("Qiime2/core-metrics-results/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

metadata <- metadata %>% 
  left_join(shannon, by= "SampleID") 

shannon_plot <- metadata %>%
  filter(Time_point %in% "day0") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Swab_site2 = case_when(
    Swab_site == "lesion" ~ "Lesion",
    Swab_site == "contralateral" ~ "Contra.Skin")) %>%
  ggplot(aes(x=Swab_site2, y=shannon_entropy, color=Swab_site2)) +
  geom_violin(adjust = .6) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="black") +
  stat_compare_means(method = "t.test", label.y = 7,
                     aes(label = ..p.signif..),paired = T,
                     label.x = 1.4, size = 6, color="black") +
  xlab("") +
  ylab("Shannon Index")

shannon_plot2 <- metadata %>%
  filter(Time_point %in% c("day0","day30","day60","day90"),
         treatment_other_drug %in% "sbv") %>% # We don't have samples enough from other time points
  mutate(Swab_site2 = case_when(
    Swab_site == "lesion" ~ "Lesion",
    Swab_site == "contralateral" ~ "Contra.Skin")) %>%
  mutate(Time_point2 = case_when(
    Time_point == "day0" ~ "Day 0",
    Time_point == "day30" ~ "Day 30",
    Time_point == "day60" ~ "Day 60",
    Time_point == "day90" ~ "Day 90")) %>%
  ggplot(aes(x=Swab_site2, y=shannon_entropy, color=Swab_site2)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none", strip.background = element_blank(), strip.text = element_text(size = 11)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", label.y = 7,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Shannon Diversity") +
  facet_wrap(. ~ Time_point2, nrow = 1)

#Treatment outcome confirmation:
metadata %>%
  filter(Time_point %in% c("day0","day30","day60","day90"),
         treatment_other_drug %in% "sbv") %>% # We don't have samples enough from other time points
  mutate(Swab_site2 = case_when(
    Swab_site == "lesion" ~ "Lesion",
    Swab_site == "contralateral" ~ "Contra.Skin")) %>%
  mutate(Time_point2 = case_when(
    Time_point == "day0" ~ "Day 0",
    Time_point == "day30" ~ "Day 30",
    Time_point == "day60" ~ "Day 60",
    Time_point == "day90" ~ "Day 90")) %>%
  ggplot(aes(x=Swab_site2, y=shannon_entropy, color=Swab_site2)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none", strip.background = element_blank(), strip.text = element_text(size = 11)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", label.y = 7.5,
                     aes(label = ..p.signif..),
                     label.x = 1.3, size = 6, color="black") +
  xlab("") +
  ylab("Shannon Diversity") +
  facet_wrap(. ~ Time_point2 + treatment_outcome, nrow = 1)

#ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf")

shannon_plot3 <- metadata %>%
  filter(Time_point == "day0" | Time_point == "day90",
         treatment_other_drug %in% "sbv",
         Swab_site %in% "lesion") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Time_point2 = case_when(
    Time_point == "day0" ~ "Day 0",
    Time_point == "day90" ~ "Day 90")) %>%
  ggplot(aes(x=Time_point2, y=shannon_entropy, color=Time_point2)) +
  geom_violin(adjust = .6) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[3],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="non") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="black") +
  stat_compare_means(method = "t.test", label.y = 6,
                     aes(label = ..p.signif..), paired = T,
                     label.x = 1.4, size = 6, color="black") +
  xlab("") +
  ylab("Shannon Index")

# Observed OTUs alpha diversity ----
obs_OTUs <- read_qza("Qiime2/core-metrics-results/observed_features_vector.qza")
obs_OTUs <- obs_OTUs$data %>% rownames_to_column("SampleID") 

metadata <- metadata %>% 
  left_join(obs_OTUs, by= "SampleID") 

# All sample day 0:
metadata %>%
  filter(Time_point %in% "day0") %>%
  ggplot(aes(x=Swab_site, y=observed_features, color=Swab_site)) +
  geom_violin(adjust = .6) +
  geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", label.y = 200,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Observed OTUs (alpha diversity)")





otu_plot <- metadata %>%
  filter(Time_point %in% "day0") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Swab_site2 = case_when(
    Swab_site == "lesion" ~ "Lesion",
    Swab_site == "contralateral" ~ "Contra.Skin")) %>%
  ggplot(aes(x=Swab_site2, y=observed_features, color=Swab_site2)) +
  geom_violin(adjust = .6) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="non") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="black") +
  stat_compare_means(method = "t.test", label.y = 280,
                     aes(label = ..p.signif..), paired = T,
                     label.x = 1.4, size = 6, color="black") +
  xlab("") +
  ylab("Number of observed OTUs")

otu_plot + shannon_plot

otu_plot2 <- metadata %>%
  filter(Time_point %in% c("day0","day30","day60","day90"),
         treatment_other_drug %in% "sbv") %>% # We don't have samples enough from other time points
  mutate(Swab_site2 = case_when(
    Swab_site == "lesion" ~ "Lesion",
    Swab_site == "contralateral" ~ "Contra.Skin")) %>%
  mutate(Time_point2 = case_when(
    Time_point == "day0" ~ "Day 0",
    Time_point == "day30" ~ "Day 30",
    Time_point == "day60" ~ "Day 60",
    Time_point == "day90" ~ "Day 90")) %>%
  ggplot(aes(x=Swab_site2, y=observed_features, color=Swab_site2)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none", strip.background = element_blank(), strip.text = element_text(size = 11)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", label.y = 280,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Number of observed OTUs") +
  facet_wrap(. ~ Time_point2, nrow = 1)

otu_plot2 / shannon_plot2

otu_plot3 <- metadata %>%
  filter(Time_point == "day0" | Time_point == "day90",
         treatment_other_drug %in% "sbv",
         Swab_site %in% "lesion") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Time_point2 = case_when(
    Time_point == "day0" ~ "Day 0",
    Time_point == "day90" ~ "Day 90")) %>%
ggplot(aes(x=Time_point2, y=observed_features, color=Time_point2)) +
  geom_violin(adjust = .6) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[3],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="non") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="black") +
  stat_compare_means(method = "t.test", label.y = 180,
                     aes(label = ..p.signif..), paired = T,
                     label.x = 1.4, size = 6, color="black") +
  xlab("") +
  ylab("Number of observed OTUs")

otu_plot3 + shannon_plot3

# All the patients with sbv:
metadata %>%
  filter(Time_point %in% c("day0","day30","day60","day90"),
         treatment_other_drug %in% "sbv") %>% # We don't have samples enough from other time points
  ggplot(aes(x=Swab_site, y=observed_features, color=Swab_site)) +
  geom_violin(adjust = .6) +
  geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", label.y = 7,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Observed OTUs (alpha diversity)") +
  facet_wrap(. ~ Time_point, nrow = 1)

#ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf")

metadata %>%
  filter(Time_point == "day0" | Time_point == "day90",
         treatment_other_drug %in% "sbv",
         Swab_site %in% "lesion",
         shannon_entropy != "NA") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  ggplot(., aes(x = Time_point, y = observed_features, group = LTCP_patient_ID,
                shape=treatment_outcome)) +
  geom_line(color=Dark24[1]) +
  geom_point(size=3) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  ylab("Observed OTUs (alpha diversity)") +
  xlab("")

# Alpha metrics vs. Clinical Outcome ----
metadata %>%
  filter(Time_point %in% "day0") %>%
  filter(treatment_other_drug == "sbv") %>%
  ggplot(aes(x=treatment_outcome, y=observed_features, color=treatment_outcome)) +
  geom_violin(adjust = .6) +
  geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", #label.y = 200,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Observed OTUs (alpha diversity)")

metadata %>%
  filter(Time_point %in% "day0") %>%
  filter(treatment_other_drug == "sbv") %>%
  ggplot(aes(x=treatment_outcome, y=shannon_entropy, color=treatment_outcome)) +
  geom_violin(adjust = .6) +
  geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", #label.y = 200,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Shannon Index (alpha diversity)")

metadata %>%
  filter(Time_point %in% "day0") %>%
  filter(treatment_other_drug == "sbv") %>%
  ggplot(aes(x=healing_time_days, y=shannon_entropy)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point() +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
  xlab("") +
  ylab("Shannon Index (alpha diversity)")

metadata %>%
  filter(Time_point %in% "day0") %>%
  filter(treatment_other_drug == "sbv") %>%
  ggplot(aes(x=healing_time_days, y=observed_features)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point() +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
  xlab("") +
  ylab("Observed OTUs (alpha diversity)")

# Specific dysbiosis samples time course ----
metadata %>%
  filter(LTCP_patient_ID %in% c("29271","29213","29256","29272","29513","29159","29231","29151")) %>% # Arcanobacterium
  filter(treatment_other_drug %in% "sbv",
         Swab_site %in% "lesion",
         shannon_entropy != "NA") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  ggplot(., aes(x = Time_point, y = observed_features, group = factor(LTCP_patient_ID,
                                                                      levels = c("29271","29213","29256","29272","29513","29159","29231","29151")),
                color= LTCP_patient_ID)) +
  geom_line(color="dark gray") +
  geom_point(size=3) +
  theme_classic() +
  scale_color_manual(values = rev(Dark24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  ylab("Observed OTUs (alpha diversity)") +
  xlab("")

metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", "")) %>%
  filter(treatment_other_drug %in% "sbv",
         Swab_site %in% "lesion",
         shannon_entropy != "NA") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  ggplot(., aes(x = Time_point, y = shannon_entropy),
                color= LTCP_patient_ID) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size=3) +
  theme_classic() +
  scale_color_manual(values = rev(Dark24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  ylab("Observed OTUs (alpha diversity)") +
  xlab("") +
  facet_wrap(. ~ treatment_outcome)

# PCoA (unweighted_unifrac) beta diversity ----
uwunifrac <- read_qza("Qiime2/core-metrics-results/unweighted_unifrac_pcoa_results.qza")

metadata <- metadata %>% 
  left_join(uwunifrac$data$Vectors %>%
              select(SampleID, PC1, PC2), by= "SampleID") %>%
  dplyr::rename(PC1_uwunifrac = PC1,
         PC2_uwunifrac = PC2)

metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion" | Swab_site == "contralateral") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  ggplot(., aes(x=PC1_uwunifrac, y=PC2_uwunifrac, color=Swab_site)) +
  geom_mark_hull(aes(label = Swab_site, fill = Swab_site), show.legend = FALSE, expand = unit(3, "mm")) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #stat_ellipse(level = 0.95) +
  theme_classic() +
  scale_color_manual(values = rev(Dark24[1:2])) + scale_fill_manual(values = rev(Dark24[1:2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(uwunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(uwunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%"))

# PCoA (weighted_unifrac) beta ----
wunifrac <- read_qza("Qiime2/core-metrics-results/weighted_unifrac_pcoa_results.qza")

metadata <- metadata %>% 
  left_join(wunifrac$data$Vectors %>%
              select(SampleID, PC1, PC2), by= "SampleID") %>%
  dplyr::rename(PC1_wunifrac = PC1,
         PC2_wunifrac = PC2)

wunifrac_plot <- metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion" | Swab_site == "contralateral") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, color=Swab_site)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = rev(Dark24[1:2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  #coord_fixed() +
  NULL

#or:
wunifrac_plot <- metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion" | Swab_site == "contralateral") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, color=Swab_site)) +
  geom_mark_hull(aes(label = Swab_site, fill = Swab_site), show.legend = FALSE, expand = unit(3, "mm")) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = rev(Dark24[1:2])) +
  scale_fill_manual(values = rev(Dark24[1:2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

wunifrac_plot2 <-  metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()
  
  metadata %>%
    filter(Swab_site == "lesion",
           treatment_other_drug %in% "sbv",
           Time_point %in% c("day0","day30","day60","day90")) %>%
    arrange(LTCP_patient_ID) %>%
    group_by(LTCP_patient_ID) %>%
    filter(n() > 2) %>%
    ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, color=Time_point)) +
    #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
    geom_point(size = 4) +
    #stat_ellipse(level = 0.95) +
    theme_classic() + scale_color_manual(values = c("#FA0011",
                                                    "#CA004A",
                                                    "#4E2A63",
                                                    "#0D2A63")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
          legend.position="none", strip.background = element_blank(),
          panel.border = element_rect(fill=NA, size=1)) +
    xlab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
    ylab(paste("PCoA, axis 2 -",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
    facet_wrap(. ~ treatment_outcome) +
  coord_fixed()
  
  
  metadata %>%
    filter(Swab_site == "lesion",
           treatment_other_drug %in% "sbv",
           Time_point %in% c("day0","day30","day60","day90")) %>%
    arrange(LTCP_patient_ID) %>%
    group_by(LTCP_patient_ID) %>%
    filter(n() > 2) %>%
  ggplot(., aes(x=PC1_uwunifrac, fill=Time_point)) +
    #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
    geom_density(stat = "density", adjust = 1) +
    #stat_ellipse(level = 0.95) +
    theme_classic() + scale_fill_manual(values = c("#FA0011",
                                                   "#CA004A",
                                                   "#4E2A63",
                                                   "#0D2A63")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
          legend.position="none", strip.background = element_blank(),
          panel.border = element_rect(fill=NA, size=1)) +
    xlab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
    facet_wrap(. ~ treatment_outcome)
  
  metadata %>%
    #filter(Swab_site == "lesion",
    #       treatment_other_drug %in% "sbv",
    #       Time_point %in% c("day0","day30","day60","day90")) %>%
    
    filter(Swab_site == "lesion" &
             treatment_other_drug %in% "sbv" &
             Time_point %in% c("day0","day30","day60",
                               "day90") | 
             Swab_site == "contralateral" &
             treatment_other_drug %in% "sbv" &
             Time_point == "day90") %>%
    mutate(Time_point3 = case_when(
      Time_point == "day90" & Swab_site == "contralateral" ~ "Contra.Skin D90",
      Time_point == "day0" & Swab_site == "lesion" ~ "D0",
      Time_point == "day30" & Swab_site == "lesion" ~ "D30",
      Time_point == "day60" & Swab_site == "lesion" ~ "D60",
      Time_point == "day90" & Swab_site == "lesion" ~ "D90")) %>%
    mutate(treatment_outcome2 = case_when(
      treatment_outcome == "cure" ~ "Cured",
      treatment_outcome == "failure" ~ "Failure")) %>%
    arrange(LTCP_patient_ID) %>%
    group_by(LTCP_patient_ID) %>%
    filter(n() > 2) %>%
    ggplot(., aes(y=PC1_wunifrac,
                  x=factor(Time_point3, levels = c("D0","D30","D60", "D90","Contra.Skin D90")),
                  #x= Time_point,
                  color=Time_point3)) +
    geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
    geom_violin(size = 0.7, trim = F) +
    geom_sina(size=2.3) +
    theme_classic() + scale_color_manual(values = c("black",
                                                    "#FA0011",
                                                    "#CA004A",
                                                    "#4E2A63",
                                                    "#0D2A63")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
          legend.position="right", strip.background = element_blank(), strip.text = element_text(size=11),
          axis.ticks.x = element_blank()) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.7, color="black") +
    stat_compare_means(#comparisons = my_comparisonsAA,
      comparisons = my_comparisonsA,
      method = "t.test", label.y = 1,
      aes(label = ..p.signif..),
      #aes(label = ..p.value..),
      #label.x = 1.5,
      size =5, color="black") +
    ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
    xlab("") + ylim(-0.8,1.2) +
    facet_wrap(. ~ treatment_outcome2)

# PCoA (Jaccard) beta ----
jaccard <- read_qza("Qiime2/core-metrics-results/jaccard_pcoa_results.qza")

metadata <- metadata %>% 
  left_join(jaccard$data$Vectors %>%
              select(SampleID, PC1, PC2), by= "SampleID") %>%
  dplyr::rename(PC1_jaccard = PC1,
         PC2_jaccard = PC2)

jaccard_plot <- metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion" | Swab_site == "contralateral") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  ggplot(., aes(x=PC1_jaccard, y=PC2_jaccard, color=Swab_site)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = rev(Dark24[1:2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1 -",round(jaccard[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(jaccard[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion" | Swab_site == "contralateral") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  ggplot(., aes(x=PC1_jaccard, y=PC2_jaccard, color=Swab_site)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = rev(Dark24[1:2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1",jaccard[["data"]][["Eigvals"]][["PC1"]],"%")) +
  ylab(paste("PCoA, axis 2",jaccard[["data"]][["Eigvals"]][["PC2"]],"%"))

jaccard_plot2 <- metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_jaccard, y=PC2_jaccard, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1 -",round(jaccard[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(jaccard[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

# PCoA (Bray-Curtis) beta ----
bray <- read_qza("Qiime2/core-metrics-results/bray_curtis_pcoa_results.qza")

metadata <- metadata %>% 
  left_join(bray$data$Vectors %>%
              select(SampleID, PC1, PC2), by= "SampleID") %>%
  dplyr::rename(PC1_bray = PC1,
         PC2_bray = PC2)

bray_plot <- metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion" | Swab_site == "contralateral") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Swab_site)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point() +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = rev(Dark24[1:2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

jaccard_plot + bray_plot + wunifrac_plot

bray_plot2 <- metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none", strip.background = element_blank(),
        panel.border = element_rect(fill=NA, size=1)) +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  facet_wrap(. ~ treatment_outcome) +
  coord_fixed()

metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_bray, fill=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_density(stat = "density", adjust = 1) +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_fill_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right", strip.background = element_blank(),
        panel.border = element_rect(fill=NA, size=1)) +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  facet_wrap(. ~ treatment_outcome)

my_comparisonsA <- list(c("D0", "D90"))
my_comparisonsAA <- list(c("D0", "D90"),
                         c("D90","Contra.Skin D90"))
metadata %>%
  #filter(Swab_site == "lesion",
  #       treatment_other_drug %in% "sbv",
  #       Time_point %in% c("day0","day30","day60","day90")) %>%
  
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point == "day90") %>%
  mutate(Time_point3 = case_when(
    Time_point == "day90" & Swab_site == "contralateral" ~ "Contra.Skin D90",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90")) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(y=PC1_bray,
                x=factor(Time_point3, levels = c("D0","D30","D60", "D90","Contra.Skin D90")),
                #x= Time_point,
                color=Time_point3)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  theme_classic() + scale_color_manual(values = c("black",
                                                  "#FA0011",
                                                 "#CA004A",
                                                 "#4E2A63",
                                                 "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right", strip.background = element_blank(), strip.text = element_text(size=11),
        axis.ticks.x = element_blank()) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(#comparisons = my_comparisonsAA,
                     comparisons = my_comparisonsA,
                     method = "t.test", label.y = 1,
                     aes(label = ..p.signif..),
                     #aes(label = ..p.value..),
                     #label.x = 1.5,
                     size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + ylim(-0.8,1.2) +
  facet_wrap(. ~ treatment_outcome2)


jaccard_plot2 + bray_plot2 + wunifrac_plot2

metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(y=PC1_bray, x=Time_point, color=Time_point)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
        legend.position="right", strip.background = element_blank()) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(comparisons = my_comparisonsA,
                     method = "wilcox.test", label.y = 1,
                     aes(label = ..p.signif..),
                     #aes(label = ..p.value..),
                     #label.x = 1.5,
                     size =6, color="black") +
  ylab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + ylim(-0.7,1.2) +
  facet_wrap(. ~ treatment_outcome)

# Compare all samples from contralateral skin (Evaluating the impact of therapy in general microbiome) ----
metadata %>%
  filter(Swab_site == "contralateral",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed() +

metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()


# Contra:
metadata %>%
  #filter(LTCP_patient_ID == "29663") %>%
  filter(Swab_site == "contralateral",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 2) +
  #geom_segment(size=1,
  #             aes(color=Time_point,
  #               xend=c(tail(PC1_bray, n=-1), NA), 
  #               yend=c(tail(PC2_bray, n=-1), NA)),
  #             arrow=arrow(length=unit(0.2,"cm"))) +
  #geom_line(size=1,color="#1CA71C",
  #          arrow=arrow(length=unit(0.2,"cm")),
  #          aes(group=LTCP_patient_ID)) +
  geom_path(size=1,#color="#1CA71C",
            arrow=arrow(length=unit(0.3,"cm")),
            aes(x=PC1_bray, y=PC2_bray, group=LTCP_patient_ID, color=Time_point)) +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  geom_text_repel(aes(label = Time_point, color=Time_point), size = 3, fontface="bold") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), strip.background = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none", panel.border = element_rect(fill=NA)) +

  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  xlim(-0.35,0.50) + ylim(-0.30,0.40) +
  #coord_fixed() +
  facet_wrap(~ LTCP_patient_ID)

# Check other metrics?!!?
# Add a lesion for comparison!?
# See gepm_pointline

#Lesion:
metadata %>%
  filter(LTCP_patient_ID == "29663") %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Time_point)) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 2) +
  #geom_segment(size=1,
  #             aes(color=Time_point,
  #               xend=c(tail(PC1_bray, n=-1), NA), 
  #               yend=c(tail(PC2_bray, n=-1), NA)),
  #             arrow=arrow(length=unit(0.2,"cm"))) +
  #geom_line(size=1,color="#1CA71C",
  #          arrow=arrow(length=unit(0.2,"cm")),
  #          aes(group=LTCP_patient_ID)) +
  geom_path(size=1,#color="#1CA71C",
            arrow=arrow(length=unit(0.3,"cm")),
            aes(x=PC1_bray, y=PC2_bray, group=LTCP_patient_ID, color=Time_point)) +
  #stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  "#CA004A",
                                                  "#4E2A63",
                                                  "#0D2A63")) +
  geom_text_repel(aes(label = Time_point, color=Time_point), size = 3, fontface="bold") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), strip.background = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none", panel.border = element_rect(fill=NA)) +
  
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  xlim(-0.35,0.50) + ylim(-0.30,0.40) +
  #coord_fixed() +
  facet_wrap(~ LTCP_patient_ID, scales = "free")

#Only Staphy dysbiosis patients (weighted unifrac):
metadata %>%
  filter(microbiome_cluster == "Staphylococcus",
         Swab_site == "contralateral",
         treatment_other_drug %in% "sbv",
         #Time_point %in% c("day0","day30","day60","day90")) %>%
         Time_point %in% c("day0","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
 # filter(n() > 2) %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, color=Time_point)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  #"#CA004A",
                                                  #"#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  #facet_wrap(~treatment_outcome) +
  coord_fixed() +
  
  metadata %>%
  filter(microbiome_cluster == "Staphylococcus",
         Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         #Time_point %in% c("day0","day30","day60","day90")) %>%
         Time_point %in% c("day0","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, color=Time_point)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  theme_classic() + scale_color_manual(values = c("#FA0011",
                                                  #"#CA004A",
                                                  #"#4E2A63",
                                                  "#0D2A63")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  #facet_wrap(~treatment_outcome) +
  coord_fixed()

# Revisit PCoAs with Microbiome clusters information (created below) ----
metadata <- metadata %>%
  left_join(phenotype_pre %>%
              select(LTCP_patient_ID, microbiome_cluster))
metadata$microbiome_cluster <- factor(metadata$microbiome_cluster,
                                      levels = c("Staphylococcus",
                                                 "Arcanobacterium",
                                                 "Corynebacterium",
                                                 "Streptococcus",
                                                 "Staphylococcus_Streptococcus",
                                                 "Heterogeneous",
                                                 "Lactobacillus",
                                                 "Finegoldia"))


# Weighted Unifrac:
metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion") %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, fill=microbiome_cluster)) +
  geom_point(size=4, shape=21) +
  geom_point(size=4) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  theme_classic() + scale_fill_manual(values = c(Dark24[4],
                                                  Dark24[2],
                                                  "gray",
                                                  Dark24[12],
                                                  Dark24[1],
                                                  Dark24[3],
                                                 Dark24[6],
                                                 Dark24[7])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion") %>%
  arrange(microbiome_cluster) %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac,
                #color=microbiome_cluster,
                NULL)) +
  #geom_hline(yintercept = 0, color="purple", linetype = 2) +
  #geom_vline(xintercept = 0, color="purple", linetype = 2) +
  geom_mark_hull(aes(fill = microbiome_cluster, label = microbiome_cluster), alpha = 0.1, show.legend = F, expand = unit(2.5, "mm")) +
  geom_point(size=3) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  theme_classic() + 
  #scale_fill_manual(values = c(Dark24[1],"dark gray",Dark24[4],Dark24[2],Dark24[3],Dark24[5],Dark24[6])) +
  #scale_color_manual(values = c(Dark24[1],"dark gray",Dark24[4],Dark24[2],Dark24[3],Dark24[5],Dark24[6])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), axis.title = element_text(size = 15),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

metadata %>%
  filter(Time_point == "day0",
         LTCP_patient_ID %in% Arca_samples_pre) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  ggplot(., aes(x=PC1_wunifrac, y=PC2_wunifrac, fill=Swab_site)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size=4, shape=21) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  theme_classic() + scale_fill_manual(values = c(Dark24[7],
                                                 Dark24[4])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

# Unweighted Unifrac
metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion") %>%
  ggplot(., aes(x=PC1_uwunifrac, y=PC2_uwunifrac, fill=microbiome_cluster)) +
  geom_point(size=4, shape=21) +
  theme_classic() + scale_fill_manual(values = c(Dark24[4],
                                                 Dark24[2],
                                                 "gray",
                                                 Dark24[12],
                                                 Dark24[1],
                                                 Dark24[3])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(uwunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(uwunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%"))

# Jaccard
metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion") %>%
  ggplot(., aes(x=PC1_jaccard, y=PC2_jaccard, fill=microbiome_cluster)) +
  geom_point(size=4, shape=21) +
  theme_classic() + scale_fill_manual(values = c(Dark24[4],
                                                 Dark24[2],
                                                 "gray",
                                                 Dark24[12],
                                                 Dark24[1],
                                                 Dark24[3])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(jaccard[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(jaccard[["data"]][["Eigvals"]][["PC2"]],2),"%"))

# Jaccard:
metadata %>%
  filter(Time_point == "day0",
         LTCP_patient_ID %in% Arca_samples_pre) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  ggplot(., aes(x=PC1_jaccard, y=PC2_jaccard, fill=Swab_site)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size=4, shape=21) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  theme_classic() + scale_fill_manual(values = c(Dark24[7],
                                                 Dark24[4])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(jaccard[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(jaccard[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

# Bray-Curtis:
metadata %>%
  filter(Time_point == "day0",
         LTCP_patient_ID %in% Arca_samples_pre) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, fill=Swab_site)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  geom_point(size=4, shape=21) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  theme_classic() + scale_fill_manual(values = c(Dark24[7],
                                                 Dark24[4])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

metadata %>%
  filter(Time_point == "day0",
         Swab_site == "lesion",
         LTCP_patient_ID %in% Arca_samples_pre) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, fill=Swab_site)) +
  geom_point(size=4, shape=21) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 5, fontface="bold",colour="black") +
  theme_classic() + scale_fill_manual(values = c(Dark24[7],
                                                 Dark24[4])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  xlab(paste("PCoA, axis 1",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  coord_fixed()

# Phyloseq ----
GlobalPatterns <- data("GlobalPatterns") # for tutorial

# Import artifacts and outputs from Qiime2 into R and converting to phyloseq object:
physeq <- qza_to_phyloseq(
  features = "Qiime2/pet-table.qza",
  tree = "Qiime2/rooted-tree.qza",
  taxonomy = "Qiime2/taxonomy.qza",
  metadata = "Qiime2/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples_wTypesCol.tsv"
)


# Exploring the phyloseq object for filtering ----
physeq
ntaxa(physeq)
nsamples(physeq)
sample_names(physeq)[1:5]
rank_names(physeq)
sample_variables(physeq)
otu_table(physeq)[1:5, 1:5]
tax_table(physeq)[1:5, 1:4]
phy_tree(physeq)

# Transform to relative abundance and filter OTUs based on the mean among samples ----
# Kinda followed https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html
prevelancedf <- apply(X = otu_table(physeq),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevelancedf <- data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(physeq),
                          tax_table(physeq))
prevelancedf[1:10,]

plyr::ddply(prevelancedf, "Kingdom", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),
             stringsAsFactors = F)
})

plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),
             stringsAsFactors = F)
})

plyr::ddply(prevelancedf, "Genus", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),
             stringsAsFactors = F)
}) %>% arrange(Genus)

View(plyr::ddply(prevelancedf, "Species", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),
             stringsAsFactors = F)
}))

#Whole Kingdom/phylum filtering
#First lets remove of the feature with ambiguous Kingdom/phylum annotation.
physeq
physeq_1 <- subset_taxa(physeq, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned"))
physeq_1

# Weird ones:
level_filters <- plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),
             stringsAsFactors = F)
}) %>%
  arrange(total_abundance) %>%
  filter(total_abundance < 100)
level_filters <- level_filters$Phylum

physeq_1 <- subset_taxa(physeq_1, !is.na(Phylum) & !Phylum %in% c("", level_filters))
physeq_1

#Individual Taxa Filtering
#Subset to the remaining phyla by prevalence.

prevelancedf1 <- subset(prevelancedf, Phylum %in% get_taxa_unique(physeq_1, taxonomic.rank = "Phylum"))
ggplot(prevelancedf1, aes(TotalAbundance, Prevalence / nsamples(physeq_1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, color="purple", linetype = 2) +
  geom_hline(yintercept = 0.05, color="green", linetype = 2) +
  geom_hline(yintercept = 0.10, color="red", linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- 0.05 * nsamples(physeq_1)
prevalenceThreshold

# Execute prevalence filter
keepTaxa <- rownames(prevelancedf1)[(prevelancedf1$Prevalence >= prevalenceThreshold)]
length(keepTaxa)

physeq_2 <- prune_taxa(keepTaxa, physeq_1)
physeq_2

#Agglomerate taxa at the Genus level (combine all with the same name) and remove all taxa without genus level assignment
length(get_taxa_unique(physeq_2, taxonomic.rank = "Genus"))
physeq_3 <- tax_glom(physeq_2, "Genus", NArm = TRUE)
physeq_3

## out of curiosity how many "reads" does this leave us at???
sum(colSums(otu_table(physeq_3)))

# Look at some trees ----
tree_frequency <- plot_tree(physeq_3, color = "Swab_site", label.tips = "Genus",
          #size = "Abundance",
          ladderize = "left", justify = "left") +
  scale_color_manual(values = Dark24) +
  theme(axis.text = element_blank())
tree_frequency
plot_bar(physeq_3, "Swab_site", "Abundance", "Family")

# Lesion and Contralateral from patients at day 0 - barplots ----
## LESION
cha1 <- subset_samples(physeq_3, Time_point=="day0")
cha1 <- subset_samples(cha1, Swab_site=="lesion")
cha1_names <- names(sort(taxa_sums(cha1), decreasing = TRUE)[1:10])
cha1_pre <- prune_taxa(cha1_names, cha1)
cha1 <- transform_sample_counts(cha1_pre, function(x) x / sum(x)) # Relative abundance

# Original plot, not ordering:
cha1_plot <- plot_bar(cha1, x="SubjectID", y="Abundance", fill="Genus") +
  scale_fill_manual(values = Dark24)

# Clustering samples from top taxa:
cha1_df <- cha1@otu_table@.Data
toptaxa <- as.data.frame(cha1@tax_table@.Data)
toptaxa <- toptaxa$Genus
rownames(cha1_df) <- toptaxa

# Just modify the arrangement of samples:
toptaxa_names <- cha1_plot[["data"]] %>%
  select(Genus,OTU) %>%
  filter(OTU %in% cha1_names) %>% 
  distinct(OTU, .keep_all = TRUE)
toptaxa_names <- toptaxa_names[match(cha1_names, toptaxa_names$OTU),] 
toptaxa_names <- toptaxa_names$Genus

## by dendrogram:
distance_che <- dist(t(cha1_df), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters_che <- hclust(distance_che, method = "ward.D2") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters_che)

save(clusters_che, file = "outputs/clusters_che")

cha1_plot$data$Genus <- factor(cha1_plot$data$Genus, levels = toptaxa_names)
cha1_plot + theme(axis.text.x = element_text(angle = 90))

## CONTRALATERAL
cha2 <- subset_samples(physeq_3, Time_point == "day0")
cha2<- subset_samples(cha2, Swab_site == "contralateral")
cha2<- subset_samples(cha2, SubjectID %in% levels(cha1_plot[["data"]][["SubjectID"]]))
cha2_names <- names(sort(taxa_sums(cha2), decreasing = TRUE)[1:10])
cha2 <- prune_taxa(cha2_names, cha2)
cha2 <- transform_sample_counts(cha2, function(x) x / sum(x)) # Relative abundance

cha2_plot <- plot_bar(cha2, x="SubjectID", y="Abundance", fill="Genus")

toptaxa_names2 <- cha2_plot[["data"]] %>%
  select(Genus,OTU) %>%
  filter(OTU %in% cha2_names) %>% 
  distinct(OTU, .keep_all = TRUE)
toptaxa_names2 <- toptaxa_names2[match(cha2_names, toptaxa_names2$OTU),] 
toptaxa_names2 <- toptaxa_names2$Genus

cha2_plot$data$Genus <- factor(cha2_plot$data$Genus, levels = toptaxa_names2)

cha1_plot /
cha2_plot
ggsave("outputs/16SRA.pdf",
       width = 10, height = 6)

cha1_plot +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15))

cha1_plot /
  cha2_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15))

# Lesion and Contralateral from patients at day 0 - heatmaps ----
myheatcol <- colorRampPalette(colors=c("white",Dark24[8]))(100)

# Microbiome modules from top taxa:
micro_module <- cutree(clusters_che, k=8)
micro_module_fac <- factor(micro_module)

#micro_colors <- function(micro_module_fac){ if (micro_module_fac==1) Dark24[3] else if (micro_module_fac==2) "gray" else if (micro_module_fac==3) Dark24[1] else if (micro_module_fac==4) Dark24[4] else if (micro_module_fac==5) Dark24[2] else if (micro_module_fac==7) Dark24[6] else Dark24[5]}
#micro_colors <- unlist(lapply(micro_module_fac, micro_colors))

micro_colors <- function(micro_module_fac){ if (micro_module_fac==1) Dark24[3] else if (micro_module_fac==2) "gray" else if (micro_module_fac==3) Dark24[12] else if (micro_module_fac==4) Dark24[1] else if (micro_module_fac==5) Dark24[4] else if (micro_module_fac==6) Dark24[2] else if (micro_module_fac==7) Dark24[5] else Dark24[6]}
micro_colors <- unlist(lapply(micro_module_fac, micro_colors))

save(micro_colors, file = "outputs/micro_colors")

microbiome_clusters <- heatmap.2(cha1_df[toptaxa_names,],
          Rowv = F,
          Colv=as.dendrogram(clusters_che),
          dendrogram = "column",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = micro_colors,
          #scale='row',
          labCol = NA,
          density.info="none", trace="none",
          cexRow=1.6, cexCol=0.5,margins = c(7,15))

heatmap.2(cha1_df[toptaxa_names,],
          Rowv = F,
          Colv=as.dendrogram(clusters_che),
          dendrogram = "column",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = micro_colors,
          #scale='row',
          labCol = NA,
          density.info="none", trace="none",
          cexRow=1, cexCol=0.5,margins = c(16,10))

heatmap.2(cha1_df[toptaxa_names,],
          Rowv = F,
          Colv=as.dendrogram(clusters_che),
          dendrogram = "column",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = micro_colors,
          #scale='row',
          density.info="none", trace="none",
          cexRow=1, cexCol=1,margins = c(20,15))



# Extract the clustering (sample ordering) from the heatmap:
colnames(cha1_df[toptaxa_names,])[microbiome_clusters$colInd]
cluster_heat_order <- as.data.frame(colnames(cha1_df[toptaxa_names,])[microbiome_clusters$colInd])
colnames(cluster_heat_order)[1] <- "SampleID"
cluster_heat_order <- cluster_heat_order %>%
  left_join(metadata %>% select(SampleID, LTCP_patient_ID), by="SampleID")
patient_order2 <- cluster_heat_order$LTCP_patient_ID

save(patient_order2, file = "outputs/patient_order2")
cha1_plot$data$SubjectID <- factor(cha1_plot$data$SubjectID, levels = patient_order2)

cha2_plot$data$SubjectID <- factor(cha2_plot$data$SubjectID, levels = patient_order2)

#Export matrix to calculate RAmean:
exportRA <- rownames_to_column(as.data.frame(t(cha1_df[toptaxa_names,cluster_heat_order$SampleID])), "SampleID")
write_tsv(exportRA, "../../../../Desktop/matrix.txt")

exportRA2 <- rownames_to_column(as.data.frame(cha1_df[toptaxa_names,cluster_heat_order$SampleID]), "SampleID")
write_tsv(exportRA2, "../../../../Desktop/matrix2.txt")


# Bar plots again AND OFFICIAL:
cha1_plot <- cha1_plot + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 11),
        legend.text = element_text(size=9),legend.key.size = unit(0.5, 'cm'))

cha2_plot <- cha2_plot + scale_fill_manual(values = c(Dark24[1], Dark24[9],
                                           Dark24[2], Dark24[13], 
                                           Dark24[19], Dark24[5],
                                           Dark24[16],Dark24[15],
                                           Dark24[10], Dark24[17])) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 11),
        legend.text = element_text(size=9), legend.key.size = unit(0.5, 'cm'), legend.title = element_blank())

cha1_plot / cha2_plot

ggsave("outputs/16SRA.pdf",
       width = 10, height = 6)

# Pie chart of number of samples per group:
micro_module_df <- rownames_to_column(as.data.frame(micro_module), "SampleID")
micro_module_df <- as.data.frame(table(factor(micro_module_df$micro_module)))
micro_module_df %>%
  mutate(Percentage=paste0(round(Freq/sum(Freq)*100,0),"%")) %>%
  ggplot(., aes(x="", y=Percentage, fill=factor(Var1))) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values = Dark24) +
  coord_polar("y", start=0) +
  geom_text(aes(label = Percentage), color = "white", size=6,
            position = position_stack(vjust = 0.5),
            x=1) 

micro_module_df %>% # Colors messed up
  mutate(Percentage=paste0(round(Freq/sum(Freq)*100,0),"%")) %>%
ggplot(., aes (x="", y = Freq, fill = factor(Var1))) + 
  theme_void() +
  geom_col(position = 'stack', width = 1) +
  geom_text(aes(label = Freq, x = 1.3),
            position = position_stack(vjust = 0.5),
            size=6, color= "white", fontface=2) +
  scale_fill_manual(values = c(Dark24[3],
                               "gray",
                               Dark24[1],
                               Dark24[4],
                               Dark24[2],
                               Dark24[5],
                               Dark24[6],
                               Dark24[12])) +
  theme(legend.position = "none") +
  coord_polar("y")

micro_module_df %>%
  mutate(Percentage=paste0(round(Freq/sum(Freq)*100,0),"%")) %>%
  ggplot(., aes (x="", y = Freq, fill = factor(Var1))) + 
  theme_void() +
  geom_col(position = 'stack', width = 1) +
  geom_text(aes(label = Percentage, x = 1.3),
            position = position_stack(vjust = 0.5),
            size=6, color= "white", fontface=2) +
  scale_fill_manual(values = c(Dark24[3],
                               "gray",
                               Dark24[1],
                               Dark24[4],
                               Dark24[2],
                               Dark24[5],
                               Dark24[6])) +
  theme(legend.position = "none") +
  coord_polar("y")

# Export this microbiome data reduction (cluster and patients) for Integrative analysis ----
# I'm sending this Robject to "../Dataset_Integration" for then export to rexposome (CHMI server)
phenotype_pre <- rownames_to_column(as.data.frame(micro_module), "SampleID")
phenotype_pre <- metadata %>%
  select(SampleID, LTCP_patient_ID,
         shannon_entropy, observed_features) %>% # Diversity metrics
  left_join(phenotype_pre) %>%
  filter(micro_module != "NA") %>%
  mutate(microbiome_cluster = case_when(
    micro_module == 3 ~ "Staphylococcus_Streptococcus",
    micro_module == 2 ~ "Heterogeneous",
    micro_module == 4 ~ "Staphylococcus",
    micro_module == 1 ~ "Streptococcus",
    micro_module == 5 ~ "Arcanobacterium",
    micro_module == 6 ~ "Corynebacterium",
    micro_module == 7 ~ "Lactobacillus",
    micro_module == 8 ~ "Finegoldia")) %>%
  select(-SampleID, -micro_module)

save(phenotype_pre, file = "outputs/phenotype_pre")

# Go to /Users/amorimc/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/Dataset_Integration for phenotype pure.

# Name of the samples from the microbiome clusters: ----
Arca_samples <- phenotype_pre %>% filter(microbiome_cluster == "Arcanobacterium")
Arca_samples_pre <- Arca_samples$LTCP_patient_ID
Arca_samples <- str_replace_all(Arca_samples$LTCP_patient_ID, "^29", "CL29")

Coryne_samples <- phenotype_pre %>% filter(microbiome_cluster == "Corynebacterium")
Coryne_samples <- str_replace_all(Coryne_samples$LTCP_patient_ID, "^29", "CL29")
  
Hetero_samples <- phenotype_pre %>% filter(microbiome_cluster == "Heterogeneous")
Hetero_samples <- str_replace_all(Hetero_samples$LTCP_patient_ID, "^29", "CL29")

Staphy_samples <- phenotype_pre %>% filter(microbiome_cluster == "Staphylococcus")
Staphy_samples <- str_replace_all(Staphy_samples$LTCP_patient_ID, "^29", "CL29")

Strep_samples <- phenotype_pre %>% filter(microbiome_cluster == "Streptococcus")
Strep_samples <- str_replace_all(Strep_samples$LTCP_patient_ID, "^29", "CL29")

StaphyStrep_samples <- phenotype_pre %>% filter(microbiome_cluster == "Staphylococcus_Streptococcus")
StaphyStrep_samples <- str_replace_all(StaphyStrep_samples$LTCP_patient_ID, "^29", "CL29")

save(Arca_samples, file = "outputs/Arca_samples")
save(Coryne_samples, file = "outputs/Coryne_samples")
save(Hetero_samples, file = "outputs/Hetero_samples")
save(Staphy_samples, file = "outputs/Staphy_samples")
save(Strep_samples, file = "outputs/Strep_samples")
save(StaphyStrep_samples, file = "outputs/StaphyStrep_samples")

# Pie charts with the RAmean for each microbiome cluster and contralateral skin ----
# Lesions:
cha1_plot[["data"]] %>%
  mutate(SampleID = Sample) %>%
  select(SampleID, Abundance, Genus) %>%
  left_join(metadata) %>%
  select(SampleID, Abundance, Genus, LTCP_patient_ID, microbiome_cluster) %>%
  group_by(microbiome_cluster, Genus) %>%
  summarise(RAmean = mean(Abundance)) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Staphylococcus_Streptococcus" ~ "M3",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1")) %>%
ggplot(., aes(x="", y=RAmean, fill=Genus)) +
  geom_bar(stat="identity", width=1, color="white", size=0.3) +
  scale_fill_manual(values = Dark24) +
  coord_polar("y", start=0) +
  facet_wrap(. ~ microbiome_cluster2, nrow = 1) +
  guides(size="none") +
  theme_void() +
  theme(strip.text = element_text(size=17))

# Contralateral skin:
cha2_plot[["data"]] %>%
  mutate(SampleID = Sample) %>%
  select(SampleID, Abundance, Genus) %>%
  left_join(metadata) %>%
  select(SampleID, Abundance, Genus, LTCP_patient_ID, microbiome_cluster) %>%
  group_by(microbiome_cluster, Genus) %>%
  summarise(RAmean = mean(Abundance)) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Staphylococcus_Streptococcus" ~ "M3",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1")) %>%
  ggplot(., aes(x="", y=RAmean, fill=Genus)) +
  geom_bar(stat="identity", width=1, color="white", size=0.3) +
  scale_fill_manual(values = c(Dark24[1], Dark24[9],
                                 Dark24[2], Dark24[13], 
                                 Dark24[19], Dark24[5],
                                 Dark24[16],Dark24[15],
                                 Dark24[10], Dark24[17])) +
  coord_polar("y", start=0) +
  facet_wrap(. ~ microbiome_cluster2, nrow = 1) +
  guides(size="none") +
  theme_void() +
  theme(strip.text = element_text(size=17))

# Pie Chart of specific patients overtime ----
cha3 <- subset_samples(physeq_3, Swab_site=="lesion")
cha3_names <- names(sort(taxa_sums(cha3), decreasing = TRUE)[1:10])
cha3_pre <- prune_taxa(cha1_names, cha3)
cha3 <- transform_sample_counts(cha3_pre, function(x) x / sum(x)) # Relative abundance

# Original plot, not ordering:
cha3_plot <- plot_bar(cha3, x="SubjectID", y="Abundance", fill="Genus") +
  scale_fill_manual(values = Dark24) +
  facet_wrap(. ~ Time_point)
cha3_plot$data$Genus <- factor(cha3_plot$data$Genus, levels = toptaxa_names)
cha3_plot


# 29414, patient 55:
cha3_plot[["data"]] %>%
  mutate(SampleID = Sample) %>%
  select(SampleID, Abundance, Genus, Time_point) %>%
  left_join(metadata) %>%
  filter(LTCP_patient_ID == "29414") %>%
  select(SampleID, Abundance, Genus, LTCP_patient_ID, microbiome_cluster, Time_point) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day0" ~ "D0",
    Time_point == "day30" ~ "D30",
    Time_point == "day60" ~ "D60",
    Time_point == "day90" ~ "D90",
    Time_point == "day120" ~ "D120",
    Time_point == "day150" ~ "D150",
    Time_point == "day180" ~ "D180",
    Time_point == "day210" ~ "D210",
    Time_point == "day240" ~ "D240",
    Time_point == "day270" ~ "D270")) %>%
ggplot(., aes(x="", y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", width=1, color="white", size=0.3) +
  scale_fill_manual(values = Dark24) +
  coord_polar("y", start=0) +
  guides(size="none") +
  facet_wrap(. ~ Time_point3, nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size=17))

cha3_plot[["data"]] %>%
  mutate(SampleID = Sample) %>%
  select(SampleID, Abundance, Genus, Time_point) %>%
  left_join(metadata) %>%
  filter(LTCP_patient_ID == "29412") %>%
  select(SampleID, Abundance, Genus, LTCP_patient_ID, microbiome_cluster, Time_point) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day0" ~ "D0",
    Time_point == "day30" ~ "D30",
    Time_point == "day60" ~ "D60",
    Time_point == "day90" ~ "D90",
    Time_point == "day120" ~ "D120",
    Time_point == "day150" ~ "D150",
    Time_point == "day180" ~ "D180",
    Time_point == "day210" ~ "D210",
    Time_point == "day240" ~ "D240",
    Time_point == "day270" ~ "D270")) %>%
  ggplot(., aes(x="", y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", width=1, color="white", size=0.3) +
  scale_fill_manual(values = Dark24) +
  coord_polar("y", start=0) +
  guides(size="none") +
  facet_wrap(. ~ factor(Time_point3, levels = c("D0", "D30","D60","D90", "D150")), nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size=17))

cha3_plot[["data"]] %>%
  mutate(SampleID = Sample) %>%
  select(SampleID, Abundance, Genus, Time_point) %>%
  left_join(metadata) %>%
  filter(LTCP_patient_ID == "29272") %>%
  select(SampleID, Abundance, Genus, LTCP_patient_ID, microbiome_cluster, Time_point) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day0" ~ "D0",
    Time_point == "day30" ~ "D30",
    Time_point == "day60" ~ "D60",
    Time_point == "day90" ~ "D90",
    Time_point == "day120" ~ "D120",
    Time_point == "day150" ~ "D150",
    Time_point == "day180" ~ "D180",
    Time_point == "day210" ~ "D210",
    Time_point == "day240" ~ "D240",
    Time_point == "day270" ~ "D270")) %>%
  ggplot(., aes(x="", y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", width=1, color="white", size=0.3) +
  scale_fill_manual(values = Dark24) +
  coord_polar("y", start=0) +
  guides(size="none") +
  facet_wrap(. ~ factor(Time_point3, levels = c("D0", "D30","D60","D90", "D120")), nrow = 1) +
  #facet_wrap(. ~ Time_point3, nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size=17))


# Pie Chart of all patients overtime ----
cha4 <- subset_samples(physeq_3, Swab_site %in% c("lesion") &
                         SubjectID == "29529")
                         #SubjectID %in% str_replace(Arca_samples,"CL",""))

cha4_names <- names(sort(taxa_sums(cha4), decreasing = TRUE)[1:20])
cha4_pre <- prune_taxa(cha4_names, cha4)
cha4 <- transform_sample_counts(cha4_pre, function(x) x / sum(x)) # Relative abundance

# Original plot, not ordering:
cha4_plot <- plot_bar(cha4, x="SubjectID", y="Abundance", fill="Genus") +
  scale_fill_manual(values = Dark24) +
  facet_wrap(. ~ Time_point)

cha4_plot

toptaxa_names4 <- cha4_plot[["data"]] %>%
  select(Genus,OTU) %>%
  filter(OTU %in% cha4_names) %>% 
  distinct(OTU, .keep_all = TRUE)
#toptaxa_names4 <- toptaxa_names4[match(cha4_names, toptaxa_names4$OTU),] 
toptaxa_names4 <- toptaxa_names4$Genus

cha4_plot$data$Genus <- factor(cha4_plot$data$Genus, levels = toptaxa_names4)

# All patients - lesions:
cha4_plot[["data"]] %>%
  select(SubjectID, Abundance, Genus, Time_point) %>%
  #left_join(metadata %>% 
  #            mutate(SubjectID = LTCP_patient_ID) %>%
  #            select(SubjectID, microbiome_cluster)) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day0" ~ "D0",
    Time_point == "day30" ~ "D30",
    Time_point == "day60" ~ "D60",
    Time_point == "day90" ~ "D90",
    Time_point == "day120" ~ "D120",
    Time_point == "day150" ~ "D150",
    Time_point == "day180" ~ "D180",
    Time_point == "day210" ~ "D210",
    Time_point == "day240" ~ "D240",
    Time_point == "day270" ~ "D270")) %>%
  ggplot(., aes(x="", y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", width=1, color="white", size=0.2) +
  scale_fill_manual(values = Dark24) +
  coord_polar("y", start=0) +
  guides(size="none") +
  facet_grid(SubjectID ~ factor(Time_point3, levels = c("D0","D30","D60","D90","D120","D150","D180","D210","D240","D270"))) +
  theme_void() +
  theme(strip.text = element_text(size=12))

# DMM (not super clear, I have to study how to choose the best number of Ks----
library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)

dat <- abundances(cha1_pre)
count <- as.matrix(t(dat))

fit <- lapply(1:2, dmn, count = count, verbose=TRUE)
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
plot(aic, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
plot(bic, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")

best <- fit[[which.min(unlist(lplc))]]
mixturewt(best)

ass <- apply(mixture(best), 1, which.max)
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}



# Association of Clinical Outcome and microbiome clusters ----


View(metadata %>% 
       filter(Time_point == "day0") %>%
       filter(Swab_site == "lesion") %>%
       #filter(treatment_outcome == "cure"| treatment_outcome == "failure") %>%
       select(treatment_outcome, LTCP_patient_ID, microbiome_cluster))


metadata %>% # TO FIX
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  filter(treatment_outcome == "cure"| treatment_outcome == "failure") %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "failure" ~ "Failure",
    treatment_outcome == "cure" ~ "Cured")) %>%
ggplot(., aes(x=treatment_outcome2, fill=microbiome_cluster)) +
  geom_bar() +
  theme_classic() + scale_fill_manual(values = c(Dark24[1],
                                                 Dark24[4],
                                                 Dark24[2],
                                                 Dark24[3],
                                                 Dark24[12],
                                                 "gray",
                                                 Dark24[5],
                                                 Dark24[6])) +
  theme(legend.position="right",legend.text = element_text(size = 12),
        axis.title = element_text(size = 12), axis.text = element_text(size = 12),
        axis.title.x = element_blank()) +
  geom_text(stat="count", aes(label=..count..), position = position_stack(vjust = 0.5), size = 4, color="white") +
  ylab("Number of patients")

my_comparisonsAAA <- list(c("Heterogeneous", "Staphylococcus"),
                         c("Heterogeneous","Arcanobacterium"))

metadata %>%  # TO FIX
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  filter(treatment_outcome == "cure"| treatment_outcome == "failure") %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "failure" ~ "Failure",
    treatment_outcome == "cure" ~ "Cured")) %>%
  ggplot(aes(x=microbiome_cluster, y=healing_time_days, color=microbiome_cluster)) +
  geom_hline(yintercept = 90, color = Dark24[8], size=1.5, alpha=0.5) +
  geom_violin(size = 0.7, trim = F, fill=NA) +
  geom_sina(size=2.3) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 3, fontface="bold",
                  max.overlaps = 50) +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[1],
                                                  Dark24[4],
                                                  Dark24[2],
                                                  Dark24[3],
                                                  Dark24[12],
                                                  "gray",
                                                  Dark24[5],
                                                  Dark24[6])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
        legend.position="none") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="black") +
  stat_compare_means(method = "wilcox.test",  comparisons = my_comparisonsAAA,
                     aes(label = ..p.signif..),
                     label.x = 1.4, size = 6, color="black") +
  xlab("") +
  ylab("Healing time (days)") +
  ylim(0,300) #+
  #coord_flip()

metadata %>% 
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  filter(treatment_outcome == "cure"| treatment_outcome == "failure") %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "failure" ~ "Failure",
    treatment_outcome == "cure" ~ "Cured")) %>%
  filter(microbiome_cluster =="Arcanobacterium") %>%
  select(LTCP_patient_ID, treatment_outcome,healing_time_days,microbiome_cluster)

my_comparisonsAAAA <- list(c("Heterogenous", "Dysbiosis"))

metadata %>% 
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  filter(treatment_outcome == "cure"| treatment_outcome == "failure") %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "failure" ~ "Failure",
    treatment_outcome == "cure" ~ "Cured")) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster %in% "Heterogeneous" ~ "Heterogeneous",
    microbiome_cluster %ni% "Heterogeneous" ~ "Dysbiosis"
    )) %>%
  ggplot(aes(x=factor(microbiome_cluster2, levels = c("Heterogeneous", "Dysbiosis")),
             y=healing_time_days, color=microbiome_cluster2)) +
  geom_hline(yintercept = 90, color = Dark24[8], size=1.5, alpha=0.5) +
  geom_violin(size = 0.7, trim = F, fill=NA) +
  geom_sina(size=2.3) +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = rev(c("gray",Dark24[15]))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
        legend.position="none") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="black") +
  #stat_compare_means(method = "wilcox.test",  comparisons = my_comparisonsAAAA,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.4, size = 6, color="black") +
  xlab("") +
  ylab("Healing time (days)") +
  ylim(0,300) +
  coord_flip()

# Add here the time to heal for all clusters:




#_____________


metadata %>% 
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  filter(treatment_outcome == "cure"| treatment_outcome == "failure") %>%
ggplot(., aes(x=microbiome_cluster, fill=treatment_outcome)) +
  geom_bar() +
  theme_classic() + scale_fill_manual(values = Dark24[1:2]) +
  theme(legend.position="right",legend.text = element_text(size = 15),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45)) +
  geom_text(stat="count", aes(label=..count..),position = position_stack(vjust = 0.5), size = 8, color="white")

# Association of Clinical METADATA and microbiome clusters (RAs Maaslin?) ----
library(Maaslin2)

metadata_for_maaslin2 <- 
metadata %>%
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  mutate(local_of_lesion2 = str_remove(local_of_lesion,"\\_.*$")) %>%
  select(SampleID, age, sex, local_of_lesion2,
         size_lesion_mm2, lymphadenopathy, DTH_mm2) %>%
  mutate(local_of_lesion2 = as.factor(local_of_lesion2)) %>%
  column_to_rownames("SampleID")

fit_data1 = Maaslin2(
  input_data = t(cha1_df[toptaxa_names,]), 
  input_metadata = metadata_for_maaslin2, 
  max_significance = 0.05,
  output = "ClinMeta_RAss")

fit_data2 = Maaslin2(
  input_data = t(cha1_df[toptaxa_names,]), 
  input_metadata = metadata_for_maaslin2, 
  max_significance = 0.746525844,
  output = "ClinMeta_RAs_lessRes")

metadata %>%
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  left_join(exportRA) %>%
  ggplot(., aes(x=size_lesion_mm2, y=Staphylococcus)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point() +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "center", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
  ylim(0,1.1) +
  ylab("Staphylococcus RA") + xlab("Size of the lesion (mm2)") +

metadata %>%
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  left_join(exportRA) %>%
  mutate(local_of_lesion2 = str_remove(local_of_lesion,"\\_.*$")) %>%
  ggplot(., aes(x=local_of_lesion2, y=Prevotella)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color=Dark24[1]) +
  geom_jitter() +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right", axis.text.x.bottom = element_text(angle = 90)) +
  
  ylab("Prevotella RA") + xlab("Body site")

metadata %>%
  filter(Time_point == "day0") %>%
  filter(Swab_site == "lesion") %>%
  left_join(exportRA) %>%
  ggplot(., aes(y=size_lesion_mm2, x=microbiome_cluster)) +
  geom_jitter() +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right", axis.text.x.bottom = element_text(angle = 90)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  xlab("Microbiome Cluster") + ylab("Size of the lesion (mm2)")



# Export metadata ----
write_tsv(metadata %>%
            filter(treatment_other_drug == "sbv") %>%
            filter(Swab_site == "lesion") %>%
            filter(Time_point == "day0") %>%
            select(LTCP_patient_ID, healing_time_days, treatment_other_drug, microbiome_cluster), "outputs/microbiome_treatment.txt")



# Only plotting Arcanobacterium samples ----
cha1_plot /
  cha2_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15))



# Evaluating the time course for Arcanobacterium samples
cha_x <- subset_samples(physeq_3, SubjectID %in% c("29271","29272","29256"))
cha_x <- subset_samples(cha_x, Swab_site == "lesion")
cha_x_names <- names(sort(taxa_sums(cha_x), decreasing = TRUE)[1:10])
cha_x_pre <- prune_taxa(cha_x_names, cha_x)
cha_x <- transform_sample_counts(cha_x_pre, function(x) x / sum(x)) # Relative abundance

cha_x@sam_data[["Time_point"]] <- factor(cha_x@sam_data[["Time_point"]], levels = c("day0","day30","day60",
                                                                                    "day90","day120","day150"))

plot_bar(cha_x, x="Time_point", y="Abundance", fill="Genus") +
  scale_fill_manual(values = Dark24) + # FIX THESE COLORS THEY ARE NOT MATCHING EVERYTHING ELSE
  theme(axis.text.x = element_text(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  facet_wrap(. ~ SubjectID, scales = "free")

# Evaluating microbiome structure of Arcanobacterium, Day0:
cha_x1 <- subset_samples(physeq_3, SubjectID %in% Arca_samples_pre)
cha_x1 <- subset_samples(cha_x1,
                         Time_point == "day0" &
                           Swab_site == "lesion")
cha_x1_names <- names(sort(taxa_sums(cha_x1), decreasing = TRUE)[1:20])
cha_x1_pre <- prune_taxa(cha_x1_names, cha_x1)
cha_x1 <- transform_sample_counts(cha_x1_pre, function(x) x / sum(x)) # Relative abundance

cha_x1_plot <- plot_bar(cha_x1, x="SubjectID", y="Abundance", fill="Genus") +
  scale_fill_manual(values = Dark24) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 11), axis.title.x = element_blank())

toptaxa_namesX <- cha_x1_plot[["data"]] %>%
  select(Genus,OTU) %>%
  filter(OTU %in% cha_x1_names) %>% 
  distinct(OTU, .keep_all = TRUE)
toptaxa_namesX <- toptaxa_namesX[match(cha_x1_names, toptaxa_namesX$OTU),] 
toptaxa_namesX <- toptaxa_namesX$Genus
cha_x1_plot$data$Genus <- factor(cha_x1_plot$data$Genus, levels = toptaxa_namesX)

cha_x1_plot[["data"]][["SubjectID"]] <- factor(cha_x1_plot[["data"]][["SubjectID"]],
                                               levels = c("29213","29256","29231",
                                                          "29159","29513","29272","29271"))
cha_x1_plot

# DESeq2 Lesion vs. contralateral ----
physeq_3
head(sample_data(physeq_3)$Swab_site, 25)
physeq_4 <- subset_samples(physeq_3,
                           Swab_site != "environmental" &
                           Time_point == "day0") # all lesions and contralateral from day0, there are unmatched samples.

library("DESeq2")
diagdds <- phyloseq_to_deseq2(physeq_4, ~ Swab_site)

# Quick step to add pseudo "1" to OTU table because I have zeros. Zeros are not handled well when calculating geometric means in DESEq2
# https://mcbl.readthedocs.io/en/master/tut-phyloseq.html
gm_mean <- function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)

# DESeq2:
diagdds <- DESeq(diagdds, test="Wald", fitType="mean")

# Results
res <- DESeq2::results(diagdds, cooksCutoff = FALSE)
alpha <- 0.05
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_4)[rownames(sigtab), ], "matrix"))
sigtab <- cbind(as(res, "data.frame"), as(tax_table(physeq_4)[rownames(res), ], "matrix"))
head(sigtab)

sigtab %>%
  arrange(log2FoldChange, padj)

# Volcano:
# Colunm of significant taxa:
sig_taxa_up <- sigtab %>%
  filter(pvalue < 0.05) %>% filter(log2FoldChange > 1)
sig_taxa_down <- sigtab %>%
  filter(pvalue < 0.05) %>% filter(log2FoldChange < -1)

sig_taxa_up <- sigtab %>%
  filter(padj < 0.05) %>% filter(log2FoldChange > 1)
sig_taxa_down <- sigtab %>%
  filter(padj < 0.05) %>% filter(log2FoldChange < -1)

sig_taxa_up <- sig_taxa_up$Genus
sig_taxa_down <- sig_taxa_down$Genus
#write_tsv(as_tibble(sig_taxa_up), "up.txt") # to be used in GO
#write_tsv(as_tibble(sig_taxa_down), "down.txt") # to be used in GO
#save(sig_taxa_up, file = "RStudio_outputs/data/sig_taxa_up")
#save(sig_taxa_down, file = "RStudio_outputs/data/sig_taxa_down")

sigtab$col0 <- sigtab$Genus
col0_s <- sigtab$col0 %in% sig_taxa_up
col0_s2 <- sigtab$col0 %in% sig_taxa_down
sigtab$col0[!col0_s & !col0_s2] <- NA

sigtab$col1 <- sigtab$Genus
col1_s <- sigtab$col1 %in% sig_taxa_up
col1_s2 <- sigtab$col1 %in% sig_taxa_down
sigtab$col1[col1_s] <- "sig_up"
sigtab$col1[col1_s2] <- "sig_down"
sigtab$col1[!col1_s & !col1_s2] <- "notsig"

sigtab$col2 <- sigtab$col0
col2_s <- sigtab$col2 %in% NA
col2_s2 <- sigtab$col2 %in% sig_taxa_up
col2_s3 <- sigtab$col2 %in% sig_taxa_down
sigtab$col2[col2_s] <- "notsig"
sigtab$col2[col2_s2 | col2_s3] <- "sig"

#sigtab$col3 <- sigtab$Genus
#col3_v <- sigtab$col3 %in% c("GZMB","GNLY","PRF1","IL1B")
#sigtab$col3[!col3_v] <- NA

ggplot(sigtab, aes(y=-log10(pvalue), x=log2FoldChange,
                      color=col1, size=col1)) +
  geom_hline(yintercept = -log10(0.05), color = Dark24[14]) +
  geom_hline(yintercept = -log10(1.037954e-01), color = Dark24[12]) +
  geom_vline(xintercept = 1, color = "gray") + geom_vline(xintercept = -1, color = "gray") +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", Dark24[2], Dark24[1])) +
  scale_size_manual(values = c(0.5,2,2)) +
  geom_text_repel(aes(label = col0, color = col1), size = 4, fontface="bold",
                  max.overlaps = 11) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  #annotate("text", x=-5, y=32, label=paste(length(sig_taxa_down)), size=6, color=Dark24[2], fontface="bold") +
  #annotate("text", x=5, y=32, label=paste(length(sig_taxa_up)), size=6, color=Dark24[1], fontface="bold") +
  xlab("logFC Lesions vs. Contra.Skin")

ggplot(sigtab, aes(y=-log10(padj), x=log2FoldChange,
                   color=col1, size=col1)) +
  geom_hline(yintercept = -log10(0.05), color = Dark24[14]) +
  geom_vline(xintercept = 1, color = "gray") + geom_vline(xintercept = -1, color = "gray") +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", Dark24[2], Dark24[1])) +
  scale_size_manual(values = c(0.5,2,2)) +
  geom_text_repel(aes(label = col0, color = col1), size = 4, fontface="bold",
                  max.overlaps = 11) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  #annotate("text", x=-5, y=32, label=paste(length(sig_taxa_down)), size=6, color=Dark24[2], fontface="bold") +
  #annotate("text", x=5, y=32, label=paste(length(sig_taxa_up)), size=6, color=Dark24[1], fontface="bold") +
  xlab("logFC Lesions vs. Contra.Skin")

# Staphylococcus relative abundance Lesions vs. Contralateral ----
#physeq_4_names <- names(sort(taxa_sums(physeq_4), decreasing = TRUE)[1:20])
#physeq_4x <- prune_taxa(physeq_4_names, physeq_4)
physeq_4x <- transform_sample_counts(physeq_4, function(x) x / sum(x)) # Relative abundance
physeq_4xx <- subset_taxa(physeq_4x, Genus=="Staphylococcus")

physeq_4xx_df <- as.data.frame(t(as.data.frame(physeq_4xx@otu_table)))
physeq_4xx_df <- rownames_to_column(physeq_4xx_df, "SampleID")
colnames(physeq_4xx_df)[2] <- "Staphylococcus_RA"

physeq_4xx_df %>%
  left_join(metadata) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Swab_site2 = case_when(
    Swab_site == "lesion" ~ "Lesion",
    Swab_site == "contralateral" ~ "Contra.Skin")) %>%
ggplot(aes(x=Swab_site2, y=Staphylococcus_RA, color=Swab_site2)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2.3) +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  #geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  stat_compare_means(method = "t.test", label.y = 0.9,
                     aes(label = ..p.signif..), paired = T,
                     label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Relative abundance") +
  ylim(0,1)

# Arcanobacterium relative abundance Lesions vs. Contralateral ----
#physeq_4_names <- names(sort(taxa_sums(physeq_4), decreasing = TRUE)[1:20])
#physeq_4x <- prune_taxa(physeq_4_names, physeq_4)
physeq_4xA <- transform_sample_counts(physeq_4, function(x) x / sum(x)) # Relative abundance
physeq_4xxA <- subset_taxa(physeq_4xA, Genus=="Arcanobacterium")

physeq_4xxA_df <- as.data.frame(t(as.data.frame(physeq_4xxA@otu_table)))
physeq_4xxA_df <- rownames_to_column(physeq_4xxA_df, "SampleID")
colnames(physeq_4xxA_df)[2] <- "Arcanobacterium_RA"

physeq_4xxA_df %>%
  left_join(metadata) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Swab_site2 = case_when(
    Swab_site == "lesion" ~ "Lesion",
    Swab_site == "contralateral" ~ "Contra.Skin")) %>%
  ggplot(aes(x=Swab_site2, y=Arcanobacterium_RA, color=Swab_site2)) +
  geom_line(aes(group = LTCP_patient_ID), color="light gray") +
  #geom_violin(size = 0.7, trim = F) +
  #geom_sina(size=2.3) +
  geom_point(size=2) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 5, fontface="bold") +
  #geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[1])) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none") +
  #stat_summary(fun = median, fun.min = median, fun.max = median,
  #             geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(method = "t.test", label.y = 0.9,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.5, size = 6, color="black") +
  xlab("") +
  ylab("Arcanobacterium relative abundance") +
  ylim(0,0.80)

# Stacked plot over time Staphylococcus Lesion vs. contralateral ----
# Option 1: show mean for all patients with time point available:
samples_timecourse <- metadata %>%
  filter(Swab_site == "lesion",
         treatment_other_drug %in% "sbv",
         Time_point %in% c("day0","day30","day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2)
samples_timecourse <- unique(samples_timecourse$LTCP_patient_ID)

# Option 2 only for Staphylococcus dysbiosis:
str_replace_all(Staphy_samples, "^CL", "")

# Making physeq:
physeq_5 <- subset_samples(physeq_3,
                           SubjectID %in% samples_timecourse) # Option 1
physeq_5 <- subset_samples(physeq_3,
                           SubjectID %in% str_replace_all(Staphy_samples, "^CL", "")) # Option 2

physeq_5 <- subset_samples(physeq_5,
                           Swab_site == "lesion")
physeq_5 <- subset_samples(physeq_5,
                           Time_point %in% c("day0","day30","day60","day90"))

physeq_5_names <- names(sort(taxa_sums(physeq_5), decreasing = TRUE)[1:20])
physeq_5 <- prune_taxa(physeq_5_names, physeq_5)
physeq_5 <- transform_sample_counts(physeq_5, function(x) x / sum(x)) # Relative abundance
physeq_5_df <- as.data.frame(physeq_5@otu_table)

label_taxa1 <- as.data.frame(physeq_5@tax_table)
label_taxa1 <- label_taxa1$Genus

rownames(physeq_5_df) <- label_taxa1
physeq_5_df <- as.data.frame(t(physeq_5_df))
physeq_5_df <- rownames_to_column(physeq_5_df, "SampleID")
physeq_5_df2 <- physeq_5_df %>%
  pivot_longer(!SampleID, names_to = "Genus", values_to = "RA") %>%
  left_join(metadata %>% select(SampleID, Time_point)) %>%
  group_by(Time_point, Genus) %>%
  summarize(RA_mean = mean(RA, na.rm = T))

order1_taxa <- physeq_5_df2 %>%
  arrange(Time_point, RA_mean) %>%
  filter(Time_point == "day0")
order1_taxa <- order1_taxa$Genus

physeq_5_df2 %>%
  arrange(Time_point, RA_mean) %>%
  mutate(Genus = factor(Genus, levels = rev(order1_taxa))) %>%
  mutate(Time_point2 = case_when(
    Time_point == "day0" ~ "Day 0",
    Time_point == "day30" ~ "Day 30",
    Time_point == "day60" ~ "Day 60",
    Time_point == "day90" ~ "Day 90")) %>%
ggplot(., aes(x = Time_point2, y = RA_mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 9), axis.title = element_text(size = 9),
        legend.position="right") +
  xlab("") +
  ylab("Relative abundance mean")

# PCoA (weighted UniFrac) - ALL samples  over time ----
metadata %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "1 Sbv",
    treatment_outcome == "failure" ~ ">1 Sbv")) %>%
  filter(microbiome_cluster != "NA") %>%
  ggplot(., aes(y=PC1_wunifrac, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                                 "D90","D120","D150",
                                                                 "D180","D210","D240",
                                                                 "D270","Contra.Skin D60")), color=microbiome_cluster)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  geom_violin(size = 0.7, trim = F) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_wrap(. ~ treatment_outcome2, nrow = 2)

# Different version:
metadata %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "1 Sbv",
    treatment_outcome == "failure" ~ ">1 Sbv")) %>%
  filter(microbiome_cluster != "NA") %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Staphylococcus_Streptococcus" ~ "M3",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1")) %>%
  #filter(microbiome_cluster2 %in% c("M0","M4","M5","M6","M7")) %>%
  filter(microbiome_cluster2 %in% c("M0","M6","M7")) %>%
  ggplot(., aes(y=PC1_wunifrac, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                                 "D90","D120","D150",
                                                                 "D180","D210","D240",
                                                                 "D270","Contra.Skin D60")),
                color=microbiome_cluster2, shape=Swab_site)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  #geom_violin(size = 0.7, trim = F) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  geom_point(size=3) +
  theme_classic() +
  #scale_color_manual(values = c("gray", Dark24[5], Dark24[6],
  #                             Dark24[12], Dark24[3], 
  #                             Dark24[2], Dark24[1],
  #                             Dark24[4])) +
  scale_color_manual(values = c("dark gray", Dark24[1], 
                                  Dark24[4])) +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_line(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  #stat_summary(fun = median, fun.min = median, fun.max = median,
  #             geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_grid(microbiome_cluster2 ~ factor(treatment_outcome2, levels = c("1 Sbv",">1 Sbv")))

#Other clusters:
metadata %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 1) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "1 Sbv",
    treatment_outcome == "failure" ~ ">1 Sbv")) %>%
  filter(microbiome_cluster != "NA") %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Staphylococcus_Streptococcus" ~ "M3",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1")) %>%
  #filter(microbiome_cluster2 %in% c("M0","M4","M5","M6","M7")) %>%
  filter(microbiome_cluster2 %in% c("M0","M1","M2","M3","M4","M5","M7")) %>%
  ggplot(., aes(y=PC1_wunifrac, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                                 "D90","D120","D150",
                                                                 "D180","D210","D240",
                                                                 "D270","Contra.Skin D60")),
                color=microbiome_cluster2, shape=Swab_site)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  #geom_violin(size = 0.7, trim = F) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  geom_point(size=3) +
  theme_classic() +
  scale_color_manual(values = c("gray",Dark24[5], Dark24[6],
                               Dark24[12], Dark24[3], 
                               Dark24[2], Dark24[4])) +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_line(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  #stat_summary(fun = median, fun.min = median, fun.max = median,
  #             geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_grid(microbiome_cluster2 ~ factor(treatment_outcome2, levels = c("1 Sbv",">1 Sbv")))
print("yep")

# PCoA (Bray-Curtis) - Staphylococcus samples  over time ----
#'%ni%' <- Negate('%in%')
metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", "")) %>%
  filter(Swab_site == "lesion" &
         treatment_other_drug %in% "sbv" &
         Time_point %in% c("day0","day30","day60",
                           "day90","day120","day150",
                           "day180","day210","day240",
                           "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point == "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_bray, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                                 "D90","D120","D150",
                                                                 "D180","D210","D240",
                                                                 "D270",
                                                             "Contra.Skin D60")), color=LTCP_patient_ID)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  geom_violin(size = 0.7, trim = F) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_wrap(. ~ treatment_outcome2, nrow = 1)

# PCoA (weighted UniFrac) - Staphylococcus samples  over time ----
metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", "")) %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_wunifrac, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                             "D90","D120","D150",
                                                             "D180","D210","D240",
                                                             "D270","Contra.Skin D60")), color=LTCP_patient_ID)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  geom_violin(size = 0.7, trim = F) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_wrap(. ~ treatment_outcome2, nrow = 2)

# PCoA (Bray-Curtis) - Staphylococcus samples  over time (MULTIDIMENTIONAL) ----
mygradcol <- colorRampPalette(colors=c(Dark24[1],Dark24[2]))(6)
#'%ni%' <- Negate('%in%')
metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", "")) %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point == "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_bray, x=PC2_bray, color=factor(Time_point3, levels = c("D0","D30","D60",
                                                                     "D90","D120","D150",
                                                                     "D180","D210","D240",
                                                                     "D270",
                                                                     "Contra.Skin D60")))) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID, alpha = 0.5)) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = c(mygradcol,Dark24[20])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  ylab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab(paste("PCoA, axis 2 -",round(bray[["data"]][["Eigvals"]][["PC2"]],2),"%")) #+
  #facet_wrap(. ~ treatment_outcome, nrow = 1)

# PCoA (weighted UniFrac) - Staphylococcus samples  over time (MULTIDIMENTIONAL) ----
metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", "")) %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point == "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_wunifrac, x=PC2_wunifrac, color=factor(Time_point3, levels = c("D0","D30","D60",
                                                                             "D90","D120","D150",
                                                                             "D180","D210","D240",
                                                                             "D270",
                                                                             "Contra.Skin D60")))) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID, alpha = 0.5)) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = c(mygradcol,Dark24[20])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab(paste("PCoA, axis 2 -",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  facet_wrap(. ~ treatment_outcome, nrow = 1)

# PCoA (weighted UniFrac) - Heterogenous samples  over time ----
metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Hetero_samples, "^CL", "")) %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_wunifrac, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                                 "D90","D120","D150",
                                                                 "D180","D210","D240",
                                                                 "D270","Contra.Skin D60")), color=LTCP_patient_ID)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  geom_violin(size = 0.7, trim = F) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_wrap(. ~ treatment_outcome2, nrow = 2)

# PCoA (Bray-Curtis) - Heterogenous samples  over time ----
#'%ni%' <- Negate('%in%')
metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Hetero_samples, "^CL", "")) %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point == "day60") %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_bray, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                             "D90","D120","D150",
                                                             "D180","D210","D240",
                                                             "D270",
                                                             "Contra.Skin D60")), color=LTCP_patient_ID)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  geom_violin(size = 0.7, trim = F) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(bray[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_wrap(. ~ treatment_outcome2, nrow = 2)

# PCoA (weighted UniFrac) - ALL micro cluster samples  over time ----
metadata %>%
  #filter(LTCP_patient_ID %in% str_replace_all(Hetero_samples, "^CL", "")) %>%
  filter(microbiome_cluster != "NA") %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day90" & Swab_site == "contralateral" ~ "Contra.Skin D90",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_wunifrac, x=factor(Time_point3, levels = c("D0","D30","D60",
                                                                 "D90","D120","D150",
                                                                 "D180","D210","D240",
                                                                 "D270","Contra.Skin D60",
                                                                 "Contra.Skin D90")),
                #shape=Time_point3,
                color=microbiome_cluster)) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  geom_line(aes(group = LTCP_patient_ID)) +
  #geom_violin(size = 0.7, trim = F) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = Dark24) +
  #scale_color_gradient() +
  #scale_shape_manual(values = c(19,19,19,19,19,19,19,19,19,22,22)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  #stat_summary(fun = median, fun.min = median, fun.max = median,
  #             geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_wrap(vars(treatment_outcome2), nrow = 2)

metadata %>%
  #filter(LTCP_patient_ID %in% str_replace_all(Hetero_samples, "^CL", "")) %>%
  filter(microbiome_cluster != "NA") %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day90" & Swab_site == "contralateral" ~ "Contra.Skin D90",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  #filter(n() > 2) %>%
  ggplot(., aes(y=PC1_wunifrac, x=PC2_wunifrac,
                #shape=Time_point3,
                color=factor(Time_point3, levels = c("D0","D30","D60",
                                                     "D90","D120","D150",
                                                     "D180","D210","D240",
                                                     "D270","Contra.Skin D60",
                                                     "Contra.Skin D90")))) +
  geom_hline(yintercept = 0, color="light gray", linetype = 1, size=2) +
  #geom_line(aes(group = LTCP_patient_ID)) +
  #geom_violin(size = 0.7, trim = F) +
  geom_point(size=3) +
  theme_classic() + 
  scale_color_manual(values = c(colorRampPalette(colors=c(Dark24[2],Dark24[1]))(9),"dark gray","black")) +
  #scale_color_gradient() +
  #scale_shape_manual(values = c(19,19,19,19,19,19,19,19,19,22,22)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none", strip.background = element_blank(), panel.border = element_rect(size = 1.5, fill = NA),
        strip.text = element_text(size = 13)) +
  #stat_summary(fun = median, fun.min = median, fun.max = median,
  #             geom = "crossbar", width = 0.7, color="black") +
  #stat_compare_means(comparisons = my_comparisonsA,
  #                   method = "t.test", label.y = 1,
  #                   aes(label = ..p.signif..),
  #                   #aes(label = ..p.value..),
  #                   #label.x = 1.5,
  #                   size =5, color="black") +
  ylab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  xlab("") + #ylim(-0.3,0.6) +
  facet_wrap(vars(microbiome_cluster), nrow = 2)

metadata %>%
  #filter(LTCP_patient_ID %in% str_replace_all(Hetero_samples, "^CL", "")) %>%
  filter(microbiome_cluster != "NA") %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150",
                             "day180","day210","day240",
                             "day270") | 
           Swab_site == "contralateral" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day60","day90")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day90" & Swab_site == "contralateral" ~ "Contra.Skin D90",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180",
    Time_point == "day210" & Swab_site == "lesion" ~ "D210",
    Time_point == "day240" & Swab_site == "lesion" ~ "D240",
    Time_point == "day270" & Swab_site == "lesion" ~ "D270",
  )) %>%
  mutate(treatment_outcome2 = case_when(
    treatment_outcome == "cure" ~ "Cured",
    treatment_outcome == "failure" ~ "Failure")) %>%
  ggplot(., aes(x=PC1_bray, y=PC2_bray, color=Time_point)) +
  geom_point(size = 2) +
  geom_path(size=1,#color="#1CA71C",
            arrow=arrow(length=unit(0.3,"cm")),
            aes(x=PC1_wunifrac, y=PC2_wunifrac, group=LTCP_patient_ID, color=Time_point3)) +
  #stat_ellipse(level = 0.95) +
  theme_classic() +
  scale_color_manual(values = c(colorRampPalette(colors=c(Dark24[2],Dark24[1]))(20))) +
  #scale_color_manual(values = c("#FA0011",
  #                              "#CA004A",
  #                              "#4E2A63",
  #                              "#0D2A63")) +
  geom_text_repel(aes(label = Time_point, color=Time_point), size = 3, fontface="bold") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), strip.background = element_blank(),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none", panel.border = element_rect(fill=NA)) +
  
  xlab(paste("PCoA, axis 1 -",round(wunifrac[["data"]][["Eigvals"]][["PC1"]],2),"%")) +
  ylab(paste("PCoA, axis 2 -",round(wunifrac[["data"]][["Eigvals"]][["PC2"]],2),"%")) +
  xlim(-0.35,0.50) + ylim(-0.30,0.40) +
  #coord_fixed() +
  facet_wrap(vars(microbiome_cluster))


# DESeq2 Day0 vs 60/90 Staphylococcus dysbiotic and Stacked plots ----
# Cured patients:
Staphy_Cured <- metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", ""),
         treatment_other_drug == "sbv",
         Time_point == "day0",
         Swab_site == "lesion",
         treatment_outcome == "cure")
Staphy_Cured <- Staphy_Cured$LTCP_patient_ID

Staphy_Failure <- metadata %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", ""),
         treatment_other_drug == "sbv",
         Time_point == "day0",
         Swab_site == "lesion",
         treatment_outcome == "failure")
Staphy_Failure <- Staphy_Failure$LTCP_patient_ID

# DESeq2 Cured patients:
physeq_5.2 <- subset_samples(physeq_3,
                           SubjectID %in% Staphy_Cured)
physeq_5.2 <- subset_samples(physeq_5.2,
                           Swab_site == "lesion")
physeq_5.2 <- subset_samples(physeq_5.2,
                           Time_point %in% c("day0","day60"))
physeq_5.2_names <- names(sort(taxa_sums(physeq_5.2), decreasing = TRUE)[1:10])
physeq_5.2 <- prune_taxa(physeq_5.2_names, physeq_5.2)
diagdds2 <- phyloseq_to_deseq2(physeq_5.2, ~ Time_point)

# Quick step to add pseudo "1" to OTU table because I have zeros. Zeros are not handled well when calculating geometric means in DESEq2
# https://mcbl.readthedocs.io/en/master/tut-phyloseq.html
geoMeans2 <- apply(counts(diagdds2), 1, gm_mean)
diagdds2 <- estimateSizeFactors(diagdds2, geoMeans = geoMeans2)

# DESeq2:
diagdds2 <- DESeq(diagdds2, test="Wald", fitType="mean")

# Results
res2 <- DESeq2::results(diagdds2, cooksCutoff = FALSE)
alpha2 <- 0.05
sigtab2 <- res2[which(res2$padj < alpha2), ]
sigtab2 <- cbind(as(sigtab2, "data.frame"), as(tax_table(physeq_5.2)[rownames(sigtab2), ], "matrix"))
sigtab2 <- cbind(as(res2, "data.frame"), as(tax_table(physeq_5.2)[rownames(res2), ], "matrix"))
head(sigtab2)

# Vector of significant taxa:
sig_taxa_up2 <- sigtab2 %>%
  filter(pvalue < 0.05) %>% filter(log2FoldChange > 1)
sig_taxa_down2 <- sigtab2 %>%
  filter(pvalue < 0.05) %>% filter(log2FoldChange < -1)

sig_taxa_up2 <- sigtab2 %>%
  filter(padj < 0.05) %>% filter(log2FoldChange > 1)
sig_taxa_down2 <- sigtab2 %>%
  filter(padj < 0.05) %>% filter(log2FoldChange < -1)

sig_taxa_up2 <- sig_taxa_up2$Genus
sig_taxa_down2 <- sig_taxa_down2$Genus

sigtab2 %>%
  select(log2FoldChange, Genus, pvalue, padj) %>%
  #filter(pvalue <0.05) %>%
  arrange(log2FoldChange)

# Stacked bar plots:
physeq_6 <- subset_samples(physeq_3, SubjectID %in% str_replace_all(Staphy_samples, "^CL", ""))
#physeq_6 <- subset_samples(physeq_6,
#                           Swab_site == "lesion")

physeq_6 <- subset_samples(physeq_6,
                           Swab_site == "lesion" |
                             Swab_site == "contralateral" & Time_point %in% c("day60","day90"))

physeq_6_names <- names(sort(taxa_sums(physeq_6), decreasing = TRUE)[1:20])
physeq_6_pre <- prune_taxa(physeq_6_names, physeq_6)
physeq_6 <- transform_sample_counts(physeq_6_pre, function(x) x / sum(x)) # Relative abundance

physeq_6@sam_data[["Time_point"]] <- factor(physeq_6@sam_data[["Time_point"]],
                                            levels = c("day0","day30","day60",
                                                       "day90","day120","day150",
                                                       "day180","day210","day240",
                                                       "day270"))


physeq_6_df <- as.data.frame(physeq_6@otu_table)

label_taxa2 <- as.data.frame(physeq_6@tax_table)
label_taxa2 <- label_taxa2$Genus

rownames(physeq_6_df) <- label_taxa2
physeq_6_df <- as.data.frame(t(physeq_6_df))
physeq_6_df <- rownames_to_column(physeq_6_df, "SampleID")

#physeq_6_df2 <- physeq_6_df %>%
#  pivot_longer(!SampleID, names_to = "Genus", values_to = "RA") %>%
#  left_join(metadata %>% select(SampleID, LTCP_patient_ID, Time_point))

physeq_6_df2 <- physeq_6_df %>%
  pivot_longer(!SampleID, names_to = "Genus", values_to = "RA") %>%
  left_join(metadata %>% select(SampleID, Time_point)) %>%
  group_by(Time_point, Genus) %>%
  summarize(RA_mean = mean(RA, na.rm = T))

order2_taxa <- physeq_6_df2 %>%
  arrange(Time_point, RA_mean) %>%
  filter(Time_point == "day0")
order2_taxa <- factor(order2_taxa$Genus)
order2_taxa <- relevel(order2_taxa, "Staphylococcus")

# Individual plots:
# Cured:
physeq_6_df %>%
  pivot_longer(!SampleID, names_to = "Genus", values_to = "RA") %>%
  left_join(metadata %>% select(SampleID, Time_point, LTCP_patient_ID,
                                treatment_outcome, Swab_site)) %>%
  filter(treatment_outcome == "cure") %>%
  #arrange(Time_point, RA_mean) %>%
  mutate(Genus = factor(Genus, levels = rev(order2_taxa))) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day60" & Swab_site == "contralateral" ~ "Contra.Skin D60",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180")) %>%
  filter(Time_point3 %in% c("D0","D30","D60","D90","Contra.Skin D60")) %>%
  ggplot(., aes(x = factor(Time_point3, levels = c("D0","D30","D60","D90","Contra.Skin D60")), y = RA, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
        legend.position="right", strip.text = element_text(size=12), panel.border = element_rect(size=1, fill = NA),
        strip.background = element_blank()) +
  xlab("") +
  ylab("Relative abundance") +
  facet_wrap(. ~ LTCP_patient_ID, nrow = 1)

# Failure:
physeq_6_df %>%
  pivot_longer(!SampleID, names_to = "Genus", values_to = "RA") %>%
  left_join(metadata %>% select(SampleID, Time_point, LTCP_patient_ID,
                                treatment_outcome, Swab_site)) %>%
  filter(treatment_outcome == "failure") %>%
  #arrange(Time_point, RA_mean) %>%
  mutate(Genus = factor(Genus, levels = rev(order2_taxa))) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day90" & Swab_site == "contralateral" ~ "Contra.Skin D90",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180")) %>%
  filter(Time_point3 %in% c("D0","D30","D60","D90","D150","D180","Contra.Skin D90")) %>%
  ggplot(., aes(x = factor(Time_point3, levels = c("D0","D30","D60","D90","D150","D180","Contra.Skin D90")), y = RA, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
        legend.position="none",
        strip.text = element_text(size=12),
        strip.background = element_blank()) +
  xlab("") +
  ylab("Relative abundance") +
  facet_wrap(. ~ LTCP_patient_ID, nrow = 3)


# Patients failure ordered by clinical pictures:
Staphy_Failure_ordered <- c(
  #"29150","29161","29168", # no follow-up
  "29434", # resolved in 30
  "29491","29519", # resolved in 60
  "29253","29539", # resolved in 90
  "29270", # resolved in 150
  "29498","29412" # lesions are worse at 90
  )

physeq_6_df %>%
  pivot_longer(!SampleID, names_to = "Genus", values_to = "RA") %>%
  left_join(metadata %>% select(SampleID, Time_point, LTCP_patient_ID,
                                treatment_outcome, Swab_site)) %>%
  filter(treatment_outcome == "failure") %>%
  filter(LTCP_patient_ID %in% Staphy_Failure_ordered) %>%
  #arrange(Time_point, RA_mean) %>%
  mutate(Genus = factor(Genus, levels = rev(order2_taxa))) %>%
  mutate(Time_point3 = case_when(
    Time_point == "day90" & Swab_site == "contralateral" ~ "Contra.Skin D90",
    Time_point == "day0" & Swab_site == "lesion" ~ "D0",
    Time_point == "day30" & Swab_site == "lesion" ~ "D30",
    Time_point == "day60" & Swab_site == "lesion" ~ "D60",
    Time_point == "day90" & Swab_site == "lesion" ~ "D90",
    Time_point == "day120" & Swab_site == "lesion" ~ "D120",
    Time_point == "day150" & Swab_site == "lesion" ~ "D150",
    Time_point == "day180" & Swab_site == "lesion" ~ "D180")) %>%
  #filter(Time_point3 != "Contra.Skin D60") %>%
  filter(Time_point3 %in% c("D0","D30","D60","D90","D150","D180","Contra.Skin D90")) %>%
  ggplot(., aes(x = factor(Time_point3, levels = c("D0","D30","D60","D90","D150","D180","Contra.Skin D90")), y = RA, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
        legend.position="none", panel.border = element_rect(size=1, fill = NA),
        strip.text = element_text(size=12),
        strip.background = element_blank()) +
  xlab("") +
  ylab("Relative abundance") +
  facet_wrap(. ~ LTCP_patient_ID, nrow = 2)

# Mean for all patients. It is NOT appropriate because Cure vs Failure are different.
# This has to be re-done because I included Contralateral skin. Not sure if RAmean is calculated correctly.
physeq_6_df %>%
  pivot_longer(!SampleID, names_to = "Genus", values_to = "RA") %>%
  left_join(metadata %>% select(SampleID, Time_point, treatment_outcome)) %>%
  group_by(Time_point, Genus) %>%
  summarize(RA_mean = mean(RA, na.rm = T)) %>%
  mutate(Genus = factor(Genus, levels = rev(order2_taxa))) %>%
  mutate(Time_point2 = case_when(
    Time_point == "day0" ~ "D0",
    Time_point == "day30" ~ "D30",
    Time_point == "day60" ~ "D60",
    Time_point == "day90" ~ "D90",
    Time_point == "day150" ~ "D150",
    Time_point == "day180" ~ "D180")) %>%
  ggplot(., aes(x = factor(Time_point2, levels = c("D0","D30","D60","D90","D150","D180")), y = RA_mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
        legend.position="right") +
  xlab("") +
  ylab("Relative abundance")

# NMIT (not appropriate for this study. Missing data is a factor, therefore the dataset should be analyzed individually, per patient basis) ----
# See: https://docs.qiime2.org/2021.4/tutorials/longitudinal/
# and: https://codeocean.com/capsule/6099971/tree
# and: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01089-8

# First, filter metadata to only calculate NMIT in Staphy dysbiosis samples:
Staphy_metadata_pre <- 
  metadata %>%
  select(SampleID, LTCP_patient_ID, Swab_site, Time_point, treatment_other_drug, microbiome_cluster) %>%
  filter(LTCP_patient_ID %in% str_replace_all(Staphy_samples, "^CL", "")) %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2)

Staphy_metadata_s <- Staphy_metadata_pre$SampleID

MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples <- read_delim("Qiime2/MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples.tsv", 
           delim = "\t", escape_double = FALSE, 
           trim_ws = TRUE)
colnames(MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples)[1] <- "SampleID"
Staphy_metadata <- MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples %>%
  mutate(StaphyDysbiosisincludenmit = case_when(
    SampleID %in% Staphy_metadata_s ~ "yes",
    SampleID %ni% Staphy_metadata_s ~ "no"))

write_tsv(Staphy_metadata, "Qiime2/StaphyDysbiosis_metadata.tsv") # add a hashtag in "#SampleID" the tsv file manually

# Filter samples with at least 3 time points collected:
samp_metadata_pre <- 
  metadata %>%
  select(SampleID, LTCP_patient_ID, Swab_site, Time_point, treatment_other_drug, microbiome_cluster) %>%
  filter(Swab_site == "lesion" &
           treatment_other_drug %in% "sbv" &
           Time_point %in% c("day0","day30","day60",
                             "day90","day120","day150","day180","day210","day240","day270")) %>%
  arrange(LTCP_patient_ID) %>%
  group_by(LTCP_patient_ID) %>%
  filter(n() > 2)

samp_metadata_s <- samp_metadata_pre$SampleID

samp_metadata <- MiSeqV1V3_38_humanLesion_metadata_removedLowReadSamples %>%
  mutate(sampDysbiosisincludenmit = case_when(
    SampleID %in% samp_metadata_s ~ "yes",
    SampleID %ni% samp_metadata_s ~ "no"))

write_tsv(samp_metadata, "Qiime2/sampDysbiosis_metadata.tsv") # add a hashtag in "#SampleID" the tsv file manually



# Import pcoa result from nmit:
# Staph samples:
nmit <- read_qza("Qiime2/core-metrics-results/StaphySamples_nmit-pc.qza")

metadata <- metadata %>% 
  left_join(nmit$data$Vectors %>%
              select(SampleID, PC1, PC2, PC3), by= "SampleID") %>%
  dplyr::rename(PC1_nmit = PC1,
                PC2_nmit = PC2,
                PC3_nmit = PC3)

metadata %>%
  ggplot(., aes(x=PC1_nmit, y=PC2_nmit,
                #shape=treatment_outcome,
                color=treatment_outcome)) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[3],
                                                  Dark24[1],Dark24[4])) +
 # scale_shape_manual(values = c(67,70,70,67)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none") +
  stat_ellipse(level = 0.95) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 2.5, fontface="bold",colour="black") +
  xlab(paste("NMIT PCoA, axis 1 -",round(nmit[["data"]][["ProportionExplained"]][["PC1"]]*100,1),"%")) +
  ylab(paste("NMIT PCoA, axis 2 -",round(nmit[["data"]][["ProportionExplained"]][["PC2"]]*100,1),"%")) +

metadata %>%
  ggplot(., aes(x=PC1_nmit, y=PC2_nmit,
                #shape=treatment_outcome,
                color=treatment_outcome)) +
  geom_point(size=4) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[3],
                                                  Dark24[1],Dark24[4])) +
  #scale_shape_manual(values = c(67,70,70,67)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none") +
  stat_ellipse(level = 0.95) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 2.5, fontface="bold",colour="black") +
  xlab(paste("NMIT PCoA, axis 1 -",round(nmit[["data"]][["ProportionExplained"]][["PC1"]]*100,1),"%")) +
  ylab(paste("NMIT PCoA, axis 2 -",round(nmit[["data"]][["ProportionExplained"]][["PC2"]]*100,1),"%")) +
  xlim(0,40) + ylim(-40,10)

# 3D plots:
#library("gg3D")
library(plotly)
metadata

#ggplot(metadata, aes(x=PC1_nmit, y=PC2_nmit, z=PC3_nmit, color=treatment_outcome)) + 
#  theme_void() +
#  axes_3D() +
#  stat_3D()

plot_ly(metadata, x = ~PC1_nmit, y = ~PC2_nmit, z = ~PC3_nmit) %>%
  add_markers(color = ~treatment_outcome,
              colors = c(rev(Dark24[1:2])))

plot_ly(x=metadata$PC1_nmit, y=metadata$PC2_nmit, z=metadata$PC3_nmit, type="scatter3d",
        mode="markers",
        color=metadata$treatment_outcome, colors = c(rev(Dark24[1:2])))

# However, 2 Failure patients had their longitudinal microbiome very distinct from all the other patients. This is actually
# in match with the whole discoveries, but it is hard to visualized in a multidimentional space, they are too far away.
# So I'm gonna recalculate the NMIT without these 2 patients.

# All pateints now (>2 time points):
nmit_samp <- read_qza("Qiime2/core-metrics-results/sampSamples_nmit-pc.qza")

metadata <- metadata %>% 
  left_join(nmit_samp$data$Vectors %>%
              select(SampleID, PC1, PC2, PC3), by= "SampleID") %>%
  dplyr::rename(PC1_nmit_samp = PC1,
                PC2_nmit_samp = PC2,
                PC3_nmit_samp = PC3)

metadata %>%
  ggplot(., aes(x=PC1_nmit_samp, y=PC2_nmit_samp,
                #shape=treatment_outcome,
                color=treatment_outcome)) +
  geom_point(size=3) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[3],
                                                  Dark24[1],Dark24[4])) +
  # scale_shape_manual(values = c(67,70,70,67)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none") +
  stat_ellipse(level = 0.95) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 2.5, fontface="bold",colour="black") +
  xlab(paste("NMIT PCoA, axis 1 -",round(nmit[["data"]][["ProportionExplained"]][["PC1"]]*100,1),"%")) +
  ylab(paste("NMIT PCoA, axis 2 -",round(nmit[["data"]][["ProportionExplained"]][["PC2"]]*100,1),"%"))
  
  metadata %>%
  ggplot(., aes(x=PC1_nmit, y=PC2_nmit,
                #shape=treatment_outcome,
                color=treatment_outcome)) +
  geom_point(size=4) +
  theme_classic() + scale_color_manual(values = c(Dark24[2],Dark24[3],
                                                  Dark24[1],Dark24[4])) +
  #scale_shape_manual(values = c(67,70,70,67)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="none") +
  stat_ellipse(level = 0.95) +
  #geom_text_repel(aes(label = LTCP_patient_ID), size = 2.5, fontface="bold",colour="black") +
  xlab(paste("NMIT PCoA, axis 1 -",round(nmit[["data"]][["ProportionExplained"]][["PC1"]]*100,1),"%")) +
  ylab(paste("NMIT PCoA, axis 2 -",round(nmit[["data"]][["ProportionExplained"]][["PC2"]]*100,1),"%")) +
  xlim(0,40) + ylim(-40,10)


# Next ----














































































