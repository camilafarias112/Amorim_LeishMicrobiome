# Here I import the results from Ba's bacteria qPCR.
# See ~Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/3rd_LesionRNAseq_Dataset/qPCRs/Notes_qPCR to see how I processed the calculations.

# Load R objects ----
load("../../RStudio_outputs/data/Dark24")
load("../../RStudio_outputs/data/targets")
#load("../../RStudio_outputs/data/targets_w_Gricedata")
load("../../RStudio_outputs/data/myDGEList.filtered.norm")
load("../../RStudio_outputs/data/CL_samples")
load("../../RStudio_outputs/data/log2.cpm.filtered.norm")

# Libraries ----
library(ggrepel)

library(ggpubr)
library(patchwork)
library(reshape2)
library(ggrepel)
library(gplots)
library(ggforce)
library(ggdist)

library(tidyverse)

# Processing and normalizing the interpolated values ----
# Prism's export:
qPCR_2 <- read_delim("Nonlin fit of 16S.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
qPCR_2 <- qPCR_2 %>%
  select(1,3) %>%
  rename(interpolated_2 = 2)

qPCR_3 <- read_delim("Nonlin fit of Saureus_removeND.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)

qPCR_3 <- qPCR_3 %>%
  rename(interpolated_3 = 4,
         detectable_Saureus_qPCR = 3) %>%
  select(1,3,4)

# Tissue weight:
tw <- read_delim("~/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/3rd_LesionRNAseq_Dataset/RStudio_inputs/TissueWeight_FirstHalfBiopsy.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)

tw <- tw %>%
  mutate(TissueWeight_mg = TissueWeight_g*1000)

# Normalize Bacteria data to tissue weight, and merge everything:
qPCR_23 <- left_join(targets, qPCR_2, by="sample_RNAseq")
qPCR_23 <- left_join(qPCR_23, qPCR_3, by="sample_RNAseq")
qPCR_23 <- left_join(qPCR_23, tw, by="sample_RNAseq")
colnames(qPCR_23)
qPCR_23 %>%
  select(sample_RNAseq, interpolated_2, interpolated_3, detectable_Saureus_qPCR,TissueWeight_mg)

qPCR_23 <- qPCR_23 %>%
  mutate(Saureus_mg = interpolated_3/TissueWeight_mg) %>%
  mutate(s16_mg = interpolated_2/TissueWeight_mg) %>%
  mutate(ratio_SaureusTotal = Saureus_mg*100/s16_mg)

qPCR_23 <- qPCR_23 %>%
  mutate(Saureus_qPCRdetect = case_when(
         Saureus_mg > 0 ~ "detect",
         Saureus_mg %in% NA ~ "non_detect"))

# Visualizing S. aureus ----

load("../../../3rd_LesionRNAseq_Dataset/mapping_Saureus_pangenome/outputs/targets_info")

# cpm vs. pg:
qPCR_23 %>%
  filter(group_RNAseq == "CL") %>%
  left_join(targets_info) %>%
ggplot(.,
       aes(x=cpm_saureus, y=Saureus_mg)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2.5, fontface="bold",colour="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
  xlab("n S. aureus (CPM/10000)") + ylab(paste("S aureus load (pg cDNA/mg biopsy)"))


orderSaureusLoad <- qPCR_23 %>%
  filter(group_RNAseq == "CL") %>%
  left_join(targets_info) %>%
  arrange(cpm_saureus)

orderSaureusLoad <- orderSaureusLoad$sample_RNAseq

plot1 <- ggplot(qPCR_23 %>%
         filter(group_RNAseq == "CL") %>%
           left_join(targets_info),
                aes(x=factor(sample_RNAseq, levels = orderSaureusLoad), y=cpm_saureus,
                    fill=treatment_outcome)) +
  geom_bar(stat="identity") + 
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text()) +
  geom_hline(yintercept = .05, color = "black", linetype=2) +
  guides(fill=guide_legend(title="Treat.Outcome", ncol = 1)) +
  ylab("S. aureus (CPM)")

plot2 <- ggplot(qPCR_23 %>%
         filter(group_RNAseq == "CL"),
       aes(x=factor(sample_RNAseq, levels = orderSaureusLoad), y=Saureus_mg,
           fill=treatment_outcome)) +
  geom_bar(stat="identity") + 
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(), legend.position = "none") +
  #geom_hline(yintercept = 0.025, color = "black", linetype=2) +
  ylab("S aureus load (pg cDNA/mg biopsy)")

plot3 <- ggplot(qPCR_23 %>%
                  filter(group_RNAseq == "CL"),
                aes(x=factor(sample_RNAseq, levels = orderSaureusLoad), y=ratio_SaureusTotal,
                    fill=treatment_outcome)) +
  geom_bar(stat="identity") + 
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(), legend.position = "none") +
  ylab("Ratio S. aureus/total 16S") + ylim(0,100)

plot1 / plot2
plot1 / plot2 / plot3


hscl <- ggplot(qPCR_23,
       aes(x=factor(group_RNAseq, levels=c("HS","CL")),
           y=Saureus_mg,
           color=treatment_outcome,
           #color=group_RNAseq,
           NULL)) +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="grey") +
  geom_jitter(position=position_jitter(0.1), size=2) +
  annotate("text", x=1.5,y=.030, label=paste("CT values < RNAse-free H2O"),
           size=2, fontface="bold",color="dark grey") +
  geom_hline(yintercept = .01, color="grey") +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  geom_text_repel(aes(label = sample_RNAseq), size = 2.5, fontface="bold",colour="black") +
  xlab("") + ylab(paste("S. aureus load (pg/ul cDNA/mg biopsy)")) #+ ylim(0,.1)
hscl

ggplot(qPCR_23 %>%
         left_join(targets_info),
       aes(x=factor(group_RNAseq, levels=c("HS","CL")),
           y=cpm_saureus,
           #color=treatment_outcome,
           color=Saureus_qPCRdetect,
           NULL)) +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="grey") +
  geom_jitter(position=position_jitter(0.1), size=2) +
  #annotate("text", x=1.5,y=.030, label=paste("CT values < RNAse-free H2O"),
  #         size=2, fontface="bold",color="dark grey") +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
 # geom_text_repel(aes(label = sample_RNAseq), size = 2.5, fontface="bold",colour="black") +
  xlab("") + ylab(paste("n S. aureus (CPM/10000)")) +
  ylim(0,0.5) # This plot need more work

# Failure vs. Cure:
treatment <- ggplot(qPCR_23 %>%
         filter(group_RNAseq == "CL") %>%
         filter(treatment_other_drug == "sbv"),
       aes(x=treatment_outcome, y=Saureus_mg,
           color=treatment_outcome)) +
  geom_jitter(position=position_jitter(0.1), size=2) +
  annotate("text", x=1.5,y=.030, label=paste("CT values < RNAse-free H2O"),
           size=2, fontface="bold",color="dark grey") +
  geom_hline(yintercept = .01, color="grey") +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, hjust = 1, angle = 45), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_compare_means(method = "t.test", #label.y = 30000,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 5, color="red") +
  xlab("") + ylab(paste("S aureus load (pg cDNA/mg biopsy)")) +
  NULL
treatment

# DGE Saureus detect NOT UPDATEDDDDDDD----
detectable_Saureus_qPCR_f <- factor(qPCR_23$detectable_Saureus_qPCR[1:51])
design10 <- model.matrix(~0 + detectable_Saureus_qPCR_f)
colnames(design10) <- levels(detectable_Saureus_qPCR_f)
v.DEGList.filtered.norm10 <- voom(myDGEList.filtered.norm[,CL_samples], design10, plot = FALSE)
fit10 <- lmFit(v.DEGList.filtered.norm10, design10)
contrast.matrix10 <- makeContrasts(Saureus_detec_qPCR = D - ND,
                                 levels=design10)
fits10 <- contrasts.fit(fit10, contrast.matrix10)
ebFit10 <- eBayes(fits10)
myTopHits10 <- topTable(ebFit10, adjust ="BH", number=40000, sort.by="logFC")
myTopHits10 <- as_tibble(myTopHits10, rownames = "geneSymbol")

sig_genes_up10 <- myTopHits10 %>%
  filter(P.Value < 0.05, logFC > 1)
sig_genes_down10 <- myTopHits10 %>%
  filter(P.Value < 0.05, logFC < -1)
sig_genes_up10 <- sig_genes_up10$geneSymbol
sig_genes_down10 <- sig_genes_down10$geneSymbol
#write_tsv(as_tibble(sig_genes_up10), "../../../../../../Desktop/sig_genes_up_pval0516Continuous.txt")
#write_tsv(as_tibble(sig_genes_down10), "../../../../../../Desktop/sig_genes_down_pval0516Continuous.txt")

myTopHits10$col0 <- myTopHits10$geneSymbol
col0_v10 <- myTopHits10$col0 %in% sig_genes_up10
col0_v1010 <- myTopHits10$col0 %in% sig_genes_down10
myTopHits10$col0[!col0_v10 & !col0_v1010] <- NA

#volcano:
ggplot(myTopHits10, aes(y=-log10(P.Value), x=logFC)) +
  geom_point(size=1, color="gray") +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  geom_text_repel(aes(label = col0), size = 2.2, fontface="bold",color="black") +
  xlab("logFC D vs. ND")

# Visualizing total 16S ----
load("../../../Dataset_Integration/outputs/phenotype")
phenotype <- rownames_to_column(phenotype, "LTCP_patient_ID")
save(qPCR_23, file = "qPCR_23")
write_tsv(qPCR_23, "qPCR_23.txt")

qPCR_23 %>%
ggplot(.,aes(x=factor(group_RNAseq, levels=c("HS","CL")),
               y=s16_mg,
               #color=treatment_outcome,
               color=group_RNAseq,
               NULL)) +
  geom_violin(size = 0.7, trim = F) +
  geom_sina(size=2) +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="grey") +
  #annotate("text", x=1.5,y=1.3, label=paste("1.07pg/mg"),
  #         size=4, fontface="bold",color="black") +
  geom_hline(yintercept = 1.07, color=Dark24[19], linetype=2, size=1) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  #geom_jitter(position=position_jitter(0.1), size=3, shape=1, color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  #stat_compare_means(method = "wilcox.test", label.y = 7.5, aes(label = ..p.signif..),
  #                   label.x = 1.3, size = 7, color="black") +
  xlab("") + ylab(paste("16S ribosomal subunit (pg/mg biopsy)")) #+ ylim(0,.1)

qPCR_23 %>%
  left_join(phenotype, by="LTCP_patient_ID") %>%
  filter(microbiome_cluster != "NA") %>%
ggplot(.,aes(x=microbiome_cluster,
                      y=s16_mg,
                      #color=treatment_outcome,
                      color=microbiome_cluster,
                      NULL)) +
  geom_violin(size = 0.7, trim = F) +
  #geom_sina(size=2) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 3, fontface="bold") +
  #annotate("text", x=1.5,y=1.3, label=paste("1.07pg/mg"),
  #         size=4, fontface="bold",color="black") +
  geom_hline(yintercept = 1.07, color=Dark24[12]) +
  theme_classic() + scale_color_manual(values = c(Dark24[4],
                                                  Dark24[2],
                                                  Dark24[6],
                                                  "gray",
                                                  Dark24[5],
                                                  Dark24[1],
                                                  Dark24[3])) + 
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  #stat_compare_means(method = "wilcox.test", label.y = 7.5, aes(label = ..p.signif..),
  #                   label.x = 1.3, size = 7, color="black") +
  xlab("") + ylab(paste("16S ribosomal subunit (pg/mg biopsy)")) #+ ylim(0,.1)

plot4 <- ggplot(qPCR_23 %>%
                  filter(group_RNAseq == "CL"),
                aes(x=factor(sample_RNAseq, levels = orderSaureusLoad), y=s16_mg,
                    fill=treatment_outcome)) +
  geom_bar(stat="identity") + 
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(), legend.position = "none") +
  ylab("16S load (pg/ul cDNA/mg biopsy)")

plot4 / plot1 / plot2 / plot3

ggplot(qPCR_23 %>%
         filter(treatment_outcome == "cure"| treatment_outcome == "failure"),
       aes(x=s16_mg_cat, fill=treatment_outcome)) +
  geom_bar() +
  theme_classic() + scale_fill_manual(values = c(Dark24[1],Dark24[3])) +
  theme(legend.position="right",legend.text = element_text(size = 15),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15)) +
  geom_text(stat="count", aes(label=..count..),position = position_stack(vjust = 0.5), size = 8, color="white")

ggplot(qPCR_23 %>%
         filter(group_RNAseq == "CL"),
       aes(x=s16_mg_cat, fill=treatment_outcome)) +
  geom_bar() +
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(legend.position="right",legend.text = element_text(size = 15),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15)) +
  geom_text(stat="count", aes(label=..count..),position = position_stack(vjust = 0.5), size = 8, color="white")

# Percentile:
qPCR_23 %>%
  filter(group_RNAseq != "HS") %>%
  ggplot(.,aes(x=group_RNAseq, y=s16_mg)) +
  geom_point(aes(fill=group_RNAseq), size=3, shape=21) + 
  geom_boxplot(fill=NA, width=0.3, color=Dark24[20], size=0.5) +
  annotate("text", x=0.63, y=0.4, label="Percentile:", color="black", size=3, fontface=2) +
  annotate("text", x=0.63, y=0.3, label=">75th", color=Dark24[1], size=3) +
  annotate("text", x=0.63, y=0.23, label="50-75th", color="#CF7487", size=3) +
  annotate("text", x=0.63, y=0.15, label="25-50th", color="dark gray", size=3) +
  annotate("text", x=0.63, y=0.05, label="0-25th", color="gray", size=3) +
  #geom_violin(size = 0.7, trim = F) +
  #geom_sina(size=2.3) +
  theme_classic() + scale_fill_manual(values = Dark24[-1]) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title = element_text(size = 13),
        legend.position="right") +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 3, fontface="bold", nudge_x = 0.5, max.overlaps = 100) +
  xlab("") + ylab(paste("Failure gene signature ES")) +
  #ylim(0,.40) +
  NULL

# Correlation Total 16S vs. S. aureus CPMs (pangenome) ----
qPCR_23 %>%
  left_join(phenotype, by="LTCP_patient_ID") %>%
  left_join(targets_info) %>%
  filter(microbiome_cluster != "NA") %>%
ggplot(.,aes(x=cpm_saureus,
             y=s16_mg,
             #color=treatment_outcome,
             color=microbiome_cluster,
             NULL)) +
  geom_text_repel(aes(label = LTCP_patient_ID), size = 3, fontface="bold") +
  theme_classic() + scale_color_manual(values = c(Dark24[4],
                                                  Dark24[2],
                                                  Dark24[6],
                                                  "gray",
                                                  Dark24[5],
                                                  Dark24[1],
                                                  Dark24[3])) + 
  geom_point(size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="right") +
  xlab("S. aureus transcripts (CPM)") + ylab(paste("16S ribosomal subunit (pg/mg biopsy)")) +
  #xlim(0,2000) + ylim(0,4) +
  #stat_cor(method = "spearman",
  #         label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
  #         output.type = "expression",
  #         geom = "text", position = "identity", na.rm = FALSE,
  #         inherit.aes = TRUE) +
  NULL


# DE qPCR bacteria biopsy, absolute (continuous) ----
qPCR_23 <- qPCR_23 %>%
  mutate(s16_mg_cat = case_when(
    s16_mg > 1.07 ~ "high_bacteria",
    s16_mg < 1.07 ~ "low_bacteria"))

# 1) Categorical
#s16_mg_factor <- factor(qPCR_23$s16_mg_cat[1:51])

#design <- model.matrix(~0 + s16_mg_factor)
#colnames(design) <- levels(s16_mg_factor)
#v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm[,1:51], design, plot = FALSE)
#fit <- lmFit(v.DEGList.filtered.norm, design)
#contrast.matrix <- makeContrasts(disease = high_bacteria - low_bacteria,
#                                 levels=design)
#fits <- contrasts.fit(fit, contrast.matrix)
#ebFit <- eBayes(fits)

# 2) Continuous
design9 <- model.matrix(~qPCR_23$s16_mg[1:51]) # pay attention here, it's not automatic! Only selecting values from patients.
colnames(design9)[2] <- "qPCR_3"
head(design9)
dim(design9)
v.DEGList.filtered.norm9 <- voom(myDGEList.filtered.norm[,CL_samples], design9)
fit9 <- lmFit(v.DEGList.filtered.norm9, design9)
fits9 <- contrasts.fit(fit9, coefficients = 2)
ebFit9 <- eBayes(fits9)

myTopHits9 <- topTable(ebFit9, adjust ="BH", number=40000, sort.by="logFC")
myTopHits9 <- as_tibble(myTopHits9, rownames = "geneSymbol")

sig_genes_up9_less <- myTopHits9 %>%
  filter(P.Value < 0.01, logFC > 0.2)
sig_genes_down9_less <- myTopHits9 %>%
  filter(P.Value < 0.01, logFC < -0.2)
sig_genes_up9_less <- sig_genes_up9_less$geneSymbol
sig_genes_down9_less <- sig_genes_down9_less$geneSymbol
write_tsv(as_tibble(sig_genes_up9_less), "../../../../../../Desktop/sig_genes_up_lessres_16Continuous.txt")
write_tsv(as_tibble(sig_genes_up9_less), "sig_genes_up_lessres_16Continuous.txt")
write_tsv(as_tibble(sig_genes_down9_less), "../../../../../../Desktop/sig_genes_down_lessres_16Continuous.txt")

sig_genes_up9 <- myTopHits9 %>%
  filter(P.Value < 0.05, logFC > 0.25)
sig_genes_down9 <- myTopHits9 %>%
  filter(P.Value < 0.05, logFC < -0.25)
sig_genes_up9 <- sig_genes_up9$geneSymbol
sig_genes_down9 <- sig_genes_down9$geneSymbol

myTopHits9$col0 <- myTopHits9$geneSymbol
col0_v9 <- myTopHits9$col0 %in% sig_genes_up9
col0_v99 <- myTopHits9$col0 %in% sig_genes_down9
myTopHits9$col0[!col0_v9 & !col0_v99] <- NA

# Genes to show:
geneshow1 <- c("IL1A","IL1B","IL1RN","IL1R2")

geneshow2 <- c("CXCL8","CXCL6","CXCL1","CXCL2","CXCL3","CXCL5",
               "CXCR2","CXCR1")



geneshow <- unique(c("CXCL6","CCL21","CCL20","IL24","CXCL1","PTGS2","IL1A","IL36A",
              "ELF3","CXCR1","IL1B","CXCR2","PROK2","S100A12","FFAR3","ITGB6",
              "CXCL6","IL11","CSF2","CCL21","CCL20","IL1R2","IL24","OSM","CXCL1",
              "IL1A","CXCR1","IL1B","CXCR2","CXCL6","CXCR1","CCL20","CXCR2",
              "PROK2","CXCL1","CMTM2","DEFB4A","CXCL6","IL1RN","CSF2","CCL21",
              "CCL20","IL1R2","IL24","OSM","CXCL1","IL1A","IL36A","IL1B","DEFB4A",
              "CXCL6","CXCR1","CCL21","CCL20","CXCR2","CXCL1","CXCL6","CXCR1",
              "CCL21","CCL20","CXCR2","CXCL1","IL1A","CXCL6","IL11","CSF2",
              "CCL20","IL1B","HCAR2","CXCL6","HCAR3","CXCR1","CCL21","CCL20",
              "CXCR2","CXCL1","CXCR1","CCL21","CXCR2","CCL21","CXCR2","CXCL1",
              "CXCL6","CXCR1","CCL21","CCL20","CXCR2","CXCL1","CXCL6","CCL21","CCL20","CXCL1")) # This list is sigup p=0.01 and .3 coeffi

geneshowA <- unique(c("CXCL6", "IL1RN", "CCL21", "CCL20", "CXCL1", "PTGS2", "IL1A", "IL36A", "ELF3", "IL1B", "CXCR2", "PROK2", "S100A12", "FFAR3", "ITGB6",
                     "CXCL6", "IL1RN", "CSF2", "CCL21", "CCL20", "IL1R2", "OSM", "HLA-G", "IL1A", "IL36A", "CXCR1", "IL1B", "CXCR2", "ITGB6", "DEFB4A", "ULBP2"))
geneshowB <- unique(c("SPRR3", "LAMB3", "LAMA3", "LAMC2", "PTHLH", "SPRR1B", "SPRR2D",
                      "SPRR3", "CDH3", "SERPINB2", "IL24", "TGFA", "NRG1", "ITGB6"))

myTopHits9$col0000 <- myTopHits9$geneSymbol
col0000_v <- myTopHits9$col0000 %in% geneshow
myTopHits9$col0000[!col0000_v] <- NA

myTopHits9$col00002 <- myTopHits9$geneSymbol
myTopHits9$col00002[col0000_v] <- "inflam"
myTopHits9$col00002[!col0000_v] <- "others"

myTopHits9$col00 <- myTopHits9$geneSymbol
col00_v9 <- myTopHits9$col00 %in% geneshow1
col00_v10 <- myTopHits9$col00 %in% geneshow2
myTopHits9$col00[!col00_v9 & !col00_v10] <- NA

myTopHits9$col000 <- myTopHits9$geneSymbol
myTopHits9$col000[col00_v9] <- "IL1s"
myTopHits9$col000[col00_v10] <- "CXCs"
myTopHits9$col000[!col00_v9 & !col00_v10] <- "others"
myTopHits9$col000 <- factor(myTopHits9$col000, levels = c("others","IL1s","CXCs"))

myTopHits9$colx <- myTopHits9$geneSymbol
col00_vx <- myTopHits9$colx %in% geneshowA
col00_vxx <- myTopHits9$colx %in% geneshowB
col000000_v <- myTopHits9$colx %in% sig_genes_up9[1:10]
myTopHits9$colx[!col000000_v & !col00_vx & !col00_vxx] <- NA

myTopHits9$colxx <- myTopHits9$geneSymbol
myTopHits9$colxx[col000000_v] <- "Top genes"
myTopHits9$colxx[col00_vx] <- "Inflammatory Immune responses"
myTopHits9$colxx[col00_vxx] <- "Wound healing and epidermis development"
myTopHits9$colxx[!col000000_v & !col00_vx & !col00_vxx] <- "others"
myTopHits9$colxx <- factor(myTopHits9$colxx, levels = c("others","Top genes","Inflammatory Immune responses","Wound healing and epidermis development"))

color1 <- c(rep(c("#AF0038"), times=10),rep(c("gray"), times=length(col00_v9)))

# Top 25 genes:
myTopHits9$col00000 <- myTopHits9$geneSymbol
col00000_v <- myTopHits9$col00000 %in% sig_genes_up9[1:10]
myTopHits9$col00000[!col00000_v] <- NA

myTopHits9$col000002 <- myTopHits9$geneSymbol
myTopHits9$col000002[col00000_v] <- "top"
myTopHits9$col000002[!col00000_v] <- "others"

#volcano:
ggplot(myTopHits9 %>%
         arrange(rev(col000)), aes(y=-log10(P.Value), x=logFC)) +
  geom_point(color="gray") + scale_size_manual(values = c(2,1)) +
  theme_classic() + scale_color_manual(values = c(Dark24[1],"gray")) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05), color=Dark24[1]) + #geom_hline(yintercept = -log10(0.002645196)) +
  #geom_vline(xintercept = 0.3) + geom_vline(xintercept = -0.3) +
  annotate("text", x=-1, y=-log10(0.05)+.15, label=paste("P=0.05"), size=4, fontface=3,color=Dark24[1]) +
  xlab("Linear Model Slope (coefficient)")

myTopHits9 %>%
  arrange(col000) %>%
  ggplot(., aes(y=-log10(P.Value), x=logFC,
                                   color=col000, size=col000)) +
  geom_point() + scale_size_manual(values = c(1,2,2)) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[1],Dark24[2])) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  #geom_hline(yintercept = -log10(0.05), color=Dark24[1]) + #geom_hline(yintercept = -log10(0.002645196)) +
  #geom_vline(xintercept = 0.3) + geom_vline(xintercept = -0.3) +
  geom_text_repel(aes(label = col00), size = 5, fontface="bold",color="black") +
  #annotate("text", x=-1, y=-log10(0.05)+.15, label=paste("P=0.05"), size=4, fontface=3,color=Dark24[1]) +
  xlab("Linear Model Slope (coefficient)")

myTopHits9 %>%
  arrange(col000002) %>%
  ggplot(., aes(y=-log10(P.Value), x=logFC,
                color=col000002, size=col000002)) +
  geom_point() + scale_size_manual(values = c(1,2,2)) +
  theme_classic() + scale_color_manual(values = c("gray",Dark24[1],Dark24[2])) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  #geom_hline(yintercept = -log10(0.05), color=Dark24[1]) + #geom_hline(yintercept = -log10(0.002645196)) +
  #geom_vline(xintercept = 0.3) + geom_vline(xintercept = -0.3) +
  geom_text_repel(aes(label = col00000), size = 3.5, fontface="bold",color="black",
                  max.overlaps = 30) +
  #annotate("text", x=-1, y=-log10(0.05)+.15, label=paste("P=0.05"), size=4, fontface=3,color=Dark24[1]) +
  xlab("Linear Model Slope (coefficient)") +
  xlim(0,0.7)

myTopHits9 %>%
  arrange(desc(col00002)) %>%
  ggplot(., aes(y=-log10(P.Value), x=logFC,
                color=col00002, size=col00002, alpha=col00002)) +
  geom_point() + scale_size_manual(values = c(2,1)) + scale_alpha_manual(values = c(1,0.5)) +
  theme_classic() + scale_color_manual(values = c(Dark24[1],"gray")) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  geom_hline(yintercept = -log10(0.05), color=Dark24[2]) + #geom_hline(yintercept = -log10(0.002645196)) +
  #geom_vline(xintercept = 0.3) + geom_vline(xintercept = -0.3) +
  geom_text_repel(aes(label = col0000), size = 4, fontface="bold",color="black",
                  max.overlaps = 100) +
  annotate("text", x=0.6, y=-log10(0.05)+.3, label=paste("P=0.05"), size=4, fontface=4,color=Dark24[2]) +
  xlim(0,0.7) +
  xlab("Linear Model Slope (coefficient)")

myTopHits9 %>%
  arrange(desc(colxx)) %>%
  ggplot(., aes(y=-log10(P.Value), x=logFC,
                color=colxx, size=colxx, alpha=colxx)) +
  geom_point() + scale_size_manual(values = c(1,2,2,2)) + scale_alpha_manual(values = c(0.5,1,1,1)) +
  theme_classic() + scale_color_manual(values = c("gray","red",Dark24[13],Dark24[4])) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  #geom_hline(yintercept = -log10(0.05), color=Dark24[17]) + geom_vline(xintercept = 0.2, color=Dark24[17]) +
  #geom_vline(xintercept = 0.3) + geom_vline(xintercept = -0.3) +
  geom_text_repel(aes(label = colx), size = 3.5, fontface=4,
                  max.overlaps = 100) +
  #annotate("text", x=0.05, y=-log10(0.05)+.3, label=paste("P=0.05"), size=3.5, fontface=4,color=Dark24[17]) +
  xlim(0,0.7) +
  xlab("Linear Model Slope (coefficient)")

# Export Supplemental table 2 ----
write_tsv(myTopHits9 %>%
            select(geneSymbol:B), "../../../../../myPapers_submissions/Microbiome/SupplementalMaterial/SupplementalTable2_DEGs_bacterial_burden.txt")

# Plotting the GO figure ----
# From DAVID output, chart (not annotation cluster):
colorheat <- colorRampPalette(colors=c("white",Dark24[1]))(50)
GOsig_genes_up_lessres_16Continuous <- read_delim("GOsig_genes_up_lessres_16Continuous.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE)

GOsig_genes_up_lessres_16Continuous %>%
  slice_head(n = 15) %>%
  mutate(Term = str_remove(Term,".*~")) %>%
  mutate(Term = replace(Term, Term == "homophilic cell adhesion via plasma membrane adhesion molecules", "cellular adhesion")) %>%
#  ggplot(., aes(x=Category, y=Term, size=Count, color=-log10(Benjamini))) +
  ggplot(., aes(x=Category, y=Term, size=Count, color=-log10(Benjamini))) +
  geom_point() +
  theme_classic() +
  scale_colour_gradient(
    low = "#DC91A9",
    high = Dark24[1]) +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank())


GOsig_genes_up_lessres_16Continuous %>%
  slice_head(n = 15) %>%
  mutate(Term = str_remove(Term,".*~")) %>%
  mutate(Term = replace(Term, Term == "homophilic cell adhesion via plasma membrane adhesion molecules", "cellular adhesion")) %>%
  arrange(desc(Benjamini)) %>%
  ggplot(., aes(x=Count, y=fct_inorder(Term), fill=-log10(Benjamini))) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradient(
    low = "#DC91A9",
    high = Dark24[1]) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size=10)) +
  xlab("Gene Count")

GOsig_genes_up_lessres_16Continuous %>%
  slice_head(n = 10) %>%
  mutate(Term = str_remove(Term,".*~")) %>%
  mutate(Term = replace(Term, Term == "homophilic cell adhesion via plasma membrane adhesion molecules", "cellular adhesion")) %>%
  arrange(desc(Benjamini)) %>%
  ggplot(., aes(x=Count, y=fct_inorder(Term), fill=-log10(Benjamini))) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_gradient(
    low = "#DC91A9",
    high = Dark24[1]) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size=10)) +
  xlab("Gene Count")





# Correlation of 16S and gene expression - individual ----
head(myTopHits9 %>% arrange(desc(logFC),P.Value))

bac_variable <- "Bacterial burden"
bacburden_gene_mat <- rbind(log2.cpm.filtered.norm[,CL_samples],
                            qPCR_23$s16_mg[1:51])
rownames(bacburden_gene_mat)[dim(bacburden_gene_mat)[1]] <- "Bacterial_burden"
tail(bacburden_gene_mat)

genes_choice <- c("Bacterial_burden",geneshowA)

bacburden_gene_df <- rownames_to_column(as.data.frame(t(bacburden_gene_mat[genes_choice,])), "sample_RNAseq")
bacburden_gene_df %>%
  pivot_longer(!sample_RNAseq & !Bacterial_burden, names_to = "Feature", values_to = "value") %>%
  left_join(qPCR_23) %>%
  ggplot(., aes(x=value, y=Bacterial_burden)) +
  geom_smooth(method=lm, color=Dark24[1]) +
  geom_point(size=2, color= "#424242") +
  theme_classic() + #scale_color_manual(values = PalC) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11),
        legend.position="right", strip.background = element_blank(), strip.text = element_text(size = 10, face="bold"),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  stat_cor(method = "pearson",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE, color=Dark24[1]) + 
  xlab("Gene expression (log2CPM)") + ylab("16S") +
  facet_wrap(. ~ Feature, scales = "free")

# Export parasite numbers data ----
#write_tsv(qPCR_2 %>%
#          select(sample_RNAseq, cpm_leishmania, n_parasites_mg),
#          "../RStudio_outputs/data/parasiteload_3rdDataset.txt")

# Individual gene expression visualization ----
gene_viz_x <- log2(myDGEList.filtered.norm$counts[c("IL1A","IL1B"),CL_samples])

melt_viz_x <- as.data.frame(melt(gene_viz_x))
melt_viz_x$Var2 <- as.character(melt_viz_x$Var2)
colnames(melt_viz_x)[2] <- "sample_RNAseq"

melt_viz_y <- left_join(melt_viz_x, qPCR_23, by="sample_RNAseq")

ggplot(melt_viz_y %>%
         filter(group_RNAseq == "CL"),
       aes(x=factor(detectable_Saureus_qPCR, levels = c("ND","D")),
                       y=value, color=detectable_Saureus_qPCR)) +
  geom_jitter(position=position_jitter(0.1), size=1) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="red") +
  stat_compare_means(method = "t.test", label.y = 4,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 5, color="purple") +
  xlab("") + ylab(paste("CPM (log2)")) + 
  #ylim(2.5,10) +
  facet_wrap(. ~ Var1, scales = "free")

# Export data for other purposes ----
write_tsv(qPCR_23 %>%
            select(sample_RNAseq,
                   group_RNAseq,
                   healing_time_days,
                   treatment_outcome,
                   treatment_other_drug,
                   cpm_saureus,
                   Saureus_mg,
                   s16_mg,
                   ratio_SaureusTotal,
                   detectable_Saureus_qPCR), "../../RStudio_outputs/data/qPCR_23_filtered.txt")








