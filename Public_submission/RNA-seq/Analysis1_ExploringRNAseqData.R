# In this script, I am starting the analyses in the 3rd lesion biopies bulk RNAseq
# dataset. 

# renv::init() ----
#renv::init()

# Libraries ----
#library(EnsDb.Hsapiens.v86)
library(corrplot)
library(reshape2)
library(biomaRt)
library(tximport)
library(edgeR)
library(patchwork)
library(ggrepel)
library(ggthemes)
library(Polychrome)
library(vegan)
library(FinCal)
library(ggExtra)
library(gplots)
library(ggpubr)
library(Hmisc)
library(plotly)
#library(JLutils)
#library(scales)
library(broom)
library(reshape2)
library(ggforce)
library(msigdbr)
library(immunedeconv)
library(GSEABase)
library(GSVA)
library(gprofiler2)
library(gt)
#library(UpSetR)

library(tidyverse)

Dark24 <- rev(as.character(dark.colors(n=24)))
Dark24 <- Dark24[c(-3,-5)]
save(Dark24, file = "RStudio_outputs/data/Dark24")

# Color Scheme:
library(scales)
show_col(Dark24)

# Importing dataset ----
# Getting sample info and factoring groups:
targets_all <- read_delim("~/Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/mappingFile_LeishMicrobiome_Project.txt", 
                      "\t", escape_double = FALSE, na = "NA", col_types = cols(RNAseq_batch = col_character(),
                                                                               LTCP_patient_ID = col_character()),
                      trim_ws = TRUE)
save(targets_all, file = "RStudio_outputs/data/targets_all")

targets <- targets_all %>%
  filter(group_RNAseq != "NA")

write_tsv(targets, "../../../../Desktop/targets.txt")

#rm(targets_all)
head(targets)

# Sample IDs
CL_samples <- targets %>% dplyr::filter(group_RNAseq == "CL") %>% dplyr::select(sample_RNAseq)
HS_samples <- targets %>% dplyr::filter(group_RNAseq == "HS") %>% dplyr::select(sample_RNAseq)
CL_samples <- CL_samples$sample_RNAseq
HS_samples <- HS_samples$sample_RNAseq

save(CL_samples, file = "RStudio_outputs/data/CL_samples")

# Factoring
groups_factor <- factor(targets$group_RNAseq, levels = c("HS","CL"))
group_order <- c("HS","CL")

# Tximports ----
path1 <- file.path("mapping_kallisto/mapping_to_host/host_kallisto_outputs/",targets$sample_RNAseq, "abundance.h5")
Hs.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
Tx <- getBM(attributes=c('ensembl_transcript_id',
                         'external_gene_name'),
            mart = Hs.anno)

Tx2 <- getBM(attributes=c('external_gene_name',
                         'description',
                         'ensembl_transcript_id'),
            mart = Hs.anno)

Tx <- as_tibble(Tx)

Txi_gene1 <- tximport(path1,
                      type = "kallisto", 
                      tx2gene = Tx, 
                      txOut = FALSE, 
                      countsFromAbundance = "lengthScaledTPM",
                      ignoreTxVersion = TRUE)
#write_tsv(as.data.frame(Txi_gene1$counts), "Txitest.txt")
save(Txi_gene1, file = "RStudio_outputs/data/Txi_gene1")
myDGEList <- DGEList(Txi_gene1$counts)
colnames(myDGEList) <- targets$sample_RNAseq

# Quick part to rename patients (xLTCP) and order them:
GSE214397_Amorim_LeishMicrobiome_SupplemmentalTable1_LeishOmics_StudyDesign <- read_delim("~/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/Amorim_LeishMicrobiome/Public_submission/RNA-seq/GSE214397_Amorim_LeishMicrobiome_SupplemmentalTable1_LeishOmics_StudyDesign.txt", 
                                                                                          delim = "\t", escape_double = FALSE, 
                                                                                          col_types = cols(LTCP_patient_ID = col_character()), 
                                                                                          trim_ws = TRUE)
submission_targets <- targets %>%
  left_join(GSE214397_Amorim_LeishMicrobiome_SupplemmentalTable1_LeishOmics_StudyDesign %>%
              select(LTCP_patient_ID, PatientID)) %>%
  mutate(PatientID = case_when(
    PatientID %in% NA ~ sample_RNAseq,
    TRUE ~ PatientID
  ))

colnames(myDGEList) <- submission_targets$PatientID
Amorim2022_raw_counts <- rownames_to_column(as.data.frame(myDGEList$counts), "geneSymbol")

write_tsv(Amorim2022_raw_counts, "../Amorim_LeishMicrobiome/Public_submission/RNA-seq/Amorim_LeishMicrobiome_raw_counts.txt")
#rm(Txi_gene1)

log2.cpm <- cpm(myDGEList, log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm)
colnames(log2.cpm.df) <- targets$sample_RNAseq
log2.cpm.df.melt <- melt(log2.cpm.df)

p1 <- ggplot(log2.cpm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) + theme_classic() + coord_flip() +
  stat_summary(fun = "median", geom = "point", shape = 124, size = 6, color = "black", show.legend = FALSE)

# Data processing ----
# Filtering
cpm <- cpm(myDGEList)

min_sample = as.vector(table(groups_factor))
num_classes = length(min_sample)
min_size = min_sample[order(min_sample,decreasing=FALSE)[1]]

keepers <- rowSums(cpm>1)>=min_size
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList)
dim(myDGEList.filtered)

# Normalization
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=FALSE)
colnames(log2.cpm.filtered.norm) <- targets$sample_RNAseq
colnames(log2.cpm.filtered.norm) <- submission_targets$PatientID
colnames(myDGEList.filtered.norm) <- targets$sample_RNAseq
write.table(log2.cpm.filtered.norm, "RStudio_outputs/data/log2.cpm.filtered.norm.txt", sep = "\t", quote = FALSE) 

write_tsv(rownames_to_column(as.data.frame(log2.cpm.filtered.norm), "geneSymbol"), "../Amorim_LeishMicrobiome/Public_submission/RNA-seq/Amorim_LeishMicrobiome_processed_norm_filt_log2cpm.txt") 

save(myDGEList.filtered.norm, file = "RStudio_outputs/data/myDGEList.filtered.norm")
save(log2.cpm.filtered.norm, file = "RStudio_outputs/data/log2.cpm.filtered.norm")
save(cpm.filtered.norm, file = "RStudio_outputs/data/cpm.filtered.norm")

log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm)
log2.cpm.filtered.norm.df.melt <- melt(log2.cpm.filtered.norm.df)

p2 <- ggplot(log2.cpm.filtered.norm.df.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(trim = FALSE, show.legend = FALSE) + theme_classic() + coord_flip() +
  stat_summary(fun = "median", geom = "point", shape = 124, size = 6, color = "black", show.legend = FALSE)

preprocessing <- p1+p2
#ggsave("RStudio_outputs/images/dataset_preprocessing.pdf", preprocessing)

# Data not filtered ----
# not filtered, normalized data, coding + non-coding:
notfilt_norm <- calcNormFactors(myDGEList, method = "TMM")
notfilt_norm_log2cpm <- cpm(notfilt_norm, log=TRUE)
notfilt_norm_cpm <- cpm(notfilt_norm, log=FALSE)
write.table(notfilt_norm_log2cpm, "RStudio_outputs/data/notfilt_norm_log2cpm.txt", sep = "\t", quote = FALSE) 
save(notfilt_norm_cpm, file = "RStudio_outputs/data/notfilt_norm_cpm")
save(notfilt_norm_log2cpm, file = "RStudio_outputs/data/notfilt_norm_log2cpm")

# PCA ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df, aes(x=PC1, y=PC2, color=targets$group_RNAseq)) +
  theme_classic() + scale_color_manual(values = Dark24) +
  geom_point(size=2.5) +
#  geom_text_repel(aes(label = targets$sample_RNAseq), size = 2.5, fontface="bold",colour="black") +
  theme(legend.position="top",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  stat_ellipse(level = 0.95) +
  xlab(paste("PC1 -",pc.per[1],"%")) +
  ylab(paste("PC2 -",pc.per[2],"%")) +
  coord_fixed()

# DE analysis CL vs HS ----
# design matrix:
design <- model.matrix(~0 + groups_factor)
colnames(design) <- levels(groups_factor)

# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design)

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix
contrast.matrix <- makeContrasts(disease = CL - HS,
                                 levels=design)

# extract the linear model fit
fits <- contrasts.fit(fit, contrast.matrix)

#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

# MyTopHits and Volcano ----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits <- as_tibble(myTopHits, rownames = "geneSymbol")

# Colunm of significant genes:
sig_genes_up <- myTopHits %>%
  filter(adj.P.Val < 0.05) %>% filter(logFC > 1)
sig_genes_down <- myTopHits %>%
  filter(adj.P.Val < 0.05) %>% filter(logFC < -1)
sig_genes_up <- sig_genes_up$geneSymbol
sig_genes_down <- sig_genes_down$geneSymbol
#write_tsv(as_tibble(sig_genes_up), "up.txt") # to be used in GO
#write_tsv(as_tibble(sig_genes_down), "down.txt") # to be used in GO
save(sig_genes_up, file = "RStudio_outputs/data/sig_genes_up")
save(sig_genes_down, file = "RStudio_outputs/data/sig_genes_down")

myTopHits$col0 <- myTopHits$geneSymbol
col0_s <- myTopHits$col0 %in% sig_genes_up
col0_s2 <- myTopHits$col0 %in% sig_genes_down
myTopHits$col0[!col0_s & !col0_s2] <- NA

myTopHits$col1 <- myTopHits$geneSymbol
col1_s <- myTopHits$col1 %in% sig_genes_up
col1_s2 <- myTopHits$col1 %in% sig_genes_down
myTopHits$col1[col1_s] <- "sig_up"
myTopHits$col1[col1_s2] <- "sig_down"
myTopHits$col1[!col1_s & !col1_s2] <- "notsig"

myTopHits$col2 <- myTopHits$col0
col2_s <- myTopHits$col2 %in% NA
col2_s2 <- myTopHits$col2 %in% sig_genes_up
col2_s3 <- myTopHits$col2 %in% sig_genes_down
myTopHits$col2[col2_s] <- "notsig"
myTopHits$col2[col2_s2 | col2_s3] <- "sig"

myTopHits$col3 <- myTopHits$geneSymbol
col3_v <- myTopHits$col3 %in% c("GZMB","GNLY","PRF1","IL1B")
myTopHits$col3[!col3_v] <- NA

ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC,
                      color=col1,size=col1)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", Dark24[2], Dark24[1])) +
  scale_size_manual(values = c(0.5,1,1)) +
  #geom_text_repel(aes(label = col3), size = 5, fontface="bold",color="black", nudge_y = 5) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  annotate("text", x=-5, y=32, label=paste(length(sig_genes_down)), size=5, color=Dark24[2], fontface="bold") +
  annotate("text", x=5, y=32, label=paste(length(sig_genes_up)), size=5, color=Dark24[1], fontface="bold") +
  xlab("logFC CL vs. HS")

ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC,
                      color=col1,size=col1)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("dark gray", "dark gray", Dark24[1])) +
  scale_size_manual(values = c(0.5,1,1)) +
  #geom_text_repel(aes(label = col3), size = 5, fontface="bold",color="black", nudge_y = 5) +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  annotate("text", x=5, y=32, label=paste(length(sig_genes_up)), size=5, color=Dark24[1], fontface="bold") +
  xlab("logFC CL vs. HS")

write_tsv(myTopHits, "RStudio_outputs/data/myTopHits.txt")
save("myTopHits", file = "RStudio_outputs/data/myTopHits")


# GSEA CL vs HS - table ----
gsea_report_for_CL_1613579408465 <- read_delim("GSEA_GO_outputs/GSEA/GSEA_KEGG_phenotype_CLvsHS.Gsea.1613579408465/gsea_report_for_CL_1613579408465.tsv",
                                               "\t", escape_double = FALSE, trim_ws = TRUE)

gsea_report_for_CL_1613579408465 %>%
  select(NAME, NES, `NOM p-val`) %>%
  slice(1:10) %>%
  mutate_at(vars(NES, `NOM p-val`), funs(round(., 2))) %>%
  rename("Top 10 enriched genesets in CL (GSEA)" = NAME) %>%
  gt()

# renv::snapshot() ----
#renv::snapshot()
