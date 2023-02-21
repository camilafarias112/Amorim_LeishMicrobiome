# In this script I import variables from the 16S-seq analysis done by me, Camila.
# Find Rscripts in: ~/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/microbiome_Lesion16S-seqDataset

# Libraries ----
library(ggforce)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggdist)
library(gplots)
library(patchwork)
library(tidyverse)

# PCA to investigate 16S classification in the RNA-seq dataset ----
# Here I have straight forward  total gene expression vs Microbiome clusters. The results won't show anything super incredible.
# That's why I chose to reduce the dimensionality of the gene expression and rexposome.
# vector with swab day0 samples:
CL_samples_16Sseq <- phenotype %>% filter(microbiome_cluster != "NA")
CL_samples_16Sseq <- CL_samples_16Sseq$sample_RNAseq

CL_samples_16Sseq_HS <- c(CL_samples_16Sseq, "HS1","HS2","HS3","HS4","HS5","HS6")

load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/log2.cpm.filtered.norm")
pca.res_2 <- prcomp(t(log2.cpm.filtered.norm[,colnames(log2.cpm.filtered.norm) %in% CL_samples_16Sseq]),
                    scale.=F, retx=T)
pca.res_2 <- prcomp(t(log2.cpm.filtered.norm[,colnames(log2.cpm.filtered.norm) %in% CL_samples_16Sseq_HS]),
                    scale.=F, retx=T)
CL_samples_16Sseq <- colnames(log2.cpm.filtered.norm[,colnames(log2.cpm.filtered.norm) %in% CL_samples_16Sseq]) # Just to capture the names

get_eigenvalue(pca.res_2)
fviz_eig(pca.res_2, addlabels = TRUE, ylim = c(0, 20))
var <- get_pca_var(pca.res_2)

pc.var_2<-pca.res_2$sdev^2
pc.per_2<-round(pc.var_2/sum(pc.var_2)*100, 1)

pca.res.df_2 <- as_tibble(pca.res_2$x)
rownames(pca.res.df_2) <- CL_samples_16Sseq_HS
pca.res.df_2 <- pca.res.df_2 %>% add_column(sample_RNAseq = CL_samples_16Sseq_HS)
pca.res.df_2 <- pca.res.df_2 %>% add_column(sample_RNAseq = CL_samples_16Sseq)

pca.res.df_2 %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M4",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  ggplot(., aes(x=PC1, y=PC2, color = microbiome_cluster2)) +
  geom_mark_hull(aes(fill = microbiome_cluster2), alpha = 0.1, show.legend = F, expand = unit(2.5, "mm")) +
  geom_point(size=2.5) +
  theme_classic() + scale_color_manual(values = c("black","gray",
                                                  Dark24[5],
                                                  Dark24[6],
                                                  Dark24[3],
                                                  Dark24[2],
                                                  Dark24[1],
                                                  Dark24[4])) +
  #geom_text_repel(aes(label = targets_w_Gricedata$LTCP_patient_ID), size = 2, fontface="bold",colour="black") +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.text = element_text(size = 15)) +
  xlab(paste("PC1 -",pc.per_2[1],"%")) + ylab(paste("PC2 -",pc.per_2[2],"%")) +
  coord_fixed()

library("gg3D")
pca.res.df_2 %>%
  left_join(phenotype) %>%
  mutate(microbiome_cluster2 = case_when(
    microbiome_cluster == "Staphylococcus" ~ "M6",
    microbiome_cluster == "Arcanobacterium" ~ "M7",
    microbiome_cluster == "Streptococcus" ~ "M3",
    microbiome_cluster == "Heterogeneous" ~ "M0",
    microbiome_cluster == "Corynebacterium" ~ "M5",
    microbiome_cluster == "Finegoldia" ~ "M2",
    microbiome_cluster == "Lactobacillus" ~ "M1",
    TRUE ~ "HS")) %>%
  ggplot(., aes(x=PC1, y=PC2, z=PC3, color = microbiome_cluster2, size=5)) +
  scale_color_manual(values = c("black",
                                  "gray",
                                  Dark24[5],
                                  Dark24[6],
                                  Dark24[3],
                                  Dark24[2],
                                  Dark24[1],
                                  Dark24[4])) +
  theme_void() +
  axes_3D() +
  stat_3D()

a <- princomp(log2.cpm.filtered.norm[,colnames(log2.cpm.filtered.norm) %in% CL_samples_16Sseq_HS])
#biplot(a)

# PCA wit Staphylococcus and Heterogeneous ----
samples_pick0 <- phenotype %>%
  filter(sample_RNAseq %in% CL_samples_16Sseq) %>%
  filter(microbiome_cluster %in% c("Staphylococcus",
                                   "Heterogeneous")) %>%
  arrange(microbiome_cluster)
staphyfactor <- factor(samples_pick0$microbiome_cluster)
samples_pick <- samples_pick0$sample_RNAseq

pca.res_3 <- prcomp(t(log2.cpm.filtered.norm[,samples_pick]), scale.=F, retx=T)
pc.var_3 <- pca.res_3$sdev^2
pc.per_3 <- round(pc.var_3/sum(pc.var_3)*100, 1)

pca.res.df_3 <- as_tibble(pca.res_3$x)
rownames(pca.res.df_3) <- samples_pick
pca.res.df_3 <- pca.res.df_3 %>% add_column(sample_RNAseq = samples_pick)

pca.res.df_3 %>%
  left_join(phenotype) %>%
  ggplot(., aes(x=PC1, y=PC2, color = microbiome_cluster)) +
  geom_mark_hull(aes(fill = microbiome_cluster), alpha = 0.1, show.legend = F, expand = unit(2.5, "mm")) +
  geom_point(size=2.5) +
  theme_classic() + scale_color_manual(values = c(Dark24[1],
                                                  "gray",
                                                  Dark24[4],
                                                  Dark24[2],
                                                  Dark24[3],
                                                  Dark24[5],
                                                  Dark24[6])) +
  geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),axis.text = element_text(size = 15)) +
  xlab(paste("PC1 -",pc.per_3[1],"%")) + ylab(paste("PC2 -",pc.per_3[2],"%")) +
  coord_fixed() #results not super useful

# DE Staph dysbiosis vs. Heterogeneous ----
table(phenotype$microbiome_cluster) # number of samples for this comparison
factor(phenotype$microbiome_cluster) # number of samples for this comparison

phenotype_bothData <- phenotype %>%
  filter(sample_RNAseq %in% CL_samples_16Sseq)
  
design1 <- model.matrix(~0 + factor(phenotype_bothData$microbiome_cluster)) 
colnames(design1) <- levels(factor(phenotype_bothData$microbiome_cluster))

library(limma)
library(edgeR)

load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/myDGEList.filtered.norm")
v.DEGList.filtered.norm1 <- voom(myDGEList.filtered.norm[,phenotype_bothData$sample_RNAseq], design1)
fit1 <- lmFit(v.DEGList.filtered.norm1, design1)
contrast.matrix1 <- makeContrasts(comparison = Staphylococcus - Heterogeneous,
                                 levels=design1)
fits1 <- contrasts.fit(fit1, contrast.matrix1)
ebFit1 <- eBayes(fits1)
myTopHits1 <- topTable(ebFit1, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits1 <- as_tibble(myTopHits1, rownames = "geneSymbol")

sig_genes_up1 <- myTopHits1 %>%
  filter(P.Value < 0.05) %>% filter(logFC > 1)
sig_genes_down1 <- myTopHits1 %>%
  filter(P.Value < 0.05) %>% filter(logFC < -1)
sig_genes_up1 <- sig_genes_up1$geneSymbol
sig_genes_down1 <- sig_genes_down1$geneSymbol

write_tsv(as_tibble(sig_genes_up1), "outputs/StaphyDysvsHetero_up.txt")
write_tsv(as_tibble(sig_genes_down1), "outputs/StaphyDysvsHetero_down.txt")

myTopHits1$col0 <- myTopHits1$geneSymbol
col0_v1 <- myTopHits1$col0 %in% sig_genes_up1
col0_v21 <- myTopHits1$col0 %in% sig_genes_down1
myTopHits1$col0[!col0_v1 & !col0_v21] <- NA

# Add cluster 1 from GO annotations:
StaphyDysvsHetero_up_GO_cluster1 <- read_csv("outputs/StaphyDysvsHetero_up_GO_cluster1.txt", 
                                             col_names = FALSE)
StaphyDysvsHetero_up_GO_cluster1 <- unique(StaphyDysvsHetero_up_GO_cluster1$X1)

StaphyDysvsHetero_up_GO_cluster2 <- read_csv("outputs/StaphyDysvsHetero_up_GO_cluster2.txt", 
                                             col_names = FALSE)
StaphyDysvsHetero_up_GO_cluster2 <- unique(StaphyDysvsHetero_up_GO_cluster2$X1)

venn_clus <- venn(list(Cluster1 = StaphyDysvsHetero_up_GO_cluster1,
                       Cluster2 = StaphyDysvsHetero_up_GO_cluster2))

myTopHits1$col1 <- myTopHits1$geneSymbol
col1_v1 <- myTopHits1$col1 %in% attr(venn_clus,"intersections")$Cluster1
col1_v21 <- myTopHits1$col1 %in% attr(venn_clus,"intersections")$Cluster2
col1_v211 <- myTopHits1$col1 %in% attr(venn_clus,"intersections")$`Cluster1:Cluster2`
myTopHits1$col1[col1_v1] <- "GO:Cluster 1"
myTopHits1$col1[col1_v21] <- "GO:Cluster 2"
myTopHits1$col1[col1_v211] <- "GO:Cluster 1 and 2"
myTopHits1$col1[!col1_v1 & !col1_v21 & !col1_v211] <- "Others"

myTopHits1$col2 <- myTopHits1$geneSymbol
col2_v1 <- myTopHits1$col2 %in% attr(venn_clus,"intersections")$Cluster1
col2_v21 <- myTopHits1$col2 %in% attr(venn_clus,"intersections")$Cluster2
col2_v211 <- myTopHits1$col2 %in% attr(venn_clus,"intersections")$`Cluster1:Cluster2`
myTopHits1$col2[col2_v1] <- "GO:Cluster 1"
myTopHits1$col2[col2_v21] <- "GO:Cluster 2"
myTopHits1$col2[col2_v211] <- "GO:Cluster 1 and 2"
myTopHits1$col2[!col2_v1 & !col2_v21 & !col2_v211] <- NA

# Pie of overrepresented:
myTopHits1_df1 <- as.data.frame(myTopHits1 %>%
                                  filter(P.Value < 0.05) %>% filter(logFC > 1))
myTopHits1_df1 <- as.data.frame(table(factor(myTopHits1_df1$col1)))

Goclus <- factor(myTopHits1_df1$Var1)

myTopHits1_df1 %>%
  ggplot(., aes (x="", y = Freq, fill = factor(Var1))) + 
  theme_void() +
  geom_col(position = 'stack', width = 1) +
  #geom_text(aes(label = Freq, x = 1.3),
  #          position = position_stack(vjust = 0.5),
  #          size=6, color= "white", fontface=2) +
  scale_fill_manual(values = c(Dark24[1],
                               Dark24[2],
                               Dark24[3],
                               Dark24[4])) +
  theme(legend.position = "right") +
  coord_polar("y")

ggplot(myTopHits1, aes(y=-log10(P.Value), x=logFC, color=col1)) +
  geom_point() +
  theme_classic() + scale_color_manual(values = rev(Dark24)) +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  #geom_text_repel(aes(label = col0), size = 3, fontface="bold",color="black",
  #                max.overlaps = 500) +
  xlab("logFC Staphylococcus vs. Heterogeneous")

ggplot(myTopHits1, aes(y=-log10(P.Value), x=logFC, color=col1)) +
  geom_point() +
  theme_classic() +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  #geom_text_repel(aes(label = col0), size = 3, fontface="bold",color="black",
  #                max.overlaps = 500) +
  xlim(1,4) + ylim(-log10(0.05), 4.5) +
  xlab("logFC Staphylococcus vs. Heterogeneous")

ggplot(myTopHits1, aes(y=-log10(P.Value), x=logFC)) +
  geom_point() +
  theme_classic() + scale_color_manual(values = rev(Dark24)) +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05)) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  geom_text_repel(aes(label = col0), size = 3, fontface="bold",color="black",
                  max.overlaps = 500) +
  xlab("logFC Staphylococcus vs. Heterogeneous")

# Export some objects for Tori ----
# Tori asked me Robjects to look for Wound healing signatures associated with Staphylococcus Dysbiosis 02/15/22.
# This chunk of script was originally written in ~/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/Dataset_Integration/Analysis3_MergingMicrobiomeProfile_RNAseq.R

# Metadata with 12 patients Staphy-dysbiosis and 6 Heterogeneous. Patients included here needed to have BOTH RNA-seq and 16S-seq samples.
metadata_StaphyFocus <- phenotype_bothData %>%
  filter(microbiome_cluster %in% c("Staphylococcus", "Heterogeneous")) %>%
  arrange(microbiome_cluster)
save(metadata_StaphyFocus, file = "../../../../Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/SharedwithTori/StaphyDysbiosisFocus/metadata_StaphyFocus")

# Not filtered, normalized gene expression
metadata_StaphyFocus_samples <- metadata_StaphyFocus$sample_RNAseq
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/notfilt_norm_log2cpm")
StaphyFocus_samples_geneExp_notFilt <- notfilt_norm_log2cpm[,metadata_StaphyFocus_samples]
save(StaphyFocus_samples_geneExp_notFilt, file = "../../../../Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/SharedwithTori/StaphyDysbiosisFocus/StaphyFocus_samples_geneExp_notFilt")

# Filtered, normalized gene expression
StaphyFocus_samples_geneExp <- myDGEList.filtered.norm[,metadata_StaphyFocus_samples]
StaphyFocus_samples_geneExp <- StaphyFocus_samples_geneExp$counts
save(StaphyFocus_samples_geneExp, file = "../../../../Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/SharedwithTori/StaphyDysbiosisFocus/StaphyFocus_samples_geneExp")

# myTopHits from DGE analysis Staphy-dysbiosis vs. Heterogeneous
save(myTopHits1, file = "../../../../Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/SharedwithTori/StaphyDysbiosisFocus/myTopHits1")

# GO Staphy vs Heterogenous table (up genes) ----
StaphyDysvsHetero_up_GO_mod <- read_delim("outputs/StaphyDysvsHetero_up_GO_mod.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
StaphyDysvsHetero_up_GO_mod$Term <- str_sub(StaphyDysvsHetero_up_GO_mod$Term, start = 12)


StaphyDysvsHetero_up_GO_mod %>%
  mutate(cluster2 = case_when(
    cluster %in% "cluster1" ~ "Cluster 1",
    cluster %in% "cluster2" ~ "Cluster 2")) %>%
  ggplot(., aes(x=cluster2, y=-log10(FDR), size=`Fold Enrichment`, color = cluster2)) +
  geom_point() +
  theme_classic() + scale_color_manual(values = c(Dark24[1:2])) +
  theme(axis.title.x = element_blank(), axis.text = element_text(size=11)) +
  geom_text_repel(aes(label = Term), size = 2.5, fontface="bold",color="black",
                                  max.overlaps = 500) +
  geom_hline(yintercept = -log10(0.05))


StaphyDysvsHetero_up_GO_mod %>%
  mutate(cluster2 = case_when(
    cluster %in% "cluster1" ~ "Cluster 1",
    cluster %in% "cluster2" ~ "Cluster 2")) %>%
ggplot(., aes(Category, factor(Term, levels = rev(StaphyDysvsHetero_up_GO_mod$Term)))) +
  geom_tile(aes(fill = Count)) +
  theme_classic() + scale_fill_gradient(low = "white", high = Dark24[15])
  





# Some gene ontology here:
library(gprofiler2)
gost.res <- gost(sig_genes_up1, organism = "hsapiens", correction_method = "fdr")
mygostplot <- gostplot(gost.res, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

# Heatmap and ggdist with DEG Staphy vs Heterogeneous - Heatmap looked too much, so I didn't continue ----
Staphy_samples2 <- phenotype_bothData %>%
  filter(microbiome_cluster == "Staphylococcus")
Staphy_samples2 <- Staphy_samples2$sample_RNAseq
Hetero_samples2 <- phenotype_bothData %>%
  filter(microbiome_cluster == "Heterogenous")
Hetero_samples2 <- Hetero_samples2$sample_RNAseq

StaphyHete_mtDEG <- log2.cpm.filtered.norm[attr(venn_clus,"intersections")$Cluster1, samples_pick0$sample_RNAseq]

# Only "GO - Inflammatory response from Annotation CLuster 1:
infla_geneClust1_sig <- c("GBP5", "TNFRSF6B", "TNFRSF18", "CCL4L2", "SLC11A1", "IL24", "ACOD1", "CXCL1", "IL27",
  "FPR2", "CXCL10", "CXCL11", "CCL8", "CCL7", "IL1B", "TNIP3", "CCL4", "CCL3", "CCL2")
StaphyHete_mtDEG <- log2.cpm.filtered.norm[infla_geneClust1_sig, samples_pick0$sample_RNAseq]
write_tsv(rownames_to_column(as.data.frame(StaphyHete_mtDEG), "geneSymbol"), "outputs/Staphy_vs_Hetero.txt")

# Cluster samples
dist_2 <- dist(t(StaphyHete_mtDEG), method = "canberra") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clus2 <- hclust(dist_2, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clus2)

# Cluster genes
dist_3 <- dist(StaphyHete_mtDEG, method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clus3 <- hclust(dist_3, method = "ward.D2") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clus3)

staphy_colors <- function(staphyfactor){ if (staphyfactor== "Staphylococcus") Dark24[1] else if (staphyfactor== "Heterogenous") Dark24[2]}
staphy_colors <- unlist(lapply(staphyfactor, staphy_colors))

heatmap.2(StaphyHete_mtDEG,
          Colv=as.dendrogram(clus2),
          Rowv = as.dendrogram(clus3),
          dendrogram = "row",
          #RowSideColors=module.color,
          col=myheatcol,
          ColSideColors = staphy_colors,
          scale='row',
          density.info="none", trace="none",
          cexRow=1, cexCol=1,margins = c(5,15))

# ggdist:
library(reshape)
melt_viz_dist <- as.data.frame(melt(StaphyHete_mtDEG))
melt_viz_dist$X2 <- as.character(melt_viz_dist$X2)
colnames(melt_viz_dist)[2] <- "sample_RNAseq"
melt_viz_dist <- left_join(melt_viz_dist, phenotype, by="sample_RNAseq")

# Individual gene
melt_viz_dist %>%
  ggplot(.,
         aes(x=factor(microbiome_cluster, levels = c("Heterogeneous", "Staphylococcus")),
                      y=value, fill=microbiome_cluster, color=microbiome_cluster)) +
  stat_gradientinterval(slab_type="pdf") + # probability density (or mass) function("pdf"), cumulative distribution function ("cdf"), complementary CDF ("ccdf"), or histogram ("histogram"
  theme_classic() + scale_fill_manual(values = c(Dark24[1],"dark gray")) + scale_color_manual(values = c(Dark24[1],"dark gray")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, angle = 90), axis.title = element_text(size = 13),
        legend.position="none", strip.background = element_blank(),
        strip.text = element_text(size = 13, angle = 90, hjust = 0)) +
  #geom_text_repel(aes(label = sample_RNAseq), size = 2, fontface="bold",colour="black") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  xlab("") + ylab(paste("CPM (log2)")) + 
  #stat_compare_means(method = "t.test",# label.y = 8.2,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.5, size = 7, color="black") +
  #  ylim(0,5) +
  # facet_grid(. ~ Var1) +
  facet_wrap(. ~ X1, scales = "free", nrow = 3)

# Fold change bar plots:
fc_desc <-  myTopHits1 %>%
  filter(geneSymbol %in% infla_geneClust1_sig) %>%
  arrange(desc(logFC))
fc_desc <- fc_desc$geneSymbol

myTopHits1 %>%
  filter(geneSymbol %in% infla_geneClust1_sig) %>%
  ggplot(., aes(x=factor(geneSymbol, levels = rev(fc_desc)), y=2^logFC)) +
  geom_bar(stat = "identity", fill= Dark24[9]) +
  theme_classic() + 
  theme(panel.grid.minor = element_line(),panel.grid.major = element_line(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 15), axis.title = element_blank(),
        legend.position="none") +
  ylim(0,13) + 
  geom_text(aes(label=geneSymbol), size=5, vjust = 0.5, hjust=-0.1) +
  coord_flip()


# Network of DEG Staphy vs Hetero (not informative) ----
library(ggraph)
library(corrr)
library(tidygraph)

res.cor <- t(StaphyHete_mtDEG[sig_genes_up1,]) %>% correlate(method = "pearson") %>%
  shave(upper = TRUE) %>%
  stretch(na.rm = TRUE) %>%
  filter(r >= 0.5 | r <= -0.55)

cor.graph <- as_tbl_graph(res.cor, directed = FALSE)

#'%ni%' <- Negate('%in%')
cor.graph <- cor.graph %>% 
  activate(nodes) %>% 
  rename(Gene = name) %>%
  #left_join(description, by = "Gene") %>%
  #mutate(impgenes = case_when(
  #  Gene %in% impgenes ~ "impgenes",
  #  Gene %ni% impgenes ~ "notimpgenes")) %>%
  activate(edges) %>%
  arrange(desc(r)) #%>%
#mutate(direction = case_when( #there aren't negative correlations with spearman < -.70
#  r > 0 ~ "positive",
#  r < 0 ~ "negative"))

cor.graph
cor.graph %>%
  ggraph(., layout = "kk") +
  geom_edge_density(aes(fill = r)) + 
  geom_edge_link(aes(width = r,edge_colour = r)) +
  #scale_edge_color_manual(values = c("blue", "red")) +
  #scale_edge_color_gradient(low = "gray", high = Salvador[6]) +
  scale_edge_width(range = c(0.2, 1.3)) +
  #geom_node_point(aes(color = family, shape = impgenes), size = 3) +
  #geom_node_point(aes(color = Family), size = 3) +
  scale_shape_manual(values = c(15,16)) +
  #scale_color_manual(values = Salvador) +
  geom_node_text(aes(label = Gene), size = 4, fontface="bold",
                 color="black", max.overlaps = 20, repel = TRUE) +
  theme_graph()











# Clinical Metadata and qPCR Staph and Strep in swabs ----
ggplot(targets_w_Gricedata, aes(x=qPCR_swab_Staph_absolute_curve, y=size_lesion_mm2)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
ggplot(targets_w_Gricedata, aes(x=qPCR_swab_Strep_absolute_curve, y=size_lesion_mm2)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE)

ggplot(targets_w_Gricedata %>%
         filter(treatment_outcome == "cure"| treatment_outcome == "failure"), aes(x=qPCR_swab_Staph_absolute_curve, y=healing_time_days)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE)
  ggplot(targets_w_Gricedata %>%
           filter(treatment_outcome == "cure"| treatment_outcome == "failure"), aes(x=qPCR_swab_Strep_absolute_curve, y=healing_time_days)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE)

ggplot(targets_w_Gricedata %>%
         filter(treatment_outcome == "cure"| treatment_outcome == "failure"), aes(x=qPCR_swab_Staph_absolute_curve, y=DTH_mm2)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
ggplot(targets_w_Gricedata, aes(x=qPCR_swab_Strep_absolute_curve, y=DTH_mm2)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE)

ggplot(targets_w_Gricedata, aes(x=qPCR_swab_Staph_absolute_curve, y=parasite_cpm)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
  ggplot(targets_w_Gricedata, aes(x=qPCR_swab_Strep_absolute_curve, y=parasite_cpm)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE)

ggplot(targets_w_Gricedata, aes(x=qPCR_swab_Staph_absolute_curve, y=age)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +
  ggplot(targets_w_Gricedata, aes(x=qPCR_swab_Strep_absolute_curve, y=age)) +
  geom_smooth(method=lm, color=Dark24[10]) +
  geom_point(size=2) +
  theme_classic() + scale_color_manual(values = Dark24) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5), 
        axis.text.y = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position="none") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = "left", label.y.npc = "top",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE)

ggplot(targets_w_Gricedata %>%
         filter(treatment_outcome == "cure"| treatment_outcome == "failure"), aes(x=treatment_outcome, y=qPCR_swab_Staph_absolute_curve,
                     color=treatment_outcome)) +
  geom_jitter(position=position_jitter(0.1), size=2) +
  theme_classic() + scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, vjust = 0.5), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="none") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", #label.y = 8.2,
                     aes(label = ..p.signif..),
                     label.x = 1.5, size = 5, color="black")

ggplot(targets_w_Gricedata, aes(x=sex, fill=staphover50)) + geom_bar() +
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(legend.position="right",legend.text = element_text(size = 15),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15)) +
  geom_text(stat="count", aes(label=..count..), position = position_stack(vjust = 0.5), size = 8, color="white") +
  ggplot(targets_w_Gricedata, aes(x=staphover50, fill=sex)) + geom_bar() +
  theme_classic() + scale_fill_manual(values = Dark24) +
  theme(legend.position="right",legend.text = element_text(size = 15),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15)) +
  geom_text(stat="count", aes(label=..count..),position = position_stack(vjust = 0.5), size = 8, color="white")

#write_tsv(mod_targets_w_Gricedata, "mod_targets_w_Gricedata.txt")

# 16S plot R object by Jordan ----
load("~/Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/FilesSentbyJordan_16Splot/Abund_Genus_Anaerobes")
load("~/Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/FilesSentbyJordan_16Splot/Desired_order_subject")
load("~/Box Sync/human_lesion_microbiome/Camila_MasterMetadata_DONOTedit/FilesSentbyJordan_16Splot/Unique_Genus_test")

# All 64 samples:
plot_genus <- ggplot(Abund_Genus_Anaerobes, aes(x=SubjectID, y=Abundance, 
                                               fill=factor(Genus, levels=c(Unique_Genus_test,"Staphylococcus")))) +
  geom_bar(stat="identity") +  
  scale_fill_manual(values=c(Dark24[2:15],Dark24[1])) + 
  theme(axis.text.x=element_text(angle = 90, hjust = 0, size = 8), 
        strip.text.x = element_text(size = 12, face="bold"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  guides(fill=guide_legend(title="Genus", ncol = 1)) +
  ggtitle("16S abundance plot") +
  ylab("16s Abundance")
plot_genus$data$SubjectID <- factor(plot_genus$data$SubjectID, levels = Desired_order_subject)
plot_genus

# Only 51 samples with matching RNAseq: 
Desired_order_subject_filt <- Desired_order_subject[Desired_order_subject %in% samples_16S_RNAseq]
plot_genus_filt <- ggplot(Abund_Genus_Anaerobes %>% filter(SubjectID %in% samples_16S_RNAseq),
                     aes(x=SubjectID, y=Abundance, 
                     fill=factor(Genus, levels=c(Unique_Genus_test,"Staphylococcus")))) +
  geom_bar(stat="identity") + 
  theme_classic() +
  scale_fill_manual(values=c(Dark24[2:14],Dark24[1])) + 
  theme(axis.text.x=element_text(angle = 90, vjust = .5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size=10)) +
  guides(fill=guide_legend(title="Genus", ncol = 1)) +
  #ggtitle("16S abundance plot") +
  ylab("16s Abundance")
plot_genus_filt$data$SubjectID <- factor(plot_genus_filt$data$SubjectID, levels = Desired_order_subject_filt)
plot_genus_filt

plot_genus / 
  plot_genus_filt

colnames(Abund_Genus_Anaerobes)
Abund_Genus_Anaerobes$Family

View(Abund_Genus_Anaerobes)
View(Abund_Genus_Anaerobes %>% filter(Genus == "Staphylococcus"))
dim(Abund_Genus_Anaerobes %>% filter(Genus == "Staphylococcus"))

# getting the raw data from the 16S plot and adding to my targets_w_Grice
collecting16Sabundance_Staphylococcus <- Abund_Genus_Anaerobes %>% filter(Genus == "Staphylococcus") %>%
  select(Abundance, SubjectID) %>%
  rename(Abundance16S_Staphylococcus = Abundance,
         LTCP_patient_ID = SubjectID) %>%
  mutate(LTCP_patient_ID = as.character(LTCP_patient_ID))

collecting16Sabundance_Streptococcus <- Abund_Genus_Anaerobes %>% filter(Genus == "Streptococcus") %>%
  select(Abundance, SubjectID) %>%
  rename(Abundance16S_Streptococcus = Abundance,
         LTCP_patient_ID = SubjectID) %>%
  mutate(LTCP_patient_ID = as.character(LTCP_patient_ID))

targets_w_Gricedata <- left_join(targets_w_Gricedata, collecting16Sabundance_Staphylococcus[,c("Abundance16S_Staphylococcus","LTCP_patient_ID")],
                                 by="LTCP_patient_ID")
targets_w_Gricedata <- left_join(targets_w_Gricedata, collecting16Sabundance_Streptococcus[,c("Abundance16S_Streptococcus","LTCP_patient_ID")],
                                 by="LTCP_patient_ID")

# I had to replace NAs for zeros to include all the samples in the plots.
# Regarding 16S abundance in 16S plot, for the samples that I replaced NAs for zeros,
# these samples simply did not pass Jordan's cutoff to be included in 16S plots.
targets_w_Gricedata <- targets_w_Gricedata %>%
  replace_na(list(Abundance16S_Staphylococcus = 0,
                  qPCR_swab_Staph_absolute_curve = 0))

targets_w_Gricedata$Abundance16S_Staphylococcus
targets_w_Gricedata$qPCR_swab_Staph_absolute_curve

# DE 16S abundance plot Staphylococcus (continuous) ---- 
design7 <- model.matrix(~targets_w_Gricedata$Abundance16S_Staphylococcus) 
colnames(design7)[2] <- "abundance16SplotStaphylococcus"
head(design7)
dim(design7)
v.DEGList.filtered.norm7 <- voom(myDGEList.filtered.norm[,CL_samples], design7)
fit7 <- lmFit(v.DEGList.filtered.norm7, design7)
fits7 <- contrasts.fit(fit7, coefficients = 2)
ebFit7 <- eBayes(fits7)
myTopHits7 <- topTable(ebFit7, adjust ="BH", number=40000, sort.by="logFC")
myTopHits7 <- as_tibble(myTopHits7, rownames = "geneSymbol")

#sig_genes_up7 <- myTopHits7 %>%
#  filter(adj.P.Val < 0.05)
#sig_genes_down7 <- myTopHits7 %>%
#  filter(adj.P.Val < 0.05)

sig_genes_up7 <- myTopHits7 %>%
  filter(P.Value < 0.05, logFC > 0)
sig_genes_down7 <- myTopHits7 %>%
  filter(P.Value < 0.05, logFC < 0)

sig_genes_up7 <- sig_genes_up7$geneSymbol
sig_genes_down7 <- sig_genes_down7$geneSymbol
#write_tsv(sig_genes_up7, "RStudio_outputs/data/sig_genes_up7.txt")
#write_tsv(sig_genes_down7, "RStudio_outputs/data/sig_genes_down7.txt")

myTopHits7$col0 <- myTopHits7$geneSymbol
col0_v7 <- myTopHits7$col0 %in% sig_genes_up7
col0_v77 <- myTopHits7$col0 %in% sig_genes_down7
myTopHits7$col0[!col0_v7 & !col0_v77] <- NA

#volcano:
ggplot(myTopHits7, aes(y=-log10(P.Value), x=logFC)) +
  geom_point(size=1, color="gray") +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  geom_hline(yintercept = -log10(0.05), color="red") +
  geom_hline(yintercept = -log10(0.01), color="blue") +
  #annotate("text", x=0, y=-log10(0.05)+0.1, label="P.value=0.05", color="red") +
  #annotate("text", x=0, y=-log10(0.01)+0.1, label="P.value=0.01", color="blue") +
  #geom_vline(xintercept = 0.035) + geom_vline(xintercept = -0.035) +
  #geom_text_repel(aes(label = col0), size = 1.8, fontface="bold",color="black") +
  xlab("logFC")

# Venn Interception between S. aureus measurements in swabs ----
venn1 <- venn(list(s16Sabundance = sig_genes_up7[1:150],
                   qPCRrelative = sig_genes_up4[1:150],
                   qPCRabsolute = sig_genes_up4.4[1:150]))

venn2 <- venn(list(s16Sabundance = sig_genes_up7[1:150],
                   qPCRrelative = sig_genes_up4[1:150],
                   qPCRabsolute = sig_genes_up4.4[1:150],
                   qPCRtotal16S = sig_genes_up6[1:150]))

venn3 <- venn(list(qPCRabsolute = sig_genes_up4.4[1:150],
                   SaureusCPM = sig_genes_up8[1:150]))

genesup16Stotal <- myTopHits6 %>%
  filter(geneSymbol %in% attr(venn2,"intersections")$qPCRtotal16S)
genesup16Stotal <- genesup16Stotal$geneSymbol



# Plotting qPCR swab Staphylococcus relative to total 16S ----
plot1 <- ggplot(targets_w_Gricedata %>% filter(LTCP_patient_ID %in% samples_16S_RNAseq),
       aes(x=factor(LTCP_patient_ID, levels=Desired_order_subject_filt), y=Abundance16S_Staphylococcus)) +
  geom_bar(stat="identity", fill=Dark24[1]) + 
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text()) +
  guides(fill=guide_legend(title="Genus", ncol = 1)) +
  #ggtitle("16S abundance plot") +
  ylab("16s Abundance")
plot2<- ggplot(targets_w_Gricedata %>% filter(LTCP_patient_ID %in% samples_16S_RNAseq),
       aes(x=factor(LTCP_patient_ID, levels=Desired_order_subject_filt), y=qPCR_swab_Staph_relative16S)) +
  geom_bar(stat="identity", fill=Dark24[1]) + 
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text()) +
  guides(fill=guide_legend(title="Genus", ncol = 1)) +
  #ggtitle("16S abundance plot") +
  ylab("Relative abundance (qPCR)")
plot3<- ggplot(targets_w_Gricedata %>% filter(LTCP_patient_ID %in% samples_16S_RNAseq),
               aes(x=factor(LTCP_patient_ID, levels=Desired_order_subject_filt), y=cpm_saureus)) +
  geom_bar(stat="identity", fill=Dark24[1]) + 
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text()) +
  guides(fill=guide_legend(title="Genus", ncol = 1)) +
  #ggtitle("16S abundance plot") +
  ylab("S. aureus (CPM)")

plot1 /plot2 /plot3











