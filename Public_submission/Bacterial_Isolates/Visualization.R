# In this script I make a table with the information of bacteria isolated from the lesions
# This table was made by Tori Lovins, and just reformated for a better table visualization.

# Libraries ----
library(gt)
library(tidyverse)
library(ggthemes)
library(ggrepel)

# Load colors ----
load("../3rd_LesionRNAseq_Dataset/RStudio_outputs/data/Dark24")


# Import table ----
MALDI_Count_prelesion_mod <- read_delim("MALDI_Count_prelesion_mod.txt", 
                                        "\t", escape_double = FALSE,
                                        col_types = cols(`Number of lesions` = col_number()), 
                                        trim_ws = TRUE)

MALDI_Count_prelesion_mod %>%
  gt()

MALDI_Count_prelesion_mod %>%
  ggplot(., aes(x=`Bacterial species`, y= `Number of lesions`, fill=`Bacterial species`)) +
  geom_bar(stat = "identity", fill="#AF0038") +
  theme_classic() +
  theme(legend.position="none",legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 75, hjust = 1),
        axis.text.y = element_text(size = 15)) +
  xlab("") + ylab("Number of lesions")
  
MALDI_Count_prelesion_mod %>%
  mutate(Percentage=paste0(round(`Number of lesions`/sum(`Number of lesions`)*100,0),"%")) %>%
  ggplot(., aes(x="", y=Percentage, fill=factor(`Bacterial species`))) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values = Dark24) +
  coord_polar("y", start=0) +
  geom_text(aes(label = Percentage), color = "white", size=6,
            position = position_stack(vjust = 0.5),
            x=1) +
  theme_classic()

MALDI_Count_prelesion_mod %>%
  mutate(csum = rev(cumsum(rev(`Number of lesions`))), 
         pos = `Number of lesions`/2 + lead(csum, 1),
         pos = if_else(is.na(pos), `Number of lesions`/2, pos)) %>%
  ggplot(., aes(x="", y=`Number of lesions`, fill=fct_inorder(`Bacterial species`))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_text_repel(aes(y = pos, label = paste0(`Number of lesions`),
                   size = 4.5, nudge_x = 5)) +
  scale_fill_manual(values = Dark24) +
  guides(size="none") +
  theme_void()

MALDI_Count_prelesion_mod_mod <- MALDI_Count_prelesion_mod %>%
  mutate(`Bacterial species` = as.factor(`Bacterial species`)) %>%
  mutate(`Bacterial species` = relevel("Staphylococcus aureus",
                                       grepl('Staphy|IFI', MALDI_Count_prelesion_mod_mod$`Bacterial species`))))
MALDI_Count_prelesion_mod_mod$`Bacterial species` <- relevel(MALDI_Count_prelesion_mod_mod$`Bacterial species`, "Staphylococcus aureus")



  ggplot(., aes(x="", y=`Number of lesions`, fill=factor(`Bacterial species`))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = Dark24) +
  guides(#fill="none",
         size="none") +
  theme_void()






