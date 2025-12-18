#PACKAGES #######
###CRAN
pkgs_cran <- c(
  "plyr", "UpSetR", "scales", "ggplot2", "GUniFrac", "stringr", "dplyr",
  "vegan", "ggdendro", "pairwiseAdonis", "randomcoloR", "readr",
  "ggpubr", "glue", "tidyverse", "cowplot", "rstatix", "Polychrome",
  "RColorBrewer", "paletteer", "Matrix", "colorspace", "lme4", "coin",
  "reshape2", "ggnewscale", "ggplotify", "writexl", "spaa", "ggtern",
  "microbiomeutilities", "janitor", "ggsignif", "ranacapa", "tiff",
  "ggh4x", "grid", "btools")

##BIOCONDUCTOR
pkgs_bioc <- c(
  "phyloseq", "ANCOMBC", "metagMisc", "metagenomeSeq",
  "MicrobiotaProcess", "microbiome"
)

# Install CRAN packages
installed <- pkgs_cran %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(pkgs_cran[!installed])
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
installed_bioc <- pkgs_bioc %in% rownames(installed.packages())
if (any(!installed_bioc)) {
  BiocManager::install(pkgs_bioc[!installed_bioc], update = FALSE, ask = FALSE)
}
# Load all libraries
invisible(lapply(c(pkgs_cran, pkgs_bioc), library, character.only = TRUE))

##loadlibraries####
library(phyloseq); library (ANCOMBC); library (plyr); library(scales); library(btools)
library(ggplot2);library(GUniFrac);library(stringr);library(dplyr);library(metagMisc);library(metagenomeSeq)
library(vegan);library(ggdendro);library(pairwiseAdonis);library(randomcoloR);library(readr);
library (ggpubr); library (glue); library(tidyverse); library (cowplot)
library(rstatix); library(Polychrome); library(RColorBrewer); library(paletteer)
library(Matrix); library(colorspace);
library(lme4); library (UpSetR); library(MicrobiotaProcess); library (microbiome) 
library(coin); library(reshape2); library(ggnewscale); library (ggplotify); library(writexl)
library (spaa); library (metagMisc); library (ggtern); library(microbiomeutilities); library (janitor)
library(ggsignif); library(ranacapa); library(ggplotify); library(tiff); library ("ggh4x"); library(grid)


source("g_unifrac.R")
source("uw_unifrac.R")
source("w_unifrac.R")
source("change16STaxaNames_KingdomSpecies.R")
source('merge_low_abundance_ra.R')
source('top_taxa_legend.R')

#import the data from qiime #####
qiimedata <- import_biom("table-with-taxonomy.biom", "tree.nwk", "dna-sequences.fasta")
qiimedata ##79334 taxa and 72 samples

#Sample names####
write.csv(sample_names(qiimedata), "sample_names.csv")

##Downloaded the file and made a couple changes to unite and make the metadata file 
map_file <- import_qiime_sample_data("metadata_updated_samples_nasalswabs.txt")

# Combining data with metadata #######
data <- merge_phyloseq(qiimedata, map_file)
data #79334 taxa and 72 samples

head(otu_table(data)) ##ASV count table
head(sample_data(data)) ##sample information
head(phyloseq::tax_table(data)) ##taxa information

### pruning data###
data <- prune_taxa(taxa_sums(data) > 0, data)
data #Still 79334 taxa and 72 samples

#PRE-PROCESSING ####
##Renaming ranks########
rank_names(data) # "Rank1" - "Rank7". Will change them
colnames(phyloseq::tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data) #properly named now

# changing the SILVA style naming (k__Bacteria, etc.)
tax.data <- data.frame(phyloseq::tax_table(data)) # extract the taxonomy table as a data frame
tax.data.names <- change16Staxa(tax.data) # this gets rid of the GG format

# now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string

# now filling in the empty slots with the highest assigned taxonomy 
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- kingdom 
  } else if (tax.data.names[i,3] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- phylum
  } else if (tax.data.names[i,4] == ""){
    class <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- class
  } else if (tax.data.names[i,5] == ""){
    order <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- order
  } else if (tax.data.names[i,6] == ""){
    family <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- family
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Species[i] <- paste("unclassified ",tax.data.names$Genus[i], sep = "")
  }
}
head(tax.data.names) # no more NAs and no more k__, they're now "unclassified 
tail(tax.data.names)# no more NAs and no more k__

#### Uncultured ones
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == "uncultured"){
    kingdom <- paste("uncultured", tax.data.names[i,1])
    tax.data.names[i, 2:7] <- kingdom
  } else if (tax.data.names[i,3] == "uncultured"){
    phylum<- paste("uncultured", tax.data.names[i,2])
    tax.data.names[i, 3:7] <- phylum
  } else if (tax.data.names[i,4] == "uncultured"){
    class <- paste("uncultured", tax.data.names[i,3])
    tax.data.names[i, 4:7] <- class
  } else if (tax.data.names[i,5] == "uncultured"){
    order <- paste("uncultured", tax.data.names[i,4])
    tax.data.names[i, 5:7] <- order
  } else if (tax.data.names[i,6] == "uncultured"){
    family <- paste("uncultured", tax.data.names[i,5])
    tax.data.names[i, 6:7] <- family
  } else if (tax.data.names[i,7] == "uncultured"){
    tax.data.names$Species[i] <- paste("uncultured",tax.data.names$Genus[i])
  }
}
tail(tax.data.names, 100) #no more blanks, no more NA's, no more only "uncultured"
head(tax.data.names)

### Re-inserting the taxonomy table into the phyloseq object##
phyloseq::tax_table(data) <- as.matrix(tax.data.names)

##Removing non Bacteria/Archaea and keeping Unassigned, to know classification percentage
data #79334 taxa, 72 samples
data2 <- subset_taxa(data, Kingdom=="Bacteria" | Kingdom=="Archaea" | Kingdom=="Unassigned")
data2 #79330 (dropped 4 taxa), 72 samples

##QC checks of the ASVs per samples########
min(sample_sums(data2)) # 1, one sample only has 1 ASV
max(sample_sums(data2)) # 5,823,716
mean(sample_sums(data2)) #558,227.2 
median(sample_sums(data2)) # 105,777.5
sort(sample_sums(data2)) # 1,6,14,73,208,227, 448 (controls and extraction blanks) 
#Most of the low ones are the individual samples

##Only Extraction blanks and NTC 
data2.blanks <- subset_samples(data2, sample_type=="eb") ##Getting only EB amples
data2.blanks <- prune_taxa(taxa_sums(data2.blanks) > 0, data2.blanks)  ##Pruning taxa with 0 total counts 
data2.blanks #2 extraction blanks, 4 taxa total
summary(sample_sums(data2.blanks))
# Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
# 1.00    2.25    3.50    3.50    4.75    6.00 
c(
  mean = mean(sample_sums(data2.blanks)),
  sd   = sd(sample_sums(data2.blanks)),
  sem  = sd(sample_sums(data2.blanks)) / sqrt(length(sample_sums(data2.blanks)))
)

data2.NTC <- subset_samples(data2, sample_type=="ntc") ##Getting only NTC samples
data2.NTC <- prune_taxa(taxa_sums(data2.NTC) > 0, data2.NTC )  ##Pruning taxa with 0 total counts
data2.NTC ##5 NTC, 96 taxa total
summary(sample_sums(data2.NTC))
# Min. 1st Qu.  Median  Mean  3rd Qu.    Max. 
# 14      73     208     194     227     448 
c(
  mean = mean(sample_sums(data2.NTC)),
  sd   = sd(sample_sums(data2.NTC)),
  sem  = sd(sample_sums(data2.NTC)) / sqrt(length(sample_sums(data2.NTC)))
)

##Dropping blanks and NTC
data4 <- prune_samples(sample_sums(data2) > 40000, data2) ##dropped those with less than 40000 reads
data4 # left with 65 samples, dropped NTCs, extraction blanks, kept all samples
data4 <- prune_taxa(taxa_sums(data4) > 0, data4) 
data4 #79290 taxa, 65 samples remain

##Recheck QC after dropping those with low reads (blanks)
min(sample_sums(data4)) #43,043
max(sample_sums(data4)) # 5,823,716 BIG DIFFERENCE 
mean(sample_sums(data4)) #618,329
median(sample_sums(data4)) # 117,202
sort(sample_sums(data4))

#COLOR PALETTES ####
pooling_palette <- c("forestgreen","cornflowerblue","#d87cec","#af49c5","#843b94","#5b146a","#33003e")

#COMPARING SEQUENCING DEPTHS#######
data4.sample.sums <- sample_sums(data4) #making a sample sums object
data4.df <- cbind(data4@sam_data, data4.sample.sums) #combining sample sums with metadata
data4.df

##ordering the pool_type variable
pool_order <- factor(data4.df$pool_type, levels = c("individual", "dna", "raw")) 
data4.df$pool_order <- factor (data4.df$pool_type, levels = c("individual", "dna", "raw"))

#Sequencing Depth (Pools DNA and Raw vs Individual)
Sequencing_depth <- ggplot(data4.df, aes(x= pool_order, y=data4.sample.sums, fill= pool_order, color= pool_order)) + 
  geom_boxplot(alpha=0.5) +
  geom_point(size=2) +
  theme_bw() +
  labs(fill= "Sample Type", color = "Sample Type",shape="Sample Type", 
       y=expression(atop("Sequencing Depth", "(Total Sum of ASVs per Sample)")),
       x = "Sample Type") +
  scale_x_discrete(labels = c("I", "D", "R"))+ 
  scale_y_continuous(expand = c(0.001, 0), limits = c(0, 6700000)) +
  scale_color_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_fill_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18, face= "bold"),
        legend.key.size = unit(1, "cm"), 
        axis.title = element_text(size=20),
        axis.text = element_text(color = "black", size= 18))+ 
  geom_pwc (method = "wilcox_test", 
            p.adjust.method = "BH", 
            label = "p = {p.adj}",
            size =0.5, 
            label.size = 5,
            tip.length = 0.02,
            step.increase = 0.05,
            hide.ns = TRUE)
Sequencing_depth
pairwise.wilcox.test(data4.df$data4.sample.sums, data4.df$pool_order, p.adjust.method = "BH")
## DNA pool vs Raw pool (p= 0.91), Individual vs DNA Pool (p=1e-10), Individual vs Raw pool (p=1e-10)
##SUPPLEMENTARY FIGURE 1 #######
sfigure1 <- Sequencing_depth
ggsave("SupplementaryFigure1.tiff", plot = sfigure1, device = "tiff", dpi = 500, width = 10, height = 8) 


# ALPHA DIVERSITY ####
##ASV LEVEL ######
## Richness, div 
alpha_div1 <- estimate_richness(data4, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) 
##Evenness
alpha_div2 <- microbiome::evenness(data4, index = "pielou", 
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0)
# Faith's pd (phylogenetic)
alpha_div3<- estimate_pd(data4)

# combine metrics with metadata
alpha_div_combined<- cbind(alpha_div1, 
                     alpha_div2,
                     alpha_div3)
alpha_div_combined
alpha_div <- alpha_div_combined%>%
  select(Observed, Shannon, pielou, PD) #Selecting only Richness, Shannons, Pielous evenness and Faiths PD 
alpha_div
alpha_div_meta <- cbind(data4@sam_data, alpha_div) 
alpha_div_meta # metadata and div metrics
#Identifier of alpha div done at the ASV level 
alpha_div_asv_meta <- alpha_div_meta %>%
  mutate(level="ASV") ##identifier of the level at which I'm doing the alpha diversity analysis

##creating factored vectors for the categorical variables
alpha_div_meta$pool_order <- factor (alpha_div_meta$pool_type, levels = c("individual", "dna", "raw"))

### RICHNESS (ASV) PLOT#######
Richness <-ggplot(alpha_div_meta, aes(x= pool_order, y= Observed, fill = pool_order, colour = pool_order)) +
  theme_bw() +
  labs(y= "Observed ASVs", title = "RICHNESS", fill= "Sample Type", colour= "Sample Type") +
  geom_boxplot(alpha = 0.5) + 
  geom_point(size=2) +
  scale_color_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_fill_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_y_continuous(limits= c(0, 27000), expand = c(0.001, 0.001)) +
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18, colour = "black"))+ 
  geom_pwc (method = "wilcox_test", 
            p.adjust.method = "BH", 
            label = "Wilcoxon, p = {p.adj}",
            size =0.5, 
            label.size = 4.2,
            tip.length = 0.02,
            step.increase = 0.05,
            hide.ns = TRUE)
Richness
pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$pool_type, p.adjust.method = "BH") 
## DNA pool vs Raw pool (p= 0.6), Individual vs DNA Pool (p=1 e-10), Individual vs Raw pool (p=1 e-1o)
count(alpha_div_meta, pool_type) # 10 DNA pools, 10 Raw pools, 45 Individuals


### SHANNONS PLOT#######
Shannon <-ggplot(alpha_div_meta, aes(x= pool_order, y= Shannon, fill = pool_order, colour = pool_order)) +
  theme_bw() +
  labs(y= "Shannon", title = "DIVERSITY", fill= "Sample Type", colour= "Sample Type") +
  geom_boxplot(alpha = 0.5) + 
  geom_point(size=2) +
  scale_color_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_fill_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_y_continuous(limits = c(0,8), expand = c(0.001,0.001)) +
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18, colour = "black"))+ 
  geom_pwc (method = "wilcox_test", 
            p.adjust.method = "BH", 
            label = "Wilcoxon, p = {p.adj}",
            size =0.5, 
            label.size = 4.2,
            tip.length = 0.02,
            y.position = c(7.5,7.8),
            hide.ns = TRUE)
Shannon
pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$pool_type, 
                     p.adjust.method = "BH") ## pool vs individual 0.079 
## DNA pool vs Raw pool (p= 0.82), Individual vs DNA Pool (p=0.77), Individual vs Raw pool (p=0.82)

### FPD PLOT #######
FPD <-ggplot(alpha_div_meta, aes(x= pool_order, y= PD, fill = pool_order, colour = pool_order)) +
  theme_bw() +
  labs(y= "Faith's PD", title = "PHYLOGENETIC DIVERSITY", fill= "Sample Type", colour= "Sample Type") +
  geom_boxplot(alpha = 0.5) + 
  geom_point(size=2) +
  scale_color_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_fill_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_y_continuous(limits= c(0, 2400), expand = c(0.0009,0.01)) +
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18, colour = "black"))+ 
  geom_pwc (method = "wilcox_test", 
            p.adjust.method = "BH", 
            label = "Wilcoxon, p = {p.adj}",
            size =0.5, 
            label.size = 4.2,
            step.increase = 0.05,
            tip.length = 0.02,
            hide.ns = TRUE)
FPD
pairwise.wilcox.test(alpha_div_meta$PD, alpha_div_meta$pool_type, p.adjust.method = "BH") 
## DNA pool vs Raw pool (p= 0.68), Individual vs DNA Pool (p= 1 e-10), Individual vs Raw pool (p=1 e-10)

##PHYLUM LEVEL ########
#Creating count table at the phylum level
data4_phylum.counts <- tax_glom(data4, taxrank = "Phylum", NArm=FALSE)
alpha_div1_phylum<- estimate_richness(data4_phylum.counts, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # richness, div
alpha_div2_phylum <- microbiome::evenness(data4_phylum.counts, index = "pielou", zeroes = FALSE, 
                                          detection = 0) ##evenness
alpha_div3_phylum <- estimate_pd(data4_phylum.counts) # faith's pd (phylogenetic)

# combine metrics with metadata
alpha_div_phylum <- cbind(alpha_div1_phylum, alpha_div2_phylum, alpha_div3_phylum)
alpha_div_phylum 
alpha_div_phylum <- alpha_div_phylum%>%
  select(Observed, Shannon, pielou, PD) #Selecting only Richness, Shannons, Pielous evenness and Faiths PD 
alpha_div_phylum ## observed, shannon, simpson, invsimpson, pielou, PD
alpha_div_phylum_meta <- cbind(data4@sam_data, alpha_div_phylum ) 
alpha_div_phylum_meta ##now with metadata
alpha_div_phylum_meta <- alpha_div_phylum_meta %>%
  mutate(level="Phylum") ##identifier of the level at which I'm doing the alpha diversity analysis


##CLASS LEVEL ########
#Creating count table at the phylum level
data4_class.counts <- tax_glom(data4, taxrank = "Class", NArm=FALSE)
# richness, div
alpha_div1_class<- estimate_richness(data4_class.counts, 
                                     measures = c("Observed", "Shannon", "Simpson","InvSimpson")) 
#Evenness
alpha_div2_class <- microbiome::evenness(data4_class.counts, index = "pielou", 
                                         zeroes = FALSE, detection = 0) 
# Faith's PD (phylogenetic)
alpha_div3_class <- estimate_pd(data4_class.counts) 

# combine metrics with metadata
alpha_div_class <- cbind(alpha_div1_class, alpha_div2_class, alpha_div3_class)
alpha_div_class 
alpha_div_class <- alpha_div_class%>%
  select(Observed, Shannon, pielou, PD)  #Selecting only Richness, Shannons, Pielous evenness and Faiths PD 
alpha_div_class ## observed, shannon, simpson, invsimpson, pielou, PD
alpha_div_class_meta <- cbind(data4@sam_data, alpha_div_class ) 
alpha_div_class_meta ##now with metadata
alpha_div_class_meta <- alpha_div_class_meta %>%
  mutate(level="Class") ##identifier of the level at which I'm doing the alpha diversity analysis
alpha_div_class_meta ##class

##ORDER LEVEL ########
#Creating count table at the order level
data4_order.counts <- tax_glom(data4, taxrank = "Order", NArm=FALSE)
# Richness, Diversity
alpha_div1_order<- estimate_richness(data4_order.counts, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) 
#Evenness
alpha_div2_order <- microbiome::evenness(data4_order.counts, index = "pielou", 
                                         zeroes = FALSE, detection = 0)
# faith's pd (phylogenetic
alpha_div3_order <- estimate_pd(data4_order.counts)

# combine metrics with metadata
alpha_div_order <- cbind(alpha_div1_order, alpha_div2_order, alpha_div3_order)
alpha_div_order 
alpha_div_order <- alpha_div_order%>%
select(Observed, Shannon, pielou, PD)  #Selecting only Richness, Shannons, Pielous evenness and Faiths PD
alpha_div_order ## observed, shannon, simpson, invsimpson, pielou, PD
alpha_div_order_meta <- cbind(data4@sam_data, alpha_div_order ) 
alpha_div_order_meta ##now with metadata
alpha_div_order_meta <- alpha_div_order_meta %>%
  mutate(level="Order") ##identifier of the level at which I'm doing the alpha diversity analysis
alpha_div_order_meta ##order

##FAMILY LEVEL ########
#Creating count table at the family level
data4_family.counts <- tax_glom(data4, taxrank = "Family", NArm=FALSE)
# richness, div
alpha_div1_family<- estimate_richness(data4_family.counts, measures = c("Observed", "Shannon", "Simpson","InvSimpson"))
#evenness
alpha_div2_family <- microbiome::evenness(data4_family.counts, index = "pielou", zeroes = FALSE, detection = 0)
# faith's pd (phylogenetic)
alpha_div3_family <- estimate_pd(data4_family.counts) 

# combine metrics with metadata
alpha_div_family <- cbind(alpha_div1_family, alpha_div2_family, alpha_div3_family)
alpha_div_family 
alpha_div_family <- alpha_div_family%>%
  select(Observed, Shannon, pielou, PD)  #Selecting only Richness, Shannons, Pielous evenness and Faiths PD
alpha_div_family ## observed, shannon, simpson, invsimpson, pielou, PD
alpha_div_family_meta <- cbind(data4@sam_data, alpha_div_family ) 
alpha_div_family_meta ##now with metadata
alpha_div_family_meta <- alpha_div_family_meta %>%
  mutate(level="Family") ##identifier of the level at which I'm doing the alpha diversity analysis
alpha_div_family_meta ##family

##GENUS LEVEL ########
#Creating count table at the genus level
data4_genus.counts <- tax_glom(data4, taxrank = "Genus", NArm=FALSE)
# Richness, diversity
alpha_div1_genus<- estimate_richness(data4_genus.counts, measures = c("Observed", "Shannon", "Simpson","InvSimpson"))
alpha_div2_genus <- microbiome::evenness(data4_genus.counts, index = "pielou", zeroes = FALSE, detection = 0) #evenness
alpha_div3_genus <- estimate_pd(data4_genus.counts) # faith's pd (phylogenetic)

# combine metrics with metadata
alpha_div_genus <- cbind(alpha_div1_genus, alpha_div2_genus, alpha_div3_genus)
alpha_div_genus 
alpha_div_genus <- alpha_div_genus%>%
  select(Observed, Shannon, pielou, PD)  #Selecting only Richness, Shannons, Pielous evenness and Faiths PD
alpha_div_genus ## observed, shannon, simpson, invsimpson, pielou, PD
alpha_div_genus_meta <- cbind(data4@sam_data, alpha_div_genus ) 
alpha_div_genus_meta ##now with metadata
alpha_div_genus_meta <- alpha_div_genus_meta %>%
  mutate(level="Genus") ##identifier of the level at which I'm doing the alpha diversity analysis
alpha_div_genus_meta ##genus

##Combining Alpha div at ASV and other taxa levels#########
alpha_diversity_taxa_meta <- rbind(alpha_div_asv_meta%>%
                                     select(!pool_order), 
                                   alpha_div_phylum_meta, 
                                   alpha_div_class_meta, 
                                   alpha_div_order_meta, 
                                   alpha_div_family_meta,
                                   alpha_div_genus_meta)

##long format, making plotting easier
alpha_div_taxa<- alpha_diversity_taxa_meta %>%
  pivot_longer(
    cols = c(Observed, Shannon, pielou, PD),
    names_to = "alpha_div_index",
    values_to = "value") 


##Factoring (order I want)
alpha_div_taxa$pool_order <- factor (alpha_div_taxa$pool_type, levels = c("individual", "dna", "raw"))
alpha_div_taxa$alphadivorder <- factor (alpha_div_taxa$alpha_div_index, levels = c("Observed", "PD", "pielou","Shannon"))
alpha_div_taxa$levelorder <- factor (alpha_div_taxa$level, levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))

##Plotting alpha div for all taxonomic levels
alpha_diversity_taxa_plot <- alpha_div_taxa %>%
  ggplot(aes(x=pool_order, y=value, color= levelorder))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Paired")+
  scale_x_discrete(labels = c("I", "D", "R"))+ 
  facet_wrap (~alphadivorder, scales = "free_y", 
              labeller = labeller(alphadivorder= c(PD = "Faith's PD", pielou = "Evenness", Observed = "Richness")),
              nrow = 1)+
  theme_bw()+
  labs(x= "Sample Type", color = "Tax Level")+
  theme(
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10, face = "bold"), 
    axis.text = element_text(size=12),
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=12, face= "bold"),
    strip.text = element_text(size = 12)) +
  guides(color = guide_legend(position = "bottom", direction = "horizontal", nrow = 1)) 
alpha_diversity_taxa_plot

###Plotting alpha div at onlyASV level 
alpha_div_taxa_asv <- alpha_div_taxa %>%
  filter(level == "ASV")

alpha_diversity_taxa_plot_asv <- alpha_div_taxa_asv %>%
  ggplot(aes(x=pool_order, y=value, color= pool_order, fill= pool_order))+
  geom_boxplot(alpha=0.5)+
  geom_point(size=2)+
  scale_color_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_fill_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_x_discrete(labels = c("I", "D", "R"))+ 
  facet_wrap (~alphadivorder, scales = "free_y", 
              labeller = labeller(alphadivorder= c(PD = "Faith's PD", pielou = "Evenness", Observed = "Richness")),
              nrow = 1)+
  theme_bw()+
  labs(x= "Sample Type", color = "Sample Type", fill= "Sample Type")+
  theme(legend.text = element_text(size=22),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size=25),
        axis.title.x =  element_text(size = 28),
        axis.title.y = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.background = element_blank(),
        strip.text = element_text(size=38),
        legend.position = "bottom")+
  geom_pwc (method = "wilcox_test", 
            p.adjust.method = "BH", 
            label = "p = {p.adj}",
            size =0.5, 
            label.size = 4.5,
            tip.length = 0.02,
            hide.ns = TRUE)
alpha_diversity_taxa_plot_asv
###FIGURE 1 #####
figure1 <- alpha_diversity_taxa_plot_asv
ggsave("Figure1pubr.tiff", plot = figure1, device = "tiff", width = 14, height = 7, dpi = 300)


###Plotting alpha div at all other taxa levels except for ASV #######
alpha_div_taxa_nonASV <- alpha_div_taxa %>%
  filter(!level == "ASV")
alpha_diversity_taxa_plot_nonASV<- alpha_div_taxa_nonASV %>%
  ggplot(aes(x=pool_order, y=value, color= pool_order, fill= pool_order))+
  geom_boxplot(alpha= 0.5)+
  geom_point()+
  scale_color_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_fill_manual(values= pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_x_discrete(labels = c("I", "D", "R"))+
  facet_nested(alphadivorder~levelorder, scales = "free", independent = "y",
               labeller = labeller(alphadivorder = c(PD = "Faith's PD", pielou = "Evenness", Observed = "Richness")))+
  theme_bw()+
  labs(x= "Sample Type", color = "Sample type", fill = "Sample type")+
  theme(legend.text = element_text(size=25),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size=27),
        axis.title.x =  element_text(size = 29),
        axis.title.y = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.background = element_blank(),
        strip.text = element_text(size=32),
        legend.position = "bottom")+
  geom_pwc (method = "wilcox_test", 
            p.adjust.method = "BH", 
            label = "p = {p.adj}",
            size =0.5, 
            label.size = 4.5,
            step.increase = 0.08,
            tip.length = 0.02,
            hide.ns = TRUE)
alpha_diversity_taxa_plot_nonASV
###SUPPLEMENTARY FIGURE 2 #######
sfigure2 <- alpha_diversity_taxa_plot_nonASV
ggsave("SupplementaryFigure2.tiff", plot = sfigure2, device = "tiff", dpi = 300, width = 20, height = 18) 



#STATS ON ALPHA DIVERSITY #####
##Checking results, to match those by geom_pwc
##PHYLUM LEVEL#####
alpha_div_taxa_Phylum <- alpha_div_taxa[which(alpha_div_taxa$level=="Phylum"),] 
#Richness
alpha_div_taxa_Phylum %>%
  filter(alpha_div_index=="Observed")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, significant between individuals and pools (both), but n.s. between DNA pools and raw pools 
#Shannon
alpha_div_taxa_Phylum %>%
  filter(alpha_div_index=="Shannon")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, n.s. between all of the groups 
#Evenness
alpha_div_taxa_Phylum %>%
  filter(alpha_div_index=="pielou")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, n.s. between all of the groups 

#Faith's PD 
alpha_div_taxa_Phylum %>%
  filter(alpha_div_index=="PD")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, significant between individuals and pools (DNA and raw), but n.s. between DNA pools and raw pools 

##CLASS LEVEL########
alpha_div_taxa_Class <- alpha_div_taxa[which(alpha_div_taxa$level=="Class"),]
#Richness
alpha_div_taxa_Class %>%
  filter(alpha_div_index=="Observed")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, significant between individuals and pools (both), but n.s. between DNA pools and raw pools 
#Shannon
alpha_div_taxa_Class %>%
  filter(alpha_div_index=="Shannon")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, n.s. between any of the groups 
#Evenness
alpha_div_taxa_Class %>%
  filter(alpha_div_index=="pielou")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, n.s. between any of the groups 

#Faith's PD 
alpha_div_taxa_Class %>%
  filter(alpha_div_index=="PD")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, significant between individuals and pools (both), but n.s. between DNA pools and raw pools 


##ORDER LEVEL#######
alpha_div_taxa_Order <- alpha_div_taxa[which(alpha_div_taxa$level=="Order"),]
#Richness
alpha_div_taxa_Order %>%
  filter(alpha_div_index=="Observed")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, significant between individuals and pools (both), but n.s. between DNA pools and raw pools 
#Shannon
alpha_div_taxa_Order %>%
  filter(alpha_div_index=="Shannon")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, n.s. between any of the groups 
#Evenness
alpha_div_taxa_Order %>%
  filter(alpha_div_index=="pielou")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, n.s. between any of the groups 

#Faith's PD 
alpha_div_taxa_Order %>%
  filter(alpha_div_index=="PD")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, significant between individuals and pools (DNA and raw), but n.s. between DNA pools and raw pools 

##FAMILY LEVEL#########
alpha_div_taxa_Family <- alpha_div_taxa[which(alpha_div_taxa$level=="Family"),]
#Richness
alpha_div_taxa_Family %>%
  filter(alpha_div_index=="Observed")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, significant between individuals and pools (DNA and raw), but n.s. between DNA pools and raw pools 
#Shannon
alpha_div_taxa_Family %>%
  filter(alpha_div_index=="Shannon")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, n.s. between any of the groups 
#Evenness
alpha_div_taxa_Family %>%
  filter(alpha_div_index=="pielou")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, n.s. between any of the groups 

#Faith's PD 
alpha_div_taxa_Family %>%
  filter(alpha_div_index=="PD")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, significant between individuals and pools (DNA and Raw), but n.s. between DNA pools and raw pools 
##GENUS LEVEL#######
alpha_div_taxa_Genus <- alpha_div_taxa[which(alpha_div_taxa$level=="Genus"),]
#Richness
alpha_div_taxa_Genus %>%
  filter(alpha_div_index=="Observed")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, significant between individuals and pools (both), but n.s. between DNA pools and raw pools 
#Shannon
alpha_div_taxa_Genus %>%
  filter(alpha_div_index=="Shannon")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, n.s. between any of the groups 
#Evenness
alpha_div_taxa_Genus %>%
  filter(alpha_div_index=="pielou")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, n.s. between any of the groups 

#Faith's PD 
alpha_div_taxa_Genus %>%
  filter(alpha_div_index=="PD")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, significant between individuals and pools (DNA and raw), but n.s. between DNA pools and raw pools 

##GENUS LEVEL########
alpha_div_taxa_asv
#Richness
alpha_div_taxa_asv %>%
  filter(alpha_div_index=="Observed")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, significant between individuals and pools (DNA and raw), but n.s. between DNA pools and raw pools 
#Shannon
alpha_div_taxa_asv %>%
  filter(alpha_div_index=="Shannon")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH") #OK, n.s. between any of the groups 
#Evenness
alpha_div_taxa_asv%>%
  filter(alpha_div_index=="pielou")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, n.s. between any of the groups 

#Faith's PD 
alpha_div_taxa_asv%>%
  filter(alpha_div_index=="PD")%>%
  rstatix::wilcox_test(value~pool_type,
                       comparisons = list(c("dna", "raw"), 
                                          c("dna", "individual"), 
                                          c("raw", "individual")), 
                       p.adjust.method = "BH")#OK, significant between individuals and pools (both), but n.s. between DNA pools and raw pools 


#BETADIVERSITY #############
## CSS TRANSFORMATION (NORAMLIZATION) AND RELATIVE ABUNDANCE ####
#normalizing to account for differences in sequencing depth 
any(sample_sums(data4)== 0) ## no samples with 0 ASVs 
data4.css <- metagMisc::phyloseq_transform_css(data4, log = F)

#Relative abundance, transforming so the counts in every sample now add up to 100%
data4.ra <- transform_sample_counts(data4.css, function(x) {x/sum(x)}*100) 

#TAX GLOMMING #############
#seeing what proportion of our ASVs are classified at each rank 
##PHYLUM####
data4_phylum.ra <- tax_glom(data4.ra, taxrank = "Phylum", NArm=FALSE) ##NA false keeping what is unclassified
write.csv(data4_phylum.ra@tax_table, "phylum_taxa.csv") #taking the tax table
write.csv(data4_phylum.ra@otu_table, "phylum_otus.csv") #taking the counts table
###Percentage Unclassified 
data4_phylum.ra%>%
  psmelt(.)%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE))%>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU across samples
  summarize(Unclassified_sum = sum(OTU_Abundance)) #1.08% Unclassified. So 100 - 1.08 = 98.92% Classified 

##CLASS####
data4_class.ra <- tax_glom(data4.ra, taxrank = "Class", NArm=FALSE) 
write.csv(data4_class.ra@tax_table, "class_taxa.csv") 
write.csv(data4_class.ra@otu_table, "class_otus.csv") 
###Percentage Unclassified 
data4_class.ra%>%
  psmelt(.)%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE))%>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU across samples
  summarize(Unclassified_sum = sum(OTU_Abundance)) #1.08% Unclassified. So 100 - 1.08 = 98.92% Classified 

##ORDER ####
data4_order.ra <- tax_glom(data4.ra, taxrank = "Order", NArm=FALSE) 
write.csv(data4_order.ra@tax_table, "order_taxa.csv") 
write.csv(data4_order.ra@otu_table, "order_otus.csv")

###Percentage Unclassified 
data4_order.ra%>%
  psmelt(.)%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE))%>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU across samples
  summarize(Unclassified_sum = sum(OTU_Abundance)) #1.13% Unclassified. So 100 - 1.13 = 98.87% Classified 

##FAMILY ####
data4_family.ra <- tax_glom(data4.ra, taxrank = "Family", NArm=FALSE) 
write.csv(data4_family.ra@tax_table, "family_taxa.csv") 
write.csv(data4_family.ra@otu_table, "family_otus.csv")

###Percentage Unclassified 
data4_family.ra%>%
  psmelt(.)%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE))%>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU across samples
  summarize(Unclassified_sum = sum(OTU_Abundance)) #1.40% Unclassified. So 100 - 1.40 = 98.6% Classified 

##GENUS#####
data4_genus.ra <- tax_glom(data4.ra, taxrank = "Genus", NArm=FALSE) 
write.csv(data4_genus.ra@tax_table, "genus_taxa.csv") 
write.csv(data4_genus.ra@otu_table, "genus_otus.csv")
###Percentage Unclassified 
data4_genus.ra%>%
  psmelt(.)%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE))%>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU across samples
  summarize(Unclassified_sum = sum(OTU_Abundance)) #6.98% Unclassified. So 100 - 6.98 = 93.02% Classified 

## DISTANCE AND ORDINATIONS - ASV LEVEL##############
###GENERALIZED UNIFRAC ######
#distance matrix
data4.css.gunifrac <- gunifrac(data4.css)
data4.css.gunifrac
##ordinate with NMDS
data4.css.gunifrac.ord <- metaMDS(comm = data4.css.gunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4.css.gunifrac.scrs <- scores(data4.css.gunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4.css.gunifrac.scrs <- cbind(as.data.frame(data4.css.gunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4.css.gunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4.css.gunifrac.scrs, FUN = mean) ## Calculating the centroids (mean method)
data4.css.gunifrac.segs <- merge(data4.css.gunifrac.scrs, setNames(data4.css.gunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4.css.gunifrac.segs <- data4.css.gunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4.css.gunifrac.segs

##ordering categorical pool type
data4.css.gunifrac.segs$pool_order <- factor (data4.css.gunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

##PLOTTING- CENTROIDS
NMDS_generalized_unifrac <-ggplot(data4.css.gunifrac.segs) +
  theme_bw() +
  labs(x= "NMDS1", y = "NMDS2", title= "GENERALIZED UNIFRAC", colour= "Sample Type", fill = "Sample Type") +
  geom_point(aes (x= MDS1, y = MDS2, colour= pool_order), size= 4, shape =18) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= pool_order, colour = pool_order), alpha = 0.32, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes (x= cMDS1, y = cMDS2, colour= pool_order), size= 12, shape =18) +
  geom_text(aes (x= cMDS1, y = cMDS2,label= pool_type.abbrv), colour= "white", size = 5) + ##text adds the abbrv. version of the pool type
  scale_fill_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_color_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 28),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title = element_text(size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.75))
NMDS_generalized_unifrac

### WEIGHTED UNIFRAC #####
#distance matrix
data4.css.wunifrac <- wunifrac(data4.css)
data4.css.wunifrac

##ordinate with NMDS
data4.css.wunifrac.ord <- metaMDS(comm = data4.css.wunifrac,k =2,  try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4.css.wunifrac.scrs <- scores(data4.css.wunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4.css.wunifrac.scrs <- cbind(as.data.frame(data4.css.wunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4.css.wunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4.css.wunifrac.scrs, FUN = mean) ## Calculating the centroids 
data4.css.wunifrac.segs <- merge(data4.css.wunifrac.scrs, setNames(data4.css.wunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4.css.wunifrac.segs <- data4.css.wunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4.css.wunifrac.segs

##ordering categorical pool type
data4.css.wunifrac.segs$pool_order <- factor (data4.css.wunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

##PLOTTING- CENTROIDS
NMDS_weighted_unifrac <-ggplot(data4.css.wunifrac.segs) +
  theme_bw() +
  labs(x= "NMDS1", y = "NMDS2", title= "WEIGHTED UNIFRAC", colour= "Sample Type", fill = "Sample Type") +
  geom_point(aes (x= MDS1, y = MDS2, colour= pool_order), size= 4, shape =18) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= pool_order, colour = pool_order), alpha = 0.32, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes (x= cMDS1, y = cMDS2, colour= pool_order), size= 12, shape =18) +
  geom_text(aes (x= cMDS1, y = cMDS2,label= pool_type.abbrv), colour= "white", size = 5) + ##text adds the abbrv. version of the pool type
  scale_fill_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_color_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 28),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title = element_text(size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.75))
NMDS_weighted_unifrac

### UNWEIGHTED UNIFRAC #####
#distance matrix
data4.css.uwunifrac <- uwunifrac(data4.css)
data4.css.uwunifrac
##ordinate with NMDS
data4.css.uwunifrac.ord <- metaMDS(comm = data4.css.uwunifrac, k =2,  try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4.css.uwunifrac.scrs <- scores(data4.css.uwunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4.css.uwunifrac.scrs <- cbind(as.data.frame(data4.css.uwunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4.css.uwunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4.css.uwunifrac.scrs, FUN = mean) ## Calculating the centroids 
data4.css.uwunifrac.segs <- merge(data4.css.uwunifrac.scrs, setNames(data4.css.uwunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4.css.uwunifrac.segs <- data4.css.uwunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4.css.uwunifrac.segs

##ordering categorical pool type
data4.css.uwunifrac.segs$pool_order <- factor (data4.css.uwunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

##PLOTTING- CENTROIDS
NMDS_unweighted_unifrac <-ggplot(data4.css.uwunifrac.segs) +
  theme_bw() +
  labs(x= "NMDS1", y = "NMDS2", title= "UNWEIGHTED UNIFRAC", colour= "Sample Type", fill = "Sample Type") +
  geom_point(aes (x= MDS1, y = MDS2, colour= pool_order), size= 4, shape =18) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= pool_order, colour = pool_order), alpha = 0.32, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes (x= cMDS1, y = cMDS2, colour= pool_order), size= 12, shape =18) +
  geom_text(aes (x= cMDS1, y = cMDS2,label= pool_type.abbrv), colour= "white", size = 5) + ##text adds the abbrv. version of the pool type
  scale_fill_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_color_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 28),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title = element_text(size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.75))
NMDS_unweighted_unifrac

###UNIFRAC-BASED ORDINATIONS TOGETHER INTO ONE PLOT####
data4asv.css.gunifrac.segs.plot <- data4.css.gunifrac.segs %>%
  mutate(Taxlevel = "ASV", UniFrac = "Generalized UniFrac")
data4asv.css.wunifrac.segs.plot <- data4.css.wunifrac.segs %>%
  mutate(Taxlevel = "ASV", UniFrac = "Weighted UniFrac")
data4asv.css.uwunifrac.segs.plot <- data4.css.uwunifrac.segs %>%
  mutate(Taxlevel = "ASV", UniFrac = "Unweighted UniFrac")

##One single dataframe with all UniFrac distances at ASV level
data4asv.css.unifrac.segs.plot <- rbind(data4asv.css.gunifrac.segs.plot,
                                        data4asv.css.wunifrac.segs.plot,
                                        data4asv.css.uwunifrac.segs.plot)
data4asv.css.unifrac.segs.plot$UniFrac <- factor (data4asv.css.unifrac.segs.plot$UniFrac, 
                                                  levels = c("Weighted UniFrac", "Generalized UniFrac","Unweighted UniFrac"))

##PLOTTING
NMDS_unifrac_ASV <-ggplot(data4asv.css.unifrac.segs.plot) +
  theme_bw() +
  facet_wrap(~UniFrac, scales = "free", ncol= 1)+
  geom_point(aes (x= MDS1, y = MDS2, colour= pool_order), size= 4, shape =18) +
  labs(x= "NMDS1", y = "NMDS2", colour= "Sample Type", fill = "Sample Type") +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= pool_order, colour = pool_order), alpha = 0.32, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes (x= cMDS1, y = cMDS2, colour= pool_order), size= 8, shape =18) +
  geom_text(aes (x= cMDS1, y = cMDS2,label= pool_type.abbrv), colour= "white", size = 3) + ##text adds the abbrv. version of the pool type
  scale_fill_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_color_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  theme(panel.spacing = unit(0, "lines"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size = 12),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.background = element_blank(),
        strip.text = element_text(size=14,margin = margin(0.01, 0.01, 0.01, 0.01, "cm")),
        legend.key.size = unit(0.5, "cm"),   
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11, face= "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.1, "cm"),
        legend.position = "bottom",
        legend.justification = "left",
        legend.location = "plot")+
  guides(
    colour = guide_legend(override.aes = list(size = 3, shape = 18)),
  )
NMDS_unifrac_ASV 

## STATS ON COMMUNITY COMPOSITION - ASV LEVEL##### 
#Need a data frame with sample metadata
data4.css.df <- as(data4.css@sam_data, "data.frame") 
###GENERALIZED UNIFRAC########
####PERMANOVA #####
set.seed(87)
gunifrac.adonis <- adonis2(data4.css.gunifrac ~ pool_type, 
                           data = data4.css.df, by = "margin",permutations = 9999)
gunifrac.adonis
## only 15.60% of the variation in composition is due to pool type. p value = 1e-04

####PAIRWISE PERMANOVA #######
set.seed(87)
gunifrac.pairwise.adonis <- pairwise.adonis2(data4.css.gunifrac ~ pool_type, 
                                             data = data4.css.df, 
                                             by = "margin", permutations = 9999, p.adjust.method = "BH")
gunifrac.pairwise.adonis 
##difference in community structure between individual vs raw pools (p= 1e-04) and individual vs dna pools(p= 1e-04)
## no significant difference between dna pool and raw pool (p=0.96)


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.gunifrac.disp <- betadisper(data4.css.gunifrac, data4.css.df$pool_type)
pooltype.gunifrac.disp


##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.gunifrac.permdisp <- permutest(pooltype.gunifrac.disp, permutations = 9999, pairwise = 1)
pooltype.gunifrac.permdisp
## dispersions of variances are the same (p=0.3)


###WEIGHTED UNIFRAC#######
####PERMANOVA #####
set.seed(87)
wunifrac.adonis <- adonis2(data4.css.wunifrac~ pool_type, data = data4.css.df, by = "margin", permutations = 9999)
wunifrac.adonis
## 6.8% of the variation in composition is due to pool type. p value = 0.019

####PAIRWISE PERMANOVA #####
set.seed(87)
wunifrac.pairwise.adonis <- pairwise.adonis2(data4.css.wunifrac ~ pool_type, data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
wunifrac.pairwise.adonis 
##difference in community structure between individual vs raw pools(p= 0.023), and between individual and dna pools (p=0.042)
## no differences between raw vs dna pools (p=0.923)


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.wunifrac.disp <- betadisper(data4.css.wunifrac, data4.css.df$pool_type)
pooltype.wunifrac.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.wunifrac.permdisp <- permutest(pooltype.wunifrac.disp, permutations = 9999, pairwise = 1)
pooltype.wunifrac.permdisp
## dispersion of variances are the same (p=0.22)

###UNWEIGHTED UNIFRAC##########
####PERMANOVA #####
set.seed(87)
uwunifrac.adonis <- adonis2(data4.css.uwunifrac~ pool_type, data = data4.css.df, by = "margin", permutations = 9999)
uwunifrac.adonis
## 23.87% of the variation in composition is due to pool type. p value = 1e-04

####PAIRWISE PERMANOVA #####
set.seed(87)
uwunifrac.pairwise.adonis <- pairwise.adonis2(data4.css.uwunifrac ~ pool_type, data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
uwunifrac.pairwise.adonis 
##difference in community structure between individual vs raw pools (p= 1e-04) and individual vs dna pools (p= 1e-04)
# no significant difference between dna pool and raw pool (p=0.90) 


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.uwunifrac.disp <- betadisper(data4.css.uwunifrac, data4.css.df$pool_type)
pooltype.uwunifrac.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.uwunifrac.permdisp <- permutest(pooltype.uwunifrac.disp, permutations = 9999, pairwise = 1)
pooltype.uwunifrac.permdisp
## Variances are NOT the same (p=0.021): DNA pools vs individual are significantly different (p = 0.015000)


#UPSET PLOT#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
upsetda <- MicrobiotaProcess::get_upset(data4, factorNames="pool_type") 

#relabeling 
names(upsetda)[names(upsetda) == "individual"] <- "Individuals"
names(upsetda)[names(upsetda) == "dna"] <- "DNA Pool"
names(upsetda)[names(upsetda) == "raw"] <- "Raw Pool"

##saving the upset plot
tiff("Upset_plot.tiff", units="in", width=7, height=6, res=300)
upset_plot <-upset(upsetda, 
                   sets = c("Individuals", "DNA Pool", "Raw Pool"),
                   sets.bar.color = c("#d87cec","cornflowerblue","forestgreen"), 
                   order.by = "freq", text.scale = c(3, 2.5, 3, 1, 3, 2), 
                   point.size = 5, line.size = 2, mainbar.y.label= "ASV count",
                   sets.x.label = "ASV count", 
                   set_size.show = F)
upset_plot
dev.off()


#PREVALENCE VS RELATIVE ABUNDANCE ASVs#####
##Calculating prevalence (number of samples with ASV present)
prevalance_asv<- apply(X = otu_table(data4),
                       MARGIN = 1,
                       FUN = function(x){sum(x > 0)}) 
##Turning into dataframe. This has prevalence as number of samples the ASV is present in
prevalencedf_asv <- as.data.frame(prevalance_asv) 
##Adding a column of ASV labels
prevalencedf_asv$ASV <-rownames(prevalencedf_asv) ##adding a column of ASV labels

prevalence_asv_percentage <- microbiome::prevalence(data4,
                                                    detection = 0) %>%
  data.frame %>%
  rownames_to_column(var="ASV")%>%
  rename("Prev_percentage" =".")
prevalence_asv_percentage ##this df has prevalence of each ASV across samples as a proportion

#RA from raw counts
data4.raw.ra <- transform_sample_counts(data4, function(x) {x/sum(x)}*100) 
data4.raw.ra@otu_table ##relative abundance values for each ASV
asv_ra_table<- otu_table(data4.raw.ra) # Extract the ASV table
asv_ra_table_df <- as.data.frame(asv_ra_table) # Convert ASV table to a dataframe
asv_ra_table_df$ASV <- rownames(asv_ra_table_df) # Add ASV names as a column

upsetda$ASV <- rownames(upsetda) ##Add ASV names as a column in upset data 

##Prevalence counts and RA percentages 
prev_ra_df <- merge(prevalencedf_asv, prevalence_asv_percentage, by = "ASV") %>%
  merge(asv_ra_table_df, by = "ASV") ##merging prevalence and RA ASV table

#Adding upset data (this tells me which sample type each ASV is present in)
prev_ra_upset_df <- merge(prev_ra_df, upsetda, by= "ASV") 

## Calculating the mean relative abundance of each ASV across samples 
mean_relative_abundance <- apply(data4.ra@otu_table, MARGIN = 1, FUN = mean) 
mean_asv_relative_abundance_df<- as.data.frame(mean_relative_abundance) ##turning into dataframe
mean_asv_relative_abundance_df$ASV <- rownames(mean_asv_relative_abundance_df) ##adding a column for ASV labels

##now data frame has ASV names, prevalence, relative abundance across samples, mean relative abundance for each ASV
prev_meanRA_plot_df  <- merge(prev_ra_upset_df , mean_asv_relative_abundance_df, by= "ASV") 

###Finding out which ASVs are only in individual samples---
##"apply" a function to the rows (indicated by the "1" argument) of the selected subset (rows dna pool, raw pool, individuals) of the data frame.
#function(row){ } defines a function to be applied to each row (ASVs)

prev_meanRA_plot_df$individual.only<- apply(prev_meanRA_plot_df [, c('DNA Pool', 'Raw Pool', 'Individuals')], 1, function(row) {
  if (row['Individuals'] == 1 && row['DNA Pool'] == 0 && row['Raw Pool'] == 0) {
    return('Y')
  } else if (row['DNA Pool'] == 1 || row['Raw Pool'] == 1) {
    return('N')
  }
})
prev_meanRA_plot_df 

###Plotting#####
prev_meanRA_plot <-ggplot (prev_meanRA_plot_df , aes(x=Prev_percentage, y=mean_relative_abundance, color= individual.only))+ 
  geom_point()+
  scale_y_log10(labels = scales::number_format(accuracy = 0.001), breaks = c(0,0.001,0.01, 0.1, 1.0, 10.0))+
  scale_x_continuous(labels = scales::percent_format(), 
                     breaks = c(0, 0.05, 0.25, 0.50, 0.75, 1))+
  theme_bw()+
  scale_color_manual(values= c("red", "black")) + 
  labs(x= "Prevalence (% of samples)", y = "Mean Relative Abundance (%)", color = "Individuals\nonly")+
  theme(
    axis.title = element_text(size=20),
    axis.text.y = element_text(size=18),
    axis.text.x = element_text(size =11),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    panel.grid.minor = element_blank()
  )
prev_meanRA_plot
ggsave("prevalence vs Abundance (ASV).tiff", plot = prev_meanRA_plot, device = "tiff", width = 12, height = 6, dpi = 300)

###Prevalence and abundance ASV calculations ######
##Which individual ASVs are low abundance (less than 0.01% RA) and low prevalence (less than 5%)
low_prv5_low_abun0.01_individual_ASV <- subset(prev_meanRA_plot_df, mean_relative_abundance <= 0.01 & Prev_percentage <= 0.05 & individual.only == "Y")
nrow(low_prv5_low_abun0.01_individual_ASV) #10931 ASVs with less than 0.01% RA and prevalence less than 5%
individual_ASV <- subset(prev_meanRA_plot_df, individual.only == "Y")
nrow(individual_ASV) #11510 ASVs only in individual samples 
percentage_low_prv5_low_abun0.01_individual_ASV <- (nrow(low_prv5_low_abun0.01_individual_ASV) / nrow(individual_ASV)) * 100
percentage_low_prv5_low_abun0.01_individual_ASV ##94.97% of ASVs

#which individual ASVs are low abundance (less than 0.001% RA) and low prevalence (less than 5%)
low_prv5_low_abun0.001_individual_ASV <- subset(prev_meanRA_plot_df, mean_relative_abundance <= 0.001 & Prev_percentage <= 0.05 & individual.only == "Y")
nrow(low_prv5_low_abun0.001_individual_ASV) #10878 ASVs with less than 0.001% RA and prevalence less than 5%
percentage_low_prv5_low_abun0.001_individual_ASV <- (nrow(low_prv5_low_abun0.001_individual_ASV) / nrow(individual_ASV)) * 100
percentage_low_prv5_low_abun0.001_individual_ASV ##94.51% of ASVs

#which individual ASVs are low abundance (less than 0.1% RA) and low prevalence (less than 5%)
low_prv5_low_abun0.1_individual_ASV <- subset(prev_meanRA_plot_df, mean_relative_abundance <= 0.1 & Prev_percentage <= 0.05 & individual.only == "Y")
nrow(low_prv5_low_abun0.1_individual_ASV) #10932 ASVs with less than 0.1% RA and prevalence less than 5%
percentage_low_prv5_low_abun0.1_individual_ASV <- (nrow(low_prv5_low_abun0.1_individual_ASV) / nrow(individual_ASV)) * 100
percentage_low_prv5_low_abun0.1_individual_ASV ##94.98% of ASVs

#which individual ASVs are low abundance (less than 0.01% RA) and low prevalence (less than 10%)
low_prv10_low_abun0.01_individual_ASV <- subset(prev_meanRA_plot_df, mean_relative_abundance <= 0.01 & Prev_percentage <= 0.1 & individual.only == "Y")
nrow(low_prv10_low_abun0.01_individual_ASV) #11212 ASVs with less than 0.01% RA and prevalence less than 10%
percentage_low_prv10_low_abun0.01_individual_ASV <- (nrow(low_prv10_low_abun0.01_individual_ASV) / nrow(individual_ASV)) * 100
percentage_low_prv10_low_abun0.01_individual_ASV ##97.41% of ASVs

#which individual ASVs are low abundance (less than 0.01% RA) and low prevalence (less than 25%)
low_prv25_low_abun0.01_individual_ASV <- subset(prev_meanRA_plot_df, mean_relative_abundance <= 0.01 & Prev_percentage <= 0.25 & individual.only == "Y")
nrow(low_prv25_low_abun0.01_individual_ASV) #11425 ASVs with less than 0.01% RA and prevalence less than 25%
percentage_low_prv25_low_abun0.01_individual_ASV <- (nrow(low_prv25_low_abun0.01_individual_ASV) / nrow(individual_ASV)) * 100
percentage_low_prv25_low_abun0.01_individual_ASV ##99.26% of ASVs

#which individual ASVs are low abundance (less than 0.01% RA) and prevalence (less than 50%)
low_prv50_low_abun0.01_individual_ASV <- subset(prev_meanRA_plot_df, mean_relative_abundance <= 0.01 & Prev_percentage <= 0.50 & individual.only == "Y")
nrow(low_prv50_low_abun0.01_individual_ASV) #11489 ASVs with less than 0.01% RA and prevalence less than 50%
percentage_low_prv50_low_abun0.01_individual_ASV <- (nrow(low_prv50_low_abun0.01_individual_ASV) / nrow(individual_ASV)) * 100
percentage_low_prv50_low_abun0.01_individual_ASV ##99.82% of ASVs


##POOL PRESENCE PLOT######
##Finding out which ASVs are present for only one sample type or shared among them
prev_meanRA_plot_df$pool_presence <- apply(prev_meanRA_plot_df[, c('DNA Pool', 'Raw Pool', 'Individuals')], 1, function(row) {
  if (row['Individuals'] == 1 && row['DNA Pool'] == 0 && row['Raw Pool'] == 0) {
    return('Only Individuals')
  } else if (row['Individuals'] == 0 && row['DNA Pool'] == 1 && row['Raw Pool'] == 0) {
    return('Only DNA Pools')
  } else if (row['Individuals'] == 0 && row['DNA Pool'] == 0 && row['Raw Pool'] == 1) {
    return('Only Raw Pools')
  } else if (row['Individuals'] == 1 && row['DNA Pool'] == 1 && row['Raw Pool'] == 1) {
    return('Shared by all sample types')
  } else if (row['Individuals'] == 0 && row['DNA Pool'] == 1 && row['Raw Pool'] == 1) {
    return('Shared by DNA and Raw pools')
  } else if (row['Individuals'] == 1 && row['DNA Pool'] == 1 && row['Raw Pool'] == 0) {
    return('Shared by Individuals and DNA pools')
  }else if (row['Individuals'] == 1 && row['DNA Pool'] == 0 && row['Raw Pool'] == 1) {
    return('Shared by Individuals and Raw pools')
  }else {
    return('Other')
  }
})
prev_meanRA_plot_df

##Plotting
prev_meanRA_plot_pool_presence <-ggplot (prev_meanRA_plot_df , aes(x=Prev_percentage, y=mean_relative_abundance, color= pool_presence))+ 
  geom_point()+
  scale_y_log10(labels = scales::number_format(accuracy = 0.001), breaks = c(0,0.001,0.01, 0.1, 1.0, 10.0))+
  scale_x_continuous(labels = scales::percent_format())+
  theme_bw()+
  scale_color_manual(values = paletteer::paletteer_d("MetBrewer::Archambault")) +
  labs(x= "Prevalence (% of samples)", y = "Mean Relative Abundance (%)", color = "Presence")+
  theme(
    axis.title = element_text(size=20),
    axis.text.y = element_text(size=18),
    axis.text.x = element_text(size =18),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    panel.grid.minor = element_blank()
  )
prev_meanRA_plot_pool_presence

####SUPPLEMENTARY FIGURE 5######
sfigure5 <- prev_meanRA_plot_pool_presence
ggsave("SupplementaryFigure5.tiff", plot = prev_meanRA_plot_pool_presence, device = "tiff", width = 10, height = 5, dpi = 300)

##Maths- which shared ASVs have higher abundance (more than 0.001% RA) and  prevalence (more than 20%)
high_prv20_high_abun_shared_ASV <- subset(prev_meanRA_plot_df, mean_relative_abundance >= 0.001 & Prev_percentage >= 0.20 & pool_presence == "Shared by all sample types")
nrow(high_prv20_high_abun_shared_ASV) #3556 ASVs with equal to or more than 0.001% RA and prevalence equal to or more than 20%
shared_ASV <- subset(prev_meanRA_plot_df, pool_presence == "Shared by all sample types")
nrow(shared_ASV) #10741 ASVs shared by all sample types
percentage_high_prv20_high_abun_shared_ASV <- (nrow(high_prv20_high_abun_shared_ASV) / nrow(shared_ASV)) * 100
percentage_high_prv20_high_abun_shared_ASV ##33.1% of ASVs


# RA BY SAMPLE TYPE ####
##PHYLUM#######
data4_phylum.ra ##54 phyla
##merging low abundance 
data4_phylum.ra_filt<- merge_low_abundance_ra(data4_phylum.ra, threshold = 0.5) 
data4_phylum.ra_filt ##8 phyla with abundance more than 0.5%
##melting
data4_phylum.ra_filt_melt <- psmelt(data4_phylum.ra_filt) ##8 phyla (including unclassified unassigned and Others)
##ordering categorical pool type
data4_phylum.ra_filt_melt$pool_order <- factor (data4_phylum.ra_filt_melt$pool_type, levels = c("individual", "dna", "raw"))

##Factoring Phylum so that "Others" are the last level (this is necessary if you want the "Others" bar at the bottom of the bar plot)
data4_phylum.ra_filt_melt <- data4_phylum.ra_filt_melt %>%
  mutate(Phylum = factor(Phylum, 
                         levels = c(setdiff(Phylum, 
                                            unique(grep("Others", Phylum, value = TRUE))), 
                                    unique(grep("Others", Phylum, value = TRUE))))) 

##Color Palette
phylum.filt.palette <- brewer.pal(n = 8, name = "Dark2") 
names(phylum.filt.palette) <- unique(data4_phylum.ra_filt_melt$Phylum) ##Assign colors from the palette to the taxa(phylum) names
others_phylum <- grepl("Others", names(phylum.filt.palette)) # Find names containing "Others"
phylum.filt.palette[others_phylum] <- "grey95" # Assign "grey95" to "Others.."

##plotting, facet grid to separate by sample type (pool or individual)
ra_phylum_sample_plot <- ggplot(data4_phylum.ra_filt_melt, aes (x=sample_ID , y=Abundance, fill= Phylum))+
  theme_minimal()+
  geom_bar(stat="identity", colour= "black")+ 
  facet_grid (~pool_order, scales = "free_x", space = "free_x",
              labeller = labeller(pool_order = c(`individual` = "Individual", `dna` = "DNA Pools", `raw` = "Raw Pools"))) + ##relabeling
  scale_y_continuous(expand= c(0.010, 0, 0.010, 0))+
  scale_x_discrete()+ 
  scale_fill_manual(values = phylum.filt.palette) +
  labs(y = "Relative Abundance (%)")+
  theme(legend.title = element_text(face = "bold", size= 18),
        legend.text = element_text(size = 18),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), # leaving some space for the legend
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(), 
        axis.text.x= element_blank(),
        strip.text = element_text(size = 24))
ra_phylum_sample_plot
##SUPPLEMENTARY FIGURE 6####
sfigure6 <- ra_phylum_sample_plot
ggsave("SupplementaryFigure6.tiff", plot = sfigure6, device = "tiff", bg = "white", width = 18, height = 10, dpi = 300)

# ##ORDER
# data4_order.ra #390 orders
# data4_order.ra_filt <- merge_low_abundance(data4_order.ra, threshold = 0.3) ##filtering the orders, to not have that many for a plot, filtering threshold 0.3, if below this it gets merged into 1 order
# data4_order.ra_filt ##33 taxa (order) with relative abundance over 0.3
# data4_order.ra_filt_melt <- psmelt(data4_order.ra_filt)## melted df
# data4_order.ra_filt_melt$pool_order <- factor (data4_order.ra_filt_melt $pool_type, levels = c("individual", "dna", "raw"))
# 
# ##plotting, facet grid to separate by sample type (pool or individual)
# ra_order_sample_plot <- ggplot(data4_order.ra_filt_melt, aes (x=sample_ID , y=Abundance, fill= Order))+
#   theme_minimal()+
#   geom_bar(stat="identity", colour= "black")+ 
#   facet_grid (~pool_order, scales = "free_x", space = "free_x",
#               labeller = labeller(pool_order = c(`individual` = "Individual", `dna` = "DNA Pools", `raw` = "Raw Pools"))) + ##relabeling
#   scale_y_continuous(expand= c(0.010, 0, 0.010, 0))+
#   scale_x_discrete()+ 
#   scale_fill_manual(values= ra_dendro_order_palette)+
#   labs(y = "Relative Abundance (%)")+
#   theme(legend.title = element_text(face = "bold", size= 16),
#         legend.text = element_text(size = 12),
#         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), # leaving some space for the legend
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         axis.line.y = element_line(linewidth = 0.7, colour = "black"),
#         axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
#         axis.title.y = element_text(size = 24),
#         axis.title.x = element_blank(), 
#         axis.text.x= element_blank(),
#         strip.text = element_text(size = 24))
# ra_order_sample_plot
# ggsave("ra_order_sample_tiff", plot = ra_order_sample_plot, device = "tiff", width = 18, height = 10, dpi = 300)
# 
# 
# ###Family 
# data4_family.ra #771 families
# data4_family.ra_filt <- merge_low_abundance(data4_family.ra, threshold = 0.5) ##filtering the families, to not have that many for a plot, filtering threshold 0.3, if below this it gets merged into 1 family
# data4_family.ra_filt ##31 taxa (family) with relative abundance over 0.5
# data4_family.ra_filt_melt <- psmelt(data4_family.ra_filt)## melted df
# data4_family.ra_filt_melt$pool_order <- factor (data4_family.ra_filt_melt $pool_type, levels = c("individual", "dna", "raw"))
# 
# ##plotting, facet grid to separate by sample type (pool or individual)
# ra_family_sample_plot <- ggplot(data4_family.ra_filt_melt, aes (x=sample_ID , y=Abundance, fill= Family))+
#   theme_minimal()+
#   geom_bar(stat="summary", colour= "black")+ 
#   facet_grid (~pool_order, scales = "free_x", space = "free_x",
#               labeller = labeller(pool_order = c(`individual` = "Individual", `dna` = "DNA Pools", `raw` = "Raw Pools"))) + ##relabeling
#   scale_y_continuous(expand= c(0.010, 0, 0.010, 0))+
#   scale_x_discrete()+ 
#   scale_fill_manual(values = distinctColorPalette(50))+
#   labs(y = "Relative Abundance (%)")+
#   theme(legend.title = element_text(face = "bold", size= 16),
#         legend.text = element_text(size = 12),
#         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), # leaving some space for the legend
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         axis.line.y = element_line(linewidth = 0.7, colour = "black"),
#         axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
#         axis.title.y = element_text(size = 24),
#         axis.title.x = element_blank(), 
#         axis.text.x= element_blank(),
#         strip.text = element_text(size = 24))
# ra_family_sample_plot
# ggsave("ra_family_sample_tiff", plot = ra_family_sample_plot, device = "tiff", width = 18, height = 10, dpi = 300)


##GENUS ##########
data4_genus.ra #2200 genera
data4_genus.ra_filt <- merge_low_abundance_ra(data4_genus.ra, threshold = 0.3) ##filtering the genera, to not have that many for a plot, filtering threshold 0.3, if below this it gets merged into 1 genus
data4_genus.ra_filt ## 48 genera with relative abundance over 0.3
data4_genus.ra_filt_melt <- psmelt(data4_genus.ra_filt) ## melted df
data4_genus.ra_filt_melt$pool_order <- factor (data4_genus.ra_filt_melt $pool_type, levels = c("individual", "dna", "raw"))

##Factoring genus so that "Others" are the last level (this is necessary if you want the "Others" bar at the bottom of the bar plot)
data4_genus.ra_filt_melt <- data4_genus.ra_filt_melt %>%
  mutate(Genus = factor(Genus, 
                        levels = c(setdiff(Genus, 
                                           unique(grep("Others", Genus, value = TRUE))), 
                                   unique(grep("Others", Genus, value = TRUE))))) 

##Color Palette
genus.filt.palette <- randomColor(count =48) 
names(genus.filt.palette) <- unique(data4_genus.ra_filt_melt$Genus) ##Assign colors from the palette to the taxa(genus) names
others_genus <- grepl("Others", names(genus.filt.palette)) # Find names containing "Others"
genus.filt.palette[others_genus] <- "grey95" # Assign "grey95" to "Others.."

##Top genera for the legend
top_genera <- top_taxa_legend(data4_genus.ra_filt_melt, taxlevel = "Genus", n = 10)
##Plotting
ra_genus_plot <- ggplot(data4_genus.ra_filt_melt, aes (x=sample_ID , y=Abundance, fill= Genus))+
  theme_minimal()+
  geom_bar(stat="identity", colour= "black")+ 
  facet_grid (~pool_order, scales = "free_x", space = "free_x",
              labeller = labeller(pool_order = c(`individual` = "Individual", `dna` = "DNA Pools", `raw` = "Raw Pools"))) + ##relabeling
  scale_y_continuous(expand= c(0.010, 0, 0.010, 0))+
  scale_fill_manual(values= genus.filt.palette, breaks = top_genera) +
  labs(y = "Relative Abundance (%)")+
  theme(legend.title = element_text(face = "bold", size= 16),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), # leaving some space for the legend
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(linewidth = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(), 
        axis.text.x= element_blank(),
        strip.text = element_text(size = 24))
ra_genus_plot
ggsave("ra_genus.png", plot = ra_genus_plot, device = png, width = 18, height = 10, dpi = 300)


# BETA DIVERSITY (AT ALL TAXA LEVELS)#######
## PAIRWISE COMPARISONS GUNIFRAC#######
### PHYLUM ######
###distance matrix
data4phylum.css <- tax_glom(data4.css, taxrank = "Phylum", NArm=FALSE) ##NA false keeping what is unclassified
data4phylum.css.gunifrac <- gunifrac(data4phylum.css)
data4phylum.css.gunifrac
##ordinate with NMDS
data4phylum.css.gunifrac.ord <- metaMDS(comm = data4phylum.css.gunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4phylum.css.gunifrac.scrs <- scores(data4phylum.css.gunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4phylum.css.gunifrac.scrs <- cbind(as.data.frame(data4phylum.css.gunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4phylum.css.gunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4phylum.css.gunifrac.scrs, FUN = mean) ## Calculating the centroids (mean method)
data4phylum.css.gunifrac.segs <- merge(data4phylum.css.gunifrac.scrs, setNames(data4phylum.css.gunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4phylum.css.gunifrac.segs <- data4phylum.css.gunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4phylum.css.gunifrac.segs

##ordering categorical pool type
data4phylum.css.gunifrac.segs$pool_order <- factor (data4phylum.css.gunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
gunifrac.phylum.adonis <- adonis2(data4phylum.css.gunifrac~ pool_type, 
                                  data = data4.css.df, by = "margin", permutations = 9999)
gunifrac.phylum.adonis
## 10.92% of the variation in composition is due to pool type. p value = 0.0021

#### PAIRWISE PERMANOVA ######
set.seed(87)
gunifrac.phylum.pairwise.adonis <- pairwise.adonis2(data4phylum.css.gunifrac ~ pool_type, 
                                                    data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
gunifrac.phylum.pairwise.adonis 
##difference in community structure between individual vs raw (p= 0.0102) and individual vs dna (p= 0.0059)
## no significant difference between dna pool and raw pool (p=0.684)

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.gunifrac.phylum.disp <- betadisper(data4phylum.css.gunifrac, data4.css.df$pool_type)
pooltype.gunifrac.phylum.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.gunifrac.phylum.permdisp <- permutest(pooltype.gunifrac.phylum.disp, permutations = 9999, pairwise = 1)
pooltype.gunifrac.phylum.permdisp
## dispersions of variances are the same (p=0.47)

###CLASS######
###distance matrix
data4class.css <- tax_glom(data4.css, taxrank = "Class", NArm=FALSE) ##NA false keeping what is unclassified
data4class.css.gunifrac <- gunifrac(data4class.css)
data4class.css.gunifrac
##ordinate with NMDS
data4class.css.gunifrac.ord <- metaMDS(comm = data4class.css.gunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4class.css.gunifrac.scrs <- scores(data4class.css.gunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4class.css.gunifrac.scrs <- cbind(as.data.frame(data4class.css.gunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4class.css.gunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4class.css.gunifrac.scrs, FUN = mean) ## Calculating the centroids (mean method)
data4class.css.gunifrac.segs <- merge(data4class.css.gunifrac.scrs, setNames(data4class.css.gunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4class.css.gunifrac.segs <- data4class.css.gunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4class.css.gunifrac.segs

##ordering categorical pool type
data4class.css.gunifrac.segs$pool_order <- factor (data4class.css.gunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
gunifrac.class.adonis <- adonis2(data4class.css.gunifrac~ pool_type, 
                                 data = data4.css.df, by = "margin", permutations = 9999)
gunifrac.class.adonis
## 9% of the variation in composition is due to pool type. p value = 0.0052

#### PAIRWISE PERMANOVA ######
set.seed(87)
gunifrac.class.pairwise.adonis <- pairwise.adonis2(data4class.css.gunifrac ~ pool_type,
                                                   data = data4.css.df, 
                                                   by = "margin", permutations = 9999, p.adjust.method = "BH")
gunifrac.class.pairwise.adonis 
##difference in community structure between individual vs raw (p= 0.0174) and individual vs dna (p= 0.009)
## no significant difference between dna pool and raw pool (p=0.6994)


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.gunifrac.class.disp <- betadisper(data4class.css.gunifrac, data4.css.df$pool_type)
pooltype.gunifrac.class.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.gunifrac.class.permdisp <- permutest(pooltype.gunifrac.class.disp, permutations = 9999, pairwise = 1)
pooltype.gunifrac.class.permdisp
## dispersion of variances are the same (p=0.265)


### ORDER ######
###distance matrix
data4order.css <- tax_glom(data4.css, taxrank = "Order", NArm=FALSE) ##NA false keeping what is unorderified
data4order.css.gunifrac <- gunifrac(data4order.css)
data4order.css.gunifrac

##ordinate with NMDS
data4order.css.gunifrac.ord <- metaMDS(comm = data4order.css.gunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4order.css.gunifrac.scrs <- scores(data4order.css.gunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4order.css.gunifrac.scrs <- cbind(as.data.frame(data4order.css.gunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4order.css.gunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4order.css.gunifrac.scrs, FUN = mean) ## Calculating the centroids (mean method)
data4order.css.gunifrac.segs <- merge(data4order.css.gunifrac.scrs, setNames(data4order.css.gunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4order.css.gunifrac.segs <- data4order.css.gunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4order.css.gunifrac.segs

##ordering categorical pool type
data4order.css.gunifrac.segs$pool_order <- factor (data4order.css.gunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

####PERMANOVA ######
set.seed(87)
gunifrac.order.adonis <- adonis2(data4order.css.gunifrac~ pool_type, data = data4.css.df, by = "margin", permutations = 9999)
gunifrac.order.adonis
## 7.3% of the variation in composition is due to pool type. p value = 0.018

#### PAIRWISE PERMANOVA ######
set.seed(87)
gunifrac.order.pairwise.adonis <- pairwise.adonis2(data4order.css.gunifrac ~ pool_type, data = data4.css.df, 
                                                   by = "margin", permutations = 9999, p.adjust.method = "BH")
gunifrac.order.pairwise.adonis 
##difference in community structure between individual vs raw (p= 0.0364) and individual vs dna (p= 0.025)
## no significant difference between dna pool and raw pool (p=0.73)


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.gunifrac.order.disp <- betadisper(data4order.css.gunifrac, data4.css.df$pool_type)
pooltype.gunifrac.order.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.gunifrac.order.permdisp <- permutest(pooltype.gunifrac.order.disp, permutations = 9999, pairwise = 1)
pooltype.gunifrac.order.permdisp
## dispersions of variances are the same (p=0.1293)

### FAMILY ######
###distance matrix
data4family.css <- tax_glom(data4.css, taxrank = "Family", NArm=FALSE) ##NA false keeping what is unfamilyified
data4family.css.gunifrac <- gunifrac(data4family.css)
data4family.css.gunifrac
##ordinate with NMDS
data4family.css.gunifrac.ord <- metaMDS(comm = data4family.css.gunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4family.css.gunifrac.scrs <- scores(data4family.css.gunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4family.css.gunifrac.scrs <- cbind(as.data.frame(data4family.css.gunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4family.css.gunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4family.css.gunifrac.scrs, FUN = mean) ## Calculating the centroids (mean method)
data4family.css.gunifrac.segs <- merge(data4family.css.gunifrac.scrs, setNames(data4family.css.gunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4family.css.gunifrac.segs <- data4family.css.gunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4family.css.gunifrac.segs

##familying categorical pool type
data4family.css.gunifrac.segs$pool_order<- factor (data4family.css.gunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

####PERMANOVA ###### 
set.seed(87)
gunifrac.family.adonis <- adonis2(data4family.css.gunifrac~ pool_type, data = data4.css.df, by = "margin", permutations = 9999)
gunifrac.family.adonis
## 6.8% of the variation in composition is due to pool type. p value = 0.02

#### PAIRWISE PERMANOVA ######
set.seed(87)
gunifrac.family.pairwise.adonis <- pairwise.adonis2(data4family.css.gunifrac ~ pool_type, 
                                                    data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
gunifrac.family.pairwise.adonis 
##difference in community structure between individual vs raw (p= 0.0468) and individual vs dna (p= 0.0334)
## no significant difference between dna pool and raw pool (p=0.7292)


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.gunifrac.family.disp <- betadisper(data4family.css.gunifrac, data4.css.df$pool_type)
pooltype.gunifrac.family.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.gunifrac.family.permdisp <- permutest(pooltype.gunifrac.family.disp, permutations = 9999, pairwise = 1)
pooltype.gunifrac.family.permdisp
## dispersions of variances are the same (p=0.1)


###GENUS ######
###distance matrix
data4genus.css <- tax_glom(data4.css, taxrank = "Genus", NArm=FALSE) ##NA false keeping what is ungenusified
#data4genis.css2 <- prune_taxa(taxa_sums(data4genus.css) > 0, data4genus.css)
data4genus.css.gunifrac <- gunifrac(data4genus.css)
data4genus.css.gunifrac
##ordinate with NMDS
data4genus.css.gunifrac.ord <- metaMDS(comm = data4genus.css.gunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4genus.css.gunifrac.scrs <- scores(data4genus.css.gunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4genus.css.gunifrac.scrs <- cbind(as.data.frame(data4genus.css.gunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4genus.css.gunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4genus.css.gunifrac.scrs, FUN = mean) ## Calculating the centroids (mean method)
data4genus.css.gunifrac.segs <- merge(data4genus.css.gunifrac.scrs, setNames(data4genus.css.gunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4genus.css.gunifrac.segs <- data4genus.css.gunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4genus.css.gunifrac.segs

##genusing categorical pool type
data4genus.css.gunifrac.segs$pool_order <- factor (data4genus.css.gunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

####PERMANOVA ######
set.seed(87)
gunifrac.genus.adonis <- adonis2(data4genus.css.gunifrac~ pool_type, data = data4.css.df, by = "margin", permutations = 9999)
gunifrac.genus.adonis
## 7.02% of the variation in composition is due to pool type. p value = 0.0173

#### PAIRWISE PERMANOVA ######
set.seed(87)
gunifrac.genus.pairwise.adonis <- pairwise.adonis2(data4genus.css.gunifrac ~ pool_type, 
                                                   data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
gunifrac.genus.pairwise.adonis 
##difference in community structure between individual vs raw (p= 0.035) and individual vs dna (p= 0.0241)
## no significant difference between dna pool and raw pool (p=0.7536)


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.gunifrac.genus.disp <- betadisper(data4genus.css.gunifrac, data4.css.df$pool_type)
pooltype.gunifrac.genus.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.gunifrac.genus.permdisp <- permutest(pooltype.gunifrac.genus.disp, permutations = 9999, pairwise = 1)
pooltype.gunifrac.genus.permdisp
## dispersions of variances are the same (p=0.0557)


###Combining PERMANOVA results at different taxonomic levels into one data frame#######
####ASV######
####Overall PERMANOVA 
gunifrac.adonis.df <- as.data.frame(gunifrac.adonis)
gunifrac.adonis.df ##dataframe of PERMANOVA result
gunifrac.adonis.df <- gunifrac.adonis.df %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "ASV") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Generalized UniFrac") ##UniFrac distance identifier
gunifrac.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
gunifrac.asv.pairwise.adonis.df<- data.frame(gunifrac.pairwise.adonis, 
                                             check.names = F)

##Making a couple changes to make it a long format df 
gunifrac.asv.pairwise.adonis.df2 <- gunifrac.asv.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Generalized UniFrac",
         `Taxonomic Level` = "ASV")
gunifrac.asv.pairwise.adonis.df2 

####PHYLUM######
####Overall PERMANOVA 
gunifrac.phylum.adonis.df <- as.data.frame(gunifrac.phylum.adonis)
gunifrac.phylum.adonis.df ##dataframe of PERMANOVA result
gunifrac.phylum.adonis.df <- gunifrac.phylum.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Phylum") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Generalized UniFrac") ##UniFrac distance identifier
gunifrac.phylum.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
gunifrac.phylum.pairwise.adonis.df<- data.frame(gunifrac.phylum.pairwise.adonis, 
                                                check.names = F)

##Making a couple changes to make it a long format df 
gunifrac.phylum.pairwise.adonis.df2 <- gunifrac.phylum.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Generalized UniFrac",
         `Taxonomic Level` = "Phylum")
gunifrac.phylum.pairwise.adonis.df2


####CLASS######
####Overall PERMANOVA 
gunifrac.class.adonis.df <- as.data.frame(gunifrac.class.adonis)
gunifrac.class.adonis.df ##dataframe of PERMANOVA result
gunifrac.class.adonis.df <- gunifrac.class.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Class") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Generalized UniFrac") ##UniFrac distance identifier
gunifrac.class.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
gunifrac.class.pairwise.adonis.df<- data.frame(gunifrac.class.pairwise.adonis, 
                                               check.names = F)

##Making a couple changes to make it a long format df 
gunifrac.class.pairwise.adonis.df2 <- gunifrac.class.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Generalized UniFrac",
         `Taxonomic Level` = "Class")
gunifrac.class.pairwise.adonis.df2


####ORDER######
####Overall PERMANOVA 
gunifrac.order.adonis.df <- as.data.frame(gunifrac.order.adonis)
gunifrac.order.adonis.df 
gunifrac.order.adonis.df <- gunifrac.order.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Order") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Generalized UniFrac") ##UniFrac distance identifier
gunifrac.order.adonis.df 

##Extracting the pairwise comparisons dataframes [2],[3],[4]
gunifrac.order.pairwise.adonis.df<- data.frame(gunifrac.order.pairwise.adonis,
                                               check.names = F)

##Making a couple changes to make it a long format df 
gunifrac.order.pairwise.adonis.df2 <- gunifrac.order.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Generalized UniFrac",
         `Taxonomic Level` = "Order")
gunifrac.order.pairwise.adonis.df2 

####FAMILY######
####Overall PERMANOVA 
gunifrac.family.adonis.df <- as.data.frame(gunifrac.family.adonis)
gunifrac.family.adonis.df ##dataframe of PERMANOVA results
gunifrac.family.adonis.df <- gunifrac.family.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Family") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Generalized UniFrac") ##UniFrac distance identifier 
gunifrac.family.adonis.df 

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
gunifrac.family.pairwise.adonis.df<- data.frame (gunifrac.family.pairwise.adonis, 
                                                 check.names = F)

##Making a couple changes to make it a long format df 
gunifrac.family.pairwise.adonis.df2 <- gunifrac.family.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Generalized UniFrac",
         `Taxonomic Level` = "Family")
gunifrac.family.pairwise.adonis.df2 

####GENUS######
####Overall PERMANOVA 
gunifrac.genus.adonis.df <- as.data.frame(gunifrac.genus.adonis)
gunifrac.genus.adonis.df ##PERMANOVA results dataframe
gunifrac.genus.adonis.df <- gunifrac.genus.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Genus")%>% ##identifier of taxonomic level 
  mutate(UniFrac= "Generalized UniFrac") ##UniFrac distance identifier
gunifrac.genus.adonis.df 

####PAIRWISE PERMANOVA
##Extracting the pairwise comparisons dataframes [2],[3],[4]
gunifrac.genus.pairwise.adonis.df<- data.frame(gunifrac.genus.pairwise.adonis, 
                                               check.names = F)

##Making a couple changes to make it a long format df 
gunifrac.genus.pairwise.adonis.df2 <- gunifrac.genus.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Generalized UniFrac",
         `Taxonomic Level` = "Genus")
gunifrac.genus.pairwise.adonis.df2 

##Putting all the PERMANOVA pairwise results from the different taxonomic levels together
gunifrac.taxalevels.pairwise.adonis.df<- rbind(gunifrac.asv.pairwise.adonis.df2,
                                               gunifrac.phylum.pairwise.adonis.df2,
                                               gunifrac.class.pairwise.adonis.df2,
                                               gunifrac.order.pairwise.adonis.df2,
                                               gunifrac.family.pairwise.adonis.df2,
                                               gunifrac.genus.pairwise.adonis.df2)

gunifrac.taxalevels.pairwise.adonis.df ##All the PERMANOVA pairwise results at the taxonomic levels are now in one dataframe

#Editing df for plotting
gunifrac.taxalevels.pairwise.adonis.df2 <- gunifrac.taxalevels.pairwise.adonis.df %>%
  rename("P_value" = `Pr(>F)`)%>%
  mutate(pval_star = case_when(P_value < 0.01 ~ "**",
                               P_value <= 0.05 ~ "*",
                               P_value > 0.05 ~ "n.s.")) %>%###making another column for significance level
  mutate(pool_type_comparison= dplyr::recode(Comparison, 
                                      "raw_vs_individual"= "R vs I", "dna_vs_raw"= "D vs R", "dna_vs_individual"= "D vs I"))

##ordering (factoring) taxa levels and comparisons          
gunifrac.taxalevels.pairwise.adonis.df2$taxa_order <- factor (gunifrac.taxalevels.pairwise.adonis.df2$`Taxonomic Level`, levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))
gunifrac.taxalevels.pairwise.adonis.df2$comparisons_order <- factor (gunifrac.taxalevels.pairwise.adonis.df2$pool_type_comparison, levels = c("D vs R", "D vs I", "R vs I"))

## PAIRWISE COMPARISONS WUNIFRAC #######
###PHYLUM ######
###distance matrix
data4phylum.css ##Already did the tax glomming at the phylum level 
data4phylum.css.wunifrac <- wunifrac(data4phylum.css) ##wunifrac distances at phylum level
data4phylum.css.wunifrac

##ordinate with NMDS
data4phylum.css.wunifrac.ord <- metaMDS(comm = data4phylum.css.wunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##now, calculate centroid coordinates
data4phylum.css.wunifrac.scrs <- scores(data4phylum.css.wunifrac.ord$points) ## getting NMDS1 and NMDS2 coordinates
data4phylum.css.wunifrac.scrs <- cbind(as.data.frame(data4phylum.css.wunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4phylum.css.wunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4phylum.css.wunifrac.scrs, FUN = mean) ## calculating centroids
data4phylum.css.wunifrac.segs <- merge(data4phylum.css.wunifrac.scrs, setNames(data4phylum.css.wunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4phylum.css.wunifrac.segs <- data4phylum.css.wunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4phylum.css.wunifrac.segs

##ordering categorical pool type
data4phylum.css.wunifrac.segs$pool_order <- factor (data4phylum.css.wunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
wunifrac.phylum.adonis <- adonis2(data4phylum.css.wunifrac~ pool_type, 
                                  data = data4.css.df, by = "margin", permutations = 9999)
wunifrac.phylum.adonis
## 6.2% of the variation in composition is due to pool type. p value = 0.0823 (n.s.) 

#### PAIRWISE PERMANOVA ######
set.seed(87)
wunifrac.phylum.pairwise.adonis <- pairwise.adonis2(data4phylum.css.wunifrac ~ pool_type, 
                                                    data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
wunifrac.phylum.pairwise.adonis 
## no significant difference between any of them. dna pool vs raw rool (p=0.7373), 
#dna pool vs individuals (p=0.0936), raw pools vs individuals (p=0.098)

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.wunifrac.phylum.disp <- betadisper(data4phylum.css.wunifrac, data4.css.df$pool_type)
pooltype.wunifrac.phylum.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.wunifrac.phylum.permdisp <- permutest(pooltype.wunifrac.phylum.disp, permutations = 9999, pairwise = 1)
pooltype.wunifrac.phylum.permdisp
## dispersion of variances are the same (p=0.378)

###CLASS######
###distance matrix
data4class.css ## Already have tax gloming at the class level 
data4class.css.wunifrac <- wunifrac(data4class.css) ##wunifrac at class level
data4class.css.wunifrac

##ordinate with NMDS
data4class.css.wunifrac.ord <- metaMDS(comm = data4class.css.wunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4class.css.wunifrac.scrs <- scores(data4class.css.wunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4class.css.wunifrac.scrs <- cbind(as.data.frame(data4class.css.wunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4class.css.wunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4class.css.wunifrac.scrs, FUN = mean) ## calculating centroids (mean method)
data4class.css.wunifrac.segs <- merge(data4class.css.wunifrac.scrs, setNames(data4class.css.wunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4class.css.wunifrac.segs <- data4class.css.wunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4class.css.wunifrac.segs

##ordering categorical pool type
data4class.css.wunifrac.segs$pool_order <- factor (data4class.css.wunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
wunifrac.class.adonis <- adonis2(data4class.css.wunifrac~ pool_type, 
                                 data = data4.css.df, by = "margin", permutations = 9999)
wunifrac.class.adonis
## 6.33% of the variation in composition is due to pool type?. p value = 0.0574 (n.s.)

#### PAIRWISE PERMANOVA ######
set.seed(87)
wunifrac.class.pairwise.adonis <- pairwise.adonis2(data4class.css.wunifrac ~ pool_type, 
                                                   data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
wunifrac.class.pairwise.adonis 
##no significant difference in community structure between raw pools and individuals (p=0.0763),  
#dna pools and individuals (p=0.0785), or dna pools and raw pools(p=0.804) 

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.wunifrac.class.disp <- betadisper(data4class.css.wunifrac, data4.css.df$pool_type)
pooltype.wunifrac.class.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.wunifrac.class.permdisp <- permutest(pooltype.wunifrac.class.disp, permutations = 9999, pairwise = 1)
pooltype.wunifrac.class.permdisp
## dispersions of variances are the same (p=0.2599)


###ORDER ######
###distance matrix
data4order.css ##already did tax glomming at the order level
data4order.css.wunifrac <- wunifrac(data4order.css) #Wunifrac at the order level
data4order.css.wunifrac
##ordinate with NMDS
data4order.css.wunifrac.ord <- metaMDS(comm = data4order.css.wunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4order.css.wunifrac.scrs <- scores(data4order.css.wunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4order.css.wunifrac.scrs <- cbind(as.data.frame(data4order.css.wunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4order.css.wunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4order.css.wunifrac.scrs, FUN = mean) ## calculating centroids (mean method)
data4order.css.wunifrac.segs <- merge(data4order.css.wunifrac.scrs, setNames(data4order.css.wunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4order.css.wunifrac.segs <- data4order.css.wunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4order.css.wunifrac.segs

##ordering categorical pool type
data4order.css.wunifrac.segs$pool_order <- factor (data4order.css.wunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
wunifrac.order.adonis <- adonis2(data4order.css.wunifrac~ pool_type, 
                                 data = data4.css.df, by = "margin", permutations = 9999)
wunifrac.order.adonis
## 6.04% of the variation in composition is due to pool type. p value = 0.0542 (n.s.)

#### PAIRWISE PERMANOVA ######
set.seed(87)
wunifrac.order.pairwise.adonis <- pairwise.adonis2(data4order.css.wunifrac ~ pool_type, 
                                                   data = data4.css.df, by = "margin", 
                                                   permutations = 9999, p.adjust.method = "BH")
wunifrac.order.pairwise.adonis 
## no difference in community structure between dna vs raw pools (p=0.9065), 
#or dna pools vs individuals (p=0.086), raw pools and individuals (p=0.054)

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.wunifrac.order.disp <- betadisper(data4order.css.wunifrac, data4.css.df$pool_type)
pooltype.wunifrac.order.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.wunifrac.order.permdisp <- permutest(pooltype.wunifrac.order.disp, permutations = 9999, pairwise = 1)
pooltype.wunifrac.order.permdisp
## Variances are the same (p=0.14)

### FAMILY ######
###distance matrix
data4family.css ##already did tax glomming at the family level 
data4family.css.wunifrac <- wunifrac(data4family.css) ##wunifrac distances at the family level 
data4family.css.wunifrac
##ordinate with NMDS
data4family.css.wunifrac.ord <- metaMDS(comm = data4family.css.wunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4family.css.wunifrac.scrs <- scores(data4family.css.wunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4family.css.wunifrac.scrs <- cbind(as.data.frame(data4family.css.wunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4family.css.wunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4family.css.wunifrac.scrs, FUN = mean) ## calculating centroids(mean method)
data4family.css.wunifrac.segs <- merge(data4family.css.wunifrac.scrs, setNames(data4family.css.wunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4family.css.wunifrac.segs <- data4family.css.wunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4family.css.wunifrac.segs

##familying categorical pool type
data4family.css.wunifrac.segs$pool_order <- factor (data4family.css.wunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
wunifrac.family.adonis <- adonis2(data4family.css.wunifrac~ pool_type,
                                  data = data4.css.df, by = "margin", permutations = 9999)
wunifrac.family.adonis
##6.15% of the variation in composition is due to pool type. p value = 0.046

#### PAIRWISE PERMANOVA ######
set.seed(87)
wunifrac.family.pairwise.adonis <- pairwise.adonis2(data4family.css.wunifrac ~ pool_type, 
                                                    data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
wunifrac.family.pairwise.adonis 
## no significant difference between individual vs raw pools (p= 0.0463), 
#dna pool and raw pool (p=0.9033), and dna pool and individuals (p=0.084)

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.wunifrac.family.disp <- betadisper(data4family.css.wunifrac, data4.css.df$pool_type)
pooltype.wunifrac.family.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.wunifrac.family.permdisp <- permutest(pooltype.wunifrac.family.disp, permutations = 9999, pairwise = 1)
pooltype.wunifrac.family.permdisp
## Variances are the same (p=0.134)


###GENUS ######
###distance matrix
data4genus.css ##already did tax glomming at the genus level
data4genus.css.wunifrac <- wunifrac(data4genus.css) ##wunifrac at the genus level
data4genus.css.wunifrac
##ordinate with NMDS
data4genus.css.wunifrac.ord <- metaMDS(comm = data4genus.css.wunifrac, k =2, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4genus.css.wunifrac.scrs <- scores(data4genus.css.wunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4genus.css.wunifrac.scrs <- cbind(as.data.frame(data4genus.css.wunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4genus.css.wunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4genus.css.wunifrac.scrs, FUN = mean) ## calculating centroids (mean method)
data4genus.css.wunifrac.segs <- merge(data4genus.css.wunifrac.scrs, setNames(data4genus.css.wunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4genus.css.wunifrac.segs <- data4genus.css.wunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4genus.css.wunifrac.segs

##ordering categorical pool type
data4genus.css.wunifrac.segs$pool_order <- factor (data4genus.css.wunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
wunifrac.genus.adonis <- adonis2(data4genus.css.wunifrac~ pool_type, 
                                 data = data4.css.df, by = "margin", permutations = 9999)
wunifrac.genus.adonis
## 6.3% of the variation in composition is due to pool type. p value = 0.04

#### PAIRWISE PERMANOVA ######
set.seed(87)
wunifrac.genus.pairwise.adonis <- pairwise.adonis2(data4genus.css.wunifrac ~ pool_type, 
                                                   data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
wunifrac.genus.pairwise.adonis 
##difference in community structure between individual vs raw pools (p= 0.0411) 
## no significant difference between dna pool and raw pool (p=0.907), and individual vs dna pool (p= 0.0776)


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.wunifrac.genus.disp <- betadisper(data4genus.css.wunifrac, data4.css.df$pool_type)
pooltype.wunifrac.genus.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.wunifrac.genus.permdisp <- permutest(pooltype.wunifrac.genus.disp, permutations = 9999, pairwise = 1)
pooltype.wunifrac.genus.permdisp
## dispersion of variances are the same (p=0.1416)

###Combining PERMANOVA results at different taxonomic levels into one data frame######
####ASV######
####Overall PERMANOVA 
wunifrac.adonis.df <- as.data.frame(wunifrac.adonis)
wunifrac.adonis.df <- wunifrac.adonis.df %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "ASV") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Weighted UniFrac") ##UniFrac distance identifier
wunifrac.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons, objects [2],[3],[4]
wunifrac.asv.pairwise.adonis.df<- data.frame(wunifrac.pairwise.adonis, 
                                             check.names = F)

##Making a couple changes to make it a long format df 
wunifrac.asv.pairwise.adonis.df2 <- wunifrac.asv.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Weighted UniFrac",
         `Taxonomic Level` = "ASV")
wunifrac.asv.pairwise.adonis.df2 

####PHYLUM######
####Overall PERMANOVA 
wunifrac.phylum.adonis.df <- as.data.frame(wunifrac.phylum.adonis)
wunifrac.phylum.adonis.df ##dataframe of PERMANOVA result
wunifrac.phylum.adonis.df <- wunifrac.phylum.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Phylum") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Weighted UniFrac") ##UniFrac distance identifier
wunifrac.phylum.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
wunifrac.phylum.pairwise.adonis.df<- data.frame(wunifrac.phylum.pairwise.adonis, 
                                                check.names = F)

##Making a couple changes to make it a long format df 
wunifrac.phylum.pairwise.adonis.df2 <- wunifrac.phylum.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Weighted UniFrac",
         `Taxonomic Level` = "Phylum")
wunifrac.phylum.pairwise.adonis.df2 


####CLASS######
####Overall PERMANOVA 
wunifrac.class.adonis.df <- as.data.frame(wunifrac.class.adonis)
wunifrac.class.adonis.df ##dataframe of PERMANOVA result
wunifrac.class.adonis.df <- wunifrac.class.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Class") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Weighted UniFrac") ##UniFrac distance identifier
wunifrac.class.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
wunifrac.class.pairwise.adonis.df<- data.frame(wunifrac.class.pairwise.adonis, 
                                               check.names = F)

##Making a couple changes to make it a long format df 
wunifrac.class.pairwise.adonis.df2 <- wunifrac.class.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Weighted UniFrac",
         `Taxonomic Level` = "Class")
wunifrac.class.pairwise.adonis.df2 


####ORDER######
####Overall PERMANOVA 
wunifrac.order.adonis.df <- as.data.frame(wunifrac.order.adonis)
wunifrac.order.adonis.df ##dataframe of PERMANOVA result
wunifrac.order.adonis.df <- wunifrac.order.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Order") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Weighted UniFrac") ##UniFrac distance identifier
wunifrac.order.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
wunifrac.order.pairwise.adonis.df<- data.frame(wunifrac.order.pairwise.adonis, 
                                               check.names = F)

##Making a couple changes to make it a long format df 
wunifrac.order.pairwise.adonis.df2 <- wunifrac.order.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Weighted UniFrac",
         `Taxonomic Level` = "Order")
wunifrac.order.pairwise.adonis.df2 ##Long format!

####FAMILY######
####Overall PERMANOVA 
wunifrac.family.adonis.df <- as.data.frame(wunifrac.family.adonis)
wunifrac.family.adonis.df ##dataframe of PERMANOVA results
wunifrac.family.adonis.df <- wunifrac.family.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Family") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Weighted UniFrac") ##UniFrac distance identifier
wunifrac.family.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
wunifrac.family.pairwise.adonis.df<- data.frame(wunifrac.family.pairwise.adonis, 
                                                check.names = F)

##Making a couple changes to make it a long format df 
wunifrac.family.pairwise.adonis.df2 <- wunifrac.family.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Weighted UniFrac",
         `Taxonomic Level` = "Family")
wunifrac.family.pairwise.adonis.df2 


####GENUS######
####Overall PERMANOVA 
wunifrac.genus.adonis.df <- as.data.frame(wunifrac.genus.adonis)
wunifrac.genus.adonis.df ##dataframe of PERMANOVA results
wunifrac.genus.adonis.df <- wunifrac.genus.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Genus") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Weighted UniFrac") ##UniFrac distance identifier
wunifrac.genus.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes [2],[3],[4]
wunifrac.genus.pairwise.adonis.df<- data.frame(wunifrac.genus.pairwise.adonis, 
                                               check.names = F)

##Making a couple changes to make it a long format df 
wunifrac.genus.pairwise.adonis.df2 <- wunifrac.genus.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Weighted UniFrac",
         `Taxonomic Level` = "Genus")
wunifrac.genus.pairwise.adonis.df2 

##Putting all the PAIRWISE PERMANOVA  results from the different taxonomic levels together
wunifrac.taxalevels.pairwise.adonis.df<- rbind(wunifrac.asv.pairwise.adonis.df2,
                                               wunifrac.phylum.pairwise.adonis.df2,
                                               wunifrac.class.pairwise.adonis.df2,
                                               wunifrac.order.pairwise.adonis.df2,
                                               wunifrac.family.pairwise.adonis.df2,
                                               wunifrac.genus.pairwise.adonis.df2)

wunifrac.taxalevels.pairwise.adonis.df ##all the PERMANOVA pairwise results at the taxonomic levels are now in one dataframe

wunifrac.taxalevels.pairwise.adonis.df2 <- wunifrac.taxalevels.pairwise.adonis.df %>%
  rename("P_value" = `Pr(>F)`)%>%
  mutate(pval_star = case_when(P_value < 0.01 ~ "**",
                               P_value <= 0.05 ~ "*",
                               P_value > 0.05 ~ "n.s.")) %>%###making another column for significance level
  mutate(pool_type_comparison= recode(Comparison, "raw_vs_individual"= "R vs I", "dna_vs_raw"= "D vs R", "dna_vs_individual"= "D vs I"))

##ordering (factoring) the taxa levels and the comparisons
wunifrac.taxalevels.pairwise.adonis.df2$taxa_order <- factor (wunifrac.taxalevels.pairwise.adonis.df2$`Taxonomic Level`, 
                                                              levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))
wunifrac.taxalevels.pairwise.adonis.df2$comparisons_order <- factor (wunifrac.taxalevels.pairwise.adonis.df2$pool_type_comparison, 
                                                                     levels = c("D vs R", "D vs I", "R vs I"))
## PAIRWISE COMPARISONS UWUNIFRAC######
### PHYLUM ######
###distance matrix
data4phylum.css ##Already did the tax glomming at the phylum level 
data4phylum.css.uwunifrac <- uwunifrac(data4phylum.css) ##uwunifrac distances at phylum level
data4phylum.css.uwunifrac

##ordinate with NMDS
data4phylum.css.uwunifrac.ord <- metaMDS(comm = data4phylum.css.uwunifrac, k = 2, try = 20, trymax = 500, autotransform = F)

##now, calculate centroid coordinates
data4phylum.css.uwunifrac.scrs <- scores(data4phylum.css.uwunifrac.ord$points) ## getting NMDS1 and NMDS2 coordinates
data4phylum.css.uwunifrac.scrs <- cbind(as.data.frame(data4phylum.css.uwunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4phylum.css.uwunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4phylum.css.uwunifrac.scrs, FUN = mean) ## calculating centroids
data4phylum.css.uwunifrac.segs <- merge(data4phylum.css.uwunifrac.scrs, setNames(data4phylum.css.uwunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4phylum.css.uwunifrac.segs <- data4phylum.css.uwunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4phylum.css.uwunifrac.segs

##ordering categorical pool type
data4phylum.css.uwunifrac.segs$pool_order <- factor (data4phylum.css.uwunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
uwunifrac.phylum.adonis <- adonis2(data4phylum.css.uwunifrac~ pool_type, 
                                   data = data4.css.df,
                                   by = "margin", permutations = 9999)
uwunifrac.phylum.adonis
## 35.58% of the variation in composition is due to pool type. p value = 1e-04 

#### PAIRWISE PERMANOVA ######
set.seed(87)
uwunifrac.phylum.pairwise.adonis <- pairwise.adonis2(data4phylum.css.uwunifrac ~ pool_type, 
                                                     data = data4.css.df,
                                                     by = "margin", permutations = 9999, p.adjust.method = "BH")
uwunifrac.phylum.pairwise.adonis 
## difference in community composition between dna pools vs individuals (p=1e-04), raw pools vs individuals (p = 1e-04)
##no difference between dna and raw pools (p=0.6997)

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.uwunifrac.phylum.disp <- betadisper(data4phylum.css.uwunifrac, data4.css.df$pool_type)
pooltype.uwunifrac.phylum.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.uwunifrac.phylum.permdisp <- permutest(pooltype.uwunifrac.phylum.disp, 
                                                permutations = 9999, pairwise = 1)
pooltype.uwunifrac.phylum.permdisp
## dispersion of variances are NOT the same (p=0.0012), dna vs individual (p=0.0006), dna vs raw (p= 0.0127)

### CLASS ######
###distance matrix
data4class.css ## Already have tax gloming at the class level 
data4class.css.uwunifrac <- uwunifrac(data4class.css) ##uwunifrac at class level
data4class.css.uwunifrac

##ordinate with NMDS - had to do with k = 4 
data4class.css.uwunifrac.ord <- metaMDS(comm = data4class.css.uwunifrac, k =4, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4class.css.uwunifrac.scrs <- scores(data4class.css.uwunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4class.css.uwunifrac.scrs <- cbind(as.data.frame(data4class.css.uwunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4class.css.uwunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4class.css.uwunifrac.scrs, FUN = mean) ## calculating centroids (mean method)
data4class.css.uwunifrac.segs <- merge(data4class.css.uwunifrac.scrs, setNames(data4class.css.uwunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4class.css.uwunifrac.segs <- data4class.css.uwunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4class.css.uwunifrac.segs

##ordering categorical pool type
data4class.css.uwunifrac.segs$pool_order <- factor (data4class.css.uwunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
uwunifrac.class.adonis <- adonis2(data4class.css.uwunifrac~ pool_type, 
                                  data = data4.css.df, by = "margin", permutations = 9999)
uwunifrac.class.adonis
## 29.93% of the variation in composition is due to pool type. p value = 1e-04

#### PAIRWISE PERMANOVA ######
set.seed(87)
uwunifrac.class.pairwise.adonis <- pairwise.adonis2(data4class.css.uwunifrac ~ pool_type, 
                                                    data = data4.css.df, by = "margin", 
                                                    permutations = 9999, p.adjust.method = "BH")
uwunifrac.class.pairwise.adonis 
##difference in community structure between raw pools and individuals (p=1e-04), 
# and between dna pools and individuals (p=1e-04)
##no difference between dna pools and raw pools(p=0.9265) 

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.uwunifrac.class.disp <- betadisper(data4class.css.uwunifrac, data4.css.df$pool_type)
pooltype.uwunifrac.class.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.uwunifrac.class.permdisp <- permutest(pooltype.uwunifrac.class.disp, permutations = 9999, pairwise = 1)
pooltype.uwunifrac.class.permdisp
## Variances are not the same (p=1e-04), dna vs individual (p = 0.0007), raw vs individual (p = 0.0099)


### ORDER ######
###distance matrix
data4order.css ##already did tax glomming at the order level
data4order.css.uwunifrac <- uwunifrac(data4order.css) #uwunifrac at the order level
data4order.css.uwunifrac
##ordinate with NMDS- had to do with k = 4 
data4order.css.uwunifrac.ord <- metaMDS(comm = data4order.css.uwunifrac, k = 5, try = 20, trymax = 1000, autotransform = F)

##first, calculate centroid coordinates
data4order.css.uwunifrac.scrs <- scores(data4order.css.uwunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4order.css.uwunifrac.scrs <- cbind(as.data.frame(data4order.css.uwunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4order.css.uwunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4order.css.uwunifrac.scrs, FUN = mean) ## calculating centroids (mean method)
data4order.css.uwunifrac.segs <- merge(data4order.css.uwunifrac.scrs, setNames(data4order.css.uwunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4order.css.uwunifrac.segs <- data4order.css.uwunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4order.css.uwunifrac.segs

##ordering categorical pool type
data4order.css.uwunifrac.segs$pool_order <- factor (data4order.css.uwunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
uwunifrac.order.adonis <- adonis2(data4order.css.uwunifrac~ pool_type, 
                                  data = data4.css.df, by = "margin", permutations = 9999)
uwunifrac.order.adonis
## 34.37% of the variation in composition is due to pool type. p value = 1e-04

#### PAIRWISE PERMANOVA ######
set.seed(87)
uwunifrac.order.pairwise.adonis <- pairwise.adonis2(data4order.css.uwunifrac ~ pool_type, 
                                                    data = data4.css.df, by = "margin", 
                                                    permutations = 9999, p.adjust.method = "BH")
uwunifrac.order.pairwise.adonis 
##difference in community structure between raw pools and individuals (p=1e-04), 
#and between dna pools and individuals (p=1e-04)
##no difference between dna pools and raw pools(p=0.7562) 

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.uwunifrac.order.disp <- betadisper(data4order.css.uwunifrac, data4.css.df$pool_type)
pooltype.uwunifrac.order.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.uwunifrac.order.permdisp <- permutest(pooltype.uwunifrac.order.disp, permutations = 9999, pairwise = 1)
pooltype.uwunifrac.order.permdisp
## dispersions of variances are NOT the same (p=0.014), dna vs individual p = 0.03, individual vs raw p = 0.04

### FAMILY ######
###distance matrix
data4family.css ##already did tax glomming at the family level 
data4family.css.uwunifrac <- uwunifrac(data4family.css) ##uwunifrac distances at the family level 
data4family.css.uwunifrac
##ordinate with NMDS
data4family.css.uwunifrac.ord <- metaMDS(comm = data4family.css.uwunifrac, k =5, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4family.css.uwunifrac.scrs <- scores(data4family.css.uwunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4family.css.uwunifrac.scrs <- cbind(as.data.frame(data4family.css.uwunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4family.css.uwunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4family.css.uwunifrac.scrs, FUN = mean) ## calculating centroids(mean method)
data4family.css.uwunifrac.segs <- merge(data4family.css.uwunifrac.scrs, setNames(data4family.css.uwunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4family.css.uwunifrac.segs <- data4family.css.uwunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4family.css.uwunifrac.segs

##familying categorical pool type
data4family.css.uwunifrac.segs$pool_order <- factor (data4family.css.uwunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
uwunifrac.family.adonis <- adonis2(data4family.css.uwunifrac~ pool_type, 
                                   data = data4.css.df, by = "margin", permutations = 9999)
uwunifrac.family.adonis
## 31.55% of the variation in composition is due to pool type. p value = 1e-04

#### PAIRWISE PERMANOVA ###### 
set.seed(87)
uwunifrac.family.pairwise.adonis <- pairwise.adonis2(data4family.css.uwunifrac ~ pool_type, 
                                                     data = data4.css.df, by = "margin", permutations = 9999, p.adjust.method = "BH")
uwunifrac.family.pairwise.adonis 
##difference in community structure between raw pools and individuals (p=1e-04),
#and between dna pools and individuals (p=1e-04)
##no difference between dna pools and raw pools(p=0.6946) 

##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.uwunifrac.family.disp <- betadisper(data4family.css.uwunifrac, data4.css.df$pool_type)
pooltype.uwunifrac.family.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.uwunifrac.family.permdisp <- permutest(pooltype.uwunifrac.family.disp, permutations = 9999, pairwise = 1)
pooltype.uwunifrac.family.permdisp
## dispersions of variances are not the same (p=0.0013), dba vs individual (p=0.004), individual vs raw (p=0.0135)


### GENUS ######
###distance matrix
data4genus.css ##already did tax glomming at the genus level
#data4genis.css2 <- prune_taxa(taxa_sums(data4genus.css) > 0, data4genus.css)
data4genus.css.uwunifrac <- uwunifrac(data4genus.css) ##uwunifrac at the genus level
data4genus.css.uwunifrac
##ordinate with NMDS
data4genus.css.uwunifrac.ord <- metaMDS(comm = data4genus.css.uwunifrac, k =5, try = 20, trymax = 500, autotransform = F)

##first, calculate centroid coordinates
data4genus.css.uwunifrac.scrs <- scores(data4genus.css.uwunifrac.ord$points) ## gettind NMDS1 and NMDS2 coordinates
data4genus.css.uwunifrac.scrs <- cbind(as.data.frame(data4genus.css.uwunifrac.scrs), pool_type = data4.css@sam_data$pool_type) ##adding meta data
data4genus.css.uwunifrac.cent <- aggregate(cbind (MDS1, MDS2) ~ pool_type, data4genus.css.uwunifrac.scrs, FUN = mean) ## calculating centroids (mean method)
data4genus.css.uwunifrac.segs <- merge(data4genus.css.uwunifrac.scrs, setNames(data4genus.css.uwunifrac.cent, c("pool_type", "cMDS1", "cMDS2")), by = "pool_type", sort = F) ##adding the centroids to the dataframe with the coordinates

#making a new column with abbreviated pool type legend
data4genus.css.uwunifrac.segs <- data4genus.css.uwunifrac.segs %>%
  mutate(pool_type.abbrv = recode(pool_type, "individual"= "I", "dna"= "D", "raw"= "R"))
data4genus.css.uwunifrac.segs

##ordering categorical pool type
data4genus.css.uwunifrac.segs$pool_order <- factor (data4genus.css.uwunifrac.segs$pool_type, levels = c("individual", "dna", "raw"))

#### PERMANOVA ######
set.seed(87)
uwunifrac.genus.adonis <- adonis2(data4genus.css.uwunifrac~ pool_type, 
                                  data = data4.css.df, by = "margin", permutations = 9999)
uwunifrac.genus.adonis
## 30.12% of the variation in composition is due to pool type. p value = 1e-04

#### PAIRWISE PERMANOVA ######
set.seed(87)
uwunifrac.genus.pairwise.adonis <- pairwise.adonis2(data4genus.css.uwunifrac ~ pool_type, 
                                                    data = data4.css.df, by = "margin", 
                                                    permutations = 9999, p.adjust.method = "BH")
uwunifrac.genus.pairwise.adonis 
##difference in community structure between raw pools and individuals (p=1e-04), 
#and between dna pools and individuals (p=1e-04)
##no difference between dna pools and raw pools(p=0.7358) 


##PERMDISP- PERMANOVA assumption that variance between the populations is the same. We'll test it 
##first we calculate the average distance to centroid
pooltype.uwunifrac.genus.disp <- betadisper(data4genus.css.uwunifrac, data4.css.df$pool_type)
pooltype.uwunifrac.genus.disp

##then we test by permuting (shuffling the groups). Calculating if we can reject null hypothesis
set.seed(87)
pooltype.uwunifrac.genus.permdisp <- permutest(pooltype.uwunifrac.genus.disp, permutations = 9999, pairwise = 1)
pooltype.uwunifrac.genus.permdisp
## Variances are NOT the same (p=2e-04); dna vs individual ( p = 0.0003), individual vs raw (p = 0.0022)


###Combining PERMANOVA results at different taxonomic levels into one data frame#######
####ASV######
####Overall PERMANOVA 
uwunifrac.adonis.df <- as.data.frame(uwunifrac.adonis) ##dataframe of PERMANOVA result at ASV level
uwunifrac.adonis.df <- uwunifrac.adonis.df %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "ASV") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Unweighted UniFrac") ##UniFrac distance identifier
uwunifrac.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons, objects
uwunifrac.asv.pairwise.adonis.df<- data.frame(uwunifrac.pairwise.adonis,
                                              check.names = F)

##Making a couple changes to make it a long format df 
uwunifrac.asv.pairwise.adonis.df2 <- uwunifrac.asv.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Unweighted UniFrac",
         `Taxonomic Level` = "ASV")
uwunifrac.asv.pairwise.adonis.df2 

####PHYLUM######
####Overall PERMANOVA 
uwunifrac.phylum.adonis.df <- as.data.frame(uwunifrac.phylum.adonis)
uwunifrac.phylum.adonis.df ##dataframe of PERMANOVA result
uwunifrac.phylum.adonis.df <- uwunifrac.phylum.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Phylum") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Unweighted UniFrac") ##UniFrac distance identifier
uwunifrac.phylum.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes 
uwunifrac.phylum.pairwise.adonis.df<- data.frame(uwunifrac.phylum.pairwise.adonis, 
                                                 check.names = F)

##Making a couple changes to make it a long format df 
uwunifrac.phylum.pairwise.adonis.df2 <- uwunifrac.phylum.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Unweighted UniFrac",
         `Taxonomic Level` = "Phylum")
uwunifrac.phylum.pairwise.adonis.df2 

####CLASS######
####Overall PERMANOVA 
uwunifrac.class.adonis.df <- as.data.frame(uwunifrac.class.adonis)
uwunifrac.class.adonis.df ##dataframe of PERMANOVA result
uwunifrac.class.adonis.df <- uwunifrac.class.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Class") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Unweighted UniFrac") ##UniFrac distance identifier
uwunifrac.class.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes 
uwunifrac.class.pairwise.adonis.df<- data.frame(uwunifrac.class.pairwise.adonis,
                                                check.names = F)

##Making a couple changes to make it a long format df 
uwunifrac.class.pairwise.adonis.df2 <- uwunifrac.class.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Unweighted UniFrac",
         `Taxonomic Level` = "Class")
uwunifrac.class.pairwise.adonis.df2 


####ORDER######
####Overall PERMANOVA 
uwunifrac.order.adonis.df <- as.data.frame(uwunifrac.order.adonis)
uwunifrac.order.adonis.df ##dataframe of PERMANOVA result
uwunifrac.order.adonis.df <- uwunifrac.order.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Order") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Unweighted UniFrac") ##UniFrac distance identifier
uwunifrac.order.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA 
##Extracting the pairwise comparisons dataframes 
uwunifrac.order.pairwise.adonis.df<- data.frame(uwunifrac.order.pairwise.adonis, 
                                                check.names = F)

##Making a couple changes to make it a long format df 
uwunifrac.order.pairwise.adonis.df2 <- uwunifrac.order.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Unweighted UniFrac",
         `Taxonomic Level` = "Order")
uwunifrac.order.pairwise.adonis.df2 

####FAMILY######
####Overall PERMANOVA 
uwunifrac.family.adonis.df <- as.data.frame(uwunifrac.family.adonis)
str(uwunifrac.family.adonis.df) ##dataframe of PERMANOVA results
uwunifrac.family.adonis.df <- uwunifrac.family.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Family") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Unweighted UniFrac") ##UniFrac distance identifier
uwunifrac.family.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA
##Extracting the pairwise comparisons dataframes 
uwunifrac.family.pairwise.adonis.df<- data.frame(uwunifrac.family.pairwise.adonis, 
                                                 check.names = F)

##Making a couple changes to make it a long format df 
uwunifrac.family.pairwise.adonis.df2 <- uwunifrac.family.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Unweighted UniFrac",
         `Taxonomic Level` = "Family")
uwunifrac.family.pairwise.adonis.df2 


####GENUS######
####Overall PERMANOVA 
uwunifrac.genus.adonis.df <- as.data.frame(uwunifrac.genus.adonis)
uwunifrac.genus.adonis.df ##dataframe of PERMANOVA results
uwunifrac.genus.adonis.df <- uwunifrac.genus.adonis.df  %>%
  rownames_to_column(var = "Result")%>%
  mutate (Taxon_level = "Genus") %>% ##identifier of taxonomic level 
  mutate(UniFrac= "Unweighted UniFrac") ##UniFrac distance identifier
uwunifrac.genus.adonis.df #now I know which taxonomic level the results are from

####PAIRWISE PERMANOVA
##Extracting the pairwise comparisons dataframes
uwunifrac.genus.pairwise.adonis.df<- data.frame(uwunifrac.genus.pairwise.adonis, 
                                                check.names = F)

##Making a couple changes to make it a long format df 
uwunifrac.genus.pairwise.adonis.df2 <- uwunifrac.genus.pairwise.adonis.df %>%
  rownames_to_column(var="Effect") %>% 
  filter(Effect=="pool_type") %>%
  pivot_longer(
    cols = matches("_vs_"),
    names_to = c("Comparison", "Statistic"),
    names_sep = "\\.",
    values_to = "value")%>%
  pivot_wider(
    id_cols = c(Effect, Comparison),
    names_from  = Statistic,
    values_from = value)%>%
  mutate(UniFrac = "Unweighted UniFrac",
         `Taxonomic Level` = "Genus")
uwunifrac.genus.pairwise.adonis.df2 

##Putting all the PAIRWISE PERMANOVA  results from the different taxonomic levels together
uwunifrac.taxalevels.pairwise.adonis.df<- rbind(uwunifrac.asv.pairwise.adonis.df2,
                                                uwunifrac.genus.pairwise.adonis.df2,
                                                uwunifrac.family.pairwise.adonis.df2,
                                                uwunifrac.order.pairwise.adonis.df2,
                                                uwunifrac.class.pairwise.adonis.df2,
                                                uwunifrac.phylum.pairwise.adonis.df2)
uwunifrac.taxalevels.pairwise.adonis.df ##all the PERMANOVA pairwise results at the taxonomic levels are now in one dataframe

uwunifrac.taxalevels.pairwise.adonis.df2 <- uwunifrac.taxalevels.pairwise.adonis.df %>%
  rename("P_value" = `Pr(>F)`)%>%
  mutate(pval_star = case_when(P_value < 0.01 ~ "**",
                               P_value <= 0.05 ~ "*",
                               P_value > 0.05 ~ "n.s.")) %>%###making another column for significance level
  mutate(pool_type_comparison= recode(Comparison, 
                                      "raw_vs_individual"= "R vs I", 
                                      "dna_vs_raw"= "D vs R", "dna_vs_individual"= "D vs I"))

##ordering (factoring) taxa levels and comparisons
uwunifrac.taxalevels.pairwise.adonis.df2$taxa_order <- factor (uwunifrac.taxalevels.pairwise.adonis.df2$`Taxonomic Level`, 
                                                               levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum"))
uwunifrac.taxalevels.pairwise.adonis.df2$comparisons_order <- factor (uwunifrac.taxalevels.pairwise.adonis.df2$pool_type_comparison, 
                                                                      levels = c("D vs R", "D vs I", "R vs I"))

#PLOTTING R2 VALUES FOR PAIRWISE COMPARISONS#####
gunifrac.taxalevels.pairwise.adonis.df2 
wunifrac.taxalevels.pairwise.adonis.df2 
uwunifrac.taxalevels.pairwise.adonis.df2 

##Putting together dataframes of Unifrac distances
combined_unifrac_pairwise_taxalevels <- rbind(gunifrac.taxalevels.pairwise.adonis.df2,
                                              wunifrac.taxalevels.pairwise.adonis.df2,
                                              uwunifrac.taxalevels.pairwise.adonis.df2)
##Factor ordering 
combined_unifrac_pairwise_taxalevels$UniFrac <- factor (combined_unifrac_pairwise_taxalevels$UniFrac, 
                                                        levels = c("Weighted UniFrac", "Generalized UniFrac","Unweighted UniFrac"))
##Plot
pairwise_unifrac_combined_plot <- 
  ggplot(combined_unifrac_pairwise_taxalevels,
         aes(x=taxa_order, y=comparisons_order, size= R2*100))+
  geom_point(aes(color=taxa_order))+
  scale_color_brewer(palette = "BuPu")+
  geom_text(aes(label=pval_star), vjust = -0.7, show.legend = FALSE, size = 3)+
  facet_wrap(~UniFrac, nrow=3)+
  theme_bw()+
  theme(
    plot.margin = unit(c(0.2,0.2,0,0.2), "cm"),
    panel.spacing = unit(0.2, "lines"),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1),
    panel.grid.major.y = element_line(linetype = 3, color = "black"),
    axis.title.y = element_text(size=12),
    axis.text = element_text(size=10),
    axis.text.x= element_text(angle= 45, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.text = element_text(size=10),
    strip.background = element_blank(),
    strip.text = element_text(size=14, margin = margin(0.1, 0.01, 0.1, 0.01, "cm")),
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(-2, 0, 0, 0),
    legend.position = "bottom")+
  labs(y="Pairwise comparisons", size= expression(R^2 * "(%)"))+
  guides(color = "none")
pairwise_unifrac_combined_plot

##Supplementary TABLE 3 of Pairwise PERMANOVA results ####
pairwise_unifrac_taxalev_adonis_df <- combined_unifrac_pairwise_taxalevels %>%
  mutate(`Taxonomic Level` = factor(`Taxonomic Level`, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")))%>%
  arrange(UniFrac, `Taxonomic Level`) %>%  
  mutate(Comparison = recode(Comparison, "dna_vs_raw"= "DNA Pools vs Raw Pools", 
                             "dna_vs_individual"= "DNA Pools vs Individuals", 
                             "raw_vs_individual"= "Raw Pools vs Individuals"))%>%
  select(`Taxonomic Level`, UniFrac, Comparison, Df, SumOfSqs, R2, `F`, P_value)

stable3 <- pairwise_unifrac_taxalev_adonis_df
write_xlsx(stable3, 
           "SupplementaryTable3.xlsx") ##Table with PERMDISP p values


##PERMANOVA Results #####
unifrac.taxalev.adonis.df <- rbind(wunifrac.adonis.df,
                                   wunifrac.genus.adonis.df,
                                   wunifrac.family.adonis.df,
                                   wunifrac.order.adonis.df,
                                   wunifrac.class.adonis.df,
                                   wunifrac.phylum.adonis.df,
                                   uwunifrac.adonis.df,
                                   uwunifrac.genus.adonis.df,
                                   uwunifrac.family.adonis.df,
                                   uwunifrac.order.adonis.df,
                                   uwunifrac.class.adonis.df,
                                   uwunifrac.phylum.adonis.df,
                                   gunifrac.adonis.df,
                                   gunifrac.genus.adonis.df,
                                   gunifrac.family.adonis.df,
                                   gunifrac.order.adonis.df,
                                   gunifrac.class.adonis.df,
                                   gunifrac.phylum.adonis.df)
unifrac.taxalev.adonis.df ##PERMANOVA results at different taxonomic levels                     

#Factor UniFrac levels
unifrac.taxalev.adonis.df$UniFrac <-factor(unifrac.taxalev.adonis.df$UniFrac, 
                                           levels = c(
                                             "Weighted UniFrac", 
                                             "Generalized UniFrac", 
                                             "Unweighted UniFrac"))

###SUPPLEMENTARY TABLE 2 of PERMANOVA results #########
unifrac.taxalev.adonis.df %>%
  filter(Result=="pool_type")%>% ##only need this pool_type result%>%
  rename(`Taxonomic Level` = "Taxon_level")%>%
  mutate(`Taxonomic Level` = factor(`Taxonomic Level`, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")))%>%
  arrange(UniFrac, `Taxonomic Level`) %>% 
  select(`Taxonomic Level`, "UniFrac", "Result", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")%>% ##selecting the columns 
  mutate(`Result` = str_replace(`Result`, "pool_type", "Sample Type"))%>%
  rename(Effect = Result)%>%
  write_xlsx(., "SupplementaryTable2.xlsx") ##excel file

##Bar plot of R2 values from PERMANOVA
unifrac.taxalev.adonis_plot <-unifrac.taxalev.adonis.df %>%
  filter(Result=="pool_type")%>%
  mutate(signif = case_when(
    `Pr(>F)` < 0.001 ~ "***",
    `Pr(>F)` < 0.01 ~ "**",
    `Pr(>F)` < 0.05 ~ "*",
    TRUE ~ "" ),
    R2= R2*100,
    Taxon_level= factor(Taxon_level, levels = c("ASV",
                                                "Genus",
                                                "Family",
                                                "Order",
                                                "Class", "Phylum"))) %>%
  ggplot(aes(x = Taxon_level, y = R2, fill = Taxon_level)) +
  facet_grid(~ UniFrac,
             labeller = labeller(UniFrac = c(`Generalized UniFrac` = "Generalized", 
                                             `Weighted UniFrac` = "Weighted", 
                                             `Unweighted UniFrac` = "Unweighted"))) +
  geom_col() +
  scale_fill_brewer(palette = "BuPu")+
  geom_text(aes(label = signif, y = R2 + 0.02), vjust = 0, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  labs(y= expression(R^2 * "(%)"), x= "Taxonomic level")+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1),
    legend.position = "none",
    strip.text = element_text(size=28),
    axis.title = element_text(size=28),
    axis.text = element_text(size=25),
    axis.text.x = element_text(angle = 90)
  )
unifrac.taxalev.adonis_plot
###SUPPLEMENTARY FIGURE 4 #####
sfigure4 <- unifrac.taxalev.adonis_plot
ggsave("SupplementaryFigure4.tiff", plot = sfigure4, device = "tiff", width = 12, height = 8, dpi = 300)

#SUPPLEMENTARY PERMDISP RESULTS TABLE########
##Extracting PERMDISP results for all Unifrac distances and taxonomic levels
##Generalized unifrac
pooltype.gunifrac.phylum.permdisp_df <- pooltype.gunifrac.phylum.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Phylum", UniFrac = "Generalized UniFrac")
pooltype.gunifrac.class.permdisp_df <- pooltype.gunifrac.class.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Class", UniFrac = "Generalized UniFrac")
pooltype.gunifrac.order.permdisp_df <- pooltype.gunifrac.order.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Order", UniFrac = "Generalized UniFrac")
pooltype.gunifrac.family.permdisp_df <- pooltype.gunifrac.family.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Family", UniFrac = "Generalized UniFrac")
pooltype.gunifrac.genus.permdisp_df <- pooltype.gunifrac.genus.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Genus", UniFrac = "Generalized UniFrac")
pooltype.gunifrac.asv.permdisp_df <- pooltype.gunifrac.permdisp$tab %>%
  mutate (`Taxonomic Level`= "ASV", UniFrac = "Generalized UniFrac")
##Weighted unifrac
pooltype.wunifrac.phylum.permdisp_df <- pooltype.wunifrac.phylum.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Phylum", UniFrac = "Weighted UniFrac")
pooltype.wunifrac.class.permdisp_df <- pooltype.wunifrac.class.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Class", UniFrac = "Weighted UniFrac")
pooltype.wunifrac.order.permdisp_df <- pooltype.wunifrac.order.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Order", UniFrac = "Weighted UniFrac")
pooltype.wunifrac.family.permdisp_df <- pooltype.wunifrac.family.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Family", UniFrac = "Weighted UniFrac")
pooltype.wunifrac.genus.permdisp_df <- pooltype.wunifrac.genus.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Genus", UniFrac = "Weighted UniFrac")
pooltype.wunifrac.asv.permdisp_df <- pooltype.wunifrac.permdisp$tab %>%
  mutate (`Taxonomic Level`= "ASV", UniFrac = "Weighted UniFrac")

##Unweighted unifrac
#Phylum
pooltype.uwunifrac.phylum.permdisp_df <- pooltype.uwunifrac.phylum.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Phylum", UniFrac = "Unweighted UniFrac")
##Pairwise 
pooltype.uwunifrac.phylum.permdisp.pairwise_df <- data.frame(pooltype.uwunifrac.phylum.permdisp$pairwise$permuted)%>%
  rename("Permuted P-value" = pooltype.uwunifrac.phylum.permdisp.pairwise.permuted)%>%
  rownames_to_column(var = 'Comparison')%>%
  mutate(Comparison = recode(Comparison, "dna-individual" = "DNA Pools vs Individuals",
         "dna-raw" = "DNA Pools vs Raw Pools", 
         "individual-raw" = "Raw Pools vs Individuals"))%>%
  mutate (`Taxonomic Level`= "Phylum", UniFrac = "Unweighted UniFrac")

#Class
pooltype.uwunifrac.class.permdisp_df <- pooltype.uwunifrac.class.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Class", UniFrac = "Unweighted UniFrac")
##Pairwise 
pooltype.uwunifrac.class.permdisp.pairwise_df <- data.frame(pooltype.uwunifrac.class.permdisp$pairwise$permuted)%>%
  rename("Permuted P-value" = pooltype.uwunifrac.class.permdisp.pairwise.permuted)%>%
  rownames_to_column(var = 'Comparison')%>%
  mutate(Comparison = recode(Comparison, "dna-individual" = "DNA Pools vs Individuals",
                             "dna-raw" = "DNA Pools vs Raw Pools", 
                             "individual-raw" = "Raw Pools vs Individuals"))%>%
  mutate (`Taxonomic Level`= "Class", UniFrac = "Unweighted UniFrac")
#Order
pooltype.uwunifrac.order.permdisp_df <- pooltype.uwunifrac.order.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Order", UniFrac = "Unweighted UniFrac")
##Pairwise 
pooltype.uwunifrac.order.permdisp.pairwise_df <- data.frame(pooltype.uwunifrac.order.permdisp$pairwise$permuted)%>%
  rename("Permuted P-value" = pooltype.uwunifrac.order.permdisp.pairwise.permuted)%>%
  rownames_to_column(var = 'Comparison')%>%
  mutate(Comparison = recode(Comparison, "dna-individual" = "DNA Pools vs Individuals",
                             "dna-raw" = "DNA Pools vs Raw Pools", 
                             "individual-raw" = "Raw Pools vs Individuals"))%>%
  mutate (`Taxonomic Level`= "Order", UniFrac = "Unweighted UniFrac")

#Family
pooltype.uwunifrac.family.permdisp_df <- pooltype.uwunifrac.family.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Family", UniFrac = "Unweighted UniFrac")
##Pairwise 
pooltype.uwunifrac.family.permdisp.pairwise_df <- data.frame(pooltype.uwunifrac.family.permdisp$pairwise$permuted)%>%
  rename("Permuted P-value" = pooltype.uwunifrac.family.permdisp.pairwise.permuted)%>%
  rownames_to_column(var = 'Comparison')%>%
  mutate(Comparison = recode(Comparison, "dna-individual" = "DNA Pools vs Individuals",
                             "dna-raw" = "DNA Pools vs Raw Pools", 
                             "individual-raw" = "Raw Pools vs Individuals"))%>%
  mutate (`Taxonomic Level`= "Family", UniFrac = "Unweighted UniFrac")

#Genus
pooltype.uwunifrac.genus.permdisp_df <- pooltype.uwunifrac.genus.permdisp$tab %>%
  mutate (`Taxonomic Level`= "Genus", UniFrac = "Unweighted UniFrac")
##Pairwise 
pooltype.uwunifrac.genus.permdisp.pairwise_df <- data.frame(pooltype.uwunifrac.genus.permdisp$pairwise$permuted)%>%
  rename("Permuted P-value" = pooltype.uwunifrac.genus.permdisp.pairwise.permuted)%>%
  rownames_to_column(var = 'Comparison')%>%
  mutate(Comparison = recode(Comparison, "dna-individual" = "DNA Pools vs Individuals",
                             "dna-raw" = "DNA Pools vs Raw Pools", 
                             "individual-raw" = "Raw Pools vs Individuals"))%>%
  mutate (`Taxonomic Level`= "Genus", UniFrac = "Unweighted UniFrac")

#ASV
pooltype.uwunifrac.asv.permdisp_df <- pooltype.uwunifrac.permdisp$tab %>%
  mutate (`Taxonomic Level`= "ASV", UniFrac = "Unweighted UniFrac")
##Pairwise 
pooltype.uwunifrac.asv.permdisp.pairwise_df <- data.frame(pooltype.uwunifrac.permdisp$pairwise$permuted)%>%
  rename("Permuted P-value" = pooltype.uwunifrac.permdisp.pairwise.permuted)%>%
  rownames_to_column(var = 'Comparison')%>%
  mutate(Comparison = recode(Comparison, "dna-individual" = "DNA Pools vs Individuals",
                             "dna-raw" = "DNA Pools vs Raw Pools", 
                             "individual-raw" = "Raw Pools vs Individuals"))%>%
  mutate (`Taxonomic Level`= "ASV", UniFrac = "Unweighted UniFrac")

#Overall PERMDISP
pooltype.unifrac.taxa.permdisp <-rbind(pooltype.wunifrac.asv.permdisp_df,
                                       pooltype.wunifrac.genus.permdisp_df,
                                       pooltype.wunifrac.family.permdisp_df,
                                       pooltype.wunifrac.order.permdisp_df,
                                       pooltype.wunifrac.class.permdisp_df,
                                       pooltype.wunifrac.phylum.permdisp_df,
                                       pooltype.uwunifrac.asv.permdisp_df,
                                       pooltype.uwunifrac.genus.permdisp_df,
                                       pooltype.uwunifrac.family.permdisp_df,
                                       pooltype.uwunifrac.order.permdisp_df,
                                       pooltype.uwunifrac.class.permdisp_df,
                                       pooltype.uwunifrac.phylum.permdisp_df,
                                       pooltype.gunifrac.asv.permdisp_df,
                                       pooltype.gunifrac.genus.permdisp_df,
                                       pooltype.gunifrac.family.permdisp_df,
                                       pooltype.gunifrac.order.permdisp_df,
                                       pooltype.gunifrac.class.permdisp_df,
                                       pooltype.gunifrac.phylum.permdisp_df
)
pooltype.unifrac.taxa.permdisp ##Now this has all the PERMDISP results 

pooltype.unifrac.taxa.permdisp_df <- pooltype.unifrac.taxa.permdisp %>%
  filter(Df == 2)%>% ##not interested in residuals, filtering them out 
  mutate(Effect = "Sample Type", 
         UniFrac = factor(UniFrac, levels = c("Weighted UniFrac", "Generalized UniFrac", "Unweighted UniFrac")))%>%
  select(`Taxonomic Level`, UniFrac, Effect, everything())
pooltype.unifrac.taxa.permdisp_df

###SUPPLEMENTARY TABLE 4 ######
stable4 <- pooltype.unifrac.taxa.permdisp_df%>%
  mutate(`Taxonomic Level` = factor(`Taxonomic Level`, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")))%>%
  arrange(UniFrac, `Taxonomic Level`) 
write_xlsx(stable4, 
           "SupplementaryTable4.xlsx") ##Table with PERMDISP p values

##Pairwise permdisps 
pooltype.uwunifrac.asv.permdisp.pairwise_df
pooltype.uwunifrac.genus.permdisp.pairwise_df
pooltype.uwunifrac.family.permdisp.pairwise_df
pooltype.uwunifrac.order.permdisp.pairwise_df
pooltype.uwunifrac.class.permdisp.pairwise_df
pooltype.uwunifrac.phylum.permdisp.pairwise_df

pooltype.unifrac.taxa.permdisp.pairwise <-rbind(pooltype.uwunifrac.asv.permdisp.pairwise_df,
                                                pooltype.uwunifrac.genus.permdisp.pairwise_df,
                                                pooltype.uwunifrac.family.permdisp.pairwise_df,
                                                pooltype.uwunifrac.order.permdisp.pairwise_df,
                                                pooltype.uwunifrac.class.permdisp.pairwise_df,
                                                pooltype.uwunifrac.phylum.permdisp.pairwise_df)


pooltype.unifrac.taxa.permdisp.pairwise_df <- pooltype.unifrac.taxa.permdisp.pairwise %>%
  select(`Taxonomic Level`, UniFrac, Comparison, `Permuted P-value`)

###SUPPLEMENTARY TABLE 5 ######
stable5 <- pooltype.unifrac.taxa.permdisp.pairwise_df%>%
  mutate(`Taxonomic Level` = factor(`Taxonomic Level`, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")))%>%
  arrange(`Taxonomic Level`)
write_xlsx(stable5, 
           "SupplementaryTable5.xlsx") ##Table with PERMDISP p values

#BOX PLOTS OF UNIFRAC DISTANCES (ASV level)----
##sample data 
data4_sampledata <- sample_data(data4)
##WUNIFRAC
wunifrac_distance_df <-metagMisc::dist2list(data4.css.wunifrac, tri=T) #This function takes a distance matrix (of class 'dist') and transforms it to a data.frame, where each row represents a single pairwise comparison.
#write.xlsx(wunifrac_distance_df, file = "wunifrac_distance_df.xlsx") ##just checking it out

### adding columns to know the types of comparisons (dna vs. individual, DNA vs. raw, raw vs. individual)
wunifrac_distance_df <- wunifrac_distance_df %>%
  mutate(pool_type_col = data4_sampledata$pool_type[match(col, data4_sampledata$sample_ID)],
         pool_type_row = data4_sampledata$pool_type[match(row, data4_sampledata$sample_ID)])%>% ##pool_type_row and pool_type_col now have the type of sample (pooling) for the comparisons
  mutate (pool_type_col.abbrv = recode(pool_type_col, "individual"= "I", "dna"= "D", "raw"= "R"),
          pool_type_row.abbrv = recode(pool_type_row, "individual"= "I", "dna"= "D", "raw"= "R"))%>% ## abbreviated the sample types
  mutate (pairwise_comparison= paste(pool_type_row.abbrv, pool_type_col.abbrv, sep= " vs "))%>%
  mutate (pairwise_comparison = ifelse(pairwise_comparison == "R vs D", "D vs R", pairwise_comparison),
          unifrac= "Weighted") ## (D vs R) and (R vs D) are the same

wunifrac_distance_df_filt <- wunifrac_distance_df %>%
  dplyr::filter(!pairwise_comparison %in% c("D vs D", "I vs I", "R vs R")) ##filtering out the within-same-group comparisons

##GUNIFRAC
gunifrac_distance_df <-metagMisc::dist2list(data4.css.gunifrac, tri=T) #This function takes a distance matrix (of class 'dist') and transforms it to a data.frame, where each row represents a single pairwise comparison.
#write.xlsx(gunifrac_distance_df, file = "gunifrac_distance_df.xlsx") ##just checking it out

### adding columns to know the types of comparisons (dna vs. individual, DNA vs. raw, raw vs. individual)
gunifrac_distance_df <- gunifrac_distance_df %>%
  mutate(pool_type_col = data4_sampledata$pool_type[match(col, data4_sampledata$sample_ID)],
         pool_type_row = data4_sampledata$pool_type[match(row, data4_sampledata$sample_ID)])%>% ##pool_type_row and pool_type_col now have the type of sample (pooling) for the comparisons
  mutate (pool_type_col.abbrv = recode(pool_type_col, "individual"= "I", "dna"= "D", "raw"= "R"),
          pool_type_row.abbrv = recode(pool_type_row, "individual"= "I", "dna"= "D", "raw"= "R"))%>% ## abbreviated the sample types
  mutate (pairwise_comparison= paste(pool_type_row.abbrv, pool_type_col.abbrv, sep= " vs "))%>%
  mutate (pairwise_comparison = ifelse(pairwise_comparison == "R vs D", "D vs R", pairwise_comparison),
          unifrac= "Generalized") ## (D vs R) and (R vs D) are the same

gunifrac_distance_df_filt <- gunifrac_distance_df %>%
  dplyr::filter(!pairwise_comparison %in% c("D vs D", "I vs I", "R vs R")) ##filtering out the within-same-group comparisons

##UWUNIFRAC
uwunifrac_distance_df <-metagMisc::dist2list(data4.css.uwunifrac, tri=T) #This function takes a distance matrix (of class 'dist') and transforms it to a data.frame, where each row represents a single pairwise comparison.
#write.xlsx(uwunifrac_distance_df, file = "uwunifrac_distance_df.xlsx") ##just checking it out

### adding columns to know the types of comparisons (dna vs. individual, DNA vs. raw, raw vs. individual)
uwunifrac_distance_df <- uwunifrac_distance_df %>%
  mutate(pool_type_col = data4_sampledata$pool_type[match(col, data4_sampledata$sample_ID)],
         pool_type_row = data4_sampledata$pool_type[match(row, data4_sampledata$sample_ID)])%>% ##pool_type_row and pool_type_col now have the type of sample (pooling) for the comparisons
  mutate (pool_type_col.abbrv = recode(pool_type_col, "individual"= "I", "dna"= "D", "raw"= "R"),
          pool_type_row.abbrv = recode(pool_type_row, "individual"= "I", "dna"= "D", "raw"= "R"))%>% ## abbreviated the sample types
  mutate (pairwise_comparison= paste(pool_type_row.abbrv, pool_type_col.abbrv, sep= " vs "))%>%
  mutate (pairwise_comparison = ifelse(pairwise_comparison == "R vs D", "D vs R", pairwise_comparison),
          unifrac= "Unweighted") ## (D vs R) and (R vs D) are the same

##Filtering out the within-same-group comparisons
uwunifrac_distance_df_filt <- uwunifrac_distance_df %>%
  dplyr::filter(!pairwise_comparison %in% c("D vs D", "I vs I", "R vs R")) 

##Binding all the Unifrac dataframes into one 
unifrac_distances_df <- rbind(gunifrac_distance_df_filt,
                              wunifrac_distance_df_filt,
                              uwunifrac_distance_df_filt)


##factor ordering for the plot 
unifrac_distances_df$pwc_ord <- factor (unifrac_distances_df$pairwise_comparison, levels = c("D vs I", "R vs I", "D vs R"))
unifrac_distances_df$unifrac <- factor (unifrac_distances_df$unifrac, levels = c("Weighted", "Generalized", "Unweighted"))

##Box plot of comparisons
unifrac_distances_boxplot <- unifrac_distances_df %>%
  ggplot(aes(x=pwc_ord, y= value, color= pwc_ord, fill= pwc_ord))+
  facet_wrap(~unifrac)+
  labs(y="UniFrac Distance")+
  geom_jitter(alpha=0.08, width = 0.2)+
  geom_boxplot(alpha=0.5)+
  scale_color_manual(values = c("#66C2A5", "#7589A9", "#9E89ED"))+
  scale_fill_manual(values = c("#66C2A5", "#7589A9", "#9E89ED"))+
  scale_y_continuous(expand = c(0.0001, 0.001), limits = c(0,1.25))+
  theme_bw()+
  theme(
    panel.border = element_rect(colour = "black", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.position = "none",
    strip.text = element_text(size = 14),
    strip.background = element_blank())+
  geom_pwc (method = "wilcox_test", 
            p.adjust.method = "BH", 
            label = "p = {p.adj}",
            size =0.5, 
            label.size = 3,
            tip.length = 0.02,
            step.increase = 0.07,
            hide.ns = TRUE)
unifrac_distances_boxplot

#BETA DIVERSITY ORDINATION PLOT FOR OTHER TAX LEVELS ####
##Generalized UniFrac
data4phylum.css.gunifrac.segs.plot <-data4phylum.css.gunifrac.segs %>%
  mutate(Taxlevel = "Phylum")
data4class.css.gunifrac.segs.plot <- data4class.css.gunifrac.segs  %>%
  mutate(Taxlevel = "Class")
data4order.css.gunifrac.segs.plot <- data4order.css.gunifrac.segs %>%
  mutate(Taxlevel = "Order")
data4family.css.gunifrac.segs.plot <- data4family.css.gunifrac.segs %>%
  mutate(Taxlevel = "Family")
data4genus.css.gunifrac.segs.plot <- data4genus.css.gunifrac.segs %>%
  mutate(Taxlevel = "Genus")
#Putting all the taxonomic levels together (generalized unifrac)
data4taxa.css.gunifrac.segs.plot <- rbind(data4phylum.css.gunifrac.segs.plot,
                                          data4class.css.gunifrac.segs.plot,
                                          data4order.css.gunifrac.segs.plot,
                                          data4family.css.gunifrac.segs.plot,
                                          data4genus.css.gunifrac.segs.plot) %>%
  mutate(UniFrac = "Generalized UniFrac")

##Weighted UniFrac
data4phylum.css.wunifrac.segs.plot <-data4phylum.css.wunifrac.segs %>%
  mutate(Taxlevel = "Phylum")
data4class.css.wunifrac.segs.plot <- data4class.css.wunifrac.segs  %>%
  mutate(Taxlevel = "Class")
data4order.css.wunifrac.segs.plot <- data4order.css.wunifrac.segs %>%
  mutate(Taxlevel = "Order")
data4family.css.wunifrac.segs.plot <- data4family.css.wunifrac.segs %>%
  mutate(Taxlevel = "Family")
data4genus.css.wunifrac.segs.plot <- data4genus.css.wunifrac.segs %>%
  mutate(Taxlevel = "Genus")

#Putting all the taxonomic levels together (weighted unifrac)
data4taxa.css.wunifrac.segs.plot <- rbind(data4phylum.css.wunifrac.segs.plot,
                                          data4class.css.wunifrac.segs.plot,
                                          data4order.css.wunifrac.segs.plot,
                                          data4family.css.wunifrac.segs.plot,
                                          data4genus.css.wunifrac.segs.plot) %>%
  mutate(UniFrac = "Weighted UniFrac")

##Unweighted UniFrac 
data4phylum.css.uwunifrac.segs.plot <-data4phylum.css.uwunifrac.segs %>%
  mutate(Taxlevel = "Phylum")
data4class.css.uwunifrac.segs.plot <- data4class.css.uwunifrac.segs  %>%
  mutate(Taxlevel = "Class")
data4order.css.uwunifrac.segs.plot <- data4order.css.uwunifrac.segs %>%
  mutate(Taxlevel = "Order")
data4family.css.uwunifrac.segs.plot <- data4family.css.uwunifrac.segs %>%
  mutate(Taxlevel = "Family")
data4genus.css.uwunifrac.segs.plot <- data4genus.css.uwunifrac.segs %>%
  mutate(Taxlevel = "Genus")

#Putting all the taxonomic levels together (unweighted unifrac)
data4taxa.css.uwunifrac.segs.plot <- bind_rows(data4phylum.css.uwunifrac.segs.plot,
                                           data4class.css.uwunifrac.segs.plot,
                                           data4order.css.uwunifrac.segs.plot,
                                           data4family.css.uwunifrac.segs.plot,
                                           data4genus.css.uwunifrac.segs.plot) %>%
  mutate(UniFrac = "Unweighted UniFrac")

##Putting together all the unifrac distance measures
data4taxa.css.unifrac.segs.plot <- bind_rows(data4taxa.css.gunifrac.segs.plot,
                                         data4taxa.css.wunifrac.segs.plot,
                                         data4taxa.css.uwunifrac.segs.plot)
data4taxa.css.unifrac.segs.plot$Taxlevel <- factor (data4taxa.css.unifrac.segs.plot$Taxlevel, 
                                                    levels = c("Genus", "Family", "Order", "Class", "Phylum"))
data4taxa.css.unifrac.segs.plot$UniFrac <- factor (data4taxa.css.unifrac.segs.plot$UniFrac, 
                                                   levels = c("Weighted UniFrac",
                                                              "Generalized UniFrac", 
                                                              "Unweighted UniFrac"))

##PLOTTING
NMDS_unifrac_taxlevels <-ggplot(data4taxa.css.unifrac.segs.plot) +
  theme_bw() +
  facet_grid(UniFrac~Taxlevel, 
             scales = "free")+
  labs(x= "NMDS1", y = "NMDS2", colour= "Sample Type", fill = "Sample Type") +
  geom_point(aes (x= MDS1, y = MDS2, colour= pool_order), size= 4, shape =18) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= pool_order, colour = pool_order), alpha = 0.32, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes (x= cMDS1, y = cMDS2, colour= pool_order), size= 12, shape =18) +
  geom_text(aes (x= cMDS1, y = cMDS2,label= pool_type.abbrv), colour= "white", size = 5) + ##text adds the abbrv. version of the pool type
  scale_fill_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  scale_color_manual(values = pooling_palette, label = c("Individual", "DNA Pool", "Raw Pool"))+
  theme(legend.text = element_text(size=22),
        legend.title = element_text(size=22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size=25),
        axis.title = element_text(size = 30),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.background = element_blank(),
        strip.text = element_text(size=30),
        legend.position = "bottom")+
  guides(
    colour = guide_legend(override.aes = list(size = 10, shape = 18)),
  )
NMDS_unifrac_taxlevels

##SUPPLEMENTARY FIGURE 3#######
sfigure3 <- NMDS_unifrac_taxlevels
ggsave("SupplementaryFigure3.tiff", plot = sfigure3, device = "tiff", width = 20, height = 15, dpi = 300)

# DIFFERENTIAL ABUNDANCE WIRH ANCOMBC####
### Ancombc uses untransformed counts 
##PHYLUM######
ancombc_counts_phylum <-data4_phylum.counts ##data4_phylum.counts has the untransformed counts glommed at the Phylum level
## Reorder sample variable, individual as "reference"
ancombc_counts_phylum@sam_data$pool_type <- factor(ancombc_counts_phylum@sam_data$pool_type, levels = c("individual", "dna", "raw"))

##running ancombc on the variable of interest (pool type)
ancombc_output_phylum <-ancombc2(data= ancombc_counts_phylum, assay_name = "counts", tax_level = "Phylum",
                                 fix_formula = "pool_type", rand_formula = NULL, 
                                 p_adj_method = "BH", prv_cut = 0.10, lib_cut = 1000, 
                                 group= "pool_type", struc_zero = T, neg_lb = T, 
                                 alpha = 0.05, n_cl = 1, verbose = T, global = T, pairwise = T)

## extract results from pairwise comparisons 
res_pair_phylum <- ancombc_output_phylum$res_pair 

#### Pivot into long form for plotting 
ancom_pool_type_phylum <- res_pair_phylum %>%
  # Pivot LFC, diff, and passed_ss together
  pivot_longer(
    cols = c(
      starts_with("lfc_pool_type"),
      starts_with("diff_pool_type"),
      starts_with("passed_ss_pool_type")
    ),
    names_to = c(".value", "comparison"),
    names_pattern = "(lfc|diff|passed_ss)_(pool_type.*)"
  ) %>%
  
  # Add label based on significance and if it passed_ss
  mutate(
    lfc_rounded = round(lfc, 2),
    significance_ss_label = case_when(
      diff == 1 & passed_ss == 1 ~ "Significant and Passed Sensitivity Test", # significant and passed SS
      diff == 1 & passed_ss != 1 ~ "Significant but Did Not Pass Sensitivity Test", # significant but failed SS
      TRUE ~ "Not Significant"  # not significant
    ),
    #Nicer comparison labels
    group = recode(
      comparison,
      pool_typeraw = "Raw vs I",
      pool_typedna = "DNA vs I",
      pool_typeraw_pool_typedna = "Raw vs DNA"
    )
  ) %>%
  arrange(taxon, comparison)
ancom_pool_type_phylum

##This ancombc was done at the genus level, changing the column name to that
ancom_pool_type_phylum<- ancom_pool_type_phylum %>%
  rename(Phylum = taxon) 

# Calculating the mean abundance across samples of each Phylum 
data4_phylum.counts ##already have this object where I did tax glomming at the Phylum level on counts (data4)
abundance_phylum_df <- data.frame(`Mean Abundance Across Samples` = rowMeans(otu_table(data4_phylum.counts)), 
                                  phyloseq::tax_table(data4_phylum.counts), 
                                  check.names = F) ##calculating mean abundance (across all samples) and adding taxonomy at phylum level
##merging ancombc results at the Phylum level with the mean abundance of each phylum
ancombc_pool_type_phylum_df <- merge (abundance_phylum_df, ancom_pool_type_phylum, by="Phylum") 

##Plotting differential abundance
dot_plot_ancombc_phylum <-ggplot(ancombc_pool_type_phylum_df, aes(x=lfc, y=Phylum, fill=significance_ss_label))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_point(aes(size=`Mean Abundance Across Samples`), 
             alpha = 1, na.rm = T, shape = 21, colour = "black") +
  scale_fill_manual(values=c("Not Significant" = "darkgrey", 
                             "Significant and Passed Sensitivity Test" = "darkred",
                             "Significant but Did Not Pass Sensitivity Test" = "black"))+
  facet_wrap(~group)+
  labs(fill = "Significance", x="LogFC")+
  theme_bw()+
  theme(
    strip.text=element_text(size=32),
    axis.text.x=element_text(size=20),
    axis.text.y = element_text(size=14, face = "bold"),
    axis.title=element_text(size=30),
    strip.background = element_blank(),
    legend.position = "bottom", 
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 16), 
    legend.title.position = "top",
    legend.direction = "vertical") +
  guides(fill = guide_legend (override.aes = list(size = 6)))
dot_plot_ancombc_phylum

##GENUS######
data4_genus.ra <- tax_glom(data4.ra, taxrank = "Genus", NArm=FALSE) 
data4_genus.counts <- tax_glom(data4, taxrank = "Genus", NArm=FALSE) 
#Preprocessing, filtering out low abundance genera (less than 0.2% RA)
data4_genus.ra_filt_0.2 <- filter_taxa(data4_genus.ra, function(x) mean(x) > 0.2, TRUE)
data4_genus.ra_filt_0.2 ##70 genera with mean RA across samples > 0.2%
# Now, subset the counts object to only include the genera >0.2% RA
data4_genus.counts_filtered <- subset_taxa(
  data4_genus.counts,
  Genus %in% tax_table(data4_genus.ra_filt_0.2)[, "Genus"]
)
data4_genus.counts_filtered ##70 genera -OK

### ANCOMBC uses untransformed counts 
ancombc_counts_genus <-data4_genus.counts_filtered ##data4_genus.counts has the untransformed counts (with RA > 0.2%) glommed at the Genus level
ancombc_counts_genus@sam_data$pool_type <- factor(ancombc_counts_genus@sam_data$pool_type, levels = c("individual", "dna", "raw"))## make a new object to reorder sample variable, individual as "reference"
ancombc_counts_genus@sam_data$pool_type

##running ancombc on the variable of interest (pool type)
ancombc_output_genus <-ancombc2(data= ancombc_counts_genus, assay_name = "counts", tax_level = "Genus",
                                fix_formula = "pool_type", rand_formula = NULL, 
                                p_adj_method = "BH", prv_cut = 0.1, lib_cut = 1000, 
                                group= "pool_type", struc_zero = T, neg_lb = T, 
                                alpha = 0.05, n_cl = 1, verbose = T, global = T, pairwise = T)

## extract results from pairwise comparisons 
res_pair_genus <- ancombc_output_genus$res_pair 

#### Pivot into long form for plotting 
ancom_pool_type_genus <- res_pair_genus %>%
  # Pivot LFC, diff, and passed_ss together
  pivot_longer(
    cols = c(
      starts_with("lfc_pool_type"),
      starts_with("diff_pool_type"),
      starts_with("passed_ss_pool_type")
    ),
    names_to = c(".value", "comparison"),
    names_pattern = "(lfc|diff|passed_ss)_(pool_type.*)"
  ) %>%
  
  # Add label based on significance and if it passed_ss
  mutate(
    lfc_rounded = round(lfc, 2),
    significance_ss_label = case_when(
      diff == 1 & passed_ss == 1 ~ "Significant and Passed Sensitivity Test", # significant and passed SS
      diff == 1 & passed_ss != 1 ~ "Significant but Did Not Pass Sensitivity Test", # significant but failed SS
      TRUE ~ "Not Significant"  # not significant
    ),
    #Nicer comparison labels
    group = recode(
      comparison,
      pool_typeraw = "Raw vs I",
      pool_typedna = "DNA vs I",
      pool_typeraw_pool_typedna = "Raw vs DNA"
    )
  ) %>%
  arrange(taxon, comparison)
ancom_pool_type_genus

##This ancombc was done at the genus level, changing the column name to that
ancom_pool_type_genus<- ancom_pool_type_genus %>%
  rename(Genus = taxon) 

# Calculating the mean abundance across samples of each genus 
data4_genus.counts_filtered ##already have this object where I did tax glomming at the genus level on counts (data4), and filtered out low abundance 
abundance_genus_df <- data.frame(`Mean Abundance Across Samples`= rowMeans(otu_table(data4_genus.counts_filtered)), 
                                 phyloseq::tax_table(data4_genus.counts_filtered), 
                                 check.names = F) ##calculating mean abundance and adding taxonomy at genus level
##Merging ancombc results at the genus level with the mean abundance of each genus
ancombc_pool_type_genus_df <- merge (abundance_genus_df, ancom_pool_type_genus, by="Genus") 

ancombc_pool_type_genus_df%>%
  filter(Phylum == "unclassified Unassigned")

##Plotting differential abundance
dot_plot_ancombc_genus <-ggplot(ancombc_pool_type_genus_df, aes(x=lfc, y=Genus, fill=significance_ss_label))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_point(aes(size=`Mean Abundance Across Samples`), 
             alpha = 1, na.rm = T, shape = 21, colour = "black") +
  scale_fill_manual(values=c("Not Significant" = "darkgrey", 
                             "Significant and Passed Sensitivity Test" = "darkred",
                             "Significant but Did Not Pass Sensitivity Test" = "black"))+
  facet_wrap(~group)+
  labs(fill = "Significance", x="LogFC")+
  theme_bw()+
  theme(
    strip.text=element_text(size=32),
    axis.text.x=element_text(size=20),
    axis.text.y = element_text(size=12, face = "bold"),
    axis.title=element_text(size=30),
    strip.background = element_blank(),
    legend.position = "bottom", 
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 16), 
    legend.title.position = "top",
    legend.direction = "vertical") +
  guides(fill = guide_legend (override.aes = list(size = 6)))
dot_plot_ancombc_genus

##Genus and Phylum together 
ANCOMBC_phylum_genus <- ggarrange(dot_plot_ancombc_phylum, 
                                  dot_plot_ancombc_genus,
                                  ncol = 1,
                                  common.legend = T,
                                  legend = "bottom",
                                  align = "v",
                                  labels = "AUTO",
                                  font.label = list(size = 25))
ANCOMBC_phylum_genus
ggsave ("ANCOMBC_phylum_genus.tiff", plot = ANCOMBC_phylum_genus, 
        bg = "white",
        device = "tiff", width = 12, height = 24, dpi = 300)

##FIGURE 4 ######
figure4 <- ggarrange(dot_plot_ancombc_phylum, 
                     dot_plot_ancombc_genus,
                     ncol = 1,
                     common.legend = T,
                     legend = "bottom",
                     align = "v",
                     labels = "AUTO",
                     font.label = list(size = 25))
figure4 
ggsave ("Figure4pubr.tiff", plot = figure4, 
        bg = 'white', device = "tiff", width = 12, height = 24, dpi = 300)



#FIGURE 2#######
figure2 <- plot_grid(
  NMDS_unifrac_ASV,
  plot_grid(pairwise_unifrac_combined_plot, unifrac_distances_boxplot, ncol = 1, labels = c("B", "C")),
  ncol = 2, 
  labels = "A"
)
figure2
ggsave("Figure2pubr.tiff", plot= figure2, device = "tiff", width= 8.5, height = 8.5, dpi = 300)

#FIGURE 3#######
gg_upset_plot <- as.ggplot(upset_plot)

##Saving Upset as TIFF
upset_plot_tiff <- readTIFF("Upset_plot.tiff")
# Then converting image to a grob
upset_plot_grob <- rasterGrob(upset_plot_tiff, interpolate = TRUE, x = 0, y = 0, hjust = 0, vjust = 0)
upset_plot_grob

figure3 <-plot_grid(upset_plot_grob, prev_meanRA_plot, labels = "AUTO", label_size = 20)
figure3 
ggsave("Figure3pubr.tiff", plot= figure3, device = "tiff", width= 14, height = 6, dpi = 300)


