#set working directory to biom files
#setwd("/mnt/md0/Benchtop_BSF_Seqs/kraken_mapped") 

library(phyloseq)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(stats)
library(rstatix)
library(ggpubr)

############################
#### Import data
############################

#import biom file from kraken-biom
data <-import_biom("./biom_files_targetcapture/bracken_genera.ITS_targetcapture.biom", parseFunction = parse_taxonomy_default)
data@tax_table@.Data <- substring(data@tax_table@.Data, 4)
#rename columns to species level
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#View(data@tax_table@.Data)

#I added a taxa to this that is likely a false positive from Kraken
data <- subset_taxa(data, (!Kingdom %in% "Viridiplantae")) # Plants
data <- subset_taxa(data, (!Kingdom %in% "Metazoa")) # Animals
data <- subset_taxa(data, (Kingdom %in% "Fungi")) # Keep fungi

#change names of samples from files
# this step is to fix column names (e.g. "CLBJ-22-1-Rh.k2_bracken_genuses" ) View(data@otu_table). 
colnames(data@otu_table@.Data) <- gsub("\\.k2_bracken_genuses", "", colnames(data@otu_table@.Data))

#change this based on what taxonomic level you imported from Bracken. e.g., if you imported genus level data, 
##change "_._report_bracken_families" to "_._report_bracken_genuses" 

#load in sample metadata
md <- read.csv(file = "./metadata/amplicon_resequencing_metadata.tsv", header = T, sep = "\t")

#create sample data for phyloseq object
samdf <- md

#make sample names the rows
row.names(samdf) <- samdf[, 1]
samdf <- samdf[, 2:ncol(samdf)]
## samdf$Bacterial.Supplement <- factor(samdf$Bacterial.Supplement, levels = c("Yes", "No"))

# adds metadata to data
# --- ALIGN METADATA TO OTU TABLE SAMPLE NAMES (critical) ---
otu_samps <- sample_names(data)  # or colnames(data@otu_table@.Data)
md_samps  <- rownames(samdf)
keep <- intersect(otu_samps, md_samps)
keep

samdf <- samdf[keep, , drop = FALSE]
samdf <- samdf[sample_names(data), , drop = FALSE] # reorder metadata to match phyloseq sample order
stopifnot(identical(rownames(samdf), sample_names(data)))

# assign sample_data
sample_data(data) <- sample_data(samdf)
data  <- prune_samples(keep, data)

# Separate by sample type
unique(data@sam_data$Type)
nodule <- subset_samples(data, Type == "nodule")
root <- subset_samples(data, Type == "root")
rhizosphere <- subset_samples(data, Type == "rhizosphere")
mock_community <- subset_samples(data, Type == "Mock community")

################################
#### Mock community taxon barplots
################################

ideal <- read.csv("mock_community_ITS.csv")

#top20 <- names(sort(taxa_sums(mock_community), decreasing = TRUE))[1:12]
#ps.top20 <- transform_sample_counts(mock_community, function(OTU) OTU/sum(OTU))
#ps.top20 <- prune_taxa(top20, ps.top20)

####
# Create plotting object
####

alltaxa <- names(sort(taxa_sums(mock_community), decreasing = TRUE))
ps.alltaxa <- transform_sample_counts(mock_community, function(OTU) OTU/sum(OTU))
ps.alltaxa <- prune_taxa(alltaxa, ps.alltaxa)

for_gg <- psmelt(ps.alltaxa)
for_gg_ideal <- ideal

for_gg %>% mutate(Genus_lumped = ifelse((! Genus %in% ideal$Genus), "Other", as.character(Genus)), Genus_lumped = factor(Genus_lumped)) -> for_gg
# Code to force all factors to show in the following ggplot even with 0 data to avoid messing up legend ordering
all_levels <- union(levels(factor(for_gg$Genus_lumped)), levels(factor(for_gg_ideal$Genus_lumped)))
for_gg_ideal %>% mutate(Genus_lumped = ifelse((! Genus %in% ideal$Genus), "Other", as.character(Genus)), Genus_lumped = factor(Genus_lumped)) -> for_gg_ideal
for_gg_ideal$Genus_lumped <- factor(for_gg_ideal$Genus_lumped, levels = all_levels)

for_gg %>% mutate(Family_lumped = ifelse((! Family %in% ideal$Family), "Other", as.character(Family)), Family_lumped = factor(Family_lumped)) -> for_gg
# Code to force all factors to show in the following ggplot even with 0 data to avoid messing up legend ordering
all_levels <- union(levels(factor(for_gg$Family_lumped)), levels(factor(for_gg_ideal$Family_lumped)))
for_gg$Family_lumped <- factor(for_gg$Family_lumped, levels = all_levels)
# Check for missing genera and add them to the ideal dataset to avoid messing up legend ordering
for_gg_ideal %>% mutate(Family_lumped = ifelse((! Family %in% ideal$Family), "Other", as.character(Family)), Family_lumped = factor(Family_lumped)) -> for_gg_ideal
for_gg_ideal$Family_lumped <- factor(for_gg_ideal$Family_lumped, levels = all_levels)

for_gg %>% mutate(Phylum_lumped = ifelse((! Phylum %in% ideal$Phylum), "Other", as.character(Phylum)), Phylum_lumped = factor(Phylum_lumped)) -> for_gg
for_gg_ideal %>% mutate(Phylum_lumped = ifelse((! Phylum %in% ideal$Phylum), "Other", as.character(Phylum)), Phylum_lumped = factor(Phylum_lumped)) -> for_gg_ideal
# Code to force all factors to show in the following ggplot even with 0 data to avoid messing up legend ordering
all_levels <- union(levels(factor(for_gg$Phylum_lumped)), levels(factor(for_gg_ideal$Phylum_lumped)))
for_gg$Phylum_lumped <- factor(for_gg$Phylum_lumped, levels = all_levels)
for_gg_ideal$Phylum_lumped <- factor(for_gg_ideal$Phylum_lumped, levels = all_levels)

#######
## Create empty entries for non-shared taxa to enforce same colors across plots
## This ended up being frail code; consider instead automatically generating shared manual colors following the resequencing code lower down in this script
#######
#
#for_gg_ideal <- bind_rows(for_gg_ideal, data.frame(Sample.ID.1 = "Mock community", Genus_lumped = "Other", Family_lumped = "Other", Phylum_lumped = "Other"))
#
#library(dplyr)
#missing_genera <- as.data.frame(setdiff(for_gg$Genus_lumped, for_gg_ideal$Genus_lumped))
#names(missing_genera)[1] <- "Genus_lumped"
#if(length(missing_genera$Genus_lumped) > 0){
#	missing_genera$Abundance <- 0
#	missing_genera$Sample.ID.1 <- "Mock community"
#	for_gg_ideal <- bind_rows(for_gg_ideal, missing_genera)
#	}
#
#missing_families <- as.data.frame(setdiff(for_gg$Family_lumped, for_gg_ideal$Family_lumped))
#names(missing_families)[1] <- "Family_lumped"
#if(length(missing_families$Family_lumped) > 0){
#	missing_families$Abundance <- 0
#	missing_families$Sample.ID.1 <- "Mock community"
#	for_gg_ideal <- bind_rows(for_gg_ideal, missing_families)
#	}
#
#missing_phyla <- as.data.frame(setdiff(for_gg$Phylum_lumped, for_gg_ideal$Phylum_lumped))
#names(missing_phyla)[1] <- "Phylum_lumped"
#if(length(missing_phyla$Phylum_lumped) > 0){
#	missing_phyla$Abundance <- 0
#	missing_phyla$Sample.ID.1 <- "Mock community"
#	for_gg_ideal <- bind_rows(for_gg_ideal, missing_phyla)
#	}


### Genus

#make pct level bar chart
pct_mockcommunity_genus <- ggplot(data = for_gg, 
    mapping = aes(Sample.ID.1, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12))+
	geom_bar(aes(color = Genus_lumped, fill = Genus_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	#scale_fill_manual(values = col_tax) + 
	#scale_color_manual(values = col_tax)+
	ggtitle("Rel. Abundance:\nMock community")

pct_mockcommunity_genus_ideal <- ggplot(data = for_gg_ideal, 
    mapping = aes(Sample.ID.1, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12)) + 
	geom_bar(aes(color = Genus_lumped, fill = Genus_lumped), stat = "identity", position = "fill") +
	scale_fill_discrete(drop = FALSE) + scale_color_discrete(drop = FALSE) +
	#facet_wrap(~Genus, scales = "free_x") + 
	#scale_fill_manual(values = col_tax) + 
	#scale_color_manual(values = col_tax)+
	ggtitle("Expected Rel. Abundance:\nMock community")
 
library(patchwork)
 
# Use guides = collect to enforce common legend and fill colors
mockcommunity_fig <- pct_mockcommunity_genus + pct_mockcommunity_genus_ideal + patchwork::plot_layout(nrow = 1, ncol = 2, guides = "collect") + patchwork::plot_annotation(tag_levels = "A")

mockcommunity_fig

### Family


#make pct level bar chart
pct_mockcommunity_family <- ggplot(data = for_gg, 
     mapping = aes(Sample.ID.1, Abundance)) + 
 theme(legend.text = element_text(face = "italic"), text = element_text(size = 12)) +
 geom_bar(aes(color = Family_lumped, fill = Family_lumped), stat = "identity", position = "fill") +
 #facet_wrap(~Genus, scales = "free_x") + 
 #scale_fill_manual(values = col_tax) + 
 #scale_color_manual(values = col_tax)+
 ggtitle("Rel. Abundance:\nMock community")

pct_mockcommunity_family_ideal <- ggplot(data = for_gg_ideal, 
     mapping = aes(Sample.ID.1, Abundance)) + 
 theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), legend.position = "none") +
 geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "fill") +
 #facet_wrap(~Genus, scales = "free_x") + 
 #scale_fill_manual(values = col_tax) + 
 #scale_color_manual(values = col_tax)+
 ggtitle("Expected Rel. Abundance:\nMock community")
 
mockcommunity_fig <- (pct_mockcommunity_family + pct_mockcommunity_family_ideal) + patchwork::plot_layout(nrow = 1, ncol = 2, guides = "collect") + patchwork::plot_annotation(tag_levels = "A")

mockcommunity_fig

### Phylum

#make pct level bar chart
pct_mockcommunity_phylum <- ggplot(data = for_gg, 
     mapping = aes(Sample.ID.1, Abundance)) + 
 theme(legend.text = element_text(face = "italic"), text = element_text(size = 12)) +
 geom_bar(aes(color = Phylum_lumped, fill = Phylum_lumped), stat = "identity", position = "fill") +
 #facet_wrap(~Genus, scales = "free_x") + 
 #scale_fill_manual(values = col_tax) + 
 #scale_color_manual(values = col_tax)+
 ggtitle("Rel. Abundance:\nMock community")

pct_mockcommunity_phylum_ideal <- ggplot(data = for_gg_ideal, 
     mapping = aes(Sample.ID.1, Abundance)) + 
 theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), legend.position = "none")+
 geom_bar(aes(color = Phylum_lumped, fill = Phylum_lumped), stat = "identity", position = "fill") +
 #facet_wrap(~Genus, scales = "free_x") + 
 #scale_fill_manual(values = col_tax) + 
 #scale_color_manual(values = col_tax)+
 ggtitle("Expected Rel. Abundance:\nMock community")

mockcommunity_fig <- (pct_mockcommunity_phylum + pct_mockcommunity_phylum_ideal) + patchwork::plot_layout(nrow = 1, ncol = 2, guides = "collect") + patchwork::plot_annotation(tag_levels = "A")

mockcommunity_fig




# No amplicon mock community



########################################
#### Import amplicon resequencing and original amplicon data
########################################

#### QIIME2 CODE TO GET COMPARABLE .BIOM FILE
#conda activate qiime2-2022.8 
#qiime tools export --input-path table.hostremoved-l6.qza --output-path collapsed-table-L6
#biom convert -i collapsed-table-L6/feature-table.biom -o table.hostremoved-l6.tsv --to-tsv --header-key taxonomy
#mv ./collapsed-table-L6/feature-table.biom ./table.hostremoved-l6.biom && rmdir ./collapsed-table-L6/
#
#qiime taxa collapse --i-table table.hostremoved.qza --i-taxonomy taxonomy.qza --p-level 5 --o-collapsed-table table.hostremoved-l5.qza
#qiime tools export --input-path table.hostremoved-l5.qza --output-path collapsed-table-L5
#biom convert -i collapsed-table-L5/feature-table.biom -o table.hostremoved-l5.tsv --to-tsv --header-key taxonomy
#mv ./collapsed-table-L5/feature-table.biom ./table.hostremoved-l5.biom && rmdir ./collapsed-table-L5/
#
#qiime tools export --input-path table.hostremoved-l2.qza --output-path collapsed-table-L2
#biom convert -i collapsed-table-L2/feature-table.biom -o table.hostremoved-l2.tsv --to-tsv --header-key taxonomy
#mv ./collapsed-table-L2/feature-table.biom ./table.hostremoved-l2.biom && rmdir ./collapsed-table-L2/

# Target sequencing
taxon_list = c("N13-2-D-No", "N13-2-D-Rh", "N13-2-D-Ro", "N34-2-D-No", "N34-2-D-Rh", "N34-2-D-Ro", "N35-4-D-No", "N35-4-D-Rh", "N35-4-D-Ro", "N36-4-D-No", "N36-4-D-Rh", "N36-4-D-Ro", "N38-2-D-No", "N38-2-D-Rh", "N38-2-D-Ro", "N39-2-D-No", "N39-2-D-Rh", "N39-2-D-Ro", "N49-3-D-No", "N49-3-D-Rh", "N49-3-D-Ro")#, "ZymoMockDNA", "ZymoMockDNA_1685")
data  <- prune_samples(taxon_list, data)
sample_names(data)

# Separate by sample type
unique(data@sam_data$Type)
nodule <- subset_samples(data, Type == "nodule")
root <- subset_samples(data, Type == "root")
rhizosphere <- subset_samples(data, Type == "rhizosphere")

amplicon <- read.table("./biom_files_amplicon/ITS/table.hostremoved-l6.tsv", sep = "\t", header = TRUE, row.names = 1, check.names=FALSE)
amplicon <- subset(amplicon, select = names(amplicon) %in% taxon_list)
OTU = otu_table(amplicon, taxa_are_rows = TRUE)

tax <- as.data.frame(t(as.data.frame(strsplit(row.names(OTU), ";"))))
# Define the taxonomic levels
tax <- lapply(tax, function(x) gsub("^[kpcofg]__", "", x, perl = TRUE))
tax <- as.data.frame(lapply(tax, function(x) gsub("__", "", x, perl = TRUE)))
tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
names(tax) <- tax_levels

TAX = tax_table(as.matrix(tax))

row.names(OTU) <- row.names(TAX)
data_amplicon = phyloseq(OTU, TAX)

# load in sample metadata
md <- read.csv(file = "./metadata/amplicon_original_metadata.tsv", header = T, sep = "\t")

# create sample data for phyloseq object
samdf <- md
row.names(samdf) <- samdf[, 1]
samdf <- samdf[, 2:ncol(samdf)]

# adds metadata to data
data_amplicon@sam_data<-sample_data(samdf)

#### Set up plot objects

# This time don't do top 20 microbial taxa, it will warp things in diverse samples
alltaxa <- names(sort(taxa_sums(data), decreasing = TRUE))
ps.alltaxa <- transform_sample_counts(data, function(OTU) OTU/sum(OTU))
ps.alltaxa <- prune_taxa(alltaxa, ps.alltaxa)

for_gg <- psmelt(ps.alltaxa)

# Merge low-abundance taxa into an "other" category
for_gg %>% mutate(Genus_lumped = ifelse(Abundance / sum(for_gg$Abundance) < 0.005, "Other", as.character(Genus)), Genus_lumped = factor(Genus_lumped)) -> for_gg
for_gg %>% mutate(Family_lumped = ifelse(Abundance / sum(for_gg$Abundance) < 0.005, "Other", as.character(Family)), Family_lumped = factor(Family_lumped)) -> for_gg
#for_gg %>% mutate(Phylum_lumped = ifelse(Abundance / sum(for_gg$Abundance) < 0.005, "Other", as.character(Phylum)), Phylum_lumped = factor(Phylum_lumped)) -> for_gg
targets <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Glomeromycota", "Rozellomycota", "Mucoromycota", "Mortierellomycota", "Neocallimastigomycota", "Olpidiomycota", "Kickxellomycota", "Zoopagomycota", "Blastocladiomycota")
for_gg$Phylum_lumped <- for_gg$Phylum
for_gg$Phylum_lumped <- replace(for_gg$Phylum_lumped, for_gg$Phylum_lumped == "GS01_phy_Incertae_sedis", "Other")
for_gg$Phylum_lumped <- replace(for_gg$Phylum_lumped, for_gg$Phylum_lumped == "Fungi_phy_Incertae_sedis", "Other")


# Reorder factor for plot order
for_gg$Type <- factor(for_gg$Type, levels = c("rhizosphere", "root", "nodule"))

alltaxa <- names(sort(taxa_sums(data_amplicon), decreasing = TRUE))
ps.alltaxa <- transform_sample_counts(data_amplicon, function(OTU) OTU/sum(OTU))
ps.alltaxa <- prune_taxa(alltaxa, ps.alltaxa)

# Merge low-abundance taxa into an "other" category -- different strategy than 16S because this behaved weirdly
for_gg_amplicon <- psmelt(ps.alltaxa)
for_gg_amplicon %>% mutate(Genus_lumped = ifelse(Abundance / sum(for_gg_amplicon$Abundance) < 0.005, "Other", as.character(Genus))) -> for_gg_amplicon
for_gg_amplicon$Genus_lumped <- replace(for_gg_amplicon$Genus_lumped, for_gg_amplicon$Genus_lumped == "", "Unassigned")
for_gg_amplicon %>% mutate(Family_lumped = ifelse(Abundance / sum(for_gg_amplicon$Abundance) < 0.005, "Other", as.character(Family))) -> for_gg_amplicon
for_gg_amplicon$Family_lumped <- replace(for_gg_amplicon$Family_lumped, for_gg_amplicon$Family_lumped == "", "Unassigned")
#for_gg_amplicon %>% mutate(Phylum_lumped = ifelse(Abundance / sum(for_gg_amplicon$Abundance) < 0.005, "Other", as.character(Phylum))) -> for_gg_amplicon
targets <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Glomeromycota", "Rozellomycota", "Mucoromycota", "Mortierellomycota", "Neocallimastigomycota", "Olpidiomycota", "Kickxellomycota", "Zoopagomycota", "Blastocladiomycota")
for_gg_amplicon$Phylum_lumped <- for_gg_amplicon$Phylum
for_gg_amplicon$Phylum_lumped <- replace(for_gg_amplicon$Phylum_lumped, for_gg_amplicon$Phylum_lumped == "unidentified", "Other")
for_gg_amplicon$Phylum_lumped <- replace(for_gg_amplicon$Phylum_lumped, for_gg_amplicon$Phylum_lumped == "", "Other")



## Remove leftover host designations
#for_gg_amplicon <- for_gg_amplicon[for_gg_amplicon$Genus != "Lupinus", ]
#for_gg_amplicon <- for_gg_amplicon[for_gg_amplicon$Family != "mitochondria", ]

## Name correction
#for_gg_amplicon$Family[for_gg_amplicon$Genus == "Sinorhizobium"] <- "Ensifer"
#for_gg_amplicon$Family[for_gg_amplicon$Family == "Bradyrhizobiaceae"] <- "Xanthobacteraceae"
#for_gg_amplicon$Family[for_gg_amplicon$Family == "Methylobacteriaceae"] <- "Beijerinckiaceae"
#for_gg_amplicon$Family[for_gg_amplicon$Family == "Rhizobium"] <- "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#for_gg_amplicon$Family[for_gg_amplicon$Family == "Neorhizobium"] <- "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#for_gg_amplicon$Family[for_gg_amplicon$Family == "Pararhizobium"] <- "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#for_gg_amplicon$Family[for_gg_amplicon$Family == "Allorhizobium"] <- "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"

# Reorder factor for plot order
for_gg_amplicon$Type <- factor(for_gg_amplicon$Type, levels = c("rhizosphere", "root", "nodule"))


### Set manual coloring (sadly needed due to ggplot shortcomings with coloring across figures for frustrating reasons)
all_genera <- unique(c(as.character(for_gg$Genus_lumped), as.character(for_gg_amplicon$Genus_lumped)))
# <- c(all_genera, "Other", "Unassigned")
library(RColorBrewer)
# Rep is to do color cycling, rev is to reverse the color palette (I didn't like it forward)
col_tax_genus <- rev(rep(RColorBrewer::brewer.pal(n = length(all_genera), name = "Set3"), length.out = length(all_genera)))
names(col_tax_genus) <- all_genera

### Set manual coloring (sadly needed due to ggplot shortcomings with coloring across figures for frustrating reasons)
all_families <- unique(c(as.character(for_gg$Family_lumped), as.character(for_gg_amplicon$Family_lumped)))
#all_families <- c(all_families, "Other", "Unassigned")
library(RColorBrewer)
# Rep is to do color cycling, rev is to reverse the color palette (I didn't like it forward)
col_tax_family <- rev(rep(RColorBrewer::brewer.pal(n = length(all_families), name = "Set3"), length.out = length(all_families)))
names(col_tax_family) <- all_families

### Set manual coloring (sadly needed due to ggplot shortcomings with coloring across figures for frustrating reasons)
all_phyla <- unique(c(as.character(for_gg$Phylum_lumped), as.character(for_gg_amplicon$Phylum_lumped)))
#all_phyla <- c(all_phyla, "Other", "Unassigned")
library(RColorBrewer)
# Rep is to do color cycling, rev is to reverse the color palette (I didn't like it forward)
col_tax_phylum <- rev(rep(RColorBrewer::brewer.pal(n = length(all_phyla), name = "Set3"), length.out = length(all_phyla)))
names(col_tax_phylum) <- all_phyla



##### 
#Genus amplicon resequencing
#####

#make pct level bar chart
ggplot(data = for_gg, 
    mapping = aes(Tribe, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Genus_lumped, fill = Genus_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_genus) + 
	scale_color_manual(values = col_tax_genus) +
	ggtitle("Rel. Abundance:\nTarget capture")
  
ggplot(data = for_gg, 
    mapping = aes(Type, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Genus_lumped, fill = Genus_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_genus) + 
	scale_color_manual(values = col_tax_genus) +
	ggtitle("Rel. Abundance:\nTarget capture")
	
	
##### 
#Family amplicon resequencing
#####

#make pct level bar chart
ggplot(data = for_gg, 
    mapping = aes(Tribe, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Family_lumped, fill = Family_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_family) + 
	scale_color_manual(values = col_tax_family)+
	ggtitle("Rel. Abundance:\nTarget capture")
  
ggplot(data = for_gg, 
    mapping = aes(Type, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Family_lumped, fill = Family_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_family) + 
	scale_color_manual(values = col_tax_family)+
	ggtitle("Rel. Abundance:\nTarget capture")
	
##### 
# Phylum amplicon resequencing
#####

ggplot(data = for_gg, 
    mapping = aes(Tribe, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Phylum_lumped, fill = Phylum_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_phylum) + 
	scale_color_manual(values = col_tax_phylum) +
	ggtitle("Rel. Abundance:\nTarget capture")

ggplot(data = for_gg, 
    mapping = aes(Type, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Phylum_lumped, fill = Phylum_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_phylum) + 
	scale_color_manual(values = col_tax_phylum)+
	ggtitle("Rel. Abundance:\nTarget capture")


##### 
# Genus amplicon original figs
#####

#make pct level bar chart
ggplot(data = for_gg_amplicon[for_gg_amplicon$Tribe != "na",], 
	mapping = aes(Tribe, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Genus_lumped, fill = Genus_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_genus) + 
	scale_color_manual(values = col_tax_genus)+
	ggtitle("Rel. Abundance:\nAmplicon sequencing")
 
ggplot(data = for_gg_amplicon, 
    mapping = aes(Sample.ID.1, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Genus_lumped, fill = Genus_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_genus) + 
	scale_color_manual(values = col_tax_genus)+
	ggtitle("Rel. Abundance:\nAmplicon sequencing")
 
##### 
# Family amplicon original
#####


#make pct level bar chart
ggplot(data = for_gg_amplicon[for_gg_amplicon$Tribe != "na",], 
	mapping = aes(Tribe, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Family_lumped, fill = Family_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_family) + 
	scale_color_manual(values = col_tax_family)+
	ggtitle("Rel. Abundance:\nAmplicon sequencing")
 
ggplot(data = for_gg_amplicon, 
    mapping = aes(Sample.ID.1, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Family_lumped, fill = Family_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_family) + 
	scale_color_manual(values = col_tax_family)+
	ggtitle("Rel. Abundance:\nAmplicon sequencing")


##### 
# Phylum amplicon original
#####


#make pct level bar chart
ggplot(data = for_gg_amplicon[for_gg_amplicon$Tribe != "na",], 
	mapping = aes(Tribe, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Phylum_lumped, fill = Phylum_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_phylum) + 
	scale_color_manual(values = col_tax_phylum)+
	ggtitle("Rel. Abundance:\nAmplicon sequencing")
 
ggplot(data = for_gg_amplicon, 
    mapping = aes(Sample.ID.1, Abundance)) + 
	theme(legend.text = element_text(face = "italic"), text = element_text(size = 12), axis.text.x = element_text(angle = 45))+
	geom_bar(aes(color = Phylum_lumped, fill = Phylum_lumped), stat = "identity", position = "fill") +
	#facet_wrap(~Genus, scales = "free_x") + 
	scale_fill_manual(values = col_tax_phylum) + 
	scale_color_manual(values = col_tax_phylum)+
	ggtitle("Rel. Abundance:\nAmplicon sequencing")


