### Code for processing 18S and 16S amplicon data with QIIME 2 input files (.qza files)
## Updated 5/27/2022
# Sean Anderson

# Load packages
library(qiime2R)
library(phyloseq)
library(Matrix)
library(vegan)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggExtra)
library(fantaxtic)
library(RColorBrewer)
library(SpiecEasi)
library(factoextra)
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
library(corrplot)

### Upload 18S count, taxonomy, and rooted tree files (.qza files)
table <- read_qza("18S-table.qza")
count_tab <- table$data %>% as.data.frame() # Convert to data frame 
taxonomy <- read_qza("18S-tax.qza")
tax_tab <- taxonomy$data %>% # Convert to data frame, tab separate and rename taxa levels, and remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom","Supergroup","Division","Class","Order","Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)
tree_18S <- read_qza("18S-rooted-tree.qza") # Upload tree and reformat
tree <- tree_18S$data
new_tre <- ape::multi2di(tree)

# Upload 18S metadata
sample_info_tab <- read.table("Sampleinfo_18S.txt", header=TRUE, row.names=1, check.names=F, sep="\t")

# Merge into phyloseq object
ps <- phyloseq(tax_table(as.matrix(tax_tab)), otu_table(count_tab, taxa_are_rows = T), sample_data(sample_info_tab),phy_tree(new_tre))

# Export 18S ASV information - Table S1
OTU = as(otu_table(ps), "matrix")
TAX = as(tax_table(ps), "matrix")
merge_18S <- cbind(OTU,TAX)
write.csv(merge_18S, file="Table_S1.csv", row.names=T)

# Rename ASVs in sequential order
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Remove 18S taxa not assigned at supergroup level or assigned to Metazoa or Streptophyta
ps <- subset_taxa(ps, Supergroup!="Unassigned", Prune = T)
ps <- subset_taxa(ps, Division!="Metazoa", Prune = T)
ps <- subset_taxa(ps, Division!="Streptophyta", Prune = T)

# Adds "Unassigned" label names for missing taxa
ps = name_taxa(ps, label = "Unassigned")

# Estimate minimum, mean, and maximum read counts for 18S 
ps_min <- min(sample_sums(ps))
ps_mean <- mean(sample_sums(ps))
ps_max <- max(sample_sums(ps))

# Rarefaction curves for 18S
rare_18S <- ggrare(ps, step = 100, plot = TRUE, parallel = FALSE, se = FALSE)
rare_18S + theme(legend.position = "none") + theme_bw()+ theme(legend.position = "none")
ggsave(filename = "18S_rare.eps", plot = last_plot(), device = eps, path = NULL, scale = 1, width = 4, height = 4, dpi = 600)

# Remove singletons
ps_filt = filter_taxa(ps, function (x) {sum(x) > 1}, prune=TRUE)

# Rarefy to even sampling depth
ps_rare <- rarefy_even_depth(ps_filt, sample.size = min(sample_sums(ps_filt)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Find out the most relatively abundant 18S groups at order level (>2% abundance)
order_18S <- tax_glom(ps_rare, taxrank = "Order")
order_18S <- transform_sample_counts(order_18S, function(x)100* x / sum(x)) # Transform to relative abundance
OTU <- otu_table(order_18S)
TAX <- tax_table(order_18S)[,"Order"]
Average <- as.data.frame(rowMeans(OTU))# Estimate average abundance at the order level
names(Average) <- c("Mean") # Rename to mean
Table <- merge(TAX, Average, by=0, all=TRUE)
Table$Row.names = NULL # View abundances at order level - we are interested in the top 10 most abundant groups

# Temporal abundance plots for the major 18S groups
abund_18S <- psmelt(order_18S) # Melt the data
abund_18S$Order <- as.character(abund_18S$Order) # Convert to character
abund_18S <- abund_18S%>%  
  mutate(Order = as.character(Order)) %>%
  mutate(Order = replace(Order, Order == 'Unassigned Dinophyceae (Class)', 'Unassigned_Dinophyceae'))
abund_18S_new <- subset(abund_18S, grepl("Unassigned_Dinophyceae|Bacillariophyta_X|Cryptomonadales|Mamiellales|Peridiniales|^Dino-Group-I$|^Dino-Group-II$|Gymnodiniales|^Strombidiida$|^Choreotrichida$", abund_18S$Order)) # Subset to major groups (>2%)
abund_18S_new$Order <- factor(abund_18S_new$Order,levels=c("Bacillariophyta_X", "Mamiellales","Cryptomonadales", "Peridiniales", "Unassigned_Dinophyceae","Dino-Group-I", "Dino-Group-II", "Gymnodiniales", "Strombidiida", "Choreotrichida" )) # Order groups in the plot
p <- ggplot(abund_18S_new, aes(factor(Date), Abundance, group=Order))
p$data$Date = factor(p$data$Date, levels = c(31617, 32217, 32917, 41217, 42017, 50317,60717, 61517, 62217, 70617, 72117, 72717, 80317, 81017, 81617, 82317, 83017, 90617, 92017, 92617, 100417, 101117, 101817, 110817, 111617, 112117, 120717, 121417, 11618, 12418, 20118, 20818, 22118)) # Set the order of sampling days on x-axis
p + geom_point(aes(fill = Order)) + 
  geom_smooth(se=TRUE, linetype="solid", size=1, level=0.95, fill = "gray45", color = "black", alpha = 0.4, method = "loess")+ 
  facet_wrap(~ Order, ncol =2, scales = "free_y") + theme_bw() + ylab("Relative Abundance (%)") + xlab("Day") + 
  theme(legend.position = "none")+ guides(fill=guide_legend(ncol=2))+ theme(axis.text.x=element_text(angle=90,vjust=0.5, size=10))
ggsave(filename = "18S_abund.eps", plot = last_plot(),  path = NULL, scale = 1, width = 5, height = 8, dpi = 600,device=cairo_ps) # Save and export ggplot figure as .eps file 

# Spearman correlations for 18S order level groups (relative abundance) and environmental factors
mean.18S = abund_18S_new %>% # Group order level abundance data and environmental factors
  group_by(Date, Order) %>%
  summarise_at(.vars = c("Temp", "Salinity", "Chlorophyll","Phosphate", "Nitrate", "Silicate", "Ammonium", "POC", "PON", "Solar", "Abundance"), .funs = mean)
factors = mean.18S[, c(3:12)] # Subset out the environmental factors
factors = factors[!duplicated(factors$Temp), ] # Remove duplicates
factors = factors[-c(23), ] # Filter out 9/6, where nutrients were not measured 
abund_split = mean.18S[, c(2,13)] # Subset the order level taxonomy and abundance
abund_split = do.call('data.frame', split(abund_split, abund_split$Order)) # Split out taxonomy and abundance columns for each order group
abund_final = abund_split[, c(2,4,6,8,10,12,14,16,18,20)] # Subset the abundance columns we want to correlate with factors
abund_final= abund_final[-c(23), ] # Filter out 9/6, where nutrients were not measured
res2 <- Hmisc::rcorr(as.matrix(factors), (as.matrix(abund_final)), type = "spearman") # Estimate Spearman correlation coefficients and p-values
coeff = res2$r # Subset the coefficients
coeff = coeff[, c(1:10)]
coeff= coeff[-c(1:10), ] 
pvalue = res2$P # Subset the p-values
pvalue[is.na(pvalue)] <- 0
pvalue = pvalue[, c(1:10)]
pvalue= pvalue[-c(1:10), ]
setEPS()
postscript("corrplot_18S.eps") # Set up .eps file to export correlation matrix
corrplot(coeff, p.mat = pvalue, outline = TRUE, type="full", insig="blank", sig.level =0.05, pch.cex = .9, tl.col = "black", tl.srt = 45, method = "color", addgrid.col = "black") # Plot correlation plot 
dev.off() 

# Create UniFrac distance matrix for clustering analysis and dbRDA
ps_rel  <- transform_sample_counts(ps_rare, function(x)100* x / sum(x))  # Normalize the data
ps_18S_subset = subset_samples(ps_rel, sample_names(ps_rel) != "090617-A-E3" & sample_names(ps_rel) != "090617-B-E4" & sample_names(ps_rel) != "090617-C-E5") # Remove samples from 9/6 that do not have nutrient data
wUF = phyloseq::distance(ps_18S_subset, method="unifrac") # Estimate UniFrac distances
wUF.table<-as.matrix(dist(wUF)) # Convert to matrix

# Cluster analysis 18S
set.seed(123)
fviz_nbclust(wUF.table, kmeans, method = "silhouette") + labs(subtitle = "Elbow method") # Use Silhouette method to determine optimal number of clusters
ggsave(filename = "18S_cluster_silhouette.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 4, dpi = 600) 

# Cluster samples with Ward's method and view dendrogram
spe.ward <- hclust(wUF, method = 'ward.D2')
fviz_dend(x = spe.ward, cex = 0.7, lwd = 0.7, k = 2, k_colors = c("#D55E00","#0072B2"),color_labels_by_k = FALSE, rect = TRUE,rect_fill = TRUE,rect_border = c("#D55E00","#0072B2"), horiz = TRUE) 
ggsave(filename = "18S_dendro.eps", plot = last_plot(), path = NULL, scale = 1, width = 5, height = 10, dpi = 600,device=cairo_ps) 

# Run dbRDA for 18S
metadata <- as(sample_data(ps_18S_subset), "data.frame") # Export metadata from phyloseq object
metadata[, 4:14] <- log1p((metadata[4:14])) # Log transform metadata variables (e.g., temperature, salinity, POC, etc.)
meta = metadata[, c(4:14)]
meta = as.matrix(meta)
merge = cbind(meta,wUF.table) # Merge metadata with distance matrix
merge = as.data.frame(merge) # Convert to data frame
fix(merge) # Visualize merged data frame
rownames(merge) <- NULL
species=merge[,12:104] # Separate species (UniFrac) and environment (variables) for dbRDA
environment=merge[,1:11] 

# Run dbRDA using UniFrac distance and log transformed metadata
dbRDA = dbrda(species ~  Temp  + Salinity + Chlorophyll + Nitrate + Ammonium + Phosphate + Silicate + Solar + POC + PON, environment, dist ="unifrac")
dbrda1 <- dbrda(species ~ ., data=environment) # Stepwise analysis to determine significant factors - DO not significant, removed from ordination 
dbrda0 <- dbrda(species ~ 1, data=environment)
step.env <- ordistep(dbrda0, scope=formula(dbrda1), direction='both',  perm.max=999)

# Extract site and species scores
scores_dbRDA=vegan::scores(dbRDA)
site_scores=scores_dbRDA$sites
fix(site_scores) # Visualize the scores
species_scores=scores_dbRDA$sites
site_scores_environment=cbind(site_scores,environment)
correlations=cor(site_scores_environment) 
fix(correlations) # Visualize the correlations
site_scores = as.data.frame(site_scores)
species_scores = as.data.frame(species_scores)

# Prepare data for 18S dbRDA plot
meta_new = metadata[, c(3,15)] # Subset from the metadata again to include categories for plotting (month and cluster)
RDAscores <- cbind(site_scores, meta_new) # Combine site scores with metadata for plotting
PCAvect <- scores(dbRDA, display = "species") %>% 
  as.data.frame()
arrowmat <- vegan::scores(dbRDA, display = "bp")  # Include factors as biplot arrows on the plot
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = dbRDA1, yend = dbRDA2, x = 0, y = 0, shape = NULL, color = NULL, label = labels) # Maps the arrows on the dbRDA plot
label_map <- aes(x = 1.3 * dbRDA1, y = 1.3 * dbRDA2, shape = NULL, color = NULL, label = labels) # Maps metadata labels
arrowhead = arrow(length = unit(0.02, "npc"))
p <- ggplot(data=RDAscores, aes(x=dbRDA1, y=dbRDA2), color = "Month", inherit.aes = FALSE) + guides(color="none")
p$data$Month <- factor(p$data$Month, levels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb")) # Set order of months
p + geom_point(aes(fill=Month, shape=Cluster),size = 5, colour = "black") + scale_fill_viridis(discrete=TRUE, option="plasma", direction = -1)+
  scale_shape_manual(values = c(21, 24)) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) + 
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) + theme_bw() + 
  geom_segment(mapping = arrow_map, size = 0.5, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4, color="black", data = arrowdf, show.legend = FALSE)
ggsave(filename = "18S_dbRDA.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 7, height = 5, dpi = 600)

# Estimate Shannon diversity and richness using phyloseq object with singletons and compare between clusters
cluster = sample_data(ps)$Cluster # Define variable from metadata (clusters)
alpha_div <- estimate_richness(ps, measures=c("Observed", "Shannon")) # Estimate richness and diversity
alpha_div = cbind(alpha_div, cluster)
shapiro.test(alpha_div$Observed) # Shapiro Wilk's normality test
shapiro.test(alpha_div$Shannon)
wilcox.test(Shannon ~ cluster, data = alpha_div) # Wilcoxon test for Shannon and paired t-test for richness
t.test(Observed ~ cluster, data = alpha_div)
compare_div <- alpha_div %>% # Compare mean and SD between clusters
  group_by(cluster) %>% summarise(across(everything(), list(mean,sd)))

# Compare metadata variables between clusters - no need to repeat for 16S (same metadata)
metadata_comp <- as(sample_data(ps_18S_subset), "data.frame") # Subset metadata 
meta_comp_new = metadata_comp[!duplicated(metadata_comp$Temp), ] # Identify repeated values - no need for triplicates here
norm_meta <-apply(meta_comp_new[,4:13], 2, function(x) shapiro.test(x)$p.value) # Perform Shapiro Wilk's test on all factors of interest (p-value > 0.05 use t-test; p-value < 0.05 use Wilcoxon test)
norm_meta <- as.data.frame(norm_meta)
ttest = lapply(meta_comp_new[,c(5,10,11)], function(x) t.test(x ~ meta_comp_new$Cluster)) # T-test on salinity, silicate, and solar
ttest$Salinity # Check significance of each factor
wilcox = lapply(meta_comp_new[,c(4,6:9,12:13)], function(x) wilcox.test(x ~ meta_comp_new$Cluster)) # Wilcox for other factors 
wilcox$Temp # Check significance 
subset = meta_comp_new[, c(4:13, 15)] #) Subset to the environmental factors and cluster category
Compare_fact <- subset %>% # Compare mean and SD of factors between clusters
  group_by(Cluster) %>% summarise(across(everything(), list(mean,sd)))

# Export metadata - Table S2
write.csv(meta_comp_new, file="Table_S2.csv", row.names=T)

### Repeat same analyses for 16S data

# Upload 16S count, taxonomy, and rooted tree files (.qza files)
table_16S <- read_qza("16S-table.qza")
count_tab_16S <- table_16S$data %>% as.data.frame() 
colnames(count_tab_16S)[1:3] <- c("011818-A-G10", "011818-B-G11", "011818-C-G12") # Some date labels were slightly off. Must match 18S samples for network analysis
colnames(count_tab_16S)[57:59] <- c("081717-A-D7", "081717-B-D8", "081717-C-D9")
colnames(count_tab_16S)[89:91] <- c("120817-A-G4", "120817-B-G5", "120817-C-G6")
taxonomy_16S <- read_qza("16S-tax.qza")
tax_tab_16S <- taxonomy_16S$data %>% 
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)
tax_tab_16S$Kingdom <- gsub("^.{0,5}", "", tax_tab_16S$Kingdom) # Clean up the taxonomy names for each level
tax_tab_16S$Phylum <- gsub("^.{0,5}", "", tax_tab_16S$Phylum)
tax_tab_16S$Class <- gsub("^.{0,5}", "", tax_tab_16S$Class)
tax_tab_16S$Order <- gsub("^.{0,5}", "", tax_tab_16S$Order)
tax_tab_16S$Family <- gsub("^.{0,5}", "", tax_tab_16S$Family)
tax_tab_16S$Genus <- gsub("^.{0,5}", "", tax_tab_16S$Genus)
tax_tab_16S$Species<- gsub("^.{0,5}", "", tax_tab_16S$Species)
tree_16S <- read_qza("16S-rooted-tree.qza") # Upload and reformat tree
tree <- tree_16S$data
new_tre <- ape::multi2di(tree)

# Load 16S metadata file
sample_info_tab_16S <- read.table("Sampleinfo_16S.txt", header=TRUE, row.names=1, check.names=F, sep="\t")

# Merge into phyloseq object
ps_16S <- phyloseq(tax_table(as.matrix(tax_tab_16S)), otu_table(count_tab_16S, taxa_are_rows = T),  sample_data(sample_info_tab_16S),phy_tree(new_tre))

# Merge 16S ASV information - Table S1
OTU = as(otu_table(ps_16S), "matrix")
TAX = as(tax_table(ps_16S), "matrix")
merge_16S <- cbind(OTU,TAX)
write.csv(merge_16S, file="Table_S1_16S.csv", row.names=T)

# Rename ASVs (bacteria ASVs labeled "bASVs" to distinguish from 18S)
taxa_names(ps_16S) <- paste0("bASV", seq(ntaxa(ps_16S)))

# Remove taxa not assigned at the kingdom level
ps_16S <- subset_taxa(ps_16S, Kingdom!="Unassigned", Prune = T)

# Remove chloroplast and mitochondria reads 
psnew_16S = subset_taxa(ps_16S, Order != "Chloroplast" |is.na(Order))
psnew_16S = subset_taxa(psnew_16S, Family !="Mitochondria" |is.na(Family))

# Adds "Unassigned" label names for missing taxa
psnew_16S = name_taxa(psnew_16S, label = "Unassigned")

# Estimate minimum, mean, and maximum read counts
psnew_16S_min <- min(sample_sums(psnew_16S))
psnew_16S_mean <- mean(sample_sums(psnew_16S))
psnew_16S_max <- max(sample_sums(psnew_16S))

# Rarefaction curve for 16S
rare_16S <- ggrare(psnew_16S, step = 100, plot = TRUE, parallel = FALSE, se = FALSE)
rare_16S + theme(legend.position = "none")  + theme_bw() + theme(legend.position = "none")
ggsave(filename = "16S_rare.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 4, height = 4, dpi = 600) 

# Remove singletons 
ps_16S_filt = filter_taxa(psnew_16S, function (x) {sum(x) > 1}, prune=TRUE) # Remove singletons

# Rarefy data to even sampling depth
ps_16S_rare <- rarefy_even_depth(ps_16S_filt, sample.size = min(sample_sums(ps_16S_filt)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Find out the most relatively abundant 16S groups (at order level)
order_16S <- tax_glom(ps_16S_rare, taxrank = "Order")
order_16S <- transform_sample_counts(order_16S, function(x)100* x / sum(x))
OTU <- otu_table(order_16S)
TAX <- tax_table(order_16S)[,"Order"]
Average <- as.data.frame(rowMeans(OTU))
names(Average) <- c("Mean")
Table <- merge(TAX, Average, by=0, all=TRUE)
Table$Row.names = NULL # View abundances at order level - we are interested in the top 10 most abundant groups
head(Table)

# Temporal abundance plots for top 16S   groups
abund_16S <- psmelt(order_16S) # Melt data
abund_16S$Order <- as.character(abund_16S$Order) # Convert to character
abund_16S_new <- subset(abund_16S, grepl("^SAR11 clade$|Flavobacteriales|Rhodobacterales|Actinomarinales|^SAR86 clade$|Betaproteobacteriales|Puniceispirillales|Oceanospirillales|Cellvibrionales|Thiomicrospirales", abund_16S$Order)) # Subset to major groups (>2%)
abund_16S_new$Order <- factor(abund_16S_new$Order,levels=c("SAR11 clade","Flavobacteriales", "Rhodobacterales", "Actinomarinales", "SAR86 clade", "Betaproteobacteriales", "Puniceispirillales","Oceanospirillales", "Cellvibrionales", "Thiomicrospirales")) # Reorder them for the plot
p <- ggplot(abund_16S_new, aes(factor(Date), Abundance, group=Order))
p$data$Date = factor(p$data$Date, levels = c(31617, 32217, 32917, 41217, 42017, 50317,60717, 61517, 62217, 70617, 72117, 72717, 80317, 81017, 81617, 82317, 83017, 90617, 92017, 92617, 100417, 101117, 101817, 110817, 111617, 112117, 120717, 121417, 11618, 12418, 20118, 20818, 22118)) # Order sampling day on x-axis
p + geom_point(aes(fill = Order)) + 
  geom_smooth(se=TRUE, linetype="solid", size=1, level=0.95, fill = "gray45", color = "black", alpha = 0.4, method = "loess")+ 
  facet_wrap(~ Order, ncol =2, scales = "free_y") + theme_bw() + ylab("Relative Abundance (%)") + xlab("Day") + 
  theme(legend.position = "none")+ guides(fill=guide_legend(ncol=2))+ theme(axis.text.x=element_text(angle=90,vjust=0.5, size=10))
ggsave(filename = "16S_abund.eps", plot = last_plot(), path = NULL, scale = 1, width = 5, height = 8, dpi = 600, device=cairo_ps) 

# Spearman correlations for 16S order level groups (relative abundance) and environmental factors
mean.16S = abund_16S_new %>% # Group order level abundance data and environmental factors
  group_by(Date, Order) %>%
  summarise_at(.vars = c("Temp", "Salinity", "Chlorophyll","Phosphate", "Nitrate", "Silicate", "Ammonium", "POC", "PON", "Solar", "Abundance"), .funs = mean)
factors = mean.16S[, c(3:12)] # Subset out the environmental factors
factors = factors[!duplicated(factors$Temp), ] # Remove duplicates
factors = factors[-c(23), ] # Filter out 9/6, where nutrients were not measured 
abund_split = mean.16S[, c(2,13)] # Subset the order level taxonomy and abundance
abund_split = do.call('data.frame', split(abund_split, abund_split$Order)) # Split out taxonomy and abundance columns for each order group
abund_final = abund_split[, c(2,4,6,8,10,12,14,16,18,20)] # Subset the abundance columns we want to correlate with factors
abund_final= abund_final[-c(23), ] # Filter out 9/6, where nutrients were not measured
res2 <- Hmisc::rcorr(as.matrix(factors), (as.matrix(abund_final)), type = "spearman") # Estimate Spearman correlation coefficients and p-values 
coeff = res2$r # Subset the coefficients
coeff = coeff[, c(1:10)]
coeff= coeff[-c(1:10), ] 
pvalue = res2$P # Subset the p-values
pvalue[is.na(pvalue)] <- 0
pvalue = pvalue[, c(1:10)]
pvalue= pvalue[-c(1:10), ]
setEPS()
postscript("corrplot_16S.eps") # Set up .eps file to export 16S correlation matrix
corrplot(coeff, p.mat = pvalue, outline = TRUE, type="full", insig="blank", sig.level =0.05, pch.cex = .9, tl.col = "black", tl.srt = 45, method = "color", addgrid.col = "black") # Plot correlation plot 
dev.off() 

# Create UniFrac distance matrix for clustering analysis and dbRDA
ps_16S_rel  <- transform_sample_counts(ps_16S_rare, function(x) x / sum(x)) # Normalize the data
ps_16S_subset <- subset_samples(ps_16S_rel, sample_names(ps_16S_rel) != "090617-A-E3" & sample_names(ps_16S_rel) != "090617-B-E4" & sample_names(ps_16S_rel) != "090617-C-E5") # Remove samples from 9/6 that do not have nutrient data
wUF = phyloseq::distance(ps_16S_subset, method="unifrac") # Estimate UniFrac distances
wUF.table<-as.matrix(dist(wUF)) # Convert to matrix

# Cluster analysis 16S
set.seed(123)
fviz_nbclust(wUF.table, kmeans, method = "silhouette") + labs(subtitle = "Elbow method") # Use Silhouette method to determine optimal number of clusters
ggsave(filename = "16S_cluster_silhouette.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 4, dpi = 600) 

# Cluster samples with Ward's method and view dendrogram
spe.ward <- hclust(wUF, method = 'ward.D2')
fviz_dend(x = spe.ward, cex = 0.7, lwd = 0.7, k = 2, k_colors = c("#D55E00","#0072B2"),color_labels_by_k = FALSE, rect = TRUE,rect_fill = TRUE,rect_border = c("#D55E00","#0072B2"), horiz = TRUE) 
ggsave(filename = "16S_dendro.eps", plot = last_plot(), path = NULL, scale = 1, width = 5, height = 10, dpi = 600,device=cairo_ps) 

# Run dbRDA for 16S
metadata <- as(sample_data(ps_16S_subset), "data.frame") # Export metadata from phyloseq object
metadata[, 4:14] <- log1p((metadata[4:14])) # Log transform metadata variables (e.g., temperature, salinity, POC, etc.)
meta = metadata[, c(4:14)]
meta = as.matrix(meta)
merge = cbind(meta,wUF.table) # Merge metadata with distance matrix
merge = as.data.frame(merge) # Convert to data frame
fix(merge) # Visualize merged data frame
rownames(merge) <- NULL
species=merge[,12:102] # Separate species (UniFrac) and environment (variables) for dbRDA
environment=merge[,1:11] 

# Run dbRDA using UniFrac distance and log transformed metadata
dbRDA = dbrda(species ~  Temp  + Salinity + Chlorophyll + Nitrate + Ammonium + Phosphate + Silicate + Solar + POC + PON, environment, dist ="unifrac")
dbrda1 <- dbrda(species ~ ., data=environment) # Stepwise analysis to determine significant factors - DO not significant, removed from ordination 
dbrda0 <- dbrda(species ~ 1, data=environment)
step.env <- ordistep(dbrda0, scope=formula(dbrda1), direction='both',  perm.max=999)

# Extract site and species scores
scores_dbRDA=vegan::scores(dbRDA)
site_scores=scores_dbRDA$sites
fix(site_scores) # Visualize the scores
species_scores=scores_dbRDA$sites
site_scores_environment=cbind(site_scores,environment)
correlations=cor(site_scores_environment) 
fix(correlations) # Visualize the correlations
site_scores = as.data.frame(site_scores)
species_scores = as.data.frame(species_scores)

# Prepare data for 16S dbRDA plot
meta_new = metadata[, c(3,15)] # Subset from the metadata again to include categories for plotting (month and cluster)
RDAscores <- cbind(site_scores, meta_new) # Combine site scores with metadata for plotting
PCAvect <- scores(dbRDA, display = "species") %>% 
  as.data.frame()
arrowmat <- vegan::scores(dbRDA, display = "bp")  # Include factors as biplot arrows on the plot
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = dbRDA1, yend = dbRDA2, x = 0, y = 0, shape = NULL, color = NULL, label = labels) # Maps the arrows on the dbRDA plot
label_map <- aes(x = 1.3 * dbRDA1, y = 1.3 * dbRDA2, shape = NULL, color = NULL, label = labels) # Maps metadata labels
arrowhead = arrow(length = unit(0.02, "npc"))
p <- ggplot(data=RDAscores, aes(x=dbRDA1, y=dbRDA2), color = "Month", inherit.aes = FALSE) + guides(color="none")
p$data$Month <- factor(p$data$Month, levels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb")) # Set order of months
p + geom_point(aes(fill=Month, shape=Cluster),size = 5, colour = "black") + scale_fill_viridis(discrete=TRUE, option="plasma", direction = -1)+
  scale_shape_manual(values = c(21, 24)) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) + 
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) + theme_bw() + 
  geom_segment(mapping = arrow_map, size = 0.5, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4, color="black", data = arrowdf, show.legend = FALSE)
ggsave(filename = "16S_dbRDA.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 7, height = 5, dpi = 600)

# Estimate 16S Shannon diversity and richness and compare between clusters
cluster = sample_data(psnew_16S)$Cluster # Define variable from metadata (clusters)
alpha_div_16S <- estimate_richness(psnew_16S, measures=c("Observed", "Shannon")) # Estimate richness and diversity
alpha_div_16S = cbind(alpha_div_16S, cluster)
shapiro.test(alpha_div_16S$Observed) # Shapiro Wilk's normality test
shapiro.test(alpha_div_16S$Shannon)
wilcox.test(Shannon ~ cluster, data = alpha_div_16S) # Wilcoxon test for Shannon and paired t-test for richness
t.test(Observed ~ cluster, data = alpha_div_16S)
compare_div <- alpha_div_16S %>% # Compare mean and SD between clusters
  group_by(cluster) %>% summarise(across(everything(), list(mean,sd)))

### Prepare 16S and 18S ASVs for SPIEC-EASI network analysis
X <- otu_table(ps) # Export the ASV count, taxonomy, and sample data from the original phyloseq object - rooted tree files are not matching
Y <- tax_table(ps)
Z <- sample_data(ps)
ps_network = merge_phyloseq(X,Y,Z) # Merge into phyloseq object
ps_network = subset_samples(ps_network, sample_names(ps_network) != "090617-A-E3" & sample_names(ps_network) != "090617-B-E4" & sample_names(ps_network) != "090617-C-E5") # Remove 9/6 - no nutrient data
ps_network <- rarefy_even_depth(ps_network, sample.size = min(sample_sums(ps_network)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) # Rarefy

x <- otu_table(psnew_16S) # Repeat for 16S data
y <- tax_table(psnew_16S)
z <- sample_data(psnew_16S)
ps_network_16S = merge_phyloseq(x,y,z)
ps_network_16S = subset_samples(ps_network_16S, sample_names(ps_network_16S) != "090617-A-E3" & sample_names(ps_network_16S) != "090617-B-E4" & sample_names(ps_network_16S) != "090617-C-E5")
ps_network_16S <- rarefy_even_depth(ps_network_16S, sample.size = min(sample_sums(ps_network_16S)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Subset 16S & 18S phyloseq objects based on clusters (Cluster 1 = March-May/November-February; Cluster 2 = June-October)
ps_18S_Clus1 <- subset_samples(ps_network, Cluster == "Cluster_1")
ps_18S_Clus2 <- subset_samples(ps_network, Cluster == "Cluster_2")
ps_16S_Clus1 <- subset_samples(ps_network_16S, Cluster == "Cluster_1")
ps_16S_Clus2 <- subset_samples(ps_network_16S, Cluster == "Cluster_2")

# Remove samples that are in the 18S but not in the 16S - sample IDs in count tables have to match 
ps_18S_Clus1  = subset_samples(ps_18S_Clus1, sample_names(ps_18S_Clus1) != "031617-B-A2")
ps_18S_Clus2 = subset_samples(ps_18S_Clus2, sample_names(ps_18S_Clus2) != "092017-C-E8")

# Trim sets to include top 150 16S and 18S ASVs - improve sparsity and avoid ambiguous relationships
filter_18S_Clus1 = prune_taxa(names(sort(taxa_sums(ps_18S_Clus1), TRUE))[1:150], ps_18S_Clus1)
filter_18S_Clus2 = prune_taxa(names(sort(taxa_sums(ps_18S_Clus2), TRUE))[1:150], ps_18S_Clus2)
filter_16S_Clus1 = prune_taxa(names(sort(taxa_sums(ps_16S_Clus1), TRUE))[1:150], ps_16S_Clus1)
filter_16S_Clus2 = prune_taxa(names(sort(taxa_sums(ps_16S_Clus2), TRUE))[1:150], ps_16S_Clus2)

# Merge phyloseq objects for each cluster and convert OTU count table to data frame
all_Clus2 <- merge_phyloseq(filter_16S_Clus2, filter_18S_Clus2)
Clus2_all.otu <- t(data.frame(phyloseq::otu_table(all_Clus2), check.names = FALSE))
all_Clus1 <- merge_phyloseq(filter_16S_Clus1, filter_18S_Clus1)
Clus1_all.otu <- t(data.frame(phyloseq::otu_table(all_Clus1), check.names = FALSE))

# Convert individual Cluster 1/2 phyloseq objects for 16S/18S for input into SPIEC-EASI
Clus1_16S.otu <- t(data.frame(phyloseq::otu_table(filter_16S_Clus1), check.names = FALSE))
Clus1_18S.otu <- t(data.frame(phyloseq::otu_table(filter_18S_Clus1), check.names = FALSE))
Clus2_16S.otu <- t(data.frame(phyloseq::otu_table(filter_16S_Clus2), check.names = FALSE))
Clus2_18S.otu <- t(data.frame(phyloseq::otu_table(filter_18S_Clus2), check.names = FALSE))

# Run SPIEC-EASI for Cluster 1 (16S + 18S)
se <- spiec.easi(list(Clus1_16S.otu, Clus1_18S.otu), method ='mb', nlambda = 20, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05))

# Function to get weights and edge data into correct format for exporting (code via biovcnet.github.io) 
convertSEToTable <- function(se,sp.names){
  sebeta <- symBeta(getOptBeta(se), mode='maxabs') 
  elist     <- summary(sebeta)
  elist[,1] <- sp.names[elist[,1]]
  elist[,2] <- sp.names[elist[,2]]
  elist[,4] <- paste(elist[,1],elist[,2])
  full_e <- expand.grid(sp.names,sp.names)
  rownames(full_e) <- paste(full_e[,1],full_e[,2])
  full_e[,"Weight"] <- 0
  full_e[elist[,4],"Weight"] <- elist[,3]
  x <- expand.grid(1:length(sp.names),1:length(sp.names))
  full_e[x[,"Var1"]>x[,"Var2"],"Weight"] <- NA
  return(as.data.frame(full_e,stringsAsFactors=F))
}

tab.se <- convertSEToTable(se, sp.names=colnames(Clus1_all.otu)) 
tab.se.filtered <- tab.se %>% filter(is.finite(Weight))
Clus1_filt <- subset(tab.se.filtered, abs(Weight) > 0) # Remove values that are 0

# Export .csv file with ASV-ASV edges and weights needed to build network (Cytoscape) - note, weights are all significant and do not have a traditional p-value
# Export taxonomy key for the Cluster 1 network to assign taxonomy to each node - add an "ASV" label to the taxonomy .txt file following export
write_csv(Clus1_filt, "Skio-Clus1-network_04292022.csv")
taxALL_Clus1 <- tax_table(all_Clus1)
taxALL_Clus1 <- as.data.frame(taxALL_Clus1, check.names=FALSE)
write.table(taxALL_Clus1,"Skio-Clus1-TAX_04292022.txt",sep="\t", quote=FALSE)

# Run SPIEC-EASI for Cluster 2 (16S + 18S)
se2 <- spiec.easi(list(Clus2_16S.otu, Clus2_18S.otu), method ='mb', nlambda = 20, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05))

convertSEToTable <- function(se2,sp.names){
  sebeta <- symBeta(getOptBeta(se2), mode='maxabs') 
  elist     <- summary(sebeta)
  elist[,1] <- sp.names[elist[,1]]
  elist[,2] <- sp.names[elist[,2]]
  elist[,4] <- paste(elist[,1],elist[,2])
  full_e <- expand.grid(sp.names,sp.names)
  rownames(full_e) <- paste(full_e[,1],full_e[,2])
  full_e[,"Weight"] <- 0
  full_e[elist[,4],"Weight"] <- elist[,3]
  x <- expand.grid(1:length(sp.names),1:length(sp.names))
  full_e[x[,"Var1"]>x[,"Var2"],"Weight"] <- NA
  return(as.data.frame(full_e,stringsAsFactors=F))
}

tab.se <- convertSEToTable(se2,sp.names=colnames(Clus2_all.otu)) 
tab.se.filtered <- tab.se %>% filter(is.finite(Weight))
Clus2_filt <- subset(tab.se.filtered, abs(Weight) > 0) # Remove values that are 0

# Export .csv file with Cluster 2 ASV pairings and associated taxonomy labels for nodes (Cytoscape)
write_csv(Clus2_filt, "Skio-Clus2-network_04292022.csv")
taxALL_Clus2 <- tax_table(all_Clus2)
taxALL_Clus2 <- as.data.frame(taxALL_Clus2, check.names=FALSE)
write.table(taxALL_Clus2,"Skio-Clus2-TAX_04292022.txt",sep="\t", quote=FALSE)

# Extract ASV-ASV pairings for Cluster 1 and 2 networks - Table S4
taxALL_Clus1 <- cbind(rownames(taxALL_Clus1), data.frame(taxALL_Clus1, row.names=NULL))
colnames(taxALL_Clus1)[1] <- "Var1"
newtab_Clus1 <- merge(Clus1_filt, taxALL_Clus1[, c("Var1","Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_Clus1)[1] <- "Var2"
newtab_final_Clus1 <- merge(newtab_Clus1, taxALL_Clus1[, c("Var2", "Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var2")
write.csv(newtab_final_Clus1, "Table_S4_Clus1.csv",row.names=F)

taxALL_Clus2 <- cbind(rownames(taxALL_Clus2), data.frame(taxALL_Clus2, row.names=NULL))
colnames(taxALL_Clus2)[1] <- "Var1"
newtab_Clus2 <- merge(Clus2_filt, taxALL_Clus2[, c("Var1","Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL_Clus2)[1] <- "Var2"
newtab_final_Clus2 <- merge(newtab_Clus2, taxALL_Clus2[, c("Var2","Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var2")
write.csv(newtab_final_Clus2, "Table_S4_Clus2.csv", row.names=F)

# Create a new column in ASV-ASV pairings table to denote positive and negative edges (or sign)
newtab_final_Clus1$Sign <- ifelse(newtab_final_Clus1$Weight > 0, "Pos", "Neg") # Positive sign if edge weight is > 0; if not, it is negative
newtab_final_Clus1$Cluster <- 'Cluster_1' # New column to denote these are Cluster 1 relationships
newtab_final_Clus2$Sign <- ifelse(newtab_final_Clus2$Weight > 0, "Pos", "Neg")
newtab_final_Clus2$Cluster <- 'Cluster_2' # New column to denote these are Cluster 2 relationships
newtab_merge <- rbind(newtab_final_Clus1, newtab_final_Clus2) # Combine ASV-ASV pairings tables
newtab_merge$Edge = paste(newtab_merge$Order.x, newtab_merge$Order.y,sep="-") # Create new column called "Edge" that has the order level ASV pairing in one cell (e.g., Bacillariophyta-Flavobacteriales)
newtab_merge$Edge2 = paste(newtab_merge$Kingdom.x, newtab_merge$Kingdom.y,sep="-") # Create similar edge column for kingdom level
newtab_merge$Edge3 = paste(newtab_merge$Class.x, newtab_merge$Class.y,sep="-") # Create similar edge column for class level

# Simplify domain level edge columns to indicate 16S, 18S, or a mixed (cross domain) relationship - new column "Type" will be created
newtab_merge[(newtab_merge$Edge2)=="Bacteria-Bacteria","Type"]<-"Prokaryote-Prokaryote"
newtab_merge[(newtab_merge$Edge2)=="Archaea-Bacteria","Type"]<-"Prokaryote-Prokaryote"
newtab_merge[(newtab_merge$Edge2)=="Bacteria-Archaea","Type"]<-"Prokaryote-Prokaryote"
newtab_merge[(newtab_merge$Edge2)=="Eukaryota-Archaea","Type"]<-"Prokaryote-Eukaryote"
newtab_merge[(newtab_merge$Edge2)=="Archaea-Eukaryota","Type"]<-"Prokaryote-Eukaryote"
newtab_merge[(newtab_merge$Edge2)=="Eukaryota-Eukaryota","Type"]<-"Eukaryote-Eukaryote"
newtab_merge[(newtab_merge$Edge2)=="Bacteria-Eukaryota","Type"]<-"Prokaryote-Eukaryote"
newtab_merge[(newtab_merge$Edge2)=="Eukaryota-Bacteria","Type"]<-"Prokaryote-Eukaryote"

# Number of domain level edges in Cluster 1 vs. 2 networks
newtab_merge_kingdom <- newtab_merge[c("Kingdom.x", "Kingdom.y")] # Subset kingdom
sign <- newtab_merge[c("Sign")] # Subset correlation sign
cluster <- newtab_merge[c("Cluster")] # Subset cluster
type <- newtab_merge[c("Type")] # Subset cluster

for (i in 1:nrow(newtab_merge_kingdom ))
{
  newtab_merge_kingdom [i, ] = sort(newtab_merge_kingdom [i, ])
}

merge_kingdom_final = cbind(newtab_merge_kingdom,cluster, sign, type) # Merge pairings, cluster, and sign of the relationship
df <- as.data.frame(table(merge_kingdom_final[c(3:5)])) # Subset columns with edge, cluster, and sign
df = df[with(df, order(-Freq)), ] # Estimate frequency of each type of kingdom level relationship (both positive and negative) - order from largest to smallest

# Move on to order level relationships - first establish the most prevalent 16S-16S, 18S-18S, and 16S-18S relationships
newtab_merge_order <- newtab_merge[c("Order.x", "Order.y")] # Subset order
type <- newtab_merge[c("Type")]

for (i in 1:nrow(newtab_merge_order ))
{
  newtab_merge_order  [i, ] = sort(newtab_merge_order [i, ])
}

newtab_merge_order <- newtab_merge_order %>%
  mutate(Order.y = as.character(Order.y)) %>%
  mutate(Order.y = replace(Order.y, Order.y == 'Unassigned Dinophyceae (Class)', 'Unassigned_Dinophyceae')) %>% # Rename Unassigned Dinophyceae group (code does not like spaces in taxa names)
  mutate(Order.x = as.character(Order.x)) %>%
  mutate(Order.x = replace(Order.x, Order.x == 'Unassigned Dinophyceae (Class)', 'Unassigned_Dinophyceae'))
  
newtab_merge_order$Edge = paste(newtab_merge_order$Order.x, newtab_merge_order$Order.y,sep="-") # Create an edge with both orders in the same cell
merge_order_final = cbind(newtab_merge_order, cluster, type) # Merge information 
df <- as.data.frame(table(merge_order_final[c(3:5)]))
df = df[with(df, order(-Freq)), ] # List of the most prevalent order level groups across both clusters

# Make a new merged data frame that includes the sign of the order level relationship
merge_order_final = cbind(newtab_merge_order,sign, cluster, type)
df <- as.data.frame(table(merge_order_final[c(3:6)])) # Subset information
df = df[with(df, order(-Freq)), ] # Order relationships based on prevalence 

# Subset to the top ten 16S-18S relationships
mixed <- df[ df[["Type"]] == "Prokaryote-Eukaryote" , ]
mixed[mixed == 0] <- NA
df_mix <- subset(mixed, grepl("Bacillariophyta_X-Flavobacteriales|Gymnodiniales-Rhodobacterales|Bacillariophyta_X-Rhodobacterales|Actinomarinales-Cryptomonadales|Flavobacteriales-Gymnodiniales|Flavobacteriales-Strombidiida|Peridiniales-Rhodobacterales|Actinomarinales-Bacillariophyta_X|Bacillariophyta_X-Betaproteobacteriales|^SAR11 clade-Unassigned_Dinophyceae$",mixed$Edge)) # Subset the top ten most prevalent cross domain relationships
p <- ggplot(df_mix, aes(x = Cluster, y = Freq))+ geom_col(aes(fill = Sign), width = 0.8, colour = "black", size = 1)
p$data$Edge <- factor(p$data$Edge, levels = c("Bacillariophyta_X-Flavobacteriales", "Gymnodiniales-Rhodobacterales", "Bacillariophyta_X-Rhodobacterales", "Actinomarinales-Cryptomonadales", "Flavobacteriales-Gymnodiniales", "Flavobacteriales-Strombidiida", "Peridiniales-Rhodobacterales", "Actinomarinales-Bacillariophyta_X", "Bacillariophyta_X-Betaproteobacteriales", "SAR11 clade-Unassigned_Dinophyceae"))#, "DGII-SAR11", "SAR11-Strombidiida", "SAR11-Unassigned_Dinophyceae", "Bacillariophyta_X-Betaproteobacteriales", "Bacillariophyta_X-Nitrosopumilales", "DGII-Puniceispirillales"))# Change to desired facet labels
edges.labs <- c("Bacill-Flavo", "Gymno-Rhodo", "Bacill-Rhodo","Actino-Crypto", "Flavo-Gymno", "Flavo-Strom", "Peridin-Rhodo", "Actino-Bacill", "Bacill-Betapro", "SAR11-Udino")
names(edges.labs) <- c("Bacillariophyta_X-Flavobacteriales", "Gymnodiniales-Rhodobacterales", "Bacillariophyta_X-Rhodobacterales", "Actinomarinales-Cryptomonadales", "Flavobacteriales-Gymnodiniales", "Flavobacteriales-Strombidiida", "Peridiniales-Rhodobacterales", "Actinomarinales-Bacillariophyta_X", "Bacillariophyta_X-Betaproteobacteriales", "SAR11 clade-Unassigned_Dinophyceae")  
p + scale_fill_manual(values=c("#AF8DC3", "#7FBF7B"))  + theme_bw() +
  theme(axis.title.x =element_blank())  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,10,2), limits=c(0, 10))+
  labs(y="Number of edges", size = 14)+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust = 1, size=14)) +  theme(axis.text.y=element_text(size=14)) + facet_wrap(~ Edge, ncol =5, labeller = labeller(Edge = edges.labs)) 
ggsave(filename = "mixed_edges.eps", last_plot(), path = NULL, scale = 1, width = 6, height = 3, dpi = 600, device=cairo_ps)

# Subset the top ten 18S-18S relationships
euk <- df[ df[["Type"]] == "Eukaryote-Eukaryote" , ]
euk[euk == 0] <- NA
df_euk <- subset(euk, grepl("Bacillariophyta_X-Bacillariophyta_X|Bacillariophyta_X-Unassigned_Dinophyceae|^Bacillariophyta_X-Dino-Group-II$|Bacillariophyta_X-Cryptomonadales|Bacillariophyta_X-Mamiellales|Bacillariophyta_X-Peridiniales|^Bacillariophyta_X-Dino-Group-I$|Bacillariophyta_X-Gonyaulacales|Bacillariophyta_X-Gymnodiniales|Bacillariophyta_X-Tintinnida",euk$Edge))
p <- ggplot(df_euk, aes(x = Cluster, y = Freq))+ geom_col(aes(fill = Sign), width = 0.8, colour = "black", size = 1)
p$data$Edge <- factor(p$data$Edge, levels = c("Bacillariophyta_X-Bacillariophyta_X", "Bacillariophyta_X-Unassigned_Dinophyceae", "Bacillariophyta_X-Dino-Group-II", "Bacillariophyta_X-Cryptomonadales", "Bacillariophyta_X-Mamiellales", "Bacillariophyta_X-Peridiniales", "Bacillariophyta_X-Dino-Group-I", "Bacillariophyta_X-Gonyaulacales", "Bacillariophyta_X-Gymnodiniales", "Bacillariophyta_X-Tintinnida"))#, "DGII-SAR11", "SAR11-Strombidiida", "SAR11-Unassigned_Dinophyceae", "Bacillariophyta_X-Betaproteobacteriales", "Bacillariophyta_X-Nitrosopumilales", "DGII-Puniceispirillales"))
edges.labs <- c("Bacill-Bacill", "Bacill-Udino", "Bacill-DGII","Bacill-Crypto", "Bacill-Mamiell", "Bacill-Peridin", "Bacill-DGI", "Bacill-Gonyau", "Bacill-Gymno", "Bacill-Tintinn")
names(edges.labs) <- c("Bacillariophyta_X-Bacillariophyta_X", "Bacillariophyta_X-Unassigned_Dinophyceae", "Bacillariophyta_X-Dino-Group-II", "Bacillariophyta_X-Cryptomonadales", "Bacillariophyta_X-Mamiellales", "Bacillariophyta_X-Peridiniales", "Bacillariophyta_X-Dino-Group-I", "Bacillariophyta_X-Gonyaulacales", "Bacillariophyta_X-Gymnodiniales", "Bacillariophyta_X-Tintinnida")  
p + scale_fill_manual(values=c("#AF8DC3", "#7FBF7B"))  + theme_bw() +
  theme(axis.title.x =element_blank())  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,50,5), limits=c(0, 50))+
  labs(y="Number of edges", size = 14)+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust = 1, size=14)) +  theme(axis.text.y=element_text(size=14)) + facet_wrap(~ Edge, ncol =5, labeller = labeller(Edge = edges.labs)) 
ggsave(filename = "euk_edges.eps", last_plot(), path = NULL, scale = 1, width = 6, height = 3, dpi = 600, device=cairo_ps)

# Subset the top ten 16S-16S relationships
bact <- df[ df[["Type"]] == "Prokaryote-Prokaryote" , ]
bact[bact == 0] <- NA
df_bact <- subset(bact, grepl("^SAR11 clade-SAR11 clade$|Flavobacteriales-Flavobacteriales|^Flavobacteriales-SAR11 clade$|Flavobacteriales-Rhodobacterales|Flavobacteriales-Puniceispirillales|^Rhodobacterales-SAR11 clade$|^Betaproteobacteriales-SAR11 clade$|Rhodobacterales-Rhodobacterales|Flavobacteriales-Micrococcales|Puniceispirillales-SAR86 clade",bact$Edge))
p <- ggplot(df_bact, aes(x = Cluster, y = Freq))+ geom_col(aes(fill = Sign), width = 0.8, colour = "black", size = 1)
p$data$Edge <- factor(p$data$Edge, levels = c("SAR11 clade-SAR11 clade", "Flavobacteriales-Flavobacteriales", "Flavobacteriales-SAR11 clade", "Flavobacteriales-Rhodobacterales", "Flavobacteriales-Puniceispirillales", "Rhodobacterales-SAR11 clade", "Betaproteobacteriales-SAR11 clade", "Rhodobacterales-Rhodobacterales", "Flavobacteriales-Micrococcales", "Puniceispirillales-SAR86 clade"))#, "DGII-SAR11", "SAR11-Strombidiida", "SAR11-Unassigned_Dinophyceae", "Bacillariophyta_X-Betaproteobacteriales", "Bacillariophyta_X-Nitrosopumilales", "DGII-Puniceispirillales"))
edges.labs <- c("SAR11-SAR11", "Flavo-Flavo", "Flavo-SAR11","Flavo-Rhodo", "Flavo-Puniceis", "Rhodo-SAR11", "Betapro-SAR11", "Rhodo-Rhodo", "Flavo-Micro", "Puniceis-SAR86")
names(edges.labs) <- c("SAR11 clade-SAR11 clade", "Flavobacteriales-Flavobacteriales", "Flavobacteriales-SAR11 clade", "Flavobacteriales-Rhodobacterales", "Flavobacteriales-Puniceispirillales", "Rhodobacterales-SAR11 clade", "Betaproteobacteriales-SAR11 clade", "Rhodobacterales-Rhodobacterales", "Flavobacteriales-Micrococcales", "Puniceispirillales-SAR86 clade")  
p + scale_fill_manual(values=c("#AF8DC3", "#7FBF7B"))  + theme_bw() +
  theme(axis.title.x =element_blank())  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,30,5), limits=c(0, 30))+
  labs(y="Number of edges", size = 14)+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust = 1, size=14)) +  theme(axis.text.y=element_text(size=14)) + facet_wrap(~ Edge, ncol =5, labeller = labeller(Edge = edges.labs)) 
ggsave(filename = "bact_edges.eps", last_plot(), path = NULL, scale = 1, width = 6, height = 3, dpi = 600, device=cairo_ps)

### Import network properties from Cytoscape and plot degree and centrality
# Plot # of edges per node (degree) for Cluster 1/2 networks - compare at domain level
network_properties <- read.csv("network_properties.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") # Information from the network on degree and closeness centrality
network_properties$Kingdom <- replace(as.character(network_properties$Kingdom), network_properties$Kingdom== "Bacteria", "Prokaryote")
network_properties$Kingdom <- replace(as.character(network_properties$Kingdom), network_properties$Kingdom== "Archaea", "Prokaryote")
shapiro.test(network_properties$Degree) # Shapiro Wilks test for network degree and centrality
shapiro.test(network_properties$ClosenessCentrality)
p <- ggplot(network_properties, aes(x=Cluster, y=Degree, fill=Cluster))
p$data$Cluster <- factor(p$data$Cluster, levels = c("Cluster_1", "Cluster_2"))
p + geom_boxplot(aes(fill=Cluster), outlier.shape = NA, alpha = 0.6, size = 0.6, width = 0.6) + stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_point(aes(fill=Cluster), size =2, shape = 21, colour = "black", position=position_jitterdodge()) + scale_fill_manual(values= c("#0072B2","#D55E00")) + theme_bw() + ylab("Degree") +
  theme(axis.title.x=element_blank()) + scale_y_continuous(expand = c(0, 0), breaks=seq(0,20,5), limits=c(0, 20))+ 
  theme(axis.text.y=element_text(size=12)) + theme(axis.title.y = element_text(size = 16))+ 
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1, size=12, color="black")) + 
  facet_wrap(~Kingdom)+ stat_compare_means(method= "wilcox.test",comparisons = list(c("Cluster_1", "Cluster_2")), label = "p.signif", exact = TRUE) # Compare means between networks (Cluster 1 vs. 2)
ggsave(filename = "Network_degree.eps", last_plot(), path = NULL, scale = 1, width = 7, height = 5, dpi = 600, device=cairo_ps)

# Repeat degree comparisons for major prokaryotic and eukaryotic groups
major_bact <- subset(network_properties, grepl("Alphaproteobacteria|Bacteroidia|Acidimicrobiia|Gammaproteobacteria", network_properties$Class)) # Subset to major groups (>2%)
major_euk <- subset(network_properties, grepl("Dinophyceae|Bacillariophyta|Cryptophyceae|Syndiniales|Spirotrichea|Mamiellophyceae", network_properties$Class)) # Subset to major groups (>2%)
p <- ggplot(major_bact, aes(x=Cluster, y=Degree, fill=Cluster)) # Repeat this with protist groups (major_euk)
p$data$Cluster <- factor(p$data$Cluster, levels = c("Cluster_1", "Cluster_2"))
p + geom_boxplot(aes(fill=Cluster), outlier.shape = NA, alpha = 0.6, size = 0.6, width = 0.6) + stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_point(aes(fill=Cluster), size =2, shape = 21, colour = "black", position=position_jitterdodge()) + scale_fill_manual(values= c("#0072B2","#D55E00")) + theme_bw() + ylab("Degree") +
  theme(axis.title.x=element_blank()) + scale_y_continuous(expand = c(0, 0), breaks=seq(0,20,5), limits=c(0, 20))+ 
  theme(axis.text.y=element_text(size=12)) + theme(axis.title.y = element_text(size = 16))+ 
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1, size=12, color="black")) + 
  facet_wrap(~Class, ncol=2, nrow = 2)+ stat_compare_means(method= "wilcox.test",comparisons = list(c("Cluster_1", "Cluster_2")), label = "p.signif", exact = TRUE)
ggsave(filename = "Network_degree_bact_class.eps", last_plot(), path = NULL, scale = 1, width = 7, height = 5, dpi = 600, device=cairo_ps)

# Plot closeness centrality for Cluster 1/2 networks - compare at domain level
p <- ggplot(network_properties, aes(x=Cluster, y=ClosenessCentrality, fill=Cluster))
p$data$Cluster <- factor(p$data$Cluster, levels = c("Cluster_1", "Cluster_2"))
p + geom_boxplot(aes(fill=Cluster), outlier.shape = NA, alpha = 0.6, size = 0.6, width = 0.6) + stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_point(aes(fill=Cluster), size =2, shape = 21, colour = "black", position=position_jitterdodge()) + scale_fill_manual(values= c("#0072B2","#D55E00")) + theme_bw() + ylab("Centrality") +
  theme(axis.title.x=element_blank()) + scale_y_continuous(expand = c(0, 0), breaks=seq(0.1,0.4,0.1), limits=c(0.1, 0.4))+ 
  theme(axis.text.y=element_text(size=12)) + theme(axis.title.y = element_text(size = 16))+ 
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1, size=12, color="black")) + facet_wrap(~Kingdom) + 
  stat_compare_means(method= "wilcox.test",comparisons = list(c("Cluster_1", "Cluster_2")), label = "p.signif", exact = TRUE)
ggsave(filename = "Network_centrality.eps", last_plot(), path = NULL, scale = 1, width = 7, height = 5, dpi = 600, device=cairo_ps)

# Repeat centrality comparisons for major prokaryotic and eukaryotic groups
p <- ggplot(major_bact, aes(x=Cluster, y=ClosenessCentrality, fill=Cluster))# Repeat this with protist groups (major_euk)
p$data$Cluster <- factor(p$data$Cluster, levels = c("Cluster_1", "Cluster_2"))
p + geom_boxplot(aes(fill=Cluster), outlier.shape = NA, alpha = 0.6, size = 0.6, width = 0.6) + stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_point(aes(fill=Cluster), size =2, shape = 21, colour = "black", position=position_jitterdodge()) + scale_fill_manual(values= c("#0072B2","#D55E00")) + theme_bw() + ylab("Centrality") +
  theme(axis.title.x=element_blank()) + scale_y_continuous(expand = c(0, 0), breaks=seq(0.1,0.4,0.1), limits=c(0.1, 0.4))+ 
  theme(axis.text.y=element_text(size=12)) + theme(axis.title.y = element_text(size = 16))+ 
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1, size=12, color="black")) + 
  facet_wrap(~Class, ncol=2, nrow = 2) + stat_compare_means(method= "wilcox.test",comparisons = list(c("Cluster_1", "Cluster_2")), label = "p.signif", exact = TRUE)
ggsave(filename = "Network_centrality_bact_class.eps", last_plot(), path = NULL, scale = 1, width = 7, height = 5, dpi = 600, device=cairo_ps)