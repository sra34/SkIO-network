### Code for processing 18S and 16S metabarcoding data with QIIME 2 input files (.qza files)

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
library(DESeq2)
library(SpiecEasi)
library(circlize)

### Upload 18S count and taxonomy .qza files
table <- read_qza("18S-table.qza")
count_tab <- table$data %>% as.data.frame() # Convert to data frame 
taxonomy <- read_qza("18S-tax.qza")
tax_tab <- taxonomy$data %>% # Convert to data frame, tab separate and rename taxa levels, and remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom","Supergroup","Division","Class","Order","Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence)
tax_tab[is.na(tax_tab)] <- "Unassigned" # Rename "NA" taxa as "Unassigned"

# Upload 18S metadata
sample_info_tab <- read.table("Sampleinfo_18S.txt", header=TRUE, row.names=1, check.names=F, sep="\t")

# Merge 18S ASV information - Table S1
merge_18S <- cbind(count_tab, tax_tab)
write.csv(merge_18S, file="Table_S1.csv", row.names=T)

# Merge into phyloseq object and rename sequential ASVs
ps <- phyloseq(tax_table(as.matrix(tax_tab)), otu_table(count_tab, taxa_are_rows = T), sample_data(sample_info_tab))
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Remove 18S taxa not assigned at supergroup level
ps <- subset_taxa(ps, Supergroup!="Unassigned", Prune = T)

# Rarefy to even sampling depth
ps_rare <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Plot diversity (richness and Shannon) for 18S community with respect to temperature (<23 vs. 23-31C) - Figure 1C
temp2 = sample_data(ps_rare)$Temp2 # Define variable from metadata with two temperature bins
alpha_div <- estimate_richness(ps_rare, measures=c("Observed", "Shannon")) # Estimate richness
alpha_div = cbind(alpha_div, temp2)
alpha_div_melt = melt(alpha_div, id.vars = "temp2", variable.name="value", value.name ="Diversity")
alpha_plot <- ggplot(data=alpha_div_melt, aes(x= temp2, y=Diversity, fill = temp2))
alpha_plot + geom_boxplot(aes(fill=temp2), outlier.shape = NA, alpha = 0.6, size = 0.6, width = 0.6) + stat_boxplot(geom ='errorbar', width = 0.6) +
geom_point(aes(fill=temp2), size =2, shape = 21, colour = "black", position=position_jitterdodge(0.4)) +
scale_fill_manual(values= c("#0072B2", "#D55E00")) + theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12)) + 
theme(axis.text.y = element_text(size=12))+ theme(legend.text = element_text(colour="black", size=12))+
theme(axis.title.x = element_blank())+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
facet_wrap(~value, scales = "free_y") +theme(strip.background = element_blank(),strip.text.x = element_blank())+
stat_compare_means(method= "t.test",comparisons = list(c("<23", "23-31")),label = "p.signif",  label.y.npc = 0.4)
ggsave(filename = "18S_alpha.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 6, height = 4, units = c("in", "cm", "mm"), dpi = 600) # Example of exporting high quality figures from R

# Shapiro Wilks tests for normality of 18S richness and Shannon diversity
shapiro.test(subset(alpha_div, temp2 == "<23")$Observed)
shapiro.test(subset(alpha_div, temp2 == "23-31")$Observed)
shapiro.test(subset(alpha_div, temp2 == "<23")$Shannon)
shapiro.test(subset(alpha_div, temp2 == "23-31")$Shannon)

# Rarefaction curves for 18S - Figure S1
rare_18S <- ggrare(ps, step = 100, plot = TRUE, parallel = FALSE, se = FALSE)
rare_18S + theme(legend.position = "none")  +theme_bw()+theme(legend.position = "none") 
ggsave(filename = "18S_rare.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 4, height = 4, units = c("in", "cm", "mm"), dpi = 600)

# Run NMDS ordination on the 18S samples - Figure 1B
ps_rare = filter_taxa(ps_rare, function (x) {sum(x) > 1}, prune=TRUE) # Remove singletons
ps_rel  <- transform_sample_counts(ps_rare,function(x)x / sum(x))  # Normalize the data
ps_ord  <- ordinate(ps_rel, "NMDS", "bray")
p = plot_ordination(ps_rel, ps_ord, color="Temp2")
p$layers <- p$layers[-1]
p + theme_bw() + theme(text = element_text(size = 16)) + geom_point(aes(fill=Temp2),size = 4, alpha = 1, shape = 21, colour = "black") + scale_fill_manual(values= c("#0072B2", "#D55E00")) + 
theme(legend.text = element_text(colour="black", size=12))+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "18S_nmds.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 6, height = 4, units = c("in", "cm", "mm"), dpi = 600) 

# ANOSIM (Analysis of Similarity) for 18S
temp2 = get_variable(ps_rel, "Temp2")
temp_ano = anosim(phyloseq::distance(ps_rel, "bray"), temp2)
temp_ano$signif
temp_ano$statistic

# Find out the most abundant 18S groups in terms of their proportion of reads across all samples - focus on groups >1% on average for subsequent plots
class_18S <- tax_glom(ps_rare, taxrank = "Class")
class_18S <- transform_sample_counts(class_18S, function(x)100* x / sum(x))
OTU <- otu_table(class_18S)
TAX <- tax_table(class_18S)[,"Class"]
Average <- as.data.frame(rowMeans(OTU))
names(Average) <- c("Mean")
Table <- merge(TAX, Average, by=0, all=TRUE)
Table$Row.names = NULL
Table <- Table[order(desc(Table$Mean)),]
head(Table)

# Abundance vs. temperature plots for the major 18S groups - Figure 3
abund_18S <- tax_glom(ps_rel, taxrank = 'Class')
abund_18S <- psmelt(abund_18S)
abund_18S$Class <- as.character(abund_18S$Class)
abund_18S$Abundance <- abund_18S$Abundance * 100
abund_18S <- subset(abund_18S, grepl("Syndiniales|Dinophyceae|Bacillariophyta|Cryptophyceae|Mamiellophyceae|Spirotrichea|Filosa-Thecofilosea|Arthropoda|Mollusca|Annelida",abund_18S$Class)) # Subset to major groups (>1%)
abund_18S$Class <- factor(abund_18S$Class,levels=c("Bacillariophyta", "Dinophyceae", "Mamiellophyceae",  "Arthropoda","Syndiniales","Cryptophyceae", "Spirotrichea", "Mollusca", "Annelida","Filosa-Thecofilosea" )) # Order groups in the plot
p <- ggplot(abund_18S, aes(Temp, Abundance, group=Class))
p + geom_point(aes(fill = Class)) + 
geom_smooth(se=TRUE, linetype="solid", size=1, level=0.95, fill = "gray45", color = "black", alpha = 0.4, method = "loess")+ 
stat_cor(cor.coef.name = "rho", method = "spearman",label.x.npc = 0.35, label.y.npc=1.0, vjust=1, size = 3) + 
facet_wrap(~ Class, ncol =2, scales = "free_y") + theme_bw() + ylab("Relative Abundance (%)") + xlab("Temperature") + 
theme(legend.position = "none")+ guides(fill=guide_legend(ncol=2))+ removeGrid(x = TRUE, y = TRUE) 
ggsave(filename = "18S_abund.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 8, units = c("in", "cm", "mm"), dpi = 600) 

# Log2 fold analysis for 18S data using DESeq (Differential abundance of 18S ASVs between temperature conditions. Again, we subset the data based on our temperature variable, Temp2)
diagdds = phyloseq_to_deseq2(ps, ~ Temp2) # Convert phyloseq to DESeq2
gm_mean = function(x, na.rm=TRUE){        # Calculate geometric means prior to estimating size factors
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType = "parametric")

# Investigate DESeq table results
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.001 # Significance cutoff
sigtab = res[which(res$padj < alpha), ] # Table of significant log2 fold values
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Export DESeq results - Table S2
write.csv(sigtab, file="Table_S2.csv", row.names=T)

# Subset to the major 18S groups and plot log2 fold values - Figure S3
sigtab <- subset(sigtab, grepl("Syndiniales|Dinophyceae|Bacillariophyta|Cryptophyceae|Mamiellophyceae|Spirotrichea|Filosa-Thecofilosea|Arthropoda|Annelida|Mollusca",sigtab$Class))
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) mean(x)) # Order them based on mean log2 fold
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
sigtab$Class <- factor(sigtab$Class,levels=c("Bacillariophyta", "Dinophyceae", "Mamiellophyceae",  "Arthropoda","Syndiniales","Cryptophyceae", "Spirotrichea", "Mollusca", "Annelida","Filosa-Thecofilosea"))
p <- ggplot(sigtab, aes(x = Order, y=log2FoldChange, fill = Class))+ geom_point(size=4, shape=21, fill = "gray45") +
theme_bw() + geom_hline(yintercept=0, linetype = "dashed", size =1)+ theme(legend.position="bottom") +facet_grid(~Class, scale="free", space = "free") + 
theme(strip.background = element_blank(),strip.text.x = element_blank()) + removeGrid(x = TRUE, y = TRUE)+ ylab("Log2 fold abundance (cold - warm)") 
p + theme(axis.text.x=element_text(angle=90,vjust=0.2, hjust=1, size=12, color="black")) + theme(axis.text.y=element_text(size=12, color = "black", angle = 90)) + theme(axis.title.x=element_blank()) 
ggsave(filename = "18S_log2fold.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 8, height = 5, units = c("in", "cm", "mm"), dpi = 300)


### Repeat analyses for 16S

# Upload 16S count and taxonomy .qza files
table_16S <- read_qza("16S-table.qza")
count_tab_16S <- table_16S$data %>% as.data.frame() 
colnames(count_tab_16S)[1:3] <- c("011818-A-G10", "011818-B-G11", "011818-C-G12") # Some date labels were slightly off, must match 18S samples for network analysis
colnames(count_tab_16S)[57:59] <- c("081717-A-D7", "081717-B-D8", "081717-C-D9")
colnames(count_tab_16S)[89:91] <- c("120817-A-G4", "120817-B-G5", "120817-C-G6")
taxonomy_16S <- read_qza("16S-tax.qza")
tax_tab_16S <- taxonomy_16S$data %>% 
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence)
tax_tab_16S$Kingdom <- gsub("^.{0,5}", "", tax_tab_16S$Kingdom) # Clean up the taxonomy names for each level
tax_tab_16S$Phylum <- gsub("^.{0,5}", "", tax_tab_16S$Phylum)
tax_tab_16S$Class <- gsub("^.{0,5}", "", tax_tab_16S$Class)
tax_tab_16S$Order <- gsub("^.{0,5}", "", tax_tab_16S$Order)
tax_tab_16S$Family <- gsub("^.{0,5}", "", tax_tab_16S$Family)
tax_tab_16S$Genus <- gsub("^.{0,5}", "", tax_tab_16S$Genus)
tax_tab_16S$Species<- gsub("^.{0,5}", "", tax_tab_16S$Species)
tax_tab_16S[is.na(tax_tab_16S)] <- "Unassigned"

# Load 16S metadata file
sample_info_tab_16S <- read.table("Sampleinfo_16S.txt", header=TRUE, row.names=1, check.names=F, sep="\t") # Load 16S metadata

# Merge 16S ASV information - Table S1
merge_16S <- cbind(count_tab_16S, tax_tab_16S)
write.csv(merge_16S, file="Table_S1_16S.csv", row.names=T)

# Merge into phyloseq object and rename ASVs (bacteria ASVs labeled "bASVs" to distinguish from 18S)
ps_16S <- phyloseq(tax_table(as.matrix(tax_tab_16S)), otu_table(count_tab_16S, taxa_are_rows = T), sample_data(sample_info_tab_16S))
taxa_names(ps_16S) <- paste0("bASV", seq(ntaxa(ps_16S)))

# Remove taxa not assigned at the kingdom level
ps_16S <- subset_taxa(ps_16S, Kingdom!="Unassigned", Prune = T)

# Remove chloroplast and mitochondria reads 
psnew_16S = subset_taxa(ps_16S, Order != "Chloroplast" |is.na(Order))
psnew_16S = subset_taxa(psnew_16S, Family !="Mitochondria" |is.na(Family))

# Rarefy data to even sampling depth
ps_16S_rare <- rarefy_even_depth(psnew_16S, sample.size = min(sample_sums(psnew_16S)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Plot diversity (richness and Shannon) for 16S community - Figure 1D
temp2 = sample_data(ps_16S_rare)$Temp2
alpha_div_16S <- estimate_richness(ps_16S_rare, measures=c("Observed", "Shannon")) 
alpha_div_16S  = cbind(alpha_div_16S , temp2)
alpha_div_16S_melt = melt(alpha_div_16S, id.vars = "temp2", variable.name="value", value.name ="Diversity")
alpha_plot_16S <- ggplot(data=alpha_div_16S_melt, aes(x= temp2, y=Diversity, fill = temp2))
alpha_plot_16S + geom_boxplot(aes(fill=temp2), outlier.shape = NA, alpha = 0.6, size = 0.6, width = 0.6) + 
stat_boxplot(geom ='errorbar', width = 0.6) +
geom_point(aes(fill=temp2), size =2, shape = 21, colour = "black", position=position_jitterdodge(0.4)) +
scale_fill_manual(values= c("#0072B2", "#D55E00")) + theme_bw() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12)) + 
theme(axis.text.y = element_text(size=12))+ theme(legend.text = element_text(colour="black", size=12))+
theme(axis.title.x = element_blank())+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
facet_wrap(~value, scales = "free_y") +theme(strip.background = element_blank(),strip.text.x = element_blank())+
stat_compare_means(method= "t.test",comparisons = list(c("<23", "23-31")),label = "p.signif",  label.y.npc = 0.4)
ggsave(filename = "16S_alpha_01222021.eps", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 4, units = c("in", "cm", "mm"), dpi = 600)

# Rarefaction curve for 16S - Figure S1
rare_16S <- ggrare(ps_16S, step = 100, plot = TRUE, parallel = FALSE, se = FALSE)
rare_16S + theme(legend.position = "none")  +theme_bw()+theme(legend.position = "none") 
ggsave(filename = "16S_rare.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 4, height = 4, units = c("in", "cm", "mm"), dpi = 600) 

# Run NMDS ordination on 16S data
ps_16S_rare = filter_taxa(ps_16S_rare, function (x) {sum(x) > 1}, prune=TRUE) # Remove singletons
ps_16S_rel  <- transform_sample_counts(ps_16S_rare, function(x) x / sum(x)) # Normalize the data
ps_16S_ord  <- ordinate(ps_16S_rel, "NMDS", "bray")
p = plot_ordination(ps_16S_rel, ps_16S_ord, color="Temp2")
p$layers <- p$layers[-1]
p + theme_bw() + theme(text = element_text(size = 16)) + geom_point(aes(fill=Temp2),size = 4, alpha = 0.9, shape = 21, colour = "black") + scale_fill_manual(values= c("#0072B2", "#D55E00")) + 
theme(legend.text = element_text(colour="black", size=12)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "16S_nmds_01222021.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 3, units = c("in", "cm", "mm"), dpi = 600) 

# Find out the most abundant 16 groups (at class level) in terms of their relative proportion of reads across all samples
class_16S <- tax_glom(ps_16S_rare, taxrank = "Class")
class_16S <- transform_sample_counts(class_16S, function(x)100* x / sum(x))
OTU_16S <- otu_table(class_16S)
TAX_16S <- tax_table(class_16S)[,"Class"]
Average_16S <- as.data.frame(rowMeans(OTU_16S))
names(Average_16S) <- c("Mean")
Table_16S <- merge(TAX_16S, Average_16S, by=0, all=TRUE)
Table_16S$Row.names = NULL
Table_16S <- Table_16S[order(dplyr::desc(Table_16S$Mean)),]
head(Table_16S)

# Abundance vs. temp for top bacteria groups - Figure 2
abund_16S <- tax_glom(ps_16S_rel, taxrank = 'Class')
abund_16S <- psmelt(abund_16S)
abund_16S$Class <- as.character(abund_16S$Class)
abund_16S$Abundance <- abund_16S$Abundance * 100
abund_16S <- subset(abund_16S, grepl("Alphaproteobacteria|Gammaproteobacteria|Actinobacteria|Acidimicrobiia|Bacteroidia|Deltaproteobacteria|Verrucomicrobiae|Thermoplasmata|Oxyphotobacteria|Nitrososphaeria", abund_16S$Class))
abund_16S$Class <- factor(abund_16S$Class,levels=c("Alphaproteobacteria","Gammaproteobacteria", "Bacteroidia", "Acidimicrobiia", "Deltaproteobacteria", "Actinobacteria", "Verrucomicrobiae","Thermoplasmata", "Oxyphotobacteria", "Nitrososphaeria"))
p <- ggplot(abund_16S, aes(Temp, Abundance, group=Class))
p + geom_point(aes(fill = Class)) + 
geom_smooth(se=TRUE, linetype="solid", size=1, fill="gray45", level=0.95, color="black", alpha = 0.4, method = "loess") + 
stat_cor(cor.coef.name = "rho", method = "spearman",label.x.npc = 0.35, label.y.npc=1.0, vjust=1, size = 3)+
facet_wrap(~ Class, ncol =2, scales = "free_y") + theme_bw() + ylab("Relative Abundance (%)") + xlab("Temperature") + theme(legend.position = "none")+ guides(fill=guide_legend(ncol=2)) + removeGrid(x = TRUE, y = TRUE) 
ggsave(filename = "16S_abund.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 8, units = c("in", "cm", "mm"), dpi = 600) 


# Log2 fold plots for 16S - differential ASV abundance between temperatures (Temp2 variable in the metadata)
diagdds_16S = phyloseq_to_deseq2(psnew_16S, ~ Temp2)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds_16S), 1, gm_mean)
diagdds_16S = estimateSizeFactors(diagdds_16S, geoMeans = geoMeans)
diagdds_16S = DESeq(diagdds_16S, test="Wald", fitType = "parametric")

# Investigate DESeq table results
res_16S = results(diagdds_16S, cooksCutoff = FALSE)
alpha = 0.001
sigtab_16S = res[which(res$padj < alpha), ]
sigtab_16S = cbind(as(sigtab_16S, "data.frame"), as(tax_table(psnew_16S)[rownames(sigtab_16S), ], "matrix"))
head(sigtab_16S)

# Export DESeq results for - Table S2
write.csv(sigtab_16S, file="Table_S2_16S.csv", row.names=T)

# Plot log2 fold changes for top 16S groups - Figure S2
sigtab <- subset(sigtab, grepl("Alphaproteobacteria|Gammaproteobacteria|Actinobacteria|Acidimicrobiia|Bacteroidia|Deltaproteobacteria|Verrucomicrobiae|Thermoplasmata|Oxyphotobacteria|Nitrososphaeria",sigtab$Class))
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) mean(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
sigtab$Class <- factor(sigtab$Class,levels=c("Alphaproteobacteria","Gammaproteobacteria", "Bacteroidia", "Acidimicrobiia", "Deltaproteobacteria", "Actinobacteria", "Verrucomicrobiae","Thermoplasmata", "Oxyphotobacteria", "Nitrososphaeria"))
p <- ggplot(sigtab, aes(x = Order, y=log2FoldChange, fill = Class))+ geom_point(size=4, shape=21, fill = "gray45") +
theme_bw() + geom_hline(yintercept=0, linetype = "dashed", size =1)+ theme(legend.position="bottom") +facet_grid(~Class, scale="free", space = "free") + 
theme(strip.background = element_blank(),strip.text.x = element_blank()) + removeGrid(x = TRUE, y = TRUE)+ ylab("Log2 fold abundance (cold - warm)") 
p + theme(axis.text.x=element_text(angle=90,vjust=0.2, hjust=1, size=12, color="black")) + theme(axis.text.y=element_text(size=12, color = "black", angle = 90)) + theme(axis.title.x=element_blank()) 
ggsave(filename = "16S_log2fold.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 8, height = 5, units = c("in", "cm", "mm"), dpi = 300) 


### Prepare 16S and 18S for SpiecEasi network analysis (Kurtz et al. 2015)

# Subset 16S & 18S phyloseq objects based on warm (23-31) and cold (<23) temperatures
ps_18S_cold <- subset_samples(ps_rare, Temp2 == "<23")
ps_18S_warm <- subset_samples(ps_rare, Temp2 == "23-31")
ps_16S_cold <- subset_samples(ps_16S_rare, Temp2 == "<23")
ps_16S_warm <- subset_samples(ps_16S_rare, Temp2 == "23-31")

# Remove taxa not assigned at the class level
ps_18S_cold <- subset_taxa(ps_18S_cold, Class!="Unassigned", Prune = T)
ps_18S_warm <- subset_taxa(ps_18S_warm, Class!="Unassigned", Prune = T)
ps_16S_cold <- subset_taxa(ps_16S_cold, Class!="Unassigned", Prune = T)
ps_16S_warm <- subset_taxa(ps_16S_warm, Class!="Unassigned", Prune = T)

# Remove samples that are in the 18S but not in the 16S - sample IDs in count table have to match 
ps_18S_cold  = subset_samples(ps_18S_cold, sample_names(ps_18S_cold) != "031617-B-A2" & sample_names(ps_18S_cold) != "092017-C-E8")
ps_18S_warm = subset_samples(ps_18S_warm, sample_names(ps_18S_warm) != "031617-B-A2" & sample_names(ps_18S_warm) != "092017-C-E8")

# Trim sets to include ASVs with >5 read counts in 50% (16S) and 10% (18S) of samples - improve sparsity and avoid ambiguous interactions
filter_18S_cold = filter_taxa(ps_18S_cold, function(x) sum(x > 5) > (0.50*length(x)), TRUE)
filter_18S_warm = filter_taxa(ps_18S_warm, function(x) sum(x > 5) > (0.50*length(x)), TRUE)
filter_16S_cold = filter_taxa(ps_16S_cold, function(x) sum(x > 5) > (0.10*length(x)), TRUE)
filter_16S_warm = filter_taxa(ps_16S_warm, function(x) sum(x > 5) > (0.10*length(x)), TRUE)

# Merge phyloseq objects for the warm and cold period and convert OTU count table to data frame
all_warm <- merge_phyloseq(filter_16S_warm, filter_18S_warm)
warm_all.otu <- t(data.frame(phyloseq::otu_table(all_warm), check.names = FALSE))
all_cold <- merge_phyloseq(filter_16S_cold, filter_18S_cold)
cold_all.otu <- t(data.frame(phyloseq::otu_table(all_cold), check.names = FALSE))

# Convert individual warm/cold phyloseq objects for the 16S/18S for input into SpiecEasi
cold_16S.otu <- t(data.frame(phyloseq::otu_table(filter_16S_cold), check.names = FALSE))
cold_18S.otu <- t(data.frame(phyloseq::otu_table(filter_18S_cold), check.names = FALSE))
warm_16S.otu <- t(data.frame(phyloseq::otu_table(filter_16S_warm), check.names = FALSE))
warm_18S.otu <- t(data.frame(phyloseq::otu_table(filter_18S_warm), check.names = FALSE))

# Run SpiecEasi for cold set (16S + 18S)
se <- spiec.easi(list(cold_16S.otu, cold_18S.otu), method ='mb', nlambda = 20, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05))

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

tab.se <- convertSEToTable(se, sp.names=colnames(cold_all.otu)) 
tab.se.filtered <- tab.se %>% filter(is.finite(Weight))
cold_filt <- subset(tab.se.filtered, abs(Weight) > 0) # Remove values that are 0

# Export .csv file with ASV-ASV edges and weights needed to build network (Cytoscape) - note, weights are all significant and do not have a traditional p-value
# Export taxonomy key for the cold network to assign taxonomy to each node - add an "ASV" label to the taxonomy .txt file following export
write_csv(cold_filt, "Skio-cold-network_021121.csv")
taxALL <- tax_table(all_cold)
taxALL <- as.data.frame(taxALL, check.names=FALSE)
write.table(taxALL,"Skio-cold-TAX.txt",sep="\t", quote=FALSE)

# Extract cold ASV-ASV pairings and taxonomic annotation per ASV - Table S3
taxALL <- cbind(rownames(taxALL), data.frame(taxALL, row.names=NULL))
colnames(taxALL)[1] <- "Var1"
newtab_cold <- merge(cold_filt, taxALL[, c("Var1","Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL)[1] <- "Var2"
newtab_final_cold <- merge(newtab_cold, taxALL[, c("Var2", "Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var2")
write.csv(newtab_final_cold, "Table_S3.csv",row.names=F)

# Prepare chord diagram for top 20 cold edges at the class level
cold_edges <- newtab_final_cold[c("Class.x", "Class.y")] # Subset class level taxonomy for each node
for (i in 1:nrow(cold_edges  ))
{
  cold_edges  [i, ] = sort(cold_edges  [i, ])
}

cold_edges <- with(cold_edges, table(Class.x, Class.y)) # Count # of pairs for each edge type, avoiding reverse duplicates
cold_edges = as.data.frame(cold_edges) # Convert to data frame
cold_edges = cold_edges[with(cold_edges, order(-Freq)), ] # Order based on # of pairs found between class level groups
cold_top20 = head(cold_edges, 20) # Subset to top 20
cold_top20 <- cold_top20%>%  # Clean-up naming for chord diagram plot
  mutate(Class.y = as.character(Class.y)) %>% 
  mutate(Class.y = replace(Class.y, Class.y == 'MAST-6', 'MAST'))%>% 
  mutate(Class.y = as.character(Class.y)) %>% 
  mutate(Class.y = replace(Class.y, Class.y == 'Filosa-Thecofilosea', 'FilosaT'))

# Plot cold chord diagram of top 20 identified pairs (circlize) - Figure 4C
circos.clear()
circos.par("gap.degree" = 5, "track.height" = 0.5)
cold_mat = cold_top20
grid.col = c(Dinophyceae = "#1B9E77", Spirotrichea = "#D95F02",            # Set the colors for each class
             Mamiellophyceae = "#A6CEE3", Cryptophyceae = "#66A61E",
             Bacillariophyta = "#E6AB02", MAST = "#E7298A", Alphaproteobacteria = "#A6761D", Gammaproteobacteria = "#666666", FilosaT = "#CAB2D6", Bacteroidia = "#4575B4")
pdf("Skio-cold_network.pdf", width = 10, height = 10)
chordDiagram(cold_mat, annotationTrack="grid", preAllocateTracks = 1, grid.col=grid.col, transparency = 0.2, link.decreasing = TRUE,  big.gap = 100, small.gap = 1, self.link = 1, order= c("Dinophyceae", "Spirotrichea", "Mamiellophyceae",  "Cryptophyceae", "Bacillariophyta", "MAST","FilosaT", "Alphaproteobacteria", "Gammaproteobacteria", "Bacteroidia")) # Set the order of the track, transparency of links, and other network properties
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# Create a new column in ASV-ASV pairings table to denote positive and negative edges (or sign)
newtab_final_cold$Sign <- ifelse(newtab_final_cold$Weight > 0, "Pos", "Neg")
cold_edges_sign <- newtab_final_cold[c("Class.x", "Class.y", "Sign")]
for (i in 1:nrow(cold_edges_sign ))
{
  cold_edges_sign [i, ] = sort(cold_edges_sign [i, ])
}
cold_edges_new <- cold_edges_sign%>% 
  mutate(Class.y= as.character(Class.y)) %>% 
  mutate(Class.y = replace(Class.y, Class.y == 'MAST-6', 'MAST_6'))%>% 
  mutate(Class.y= as.character(Class.y)) %>% 
  mutate(Class.y = replace(Class.y, Class.y == 'Filosa-Thecofilosea', 'FilosaT'))

cold_edges_new$Edge <- paste(cold_edges_new$Class.x, cold_edges_new$Class.y,sep="-") # Add separator to distinguish edge type (e.g., Bacillariophyta-Dinophyceae)
df <- as.data.frame(table(cold_edges_new[,3:4]))
df = df[with(df, order(-Freq)), ] # Order based on frequency and subset to the top 20
df<- subset(df, grepl("Bacillariophyta-Bacillariophyta|Alphaproteobacteria-Gammaproteobacteria|Alphaproteobacteria-Alphaproteobacteria|Bacillariophyta-Dinophyceae|Bacillariophyta-Spirotrichea|Alphaproteobacteria-Bacteroidia|Gammaproteobacteria-Gammaproteobacteria|Dinophyceae-Spirotrichea|Bacillariophyta-Mamiellophyceae|Bacteroidia-Gammaproteobacteria|Alphaproteobacteria-Dinophyceae|Bacillariophyta-Cryptophyceae|Bacillariophyta-Gammaproteobacteria|Bacteroidia-Bacteroidia|Dinophyceae-Dinophyceae|Alphaproteobacteria-Bacillariophyta|Spirotrichea-Spirotrichea|Bacillariophyta-MAST_6|Bacillariophyta-FilosaT|Bacillariophyta-Bacteroidia",df$Edge))

# Plot top 20 interactions in cold network with positive and negative indicated - Figure 4D
p <- ggplot(df, aes(x = Edge, y = Freq))+ geom_col(aes(fill = Sign), width = 0.8, colour = "black", size = 1)
p$data$Edge <- factor(p$data$Edge, levels = c("Bacillariophyta-Bacillariophyta","Alphaproteobacteria-Alphaproteobacteria","Bacillariophyta-Dinophyceae","Alphaproteobacteria-Bacteroidia","Alphaproteobacteria-Gammaproteobacteria", "Bacillariophyta-Gammaproteobacteria", "Dinophyceae-Spirotrichea","Bacillariophyta-Spirotrichea", "Alphaproteobacteria-Bacillariophyta","Bacteroidia-Gammaproteobacteria","Bacillariophyta-Cryptophyceae","Gammaproteobacteria-Gammaproteobacteria","Spirotrichea-Spirotrichea", "Bacillariophyta-Mamiellophyceae","Bacillariophyta-FilosaT", "Dinophyceae-Dinophyceae", "Bacillariophyta-MAST_6","Bacillariophyta-Bacteroidia", "Bacteroidia-Bacteroidia","Alphaproteobacteria-Dinophyceae"))
labels <- c("Bacill-Bacill", "Alpha-Alpha", "Bacill-Dino","Alpha-Bacter", "Alpha-Gamma", "Bacill-Gamma", "Dino-Spiro", "Bacill-Spiro", "Alpha-Bacill", "Bacter-Gamma", "Bacill-Crypto", "Gamma-Gamma", "Spiro-Spiro", "Bacill-Mamiell", "Bacill-FilosaT", "Dino-Dino", "Bacill-MAST", "Bacill-Bacter", "Bacter-Bacter", "Alpha-Dino")
p + scale_fill_manual(values=c("#AF8DC3", "#7FBF7B"))  + theme_bw() +scale_x_discrete(labels= labels)+
theme(axis.title.x =element_blank())  + removeGrid(x = TRUE, y = TRUE) +  scale_y_continuous(expand = c(0, 0), breaks=seq(0,65,5), limits=c(0, 65))+
labs(y="Number of edges", size = 14)+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust = 1, size=14)) +  theme(axis.text.y=element_text(size=14))
ggsave(filename = "Network_edges_top20_cold.eps", last_plot(), device = "eps", path = NULL, scale = 1, width = 8, height = 5, dpi = 600)

# Run SpiecEasi mb method with warm samples (23-31 C)
se2 <- spiec.easi(list(warm_16S.otu, warm_18S.otu), method ='mb', nlambda = 20, lambda.min.ratio=1e-2, pulsar.params=list(rep.num=100, thresh=0.05))

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

tab.se <- convertSEToTable(se2,sp.names=colnames(warm_all.otu)) 
tab.se.filtered <- tab.se %>% filter(is.finite(Weight))
warm_filt <- subset(tab.se.filtered, abs(Weight) > 0) # Remove values that are 0

# Export .csv file with warm ASV pairings and associated taxonomy labels for nodes (Cytoscape)
write_csv(warm_filt, "Skio-warm-network_021121.csv")
taxALL <- tax_table(all_warm)
taxALL <- as.data.frame(taxALL, check.names=FALSE)
write.table(taxALL,"Skio-warm-TAX.txt",sep="\t", quote=FALSE)

# Extract warm ASV-ASV pairings with taxonomic annotation - Table S3
taxALL <- cbind(rownames(taxALL), data.frame(taxALL, row.names=NULL))
colnames(taxALL)[1] <- "Var1"
newtab <- merge(warm_filt, taxALL[, c("Var1","Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var1")
colnames(taxALL)[1] <- "Var2"
newtab_final_warm <- merge(newtab, taxALL[, c("Var2","Kingdom","Class", "Order", "Family", "Genus", "Species")], by="Var2")
write.csv(newtab_final_warm, "Table_S3_warm.csv", row.names=F)

# Prepare chord diagram of top 20 warm edges - Figure 4A
warm_edges <- newtab_final_warm[c("Class.x", "Class.y")]
for (i in 1:nrow(warm_edges  ))
{
  warm_edges  [i, ] = sort(warm_edges  [i, ])
}

warm_edges <- with(warm_edges, table(Class.x, Class.y)) # Count # of pairs for each edge type, avoiding reverse duplicates
warm_edges = as.data.frame(warm_edges)
warm_edges =  warm_edges[with(warm_edges, order(-Freq)), ]
warm_top20 = head(warm_edges, 20)
circos.clear()
circos.par("gap.degree" = 5, "track.height" = 0.5)
warm_mat = warm_top20
grid.col = c(Syndiniales = "#532688", Dinophyceae = "#1B9E77", Spirotrichea = "#D95F02", CONThreeP = "#FDBF6F",
             Cryptophyceae = "#66A61E", Mamiellophyceae = "#A6CEE3", Bacillariophyta = "#E6AB02",
             Alphaproteobacteria = "#A6761D", Gammaproteobacteria = "#666666",  Bacteroidia = "#4575B4", 
             Acidimicrobiia = "#D73027")
pdf("Skio-warm_network.pdf", width = 10, height = 10)
chordDiagram(warm_mat, annotationTrack="grid", preAllocateTracks = 1, grid.col=grid.col, transparency = 0.2, link.decreasing = TRUE,  big.gap = 100, small.gap = 1, self.link = 1, order= c("Syndiniales", "Dinophyceae", "Spirotrichea", "CONThreeP", "Cryptophyceae", "Mamiellophyceae", "Bacillariophyta", "Alphaproteobacteria", "Gammaproteobacteria","Bacteroidia", "Acidimicrobiia"))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# Create a new column in warm ASV-ASV pairings table to denote positive and negative edges (or sign)
newtab_final_warm$Sign <- ifelse(newtab_final_warm$Weight > 0, "Pos", "Neg") # Create new "sign" column
warm_edges_sign <- newtab_final_warm[c("Class.x", "Class.y", "Sign")]

for (i in 1:nrow(warm_edges_sign))
{
  warm_edges_sign [i, ] = sort(warm_edges_sign [i, ])
}
warm_edges_new$Edge <- paste(warm_edges_sign$Class.x, warm_edges_sign$Class.y,sep="-")
df_warm <- as.data.frame(table(warm_edges_new[,3:4]))
df_warm = df_warm[with(df_warm, order(-Freq)), ]
df_warm<- subset(df_warm, grepl("Bacillariophyta-Bacillariophyta|Bacillariophyta-Dinophyceae|Alphaproteobacteria-Alphaproteobacteria|Dinophyceae-Syndiniales|Dinophyceae-Spirotrichea|Dinophyceae-Dinophyceae|Bacillariophyta-Cryptophyceae|Alphaproteobacteria-Bacillariophyta|Bacillariophyta-Syndiniales|Cryptophyceae-Dinophyceae|Spirotrichea-Syndiniales|Cryptophyceae-Spirotrichea|Gammaproteobacteria-Gammaproteobacteria|Bacillariophyta-Mamiellophyceae|Bacillariophyta-Spirotrichea|Acidimicrobiia-Alphaproteobacteria|Alphaproteobacteria-Gammaproteobacteria|Alphaproteobacteria-Bacteroidia|Bacteroidia-Gammaproteobacteria|Bacillariophyta-CONThreeP",df_warm$Edge))

# Plot top 20 interactions (positive and negative) from warm network
p <- ggplot(df_warm, aes(x = Edge, y = Freq))+ geom_col(aes(fill = Sign), width = 0.8, colour = "black", size = 1)
p$data$Edge <- factor(p$data$Edge, levels = c("Bacillariophyta-Bacillariophyta","Bacillariophyta-Dinophyceae", "Alphaproteobacteria-Alphaproteobacteria","Bacillariophyta-Spirotrichea","Dinophyceae-Spirotrichea", "Dinophyceae-Syndiniales","Bacillariophyta-Syndiniales", "Dinophyceae-Dinophyceae", "Bacillariophyta-Cryptophyceae","Spirotrichea-Syndiniales","Alphaproteobacteria-Gammaproteobacteria","Alphaproteobacteria-Bacillariophyta", "Alphaproteobacteria-Bacteroidia", "Gammaproteobacteria-Gammaproteobacteria", "Cryptophyceae-Dinophyceae","Bacteroidia-Gammaproteobacteria","Acidimicrobiia-Alphaproteobacteria", "Bacillariophyta-Mamiellophyceae","Cryptophyceae-Spirotrichea","Bacillariophyta-CONThreeP"))
labels <- c("Bacill-Bacill", "Bacill-Dino","Alpha-Alpha","Bacill-Spiro", "Dino-Spiro", "Dino-Syn", "Bacill-Syn", "Dino-Dino", "Bacill-Crypto", "Spiro-Syn", "Alpha-Gamma", "Alpha-Bacill", "Alpha-Bacter", "Gamma-Gamma", "Crypto-Dino", "Bacter-Gamma", "Acidi-Alpha", "Bacill-Mamiell", "Crypto-Spiro", "Bacill-CON3")  
p + scale_fill_manual(values=c("#AF8DC3", "#7FBF7B"))  + theme_bw() +scale_x_discrete(labels= labels)+
theme(axis.title.x =element_blank())  + removeGrid(x = TRUE, y = TRUE) +  scale_y_continuous(expand = c(0, 0), breaks=seq(0,80,5), limits=c(0, 80))+
labs(y="Number of edges", size = 14)+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust = 1, size=14)) +  theme(axis.text.y=element_text(size=14))
ggsave(filename = "Network_edges_top20_warm.eps", last_plot(), device = "eps", path = NULL, scale = 1, width = 8, height = 5, dpi = 600)

# Plot # of edges per node for 16S/18S in warm vs. cold networks - Figure 5A
network_properties <- read.csv("network_properties.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") # Information from the network on degree and closeness centrality
p <- ggplot(network_properties, aes(x=Taxa, y=Degree, fill=Temp)) 
p + geom_boxplot(aes(fill=Temp), outlier.shape = NA, alpha = 0.6, size = 0.6, width = 0.6) + stat_boxplot(geom ='errorbar', width = 0.6) +
geom_point(aes(fill=Temp), size =2, shape = 21, colour = "black", position=position_jitterdodge()) + scale_fill_manual(values= c("#0072B2", "#D55E00")) + theme_bw() + ylab("Degree") +theme(axis.title.x=element_blank()) + 
scale_y_continuous(expand = c(0, 0), breaks=seq(0,25,5), limits=c(0, 25))+ theme(axis.text.y=element_text(size=12)) + theme(axis.title.y = element_text(size = 16))+ theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1, size=12, color="black")) +removeGrid(x = TRUE, y = TRUE) +
stat_compare_means(method= "wilcox.test",comparisons = list(c("16S_cold", "16S_warm"), c("16S_cold", "18S_cold"), c("18S_cold", "18S_warm"), c("16S_warm", "18S_warm")), label = "p.signif") # Compare means between loci (16S vs. 18S) and networks (18S cold vs. 18S warm)
ggsave(filename = "edges_warmvscold.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 6, height = 4, units = c("in", "cm", "mm"), dpi = 600) 

# Prepare for plot of total edge numbers between kingdom level (e.g., 18S and 16S) and temperatures - start with cold
cold_edges_king <- newtab_final_cold[c("Kingdom.x", "Kingdom.y", "Sign")] # Subset to kingdom and sign
for (i in 1:nrow(cold_edges_king ))
{
  cold_edges_king [i, ] = sort(cold_edges_king [i, ])
}
cold_edges_king$Kingdom.x <- replace(as.character(cold_edges_king$Kingdom.x), cold_edges_king$Kingdom.x== "Archaea", "Bacteria") # Rename the few Archaea as Bacteria
cold_edges_king$Edge <- paste(cold_edges_king$Kingdom.x, cold_edges_king$Kingdom.y,sep="-") # Create new column with separator (Bacteria-Eukaryota)
df <- as.data.frame(table(cold_edges_king[,3:4]))
df = df[with(df, order(-Freq)), ] # Data frame with the frequency (positive and negative) of each kingdowm interaction type for cold network

# Add the warm pairings - create new column to distinguish warm from cold
warm_edges_king <- newtab_final_warm[c("Kingdom.x", "Kingdom.y")]
for (i in 1:nrow(warm_edges_king ))
{
  warm_edges_king [i, ] = sort(warm_edges_king [i, ])
}
warm_edges_king$Kingdom.x <- replace(as.character(warm_edges_king$Kingdom.x), warm_edges_king$Kingdom.x== "Archaea", "Bacteria")
warm_edges_king$Edge <- paste(warm_edges_king$Kingdom.x, warm_edges_king$Kingdom.y,sep="-")
df_warm <- as.data.frame(table(warm_edges_king2[,3:4]))
df_warm = df_warm[with(df_warm, order(-Freq)), ]
df_total <- bind_rows(df, df_warm) # Bind cold and warm kingdom level edges
df_total$Temp <- c("cold", "cold", "cold", "cold", "cold", "cold", "warm", "warm", "warm", "warm", "warm", "warm") # Add new column "Temp"
df_total$Edge <- replace(as.character(df_total$Edge), df_total$Edge== "Eukaryota-Eukaryota", "18S-18S") # Change to 18S and 16S labels
df_total$Edge <- replace(as.character(df_total$Edge), df_total$Edge== "Bacteria-Eukaryota", "16S-18S")
df_total$Edge <- replace(as.character(df_total$Edge), df_total$Edge== "Bacteria-Bacteria", "16S-16S")

# Plot kingdom level temperature network pairings - Figure 5B
p <- ggplot(df_total, aes(x=Temp, y=Freq, fill=Sign))
p$data$Edge <- factor(p$data$Edge, levels = c("18S-18S", "16S-18S", "16S-16S"))
p + geom_col(aes(fill = Sign), width = 0.8, colour = "black", size = 1)+ scale_fill_manual(values=c("#AF8DC3", "#7FBF7B"))  + theme_bw()+
scale_y_continuous(name="# of Edges", limits=c(0, 700), expand = c(0, 0)) + theme(axis.title.x=element_blank()) + facet_wrap(~Edge)+
theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1, size=12, color="black"))+ theme(axis.text.y=element_text(size=12)) + theme(axis.title.y = element_text(size = 16))+ removeGrid(x = TRUE, y = TRUE) 
ggsave(filename = "pairings_warmvscold.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 6, height = 4, units = c("in", "cm", "mm"), dpi = 600) 

