## Microbial relationships in the Skidaway River Estuary

This repository has QIIME 2 compatible files (.qza), metadata tables, network files (Cytoscape) and R code associated with the following manuscript:

**Estuarine microbial networks and relationships vary between environmentally distinct communities**<br/>
*Sean R. Anderson and Elizabeth L. Harvey,  (2022)*<br/>

The manuscript is currently under review.

### Project workflow

#### Summary

Microbial relationships remain enigmatic in the ocean, despite their importance for food web dynamics and biogeochemistry. Understanding how microbial associations vary over space and time, and with respect to changing environmental conditions, offers important insight into how associations may change in future ecosystems. We explored microbial relationships (prokaryotes and protists) and networks between two contrasting periods in the Skidaway River Estuary (GA, USA), performing co-occurrence network analysis with 16S and 18S metabarcoding data. Raw sequence data are available in NCBI SRA under BioProject ID's PRJNA575563 (18S) and PRJNA680039 (16S). 

#### Bioinformatics and data analysis

* 16S and 18S sequences processed separately in QIIME 2
* ASVs inferred with paired-end DADA2
* 16S taxonomy assigned with SILVA (v.138.1); 18S taxonomy assigned with PR2 (v.4.12)
* QIIME 2 taxonomy, ASV count table, and phylogenetic tree files uploaded to R using [qiime2R](https://github.com/jbisanz/qiime2R)
* Hierarchical clustering (Ward's method) of beta diversity
* Relative abundance plots of major order level groups
* Spearman correlation between group abundance and environmental factors
* [SPIEC-EASI](https://github.com/zdk123/SpiecEasi) for network analysis of distinct clusters 
* Comparison of network properties (e.g., degree, centrality, and edge density) and microbial relationships between clusters
