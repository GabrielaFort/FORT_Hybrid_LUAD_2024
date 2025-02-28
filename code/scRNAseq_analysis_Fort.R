### scRNAseq analysis for Fort et al HNF4α Manuscript
# Related to Figures 2B-F, S2B-F, S2H, and 3B-C 

# This R script was used to perform plotting and downstream analyses of scRNAseq data
# This script utilizes the preprocessed tumor cell seurat object generated by 'scRNAseq_celltype_identification_Fort.R'

# This code follows standard seurat workflows for Seurat V4 (https://satijalab.org/seurat/articles/pbmc3k_tutorial)


# Set working directory
setwd("/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/scRNAseq_analysis")


# Load (and install) required packages 
library(ggplot2)
#remotes::install_version("Seurat", version = "4.4.0")
library(Seurat)
#install.packages('ggpubr')
library(ggpubr)
#install.packages("R.utils")
library(R.utils)
#remotes::install_github('satijalab/seurat-wrappers@community-vignette')
library(SeuratWrappers)
#BiocManager::install("biomaRt")
library(biomaRt)
#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
library(CytoTRACE2)


### Plot UMAP to visualize tumor cell populations
### RELATED TO FIGURE 2B-C ###
# Read in tumor cell only seurat object generated with 'scRNAseq_celltype_identification_Fort.R'
KPvsKPH_tumor <- readRDS(file="tumor_cells_seuratobj.rds")

# Plot by seurat cluster
cluster_plot <- DimPlot(KPvsKPH_tumor, reduction = "umap",label=TRUE,pt.size=.4)
cluster_plot + NoLegend()

# Plot by genotype (KP vs KPH)
geno_plot <- DimPlot(KPvsKPH_tumor, reduction = "umap", group.by = "Geno_Type",label=FALSE,pt.size=0.4,
                     order=FALSE, cols=c("KP"="darkcyan","KPH"="maroon"))
geno_plot + NoLegend()


# Bar plot showing percent of each genotype in each cluster
### RELATED TO FIGURE 2D ###
levels(KPvsKPH_tumor) <- c('4','1','6','9','2','5','8','0','10','7','3')
barplot_table <- table(KPvsKPH_tumor$Geno_Type, Idents(KPvsKPH_tumor))
barplot_table <- as.data.frame(barplot_table)
barplot_table$Var1 <- as.character(barplot_table$Var1)

ggplot(barplot_table, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values=c("darkcyan","maroon")) +
  coord_cartesian(ylim=c(0.00,1.00))


###### Integrate ALRA imputation ######
KPvsKPH_tumor <- RunALRA(KPvsKPH_tumor) # Run imputation using ALRA model


# Use to toggle between assays for plotting individual genes
DefaultAssay(KPvsKPH_tumor) <- "RNA"
DefaultAssay(KPvsKPH_tumor) <- "alra"


### Plot genes of interest
# visualize original and imputed values for all genes
### RELATED TO FIGURES 2E, S2F, 3B ###

genes_oi = c('Hnf1a','Lgals4','Hnf4g','Pklr','Nkx2-1','Spock2','Zeb2','Cldn4','Vil1','Apln','Mcm2','Lyz1',
             'Hopx','Ager','Cxcl15','Sftpc')

for (gene in genes_oi) {
pdf(file = paste0("/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/scRNAseq_analysis/imputation_plots_010625/",gene,"_alra.pdf"), width = 5, height = 4.5)
print(FeaturePlot(KPvsKPH_tumor, features = gene, pt.size=0.3, 
            cols=c("peachpuff1","indianred","red4"), order=TRUE))
dev.off()
}


### Gene signature scoring and plotting 
### RELATED TO FIGURES 2E, 3B ###
gene_sets <-read.csv(file="gene_sets_for_scrnaseq_GF.csv")
gene_sets <- as.data.frame(gene_sets)

for (geneset in colnames(gene_sets)[1:13]) {
  KPvsKPH_tumor<-AddModuleScore(object = KPvsKPH_tumor, features = list(gene_sets[[geneset]]), name = geneset)
  print(FeaturePlot(KPvsKPH_tumor, features = paste0(geneset,"1"), pt.size=0.3, 
        cols=c("dodgerblue3","lightskyblue","wheat1","orange1","firebrick2")))
}



### Do DEG analysis (only on UNIMPUTED RNA counts)
# DEGs between clusters 3 and 4 were used as input for GSEA analysis found in Figures 2F and 3C

# Load human and mouse gene names from biomaRt
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")  

### Find differential genes between clusters 3 and 4 (C3 = KPH, C4 = KP)
cluster3v4 <- FindMarkers(object = KPvsKPH_tumor, ident.1 = 3, ident.2 = 4,min.pct = 0,only.pos = FALSE)
cluster3v4$Gene_Name = rownames(cluster3v4)
mouse_list = cluster3v4$Gene_Name

# Get corresponding attributes
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_list, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
colnames(genesV2)[1] = "Gene_Name"
# Delete duplicated mouse IDs
genesV2 = genesV2[!duplicated(genesV2$Gene_Name),]
cluster3v4 = merge(genesV2, cluster3v4, by='Gene_Name',all = TRUE)
write.csv(cluster3v4, file="/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/scRNAseq_analysis/cluster3vs4.csv",row.names = FALSE)


### Cytotrace analysis 
### RELATED TO FIGURE 2E

# Cytotrace analysis was run on scRNAseq dataset according to the following workflow:
# https://github.com/digitalcytometry/cytotrace2

# Run CytoTRACE2 main function
cytotrace2_result <- cytotrace2(KPvsKPH_tumor, is_seurat = TRUE, slot_type = "counts")


# The following code was adapted from CytoTRACE2 source code to plot relative CytoTRACE2 scores
rel_order_umap <- FeaturePlot(KPvsKPH_tumor, "CytoTRACE2_Relative") +
  scale_colour_gradientn(colours = (c("dodgerblue3","lightskyblue","wheat1",'khaki2',"orange1","firebrick2")),
                         na.value = "transparent",
                         limits=c(0,1),
                         breaks = seq(0,1, by = 0.2),
                         labels=c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"),
                         name = "Relative\norder \n" ,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  ggtitle("CytoTRACE 2") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),     
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20))) +
  theme(aspect.ratio = 1) 
rel_order_umap





