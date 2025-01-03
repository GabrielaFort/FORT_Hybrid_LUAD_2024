### scRNAseq analysis for Fort et al HNF4α Manuscript
# Related to Figures 2B-F, S2B-F, S2H, and 3B-C 

# This R script was used to download and perform preprocessing steps on the KP/KPH scRNAseq dataset
# including quality control, normalization, feature selection, scaling, clustering and dimensional reduction 

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

#####################################################################
##################### Read in MEX files  ############################
#####################################################################

# 1 : KP Pre Immune Depletion 
data_dir <- "MEX/20463X1"
KPpre <- Read10X(data_dir)
KPpre <- CreateSeuratObject(counts = KPpre,project="KPpre")
sample<-"KPpre"
KPpre@meta.data$Tumor_Type<-sample
geno<-"KP"                                
KPpre@meta.data$Geno_Type<-geno


# 2 : KP Post Immune Depletion 
data_dir <- "MEX/20463X2"
KPpost <- Read10X(data_dir)
KPpost <- CreateSeuratObject(counts = KPpost,project="KPpost")
sample<-"KPpost"
KPpost@meta.data$Tumor_Type<-sample
geno<-"KP"                                
KPpost@meta.data$Geno_Type<-geno


# 3 : KPH Pre Immune Depletion 
data_dir <- "MEX/20463X3"
KPHpre <- Read10X(data_dir)
KPHpre <- CreateSeuratObject(counts = KPHpre,project="KPHpre")
sample<-"KPHpre"
KPHpre@meta.data$Tumor_Type<-sample
geno<-"KPH"                                
KPHpre@meta.data$Geno_Type<-geno


# 4 : KPH Post Immune Depletion 
data_dir <- "MEX/20463X4"
KPHpost <- Read10X(data_dir)
KPHpost <- CreateSeuratObject(counts = KPHpost,project="KPHpost")
sample<-"KPHpost"
KPHpost@meta.data$Tumor_Type<-sample
geno<-"KPH"                                
KPHpost@meta.data$Geno_Type<-geno


#####################################################################
############## Merge into single Seurat object ######################
#####################################################################

### Combine into one seurat object ### 
KPvsKPH<-merge(KPpre,c(KPpost, KPHpre, KPHpost))


# Number of cells captured per sample
table(KPvsKPH@meta.data$Tumor_Type)
#KPHpost  KPHpre  KPpost   KPpre 
##7867    6558    7284    7658 

#####################################################################
######################### QC Filtering ##############################
#####################################################################

KPvsKPH[["percent.mt"]] <- PercentageFeatureSet(KPvsKPH, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(KPvsKPH, features = c("percent.mt", "nFeature_RNA","nCount_RNA"), ncol=3, group.by="Tumor_Type", pt.size=0.2)
# Choosing <12 percent mitochondrial content as cutoff

FeatureScatter(KPvsKPH, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Tumor_Type")

FeatureScatter(KPvsKPH, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Tumor_Type")
# Setting 500 as lower cutoff and 7500 as upper cutoff

# Subsetting quality cells 
KPvsKPH<- subset(KPvsKPH, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 12)

#Cell counts post QC filtering
table(KPvsKPH@meta.data$Tumor_Type) 
#KPHpost  KPHpre  KPpost   KPpre 
#7064    5771    6562    6783 


#####################################################################
################### Normalize, Scale, Cluster #######################
#####################################################################

# Normalize data
KPvsKPH <- NormalizeData(KPvsKPH, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
KPvsKPH <- FindVariableFeatures(KPvsKPH, selection.method = "vst", nfeatures = 2000)

# Scale Data
all.genes <- rownames(KPvsKPH)
KPvsKPH <- ScaleData(KPvsKPH, features = all.genes)

# Perform linear dimension reduction/PCA
KPvsKPH <- RunPCA(KPvsKPH, features = VariableFeatures(object = KPvsKPH))

# Use elbow plot to visualize dimensions that encompass most of variance
ElbowPlot(KPvsKPH, ndims=50)
# shows ~30 dims captures most of the variance

# Choosing 30 PCs for dimensionality
KPvsKPH <- FindNeighbors(KPvsKPH, dims = 1:30)
KPvsKPH <- FindClusters(KPvsKPH)

# Run Non-Linear Dimension reduction (UMAP)
KPvsKPH <- RunUMAP(KPvsKPH, dims = 1:30)

# Save combined seurat file
saveRDS(KPvsKPH, file="all_cells_seuratobj.rds")
#KPvsKPH <- readRDS(file="all_cells_seuratobj.rds")


