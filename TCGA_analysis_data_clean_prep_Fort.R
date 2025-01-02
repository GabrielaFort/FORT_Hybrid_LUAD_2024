### TCGA analysis for Fort et al HNF4α Manuscript
# Related to Fig 1A, S1A, S2A, S2J

# This R script was used to download and clean clinical and gene expression data
# from TCGA and to run single sample GSEA for downstream plotting/analysis 

####################### TCGA analysis for KPH manuscript #######################
######################### Gabriela Fort 08.03.2024 #############################

# Load required packages
library(ggplot2)
library(data.table)
library(dplyr)
library(survival)
#install.packages('survminer')
library(survminer)
#install.packages('tab')
library(tab)
#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
#BiocManager::install("GSVA")
library(GSVA)

# Set working directory 
setwd('/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis')

# Download data from cBioPortal query of only patients of interest:
# NKX2-1 z score >-0.2
# HNF4a z score > 0.8 and not mucinous = HNF4A-high (n=47)
# HNF4a z score < 0.8 = HNF4A-low (n=339)

# Download clinical data from cBioPortal for patients of interest
# Query link: 
# https://www.cbioportal.org/study/clinicalData?id=luad_tcga_pan_can_atlas_2018
# only patients that have been curated by E.Snyder and G.Fort 
tcga_data = fread('/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/luad_tcga_pan_can_atlas_2018_clinical_data.tsv',
                       header = T)
tcga_data = as.data.frame(tcga_data)

# Replace spaces in column names with underscores
colnames(tcga_data) = gsub(" ", "_", colnames(tcga_data))
# Remove parentheses in Survival Time column
colnames(tcga_data)[35] = "Overall_Survival_Months"

# Read in table containing curated HNF4A-high and HNF4A-low patient IDs 
samples_oi = read.csv('/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/patient_samples_oi.csv', header = T)

# Make new column with HNF4A-high or low designation 
tcga_data = tcga_data %>% dplyr::mutate(HNF_status = dplyr::case_when(Patient_ID %in% samples_oi$HNF_high ~ "HNF4a-High",
                                                                      Patient_ID %in% samples_oi$HNF_low ~ "HNF4a-Low"))
table(tcga_data$HNF_status)
# n=47 hnf-high and n=339 hnf-low

# Clean up overall survival column notation 
tcga_data = tcga_data %>% dplyr::mutate(OS_status = dplyr::case_when(Overall_Survival_Status == "0:LIVING" ~ 0,
                                                                     Overall_Survival_Status == "1:DECEASED" ~ 1))

# Save as a csv file
write.csv(tcga_data,"/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/TCGA_KPH_survival_HNFstatus.csv", row.names = FALSE)

#tcga_data = read.csv("/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/TCGA_KPH_survival_HNFstatus.csv")


### Download gene expression data from TCGA curated cohort

# Retrieve sample IDs of interest
samples_oi = tcga_data$Sample_ID

# Generate query to download all TPM counts from TCGA-LUAD patients within tcga_data table
query = GDCquery(project = 'TCGA-LUAD',
                 data.category = 'Transcriptome Profiling',
                 experimental.strategy = 'RNA-Seq',
                 workflow.type = 'STAR - Counts',
                 access = 'open',
                 barcode = samples_oi)

# Download RNA data associated with samples of interest
GDCdownload(query)

# Create expression matrix with TPM values 
luad_tcga_data = GDCprepare(query, summarizedExperiment = T)
expression_matrix = assay(luad_tcga_data, 'tpm_unstrand')

# Filter expression matrix for unexpressed or lowly expressed genes
expression_matrix = expression_matrix[rowSums(expression_matrix[])>10,]

# Remove version information from ensembl IDs
rownames(expression_matrix) = gsub("\\..*","",rownames(expression_matrix))

# Convert ensembl ID to human gene names and retain only protein-coding genes
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mrna_attributes <- getBM(attributes=c("external_gene_name",
                                      "ensembl_gene_id",
                                      "gene_biotype"),
                         filters = c("ensembl_gene_id"),
                         values = rownames(expression_matrix),
                         mart = mart)
mrna_attributes <- mrna_attributes[which(mrna_attributes$gene_biotype == "protein_coding"),][1:2]
expression_matrix <- expression_matrix[which(rownames(expression_matrix) %in% mrna_attributes$ensembl_gene_id),]
expression_matrix = as.data.frame(expression_matrix)

# Include gene names in expression matrix
expression_matrix$ensembl_gene_id = rownames(expression_matrix)
expression_matrix = merge(expression_matrix, mrna_attributes, by='ensembl_gene_id',all = F)

# Remove duplicate gene IDs and merge
expression_matrix = expression_matrix[!duplicated(expression_matrix$external_gene_name), ]
rownames(expression_matrix) = expression_matrix$external_gene_name
expression_matrix = subset(expression_matrix, select = -c(ensembl_gene_id, external_gene_name))

# Fix sample names
colnames(expression_matrix) = substr(colnames(expression_matrix), 1, 15)
expression_matrix <- expression_matrix[, !duplicated(colnames(expression_matrix))]

# Extract transformed log2(TPM+1) counts from matrix for genes of interest
# Calculate log2(TPM+1)
logTPM <- function(x) {return(log2(x+1))}
expression_matrix_transformed = t(expression_matrix %>% dplyr::mutate_if(is.numeric, logTPM))
expression_matrix_transformed = as.data.frame(expression_matrix_transformed)

genes_oi = c('IDH1','PKLR',"HNF4A","NKX2-1","ALDOB","SORD","EPS8L3","LGALS4","HNF1A")

expression_matrix_subset = expression_matrix_transformed[genes_oi]
expression_matrix_subset$Sample_ID = rownames(expression_matrix_subset)

# Merge clinical data with gene expression data into a single dataframe
tcga_data = merge(tcga_data, expression_matrix_subset, by='Sample_ID')

# Write to CSV file
write.csv(tcga_data,"/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/TCGA_KPH_survival_HNFstatus.csv", row.names = FALSE)
#tcga_data = read.csv("/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/TCGA_KPH_survival_HNFstatus.csv")


############## Run ssGSEA ###############
# Want to run ssGSEA on gene sets identified from bulk RNAseq data
# Add these ssGSEA scores to existing dataframe with clinical and gene expression information 

# Table with gene sets derived from in house bulk RNA-seq data
ssgsea_sets = read.csv("/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/ssgsea_hnf4a_genesets.csv",na.strings=c(""," ","NA"))
ssgsea_list = as.list(ssgsea_sets)
ssgsea_list = lapply(ssgsea_list, function(x) x[!is.na(x)])

expression_matrix = as.matrix(expression_matrix)

# Run ssGSEA on gene expression matrix using custom gene sets
ssgsea_par = ssgseaParam(expression_matrix,ssgsea_sets)
ssgsea_res = gsva(ssgsea_par)

ssgsea_res = as.data.frame(t(ssgsea_res))
ssgsea_res$Sample_ID = rownames(ssgsea_res)

# Merge results with existing dataframe 
tcga_data = merge(tcga_data, ssgsea_res, by='Sample_ID')

# Write to CSV file
write.csv(tcga_data,"/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/TCGA_KPH_survival_HNFstatus.csv", row.names = FALSE)
#tcga_data = read.csv("/uufs/chpc.utah.edu/common/home/snydere-group2/Gabby/KPH_TCGA_analysis/TCGA_KPH_survival_HNFstatus.csv")










