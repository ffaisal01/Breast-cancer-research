

library(dplyr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(readxl)

#TCGA, Firehose Legacy
filepath <- "C:/Users/ffais/OneDrive/Documents/rbif_120/week1/gdac.broadinstitute.org_BRCA.mRNAseq_Preprocess.Level_3.2016012800.0.0.tar.gz"

untar(filepath, exdir = "C:/Users/ffais/OneDrive/Documents/rbif_120/week1/extracted_data/")

rc_data <- read.table("C:/Users/ffais/OneDrive/Documents/rbif_120/week1/extracted_data/gdac.broadinstitute.org_BRCA.mRNAseq_Preprocess.Level_3.2016012800.0.0/BRCA.mRNAseq_raw_counts.txt", header = TRUE, sep = '\t')

head(rc_data)

newdata <- rc_data

num_characters_before = 10  # Specify the number of characters before "01" or "11"
cancer_pattern <- paste0("^(.{", num_characters_before, "})(.*.(01|06).*)")
control_pattern <- paste0("^(.{", num_characters_before, "})(.*.11.*)")


DEcondition <- sapply(colnames(newdata), function(name) {
  if (grepl(cancer_pattern, name)) {
    gsub(cancer_pattern, "Cancer", name)
  } else if (grepl(control_pattern, name)) {
    gsub(control_pattern, "Control", name)
  } else {
    name  # Keep the original name if no pattern matches
  }
})
DEcondition <- as.vector(DEcondition)
DEcondition <- DEcondition[-1]
DEcondition <- factor(DEcondition, levels = c("Control", "Cancer"), labels = c(0,1)) 
DEcondition


rownames(newdata) <- newdata$HYBRIDIZATION.R
newdata <- newdata[,-1] 

# test <- newdata
# colnames(test) <- DEcondition


###################DIFFERENTIAL ANALYSIS###############

metadata <- data.frame(
  sample = colnames(newdata),
  condition = DEcondition
)

dds <- DESeqDataSetFromMatrix(countData = newdata, colData = metadata, design = ~ condition)

dds <- DESeq(dds)
results <- results(dds)
head(results)

sig_genes <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1.5)
head(sig_genes)

#################################################################################33

# Volcano plot
results_df <- as.data.frame(results)
results_df$color <- ifelse(abs(results_df$log2FoldChange) > 1.5 & results_df$padj<0.05, "red", "black")
ggplot(as.data.frame(results_df), aes(x = log2FoldChange, y = -log10(pvalue), color=color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() + 
  theme_minimal() +
  xlim(c(-5, 5)) +
  ylim(c(0, 10))

pdata <- rownames(newdata)
gene.sym <- sapply(pdata, function(x) strsplit(x, "\\|")[[1]][1])

TCGA_rank <- results_df$log2FoldChange
TCGA_pval <- results_df$padj

ranked_output <- data.frame(gene.sym, TCGA_rank, TCGA_pval)
ranked_output <- na.omit(ranked_output)

TCGA_ranked_list <- data.frame(gene.sym, TCGA_rank)
# write.table(TCGA_ranked_list, "ranked_results", quote = FALSE, row.names = FALSE, sep = '\t')

sig.gene.names.TCGA <- ranked_output %>% 
  filter(TCGA_rank>1.5) %>% 
  filter(TCGA_pval<0.05) %>%
  select(gene.sym)
#sig.gene.names.TCGA <- as.vector(sig.gene.names.TCGA)
head(sig.gene.names.TCGA)


############################
# data from GSEA analysis of TCGA ranked list
# up_table <- read.table("C:/Users/ffais/OneDrive/Documents/rbif_120/week1/Book1.txt")
# col_names <- c("MSigDB", "SIZE",    "ES",   "NES",  "NOM p-val", "FDR", "FWER p-val", "RANK AT MAX")
# colnames(up_table) <- col_names
# up_table <- data.frame(up_table)
# 
# ggplot(up_table, aes(x=NES, y=MSigDB))+
#   geom_point(aes(color=as.numeric(FDR), size = as.numeric(SIZE)))+
#   scale_color_gradient(low = "blue", high = "red") + 
#   scale_size_continuous(range = c(2, 6)) + 
#   ggtitle("Dot Plot for upregulated pathways showing NES, FDR and size") +
#   theme_minimal() 
# +
#  facet_wrap(~ as.factor(Regulation), scales = "free")

##################################################


# # GSE158085
# # expression/count data
# # load counts table from GEO
# urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# path <- paste(urld, "acc=GSE158085", "file=GSE158085_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
# tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
# 
# # load gene annotations 
# apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
# annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
# rownames(annot) <- annot$GeneID
# 
# # sample selection
# gsms <- "111XXX000XXX000XXX000XXX"
# sml <- strsplit(gsms, split="")[[1]]
# 
# # filter out excluded samples (marked as "X")
# sel <- which(sml != "X")
# sml <- sml[sel]
# tbl <- tbl[ ,sel]
# 
# 
# # group membership for samples
# gs <- factor(sml)
# groups <- make.names(c("cancer","control"))
# levels(gs) <- groups
# sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
# 
# # pre-filter low count genes
# # keep genes with at least N counts > 10, where N = size of smallest group
# keep <- rowSums( tbl >= 10 ) >= min(table(gs))
# tbl <- tbl[keep, ]
# 
# head(tbl)
# GEO.exp.table <- data.frame(tbl)
# GEO.exp.table <- subset(GEO.exp.table, rownames(GEO.exp.table) %in% GEOdata$GeneID)
# subset.GEO.for.ID <- subset(GEOdata, GEOdata$GeneID %in% rownames(GEO.exp.table))
# rownames(GEO.exp.table) <- subset.GEO.for.ID$Symbol
# head(GEO.exp.table)


###################################
####pval and LogFC analysis
GEOdata <- read.table("C:/Users/ffais/OneDrive/Documents/rbif_120/week1/GSE158085.top.table.tsv", 
                      sep = "\t", 
                      header = TRUE, 
                      na.strings = c("", "NA", " "))

# "C:\Users\ffais\OneDrive\Documents\rbif_120\week1\GSE158085.top.table.tsv"


GEOdata <- as.data.frame(GEOdata)
str(GEOdata)
GEOdata <- na.omit(GEOdata)
colnames(GEOdata)

# GEOrank <- GEOdata$LogFC
GEOrank <- GEOdata$log2FoldChange

#GEOsymbol <- GEOdata$Gene.Symbol
GEOsymbol <- GEOdata$Symbol
# GEOsymbol <- gsub("///.*", "", GEOsymbol)

# GEOpval <- GEOdata$adj.P.Val
GEOpval <- GEOdata$padj

GEO.ranked.db <- data.frame(GEOsymbol, GEOrank, GEOpval)
GEO.ranked.db <- na.omit(GEO.ranked.db)

GEO.ranked.list <- data.frame(GEOsymbol, GEOrank)
GEO.ranked.list <- na.omit(GEO.ranked)
GEO.ranked.list <- GEO.ranked.list[order(GEO.ranked.list$GEOrank, decreasing = TRUE),]

# TABLE
# write.table(GEO.ranked.list, "GEO_ranked_new", sep = '\t', row.names = FALSE, quote = FALSE)


# "C:\Users\ffais\GEO_ranked_new"    bghfuyt7o89p[0-\
sig.gene.names.GEO <- GEO.ranked.db %>% 
  filter(GEOrank>1.5 & GEOpval<0.05) %>% 
  select(GEOsymbol)
#sig.gene.names.GEO <- as.vector(sig.gene.names.GEO)

head(sig.gene.names.GEO)
nrow(sig.gene.names.GEO)
####################################################################################


## PC highly expressed genes
##"C:\Users\ffais\OneDrive\Documents\rbif_120\week1\crc-22-0491-s02-ProstateCancer-final-list.xlsx"
file.path.pc <- ("C:/Users/ffais/OneDrive/Documents/rbif_120/week1/crc-22-0491-s02-ProstateCancer-final-list.xlsx")
PC_data <- read_excel(file.path.pc, , sheet="Gene-Level")

load("C:/Users/ffais/Downloads/surface_protein_lists.RData")
combined_surface_proteins_list <- c(CD.surface.proteins, CSPA_hiconf, CSPA_hiconf1, CSPA_putative, ion.goids, transmembrane.goids, GO0009897)

#######################################################################################


#Overlap of genes between TCGA and PC paper combined surface proteins

overlap.TCGA <- intersect(sig.gene.names.TCGA$gene.sym, combined_surface_proteins_list)
print(overlap.TCGA)

DE.upreg.TCGA <- ranked_output[ranked_output$gene.sym %in% overlap, ]


#Overlap of genes between GEO and Combined surface proteins

overlapGEO <- intersect(sig.gene.names.GEO$GEOsymbol, combined_surface_proteins_list)
print(overlapGEO)

DE.upreg.GEO <- GEO.ranked.db[GEO.ranked.db$GEOsymbol %in% overlapGEO, ]



#overlap of TCGA sig proteins with PC significant expressed surface proteins
sig.surface.proteins.overlap <- intersect(sig.gene.names.TCGA$gene.sym, PC_data$Gene)
sig.surface.proteins.overlap


# Overlap between sigly expressed TCGA and GEO surface proteins
intersect(overlap.TCGA, overlapGEO)

TCGA.GEO.upreg.genes <- as.vector(intersect(overlap.TCGA, overlapGEO))
##################################################################################


#subset newdata dataframe with the genes in the benchmark for highly expressed genes in humors
head(rownames(newdata))

newdata$gene.symbol <- gene.sym

# HKG from PLOS paper
#benchmark.genes <- c("ACTB", "B2M", "TUBB", "GAPDH", "TBP", "18S rRNA", "RPL13a", "G6PD", "SDHA", "HMBS", "HPRT1", "YWHAZ", "PPIA", "GUSB")

# Human cell selective ref genes
# https://housekeeping.unicamp.br/?download
H.ref.genes <- read.csv("C:/Users/ffais/OneDrive/Documents/rbif_120/week1/MostStable.csv")
H.ref.genes$Gene.names <- gsub(".*;", "", H.ref.genes$Ensembl.ID.Gene.name)

benchmark.df <- subset(newdata, newdata$gene.symbol %in% H.ref.genes$Gene.names)
rownames(benchmark.df) <- benchmark.df$gene.symbol
benchmark.df$gene.symbol <- NULL


# normalize and add Gene column based on rownames
benchmark.normalized <- log2(benchmark.df+1)
# Reshape data to long format for ggplot
benchmark_long <- melt(benchmark.normalized, variable.name = "Sample", value.name = "Expression")

benchmark_long$Gene <- rownames(benchmark.normalized)

# Create box plot for each gene
ggplot(benchmark_long, aes(x = Gene, y = Expression)) +
  geom_boxplot() +
  labs(title = "Expression Data for Housekeeping Genes",
       x = "Gene",
       y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################################3

#subset the expression data for TCGA sig genes identified
TCGA.exp.df <- subset(newdata, newdata$gene.symbol %in% sig.gene.names.TCGA$gene.sym)

rownames(TCGA.exp.df) <- TCGA.exp.df$gene.symbol
TCGA.exp.df$gene.symbol <- NULL


# normalize and add Gene column based on rownames
TCGA.normalized <- log2(TCGA.exp.df+1)
# Reshape data to long format for ggplot
TCGA_long <- melt(TCGA.normalized, variable.name = "Sample", value.name = "Expression")

TCGA_long$Gene <- rownames(TCGA.normalized)

# Create box plot for each gene
# ggplot(TCGA_long, aes(x = Gene, y = Expression)) +
#   geom_boxplot() +
#   labs(title = "Expression Data for TCGA Genes",
#        x = "Gene",
#        y = "Expression Level") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

#####################################################################3


#subset the expression data for Prostate cancer sig genes identified
PC.exp.df <- subset(newdata, newdata$gene.symbol %in% PC_data$Gene)

rownames(PC.exp.df) <- PC.exp.df$gene.symbol
PC.exp.df$gene.symbol <- NULL


# normalize and add Gene column based on rownames
PC.normalized <- log2(PC.exp.df+1)
# Reshape data to long format for ggplot
PC_long <- melt(PC.normalized, variable.name = "Sample", value.name = "Expression")

PC_long$Gene <- rownames(PC.normalized)

# Create box plot for each gene
ggplot(PC_long, aes(x = Gene, y = Expression)) +
  geom_boxplot() +
  labs(title = "Expression Data for Prostate Cancer Genes",
       x = "Gene",
       y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################


# GEO.hkg.df <- subset(GEO.exp.table, rownames(GEO.exp.table) %in% H.ref.genes$Gene.names)
# #rownames(GEO.hkg.df) <- GEO.hkg.df$gene.symbol
# #GEO.hkg.df$gene.symbol <- NULL
# # normalize and add Gene column based on rownames
# GEO.normalized <- log2(GEO.hkg.df+1)
# # Reshape data to long format for ggplot
# GEO_long <- melt(GEO.normalized, variable.name = "Sample", value.name = "Expression")
# 
# GEO_long$Gene <- rownames(GEO.normalized)
# 
# # Create box plot for each gene
# ggplot(GEO_long, aes(x = Gene, y = Expression)) +
#   geom_boxplot() +
#   labs(title = "Expression Data for GEO Genes",
#        x = "Gene",
#        y = "Expression Level") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

##############################################################

# combined_df <- merge(benchmark_long, TCGA_long, by="Sample", all.x = TRUE)


GEO.surface.markers <- intersect(sig.gene.names.GEO$GEOsymbol, combined_surface_proteins_list)

intersect(GEO.surface.markers, sig.gene.names.TCGA$gene.sym)

###################################################################

TCGA_enrichment <- read_excel("C:/Users/ffais/OneDrive/Documents/rbif_120/week2/TCGA_enrichment_NES.xlsx")
GEO_enrichment <- read_excel("C:/Users/ffais/OneDrive/Documents/rbif_120/week1/GEO_enrichment_NES_only.xlsx")
enrichment <- merge(TCGA_enrichment, GEO_enrichment, by="NAME", all=TRUE)
enrichment[is.na(enrichment)] <- 0
# enrichment <- as.data.frame((enrichment))

enrichment$NAME <- gsub("^HALLMARK_", "", enrichment$NAME)
rownames(enrichment) <- enrichment$NAME
enrichment$NAME <- NULL
sub_enrichment <- enrichment[, c("TCGA-NES", "GEO-NES")]


# Determine overlapping pathways 
overlapping_upregulated <- rownames(sub_enrichment)[sub_enrichment$`TCGA-NES` > 0 & sub_enrichment$`GEO-NES` > 0] 
overlapping_downregulated <- rownames(sub_enrichment)[sub_enrichment$`TCGA-NES` < 0 & sub_enrichment$`GEO-NES` < 0] 
#order non-overlapping pathways
non_overlapping <- setdiff(rownames(sub_enrichment), c(overlapping_upregulated, overlapping_downregulated)) 
non_overlapping <- non_overlapping[order(rowSums(sub_enrichment[non_overlapping, ]), decreasing = TRUE)]

# Order pathways: upregulated at the top, downregulated at the bottom 
ordered_pathways <- c(overlapping_upregulated, non_overlapping, overlapping_downregulated) 
sub_enrichment <- sub_enrichment[ordered_pathways, ]                                                                  
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

pheatmap(
  sub_enrichment,
  cluster_rows = FALSE,    # Cluster rows if needed
  cluster_cols = TRUE,    # Cluster columns if needed
  color = colorRampPalette(c("blue", "white", "red"))(50), # Color scale (low to high NES)
  main = "Heatmap of NES Scores"
)


#################################################################
# Venn Diagram for significantly upregulated genes from GEO and TCGA datasets

#install.packages("VennDiagram")
library(VennDiagram)

# gene lists from two datasets
GEO.sig.upreg <- sig.gene.names.GEO$GEOsymbol
TCGA.sig.upreg <- sig.gene.names.TCGA$gene.sym

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = list(Dataset1 = GEO.sig.upreg, Dataset2 = TCGA.sig.upreg),
  category.names = c("GEO Upreg genes", "TCGA upreg genes"),
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = c("blue", "orange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.col = c("blue", "orange"),
  main = "Venn Diagram of Upregulated Genes for GEO and TCGA datasets"
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Venn diagram for surface markers overlap between TCGA an dGEO datasets
GEO.surface.markers <- overlapGEO
TCGA.surface.markers <- overlap.TCGA

Venn.surface.markers <- venn.diagram(
  x = list(Dataset1 = GEO.surface.markers, Dataset2 = TCGA.surface.markers),
  category.names = c("GEO Upreg surface markers", "TCGA upreg surface markers"),
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.col = c("blue", "red"),
  main = "Venn Diagram of Upregulated Genes for GEO and TCGA datasets"
)

grid.draw(Venn.surface.markers)

Overlap.surface.markers <- intersect(GEO.surface.markers, TCGA.surface.markers)

#############################################################################

# if (!requireNamespace("Seurat", quietly = TRUE)) {
#   install.packages("Seurat")
# }

# Load the Seurat package
# library(Seurat)
# 
# # Read the files
# counts <- read.table("C:/Users/ffais/OneDrive/Documents/rbif_120/scRNA_analysis/Allcells_raw_count_out/BrCa_Atlas_Count_out/matrix.mtx.gz", sep = "\t", header = FALSE)
# features <- read.table("C:/Users/ffais/OneDrive/Documents/rbif_120/scRNA_analysis/Allcells_raw_count_out/BrCa_Atlas_Count_out/features.tsv.gz", sep = "\t", header = FALSE)
# barcodes <- read.table("C:/Users/ffais/OneDrive/Documents/rbif_120/scRNA_analysis/Allcells_raw_count_out/BrCa_Atlas_Count_out/barcodes.tsv.gz", sep = "\t", header = FALSE)
# 
# # Create a Seurat object
# seurat_object <- CreateSeuratObject(counts = counts, features = features, barcodes = barcodes)

#################################
# install.packages("Seurat")
# install.packages("Matrix")
# install.packages("data.table")
library(Seurat)
library(Matrix)
library(data.table)
library(R.utils)
library(ggplot2)

# Define the file paths
matrix_file <- "C:/Users/ffais/OneDrive/Documents/rbif_120/scRNA_analysis/Allcells_raw_count_out/BrCa_Atlas_Count_out/matrix.mtx.gz"
features_file <- "C:/Users/ffais/OneDrive/Documents/rbif_120/scRNA_analysis/Allcells_raw_count_out/BrCa_Atlas_Count_out/features.tsv.gz"
barcodes_file <- "C:/Users/ffais/OneDrive/Documents/rbif_120/scRNA_analysis/Allcells_raw_count_out/BrCa_Atlas_Count_out/barcodes.tsv.gz"

# Read the matrix file
expression_matrix <- readMM(file = gzfile(matrix_file))

# Read the features and barcodes
features <- fread(features_file, header = FALSE)
barcodes <- fread(barcodes_file, header = FALSE)

# Assign row and column names
rownames(expression_matrix) <- features$V1
colnames(expression_matrix) <- barcodes$V1

dim(expression_matrix)

# Convert to sparse matrix
counts_sparse <- as(expression_matrix, "sparseMatrix")
seurat_object <- CreateSeuratObject(counts = counts_sparse)
seurat_object <- NormalizeData(seurat_object, assay = "RNA")
                               
                               
rownames(seurat_object@assays$RNA@layers$counts) <- features$V1
colnames(seurat_object@assays$RNA@layers$counts) <- barcodes$V1
seurat_object@assays$RNA@layers$counts[1:10,1:10]

head(x = rownames(x = seurat_object))
head(x = colnames(x = seurat_object))
str(seurat_object)

GetAssayData(object = seurat_object, layer = 'counts')[1:10, 1:10]
head(x = seurat_object@meta.data[c('nCount_RNA', 'nFeature_RNA')])

VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(object = seurat_object, features = c("nFeature_RNA"), group.by = c('orig.ident'))

# signifacntly upregulated cell surface genes from TCCGA and GEO analysis above
TCGA.GEO.upreg.genes <- c("ABCC12","ATP1A3","ATP2A3","CACNA1H","CACNA1I","CACNG4","CACNG6","CDH15","CEACAM5", "CLIC3","CLIC6","COL5A1","COX7B2","DPEP1","ERBB2","FABP6","FGFR4","GRIN2D","HMMR","HPX","IL9R","KCNE4","KCNG3","KCNH2","KCNJ3","LCN12","LRRC15", "NCAM2","NUP210","P2RX2","PCSK9","PODXL2","RET","SEZ6L2","SLC12A5","SLC17A9", "SLC34A3","SLC39A4","SLC5A8","SLC8A2","TFR2","TNFRSF18","UNC5A")

upreg.genes.found <- which(toupper(rownames(seurat_object@assays$RNA@layers$counts)) %in% TCGA.GEO.upreg.genes)

n.expressed.genes <- Matrix::colSums(seurat_object@assays$RNA@layers$counts[upreg.genes.found, ] > 0)
seurat <- AddMetaData(object = seurat_object, metadata = n.expressed.genes, col.name = "n.expressed.genes")
VlnPlot(object = seurat, features = c("n.expressed.genes"), ncol = 1)













#########subsetting for a GENE of INTEREST###########
# Define the gene of interest
gene.of.interest <- "ABCC12"

# Subset the Seurat object by the gene of interest
# gene_subset <- seurat_object@assays$RNA$counts[1:10,1:10]

# Plot for each patient using VlnPlot
VlnPlot(object = gene_subset, features = c("nFeature_RNA"), group.by = 'orig.ident')
##############################################################33
# #OR#
# 
# set.seed(42) # For reproducibility 
# subset_indices <- sample(ncol(expression_matrix), 10000)
# expression_matrix_subset <- expression_matrix[, subset_indices]


# Create Seurat object
# seurat_object <- CreateSeuratObject(counts = expression_matrix)


# meta file
meta.sc <- read.csv("C:/Users/ffais/OneDrive/Documents/rbif_120/scRNA_analysis/Whole_miniatlas_meta.csv", row.names = 1)
meta.sc <- meta.sc[-1, ]
head(meta.sc)
# 
seurat_object <- AddMetaData(seurat_object, metadata = meta.sc)

#########VIOLIN PLOT by CELL-TYPE##################333

# Define the gene of interest and identity
genes.of.interest <- as.vector(TCGA.GEO.upreg.genes)
Idents(seurat_object) <- seurat_object@meta.data$celltype_major


VlnPlot(object = seurat_object, features = c("ABCC12")) + ggtitle(paste("Violin plot for", gene, "by Cell Type"))

# fixed_labels <- unique(seurat_object@meta.data$celltype_major)
                       
# Ensure the RNA assay contains the gene of interest
for (gene in genes.of.interest){
if(gene %in% rownames(seurat_object@assays$RNA@layers$counts)) {
  #using only gene data with non-zero values
  gene_data <- FetchData(seurat_object, vars = gene)
  non_zero_cells <- rownames(gene_data)[gene_data[, gene] > 0]
  gene_subset <- subset(seurat_object, cells = non_zero_cells)

  # Set identities using the celltype_major column
  Idents(gene_subset) <- gene_subset@meta.data$celltype_major
  
  print(paste("Generating plot for", gene))
  # Plot for each patient using VlnPlot, ensuring cell type is on the x-axis
  plot <- VlnPlot(object = gene_subset, features = c(gene)) + ggtitle(paste("Violin plot for", gene, "by Cell Type"))
  print(plot)
} else {
  message(paste("Gene", gene, "not found in the RNA assay."))
}
}

# taking out fixed labels : caused legend misalignment


patient.id <- "CID3586"
patient.data <- subset(seurat_object, subset = orig.ident == patient.id)
#Idents(seurat_object) <- seurat_object@meta.data$celltype_major
Idents(patient.data) <- patient.data@meta.data$celltype_major
head(patient.data)
# patient.data <- NormalizeData(patient.data)
# patient.data <- FindVariableFeatures(patient.data)
# patient.data <- ScaleData(patient.data)
#layers <- slotNames(seurat_object@assays$RNA) 
#print(layers)


# plot
VlnPlot(object = patient.data, features = c(genes.of.interest[1]))

patient.data@assays$RNA@layers





for (gene in TCGA.GEO.upreg.genes) {
  tryCatch({
    plot <- VlnPlot(patient.data, features = gene, assay = "RNA", slot = "counts") + ggtitle(paste("Violin Plot of", gene))
    print(plot)
  }, error = function(e) {
    message("Error in plotting gene: ", gene, "\n", e)
  })
}



patient.data <- FindVariableFeatures(patient.data)
patient.data <- ScaleData(patient.data)
patient.data <- RunPCA(patient.data)
patient.data <- FindNeighbors(patient.data, dims = 1:10)
patient.data <- FindClusters(patient.data)
patient.data <- RunUMAP(patient.data, dims = 1:10)
DimPlot(patient.data, reduction = "umap", group.by = "celltype_minor")



# normalization and analysis
# seurat_object <- NormalizeData(seurat_object)
# seurat_object <- FindVariableFeatures(seurat_object)
# seurat_object <- ScaleData(seurat_object)
# seurat_object <- RunPCA(seurat_object)
# seurat_object <- FindNeighbors(seurat_object)
# seurat_object <- FindClusters(seurat_object)
# seurat_object <- RunUMAP(seurat_object)

for (gene in TCGA.GEO.upreg.genes) { VlnPlot(seurat_object, features = gene) + ggtitle(paste("Violin Plot of", gene)) }





