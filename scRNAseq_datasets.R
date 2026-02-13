#ANALYSES FROM VISVADER DATASET
library(Seurat) 
library(dplyr) 
library(cowplot) 
library(ggrepel)


####################################
### Figure 8A
####################################
#LOAD DATA FROM TUMOR CELLS FROM TRIPLE NEGATIVE CASES
setwd("/projects/eshenderov-hpc/Shivang/projects/cutntag_RN_VW/Shivang/Shivang/Shivang/scRNAseq_datasets/scRNAseq_Triple_Negative")
load("TripleNeg_TumorCells.RData")
#str(Tumor)

#fetch genes of interest in the order of interest
#Tumor_dat <- FetchData(Tumor, vars = c("sample", "NCOR2", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-G", "B2M"))
Tumor_dat <- FetchData(Tumor, vars = c("sample", "NCOR2", "HLA-A", "HLA-B", "B2M"))

#any(is.na(Tumor_dat))
#prepare matrix
mat <- Tumor_dat
Metadata <- mat[,c("sample")] #Metadata for heatmap annotation, must be done here to keep the order of the rows.

mat$sample <- NULL # eliminate column sample before scale data
#any(is.na(mat))
#dim(mat)
#View(mat)


#scale data
mat <- scale(mat) # scale and center columns
mat <- t(scale(t(mat))) # scale and center rows
mat <- t(mat)

mat_df <- as.data.frame(t(mat))

cor_HLAA <- cor.test(mat_df$NCOR2, mat_df$`HLA-A`, method = "spearman")
cor_HLAA$estimate
cor_HLAA$p.value


#head(metadata)
library(ComplexHeatmap)
library(circlize)

annotation_colors <- list("Sample"= c("TN1" = "aquamarine", 
                                      "TN2" = "darkgreen",
                                      "TN3" = "brown2",
                                      "TN4" = "bisque1",
                                      "TN5" = "lightgreen",
                                      "TN6" = "blueviolet",
                                      "TN7" = "orange",
                                      "TN8" = "yellow3"))


ha = HeatmapAnnotation(Sample = Metadata, col = annotation_colors)

p <- Heatmap(mat,
             top_annotation = ha, 
             cluster_rows = F,
             cluster_columns = T,
             show_column_names = FALSE,
             show_row_names = T,
             column_title = NULL,
             name = "Z-score",
             col=colorRamp2(c(-2, 0, 2), c("dodgerblue", "white", "red3")),border = T)

p

pdf("Figure 8A.pdf", width = 8, height = 6)
draw(p)
dev.off()

tiff("Figure 8A.tiff",
     width = 8, height = 6,
     units = "in",
     res = 600,
     compression = "lzw")

draw(p)
dev.off()


#SAVE AS: PNG 781/267
dev.off()
#save 838 269
save(ha, file = "Annotation.RData")
save(mat, file = "MAT.RData")

Totalcells <- table(Tumor$sample)
#  TN1   TN2   TN3   TN4   TN5   TN6   TN7   TN8 
# 1188 13538    33   116  5135  1115   132   213 = 21420 cells
View(Totalcells)
Totalcells <- as.data.frame(Totalcells)
sum(Totalcells$Freq)

####################################
### Figure 8B
####################################

setwd("/projects/eshenderov-hpc/Shivang/projects/cutntag_RN_VW/Shivang/Shivang/Shivang/scRNAseq_datasets/scRNAseq_HER2positive")
load(file = "HER2posTumorCells.RData")
Tumor_dat <- FetchData(Tumor, vars = c("sample", "NCOR2", "HLA-A", "HLA-B", "B2M"))
View(Tumor_dat)
unique(Tumor_dat$sample)

mat <- Tumor_dat
Metadata <- mat[,c("sample")] #Metadata for heatmap annotation, must be done here to keep the order of the rows.
mat$sample <- NULL # eliminate column sample before scale data

mat <- scale(mat) # scale and center columns
mat <- t(scale(t(mat))) # scale and center rows
mat <- t(mat)

#head(metadata)
library(ComplexHeatmap)
library(circlize)

annotation_colors <- list("Sample"= c("HER1" = "aquamarine", 
                                      "HER2" = "darkgreen",
                                      "HER3" = "brown2",
                                      "HER4" = "bisque1",
                                      "HER5" = "lightgreen",
                                      "HER6" = "blueviolet",
                                      "HER7" = "orange",
                                      "HER8" = "yellow3"))


ha = HeatmapAnnotation(Sample = Metadata, col = annotation_colors)

p <- Heatmap(mat,
             top_annotation = ha, 
             cluster_rows = F,
             cluster_columns = T,
             show_column_names = FALSE,
             show_row_names = T,
             column_title = NULL,
             name = "Z-score",
             col=colorRamp2(c(-2, 0, 2), c("dodgerblue", "white", "red3")),border = T)

p

pdf("Figure 8B.pdf", width = 8, height = 6)
draw(p)
dev.off()

tiff("Figure 8B.tiff",
     width = 8, height = 6,
     units = "in",
     res = 600,
     compression = "lzw")

draw(p)
dev.off()

#SAVE AS: PNG 781/267
dim(mat)




####################################
### Figure 8C
####################################
#LOAD DATA FROM TUMOR CELLS FROM TRIPLE NEGATIVE CASES
setwd("/projects/eshenderov-hpc/Shivang/projects/cutntag_RN_VW/Shivang/Shivang/Shivang/scRNAseq_datasets/scRNAseq_Triple_Negative")
load("TripleNeg_TumorCells.RData")
#str(Tumor)

#fetch genes of interest in the order of interest
#Tumor_dat <- FetchData(Tumor, vars = c("sample", "NCOR2", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-G", "B2M"))
Tumor_dat <- FetchData(Tumor, vars = c("sample", "NCOR2", "HLA-A", "HLA-B", "B2M"))

#any(is.na(Tumor_dat))
#prepare matrix
mat <- Tumor_dat
Metadata <- mat[,c("sample")] #Metadata for heatmap annotation, must be done here to keep the order of the rows.

mat$sample <- NULL # eliminate column sample before scale data
#any(is.na(mat))
#dim(mat)
#View(mat)


#scale data
mat <- scale(mat) # scale and center columns
mat <- t(scale(t(mat))) # scale and center rows
mat <- t(mat)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

# -----------------------------
# 1. Compute Spearman correlation + p-values
# -----------------------------
mat <- as.matrix(t(mat))

genes <- colnames(mat)
n <- length(genes)

cor_mat <- matrix(NA, n, n, dimnames = list(genes, genes))
p_mat   <- matrix(NA, n, n, dimnames = list(genes, genes))

for (i in 1:n) {
  for (j in 1:n) {
    ct <- suppressWarnings(cor.test(mat[, i], mat[, j], method = "spearman"))
    cor_mat[i, j] <- ct$estimate
    p_mat[i, j]   <- ct$p.value
  }
}

# Format values for printing (2 decimals)
cor_lab <- matrix(sprintf("%.2f", cor_mat), n, n)
p_lab   <- matrix(sprintf("%.2f", p_mat), n, n)

# -----------------------------
# 2. Color scales
# -----------------------------

# Correlation: blue → white → red
cor_col_fun <- colorRamp2(c(-1, 0, 1), c("#313695", "white", "#A50026"))

# P-values: cyan → yellow
p_col_fun <- colorRamp2(c(0, 1), c("#00FFFF", "#FFFF00"))

# -----------------------------
# 3. Heatmap annotations
# -----------------------------
ha_top <- HeatmapAnnotation(
  Gene = colnames(mat),
  show_annotation_name = TRUE,
  show_legend = FALSE  # ← hide annotation legend
)

ha_left <- rowAnnotation(
  Gene = rownames(cor_mat),  # = gene names
  show_annotation_name = TRUE,
  show_legend = FALSE        # ← hide annotation legend
)

# -----------------------------
# 4. Plot: Correlation heatmap
# -----------------------------
ht1 <- Heatmap(
  cor_mat, name = "Spearman\nCorrelation",
  col = cor_col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(cor_lab[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
  },
  top_annotation = ha_top,
  left_annotation = ha_left
)

pdf("Figure 8C.pdf", width = 8, height = 6)
draw(ht1)
dev.off()

tiff("Figure 8C.tiff",
     width = 8, height = 6,
     units = "in",
     res = 600,
     compression = "lzw")

draw(ht1)
dev.off()


####################################
### Figure 8D
####################################
# -----------------------------
# 5. Plot: P-value heatmap
# -----------------------------
ht2 <- Heatmap(
  p_mat, name = "P-value",
  col = p_col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(p_lab[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
  },
  top_annotation = ha_top,
  left_annotation = ha_left
)

pdf("Figure 8D.pdf", width = 8, height = 6)
draw(ht2)
dev.off()

tiff("Figure 8D.tiff",
     width = 8, height = 6,
     units = "in",
     res = 600,
     compression = "lzw")

draw(ht2)
dev.off()

# -----------------------------
# 6. Draw together
# -----------------------------
ComplexHeatmap::draw(ht1, heatmap_legend_side = "right")
ComplexHeatmap::draw(ht2, heatmap_legend_side = "right")

