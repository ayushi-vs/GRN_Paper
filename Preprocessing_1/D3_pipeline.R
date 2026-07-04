#D3 (GSE138121)
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
  library(harmony); library(presto); library(Matrix)
})

setwd("~/Desktop/GRN_Paper_Deposit")
options(timeout = 600, future.globals.maxSize = 8 * 1024^3)
mem.maxVSize(vsize = 32000)
set.seed(42)

for (d in c("data/raw/GSE138121",
            "results/objects", "results/tables",
            "results/figures/GSE138121/QC",
            "results/figures/GSE138121/dimred",
            "results/figures/GSE138121/annotation"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


#Download

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138121/suppl"
dl_dir   <- "data/raw/GSE138121"

batch_files <- c(
  BatchA  = "GSE138121_BatchA_ddSEQ_UMIcounts.csv.gz",
  BatchB1 = "GSE138121_BatchB1_ddSEQ_UMIcounts.csv.gz",
  BatchB2 = "GSE138121_BatchB2_ddSEQ_UMIcounts.csv.gz",
  BatchC  = "GSE138121_BatchC_ddSEQ_UMIcounts.csv.gz",
  BatchD  = "GSE138121_BatchD_ddSEQ_UMIcounts.csv.gz",
  BatchE  = "GSE138121_BatchE_10X_UMIcounts.csv.gz",
  BatchF  = "GSE138121_BatchF_10X_UMIcounts.csv.gz"
)

for (f in c(batch_files, "GSE138121_CellBarcodes.csv.gz")) {
  dest <- file.path(dl_dir, f)
  if (!file.exists(dest))
    download.file(file.path(base_url, f), dest, mode = "wb", quiet = TRUE)
}


#Perbatch loader: controls only

load_batch_ctr <- function(path, batch_name) {
  mat_df <- read.csv(path, check.names = FALSE)
  genes  <- gsub("_", "-", mat_df[[1]])
  cn     <- colnames(mat_df)[-1]
  parts  <- strsplit(cn, "_")
  condition <- sapply(parts, `[`, 2)
  condition <- ifelse(condition == "C9TR", "CTR", condition)
  keepcols <- which(condition == "CTR")
  if (length(keepcols) == 0) { rm(mat_df); gc(); return(NULL) }
  mat <- as.matrix(mat_df[, keepcols + 1, drop = FALSE])
  rm(mat_df)
  if (any(duplicated(genes)))
    mat <- rowsum(mat, group = genes)
  else
    rownames(mat) <- genes
  mat <- as(mat, "dgCMatrix")
  
  kp <- parts[keepcols]
  f6 <- sapply(kp, `[`, 6)
  meta <- data.frame(
    cell_id   = colnames(mat),
    condition = "CTR",
    cell_line = sapply(kp, `[`, 3),
    platform  = sapply(kp, `[`, 5),
    timepoint = ifelse(grepl("^d[0-9]+$", f6), f6, "d18"),
    batch     = batch_name,
    row.names = colnames(mat),
    stringsAsFactors = FALSE
  )
  
  obj <- CreateSeuratObject(counts = mat, meta.data = meta,
                            project = batch_name,
                            min.cells = 3, min.features = 200)
  rm(mat, meta); gc()
  obj
}

seurat_list <- list()
for (b in names(batch_files)) {
  obj <- load_batch_ctr(file.path(dl_dir, batch_files[[b]]), b)
  if (!is.null(obj)) seurat_list[[b]] <- obj
}

gse138121 <- merge(seurat_list[[1]], y = seurat_list[-1],
                   add.cell.ids = names(seurat_list),
                   project = "GSE138121")
rm(seurat_list)
saveRDS(gse138121, "results/objects/GSE138121_raw.rds")


#QC

gse138121[["percent.mt"]] <- PercentageFeatureSet(gse138121, pattern = "^MT-")
gse138121$timepoint <- factor(gse138121$timepoint,
                              levels = c("d00", "d06", "d12", "d18"))

ggsave("results/figures/GSE138121/QC/violin_by_timepoint.png",
       VlnPlot(gse138121,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "timepoint", ncol = 3, pt.size = 0),
       width = 14, height = 5, dpi = 300, bg = "white")

ggsave("results/figures/GSE138121/QC/violin_by_platform.png",
       VlnPlot(gse138121,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "platform", ncol = 3, pt.size = 0),
       width = 14, height = 5, dpi = 300, bg = "white")

# ddSEQ has ~0 percent.mt by design, so the <15% cap only applies to 10X (TENEX)
keep_mt <- (gse138121$platform == "DDSEQ") | (gse138121$percent.mt < 15)
gse138121 <- subset(gse138121,
                    subset = nFeature_RNA >= 500 & nFeature_RNA <= 8000 &
                      nCount_RNA   >= 1000)
gse138121 <- gse138121[, keep_mt[colnames(gse138121)]]

gse138121 <- JoinLayers(gse138121)
saveRDS(gse138121, "results/objects/GSE138121_filtered.rds")


#Normalisation & PCA

gse138121 <- gse138121 |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData(vars.to.regress = c("nCount_RNA")) |>   
  RunPCA(npcs = 50, verbose = FALSE)

ggsave("results/figures/GSE138121/dimred/elbow.png",
       ElbowPlot(gse138121, ndims = 50),
       width = 8, height = 5, dpi = 300, bg = "white")

n_pcs <- 20


#UMAP

gse138121 <- RunUMAP(gse138121, dims = 1:n_pcs, reduction = "pca",
                     reduction.name = "umap_unintegrated", verbose = FALSE)

ggsave("results/figures/GSE138121/dimred/preintegration.png",
       DimPlot(gse138121, reduction = "umap_unintegrated",
               group.by = "platform", pt.size = 0.2) +
         ggtitle("Preintegration-by platform") +
         DimPlot(gse138121, reduction = "umap_unintegrated",
                 group.by = "batch", pt.size = 0.2) +
         ggtitle("Preintegration-by batch") +
         DimPlot(gse138121, reduction = "umap_unintegrated",
                 group.by = "timepoint", pt.size = 0.2) +
         ggtitle("Preintegration-by timepoint"),
       width = 24, height = 7, dpi = 300, bg = "white")


#Harmony

gse138121 <- RunHarmony(gse138121,
                        group.by.vars = c("platform", "batch"),
                        theta = c(4, 2),
                        reduction = "pca",
                        reduction.save = "harmony",
                        max_iter = 20,
                        verbose = FALSE)

gse138121 <- gse138121 |>
  RunUMAP(dims = 1:n_pcs, reduction = "harmony",
          reduction.name = "umap", verbose = FALSE) |>
  FindNeighbors(dims = 1:n_pcs, reduction = "harmony") |>
  FindClusters(resolution = 0.3)

ggsave("results/figures/GSE138121/dimred/postharmony.png",
       DimPlot(gse138121, reduction = "umap",
               group.by = "timepoint", pt.size = 0.2) +
         ggtitle("Postharmony: by timepoint") +
         DimPlot(gse138121, reduction = "umap",
                 group.by = "platform", pt.size = 0.2) +
         ggtitle("Post harmony: by platform") +
         DimPlot(gse138121, reduction = "umap",
                 group.by = "batch", pt.size = 0.2) +
         ggtitle("Post harmony: by batch"),
       width = 24, height = 7, dpi = 300, bg = "white")

saveRDS(gse138121, "results/objects/GSE138121_integrated.rds")


#Celltype annotation 

DefaultAssay(gse138121) <- "RNA"

traj_markers <- c("POU5F1", "NANOG", "LIN28A",
                  "SOX2", "NES", "PAX6", "HES5", "VIM",
                  "OLIG2", "NKX6-1",
                  "DCX", "STMN2", "MAP2", "TUBB3", "SYT1",
                  "ISL1", "MNX1", "LHX3", "CHAT", "SLC18A3")
traj_markers <- traj_markers[traj_markers %in% rownames(gse138121)]

ggsave("results/figures/GSE138121/dimred/trajectory_markers.png",
       DotPlot(gse138121, features = traj_markers,
               group.by = "timepoint") + RotatedAxis(),
       width = 12, height = 4, dpi = 300, bg = "white")

Idents(gse138121) <- "seurat_clusters"
cell_type_map <- c(
  "0"  = "Motor_Neuron",     "1"  = "d18_Progenitor",
  "2"  = "Early_Progenitor", "3"  = "Neural_Progenitor",
  "4"  = "iPSC",             "5"  = "Motor_Neuron",
  "6"  = "Neuron",           "7"  = "MN_Progenitor",
  "8"  = "Neuron",           "9"  = "Neuron",
  "10" = "Neuron"
)
gse138121$cell_type <- unname(cell_type_map[as.character(Idents(gse138121))])
Idents(gse138121) <- "cell_type"

ggsave("results/figures/GSE138121/annotation/annotated_umap.png",
       DimPlot(gse138121, reduction = "umap", group.by = "cell_type",
               label = TRUE, label.size = 4, pt.size = 0.2) +
         ggtitle("GSE138121: annotated cell types"),
       width = 11, height = 8, dpi = 300, bg = "white")

saveRDS(gse138121, "results/objects/GSE138121_annotated.rds")


#Differential expression: d18 vs d00

exclude <- grepl("^MT-",    rownames(gse138121)) |
  grepl("^RP[SL]", rownames(gse138121)) |
  grepl("^MRP[SL]", rownames(gse138121))
gse_clean <- subset(gse138121, features = rownames(gse138121)[!exclude])

Idents(gse_clean) <- "timepoint"
deg <- FindMarkers(gse_clean,
                   ident.1 = "d18", ident.2 = "d00",
                   test.use = "wilcox",
                   min.pct = 0.1, logfc.threshold = 0.25) |>
  rownames_to_column("gene") |>
  arrange(desc(avg_log2FC))

rm(gse_clean); gc()

write.csv(deg,
          "results/tables/GSE138121_DEG_d18_vs_d00.csv",
          row.names = FALSE)

