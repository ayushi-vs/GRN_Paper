#D2 (GSE185275)
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
  library(harmony); library(presto)
})

setwd("~/Desktop/GRN_Paper_Deposit")
options(timeout = 600, future.globals.maxSize = 8 * 1024^3)
mem.maxVSize(vsize = 32000)
set.seed(42)

for (d in c("data/raw/GSE185275", "data/processed/GSE185275",
            "results/objects", "results/tables",
            "results/figures/GSE185275/QC",
            "results/figures/GSE185275/dimred"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


#from GEO

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185275/suppl"
tar_dest <- "data/raw/GSE185275/GSE185275_RAW.tar"
if (!file.exists(tar_dest))
  download.file(file.path(base_url, "GSE185275_RAW.tar"), tar_dest,
                mode = "wb", quiet = TRUE)
untar(tar_dest, exdir = "data/raw/GSE185275")

samples <- c(
  Coculture1 = "GSM5609927_Coculture1",
  Coculture2 = "GSM5609928_Coculture2",
  Coculture3 = "GSM5609929_Coculture3",
  Purified1  = "GSM5609930_Purified1",
  Purified2  = "GSM5609931_Purified2",
  Purified3  = "GSM5609932_Purified3"
)

for (s in names(samples)) {
  sd <- file.path("data/processed/GSE185275", s)
  dir.create(sd, showWarnings = FALSE, recursive = TRUE)
  for (kind in c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"))
    file.copy(file.path("data/raw/GSE185275",
                        paste0(samples[[s]], "_", kind)),
              file.path(sd, kind), overwrite = TRUE)
}


#persample seurat object build

seurat_list <- lapply(names(samples), function(s) {
  obj <- CreateSeuratObject(Read10X(file.path("data/processed/GSE185275", s)),
                            project = s,
                            min.cells = 3, min.features = 200)
  obj$sample_id <- s
  obj$system    <- ifelse(grepl("Coculture", s), "Coculture", "Purified")
  obj$replicate <- gsub("[A-Za-z]", "", s)
  obj
})
names(seurat_list) <- names(samples)

gse185275 <- merge(seurat_list[[1]], y = seurat_list[-1],
                   add.cell.ids = names(samples),
                   project = "GSE185275")
rm(seurat_list)
saveRDS(gse185275, "results/objects/GSE185275_raw.rds")


#QC

gse185275[["percent.mt"]] <- PercentageFeatureSet(gse185275, pattern = "^MT-")

ggsave("results/figures/GSE185275/QC/violin_pre.png",
       VlnPlot(gse185275,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "system", ncol = 3, pt.size = 0),
       width = 14, height = 5, dpi = 300, bg = "white")

gse185275 <- subset(gse185275,
                    subset = nFeature_RNA >= 500 & nFeature_RNA <= 8000 &
                      nCount_RNA   >= 1000 & percent.mt   <  15)

ggsave("results/figures/GSE185275/QC/violin_post.png",
       VlnPlot(gse185275,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "system", ncol = 3, pt.size = 0),
       width = 14, height = 5, dpi = 300, bg = "white")

gse185275 <- JoinLayers(gse185275)
saveRDS(gse185275, "results/objects/GSE185275_filtered.rds")


#Normalisation and pca

gse185275 <- gse185275 |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt")) |>
  RunPCA(npcs = 50, verbose = FALSE)

ggsave("results/figures/GSE185275/dimred/elbow.png",
       ElbowPlot(gse185275, ndims = 50),
       width = 8, height = 5, dpi = 300, bg = "white")

n_pcs <- 25


#UMAP

gse185275 <- RunUMAP(gse185275, dims = 1:n_pcs, reduction = "pca",
                     reduction.name = "umap_unintegrated", verbose = FALSE)

p_sample <- DimPlot(gse185275, reduction = "umap_unintegrated",
                    group.by = "sample_id", pt.size = 0.2) +
  ggtitle("Preintegration - by sample")
p_system <- DimPlot(gse185275, reduction = "umap_unintegrated",
                    group.by = "system", pt.size = 0.2) +
  ggtitle("Preintegration - by system")

ggsave("results/figures/GSE185275/dimred/preintegration.png",
       p_sample + p_system,
       width = 18, height = 7, dpi = 300, bg = "white")


#Harmony integration applied on sample_id

gse185275 <- RunHarmony(gse185275,
                        group.by.vars = "sample_id",
                        reduction = "pca",
                        reduction.save = "harmony",
                        verbose = FALSE)

gse185275 <- gse185275 |>
  RunUMAP(dims = 1:n_pcs, reduction = "harmony",
          reduction.name = "umap", verbose = FALSE) |>
  FindNeighbors(dims = 1:n_pcs, reduction = "harmony") |>
  FindClusters(resolution = 0.3)

p_after_sample <- DimPlot(gse185275, reduction = "umap",
                          group.by = "sample_id", pt.size = 0.2) +
  ggtitle("Post Harmony-by sample")
p_after_system <- DimPlot(gse185275, reduction = "umap",
                          group.by = "system", pt.size = 0.2) +
  ggtitle("Post Harmony -by system")
p_after_clust  <- DimPlot(gse185275, reduction = "umap",
                          label = TRUE, pt.size = 0.2) +
  ggtitle("Post Harmony - clusters")

ggsave("results/figures/GSE185275/dimred/postharmony.png",
       p_after_sample + p_after_system + p_after_clust,
       width = 24, height = 7, dpi = 300, bg = "white")

saveRDS(gse185275, "results/objects/GSE185275_clustered.rds")


#Celltype annotation

DefaultAssay(gse185275) <- "RNA"

markers_check <- list(
  Progenitor   = c("SOX2", "NES", "PAX6", "HES5", "VIM"),
  Cycling      = c("MKI67", "TOP2A"),
  Neuron       = c("DCX", "MAP2", "STMN2", "TUBB3", "SYT1"),
  Dopaminergic = c("TH", "NR4A2", "LMX1A", "FOXA2", "KCNJ6"),
  Proneural    = c("ASCL1", "NEUROG2", "NEUROD1")
)
markers_check <- lapply(markers_check, function(g) g[g %in% rownames(gse185275)])

ggsave("results/figures/GSE185275/dimred/marker_check.png",
       DotPlot(gse185275, features = markers_check,
               cols = c("lightgrey", "darkblue")) + RotatedAxis() +
         theme(axis.text.x = element_text(size = 8)),
       width = 18, height = 7, dpi = 300, bg = "white")

#Cluster to celltype map
cell_type_map <- c(
  "0"  = "Progenitor",         "1"  = "Neuron",
  "2"  = "DA_Neuron",          "3"  = "Progenitor",
  "4"  = "Transitional",       "5"  = "Cycling_Progenitor",
  "6"  = "Proneural",          "7"  = "Neuron",
  "8"  = "Neuron",             "9"  = "Neuron",
  "10" = "Neuron",             "11" = "Proneural",
  "12" = "Progenitor"
)
gse185275$cell_type <- unname(cell_type_map[as.character(Idents(gse185275))])
Idents(gse185275) <- "cell_type"

ggsave("results/figures/GSE185275/dimred/annotated.png",
       DimPlot(gse185275, reduction = "umap", group.by = "cell_type",
               label = TRUE, label.size = 4, pt.size = 0.2) +
         ggtitle("GSE185275_annotated cell types"),
       width = 11, height = 8, dpi = 300, bg = "white")

saveRDS(gse185275, "results/objects/GSE185275_annotated.rds")


#Differential expression: Neuron vs Progenitor

gse185275$lineage_stage <- case_when(
  gse185275$cell_type %in% c("Progenitor", "Cycling_Progenitor",
                             "Transitional") ~ "Progenitor",
  gse185275$cell_type %in% c("Neuron", "DA_Neuron",
                             "Proneural") ~ "Neuron",
  TRUE ~ "Other"
)
Idents(gse185275) <- "lineage_stage"

deg <- FindMarkers(gse185275,
                   ident.1 = "Neuron", ident.2 = "Progenitor",
                   test.use = "wilcox",
                   min.pct = 0.1, logfc.threshold = 0.25) |>
  rownames_to_column("gene") |>
  arrange(desc(avg_log2FC))
write.csv(deg,
          "results/tables/GSE185275_DEG_Neuron_vs_Progenitor.csv",
          row.names = FALSE)


