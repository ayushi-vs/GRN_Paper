#D1 

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(tidyverse)
  library(patchwork); library(ggrepel); library(scDblFinder)
})


setwd("~/Desktop/GRN_Paper_Deposit")
options(timeout = 600, future.globals.maxSize = 8 * 1024^3)
mem.maxVSize(vsize = 32000)
set.seed(42)

for (d in c("data/raw/GSE208625", "data/processed",
            "results/objects", "results/tables",
            "results/figures/QC", "results/figures/dimred",
            "results/figures/annotation", "results/figures/DEG"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


#Download from GEO

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208625/suppl"
samples  <- c("N_I", "N_N")   # iPSC, NSC

files <- unlist(lapply(samples, function(s)
  paste0("GSE208625_", s, c("_filtered_barcodes.tsv.gz",
                            "_filtered_features.tsv.gz",
                            "_filtered_matrix.mtx.gz"))))

for (f in files) {
  dest <- file.path("data/raw/GSE208625", f)
  if (!file.exists(dest) || file.info(dest)$size < 1e4)
    download.file(file.path(base_url, f), dest, mode = "wb", quiet = TRUE)
}

for (s in samples) {
  sd <- file.path("data/processed", paste0("GSE208625_", s))
  dir.create(sd, showWarnings = FALSE, recursive = TRUE)
  for (kind in c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"))
    file.copy(file.path("data/raw/GSE208625",
                        paste0("GSE208625_", s, "_filtered_", kind)),
              file.path(sd, kind), overwrite = TRUE)
}


#Build and merge Seurat objects

ipsc <- CreateSeuratObject(Read10X("data/processed/GSE208625_N_I"),
                           project = "GSE208625_iPSC",
                           min.cells = 3, min.features = 200)
ipsc$condition <- "iPSC"; ipsc$dataset <- "GSE208625"

nsc <- CreateSeuratObject(Read10X("data/processed/GSE208625_N_N"),
                          project = "GSE208625_NSC",
                          min.cells = 3, min.features = 200)
nsc$condition <- "NSC"; nsc$dataset <- "GSE208625"

gse208625 <- merge(ipsc, y = nsc, add.cell.ids = c("iPSC", "NSC"),
                   project = "GSE208625")
saveRDS(gse208625, "results/objects/GSE208625_raw.rds")


#QC

gse208625[["percent.mt"]] <- PercentageFeatureSet(gse208625, pattern = "^MT-")

ggsave("results/figures/QC/violin_pre_filter.png",
       VlnPlot(gse208625,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "condition", ncol = 3, pt.size = 0),
       width = 14, height = 5, dpi = 300, bg = "white")

gse208625 <- subset(gse208625,
                    subset = nFeature_RNA >= 500 & nFeature_RNA <= 8000 &
                      nCount_RNA   >= 1000 & percent.mt   <  20)

ggsave("results/figures/QC/violin_post_filter.png",
       VlnPlot(gse208625,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "condition", ncol = 3, pt.size = 0),
       width = 14, height = 5, dpi = 300, bg = "white")

gse208625 <- JoinLayers(gse208625)
saveRDS(gse208625, "results/objects/GSE208625_filtered.rds")


#Normalise, PCA, cluster

gse208625 <- gse208625 |>
  NormalizeData() |>
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
  ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt")) |>
  RunPCA(npcs = 50, verbose = FALSE)

ggsave("results/figures/dimred/elbow.png",
       ElbowPlot(gse208625, ndims = 50),
       width = 8, height = 5, dpi = 300, bg = "white")

ggsave("results/figures/dimred/variable_features.png",
       LabelPoints(VariableFeaturePlot(gse208625),
                   points = head(VariableFeatures(gse208625), 20),
                   repel = TRUE),
       width = 10, height = 7, dpi = 300, bg = "white")

n_pcs <- 15
gse208625 <- gse208625 |>
  FindNeighbors(dims = 1:n_pcs) |>
  FindClusters(resolution = 0.1) |>
  RunUMAP(dims = 1:n_pcs, verbose = FALSE)


#Doublet detection

sce <- as.SingleCellExperiment(gse208625)
sce <- scDblFinder(sce, clusters = gse208625$seurat_clusters)
gse208625$doublet_class <- sce$scDblFinder.class
rm(sce)

gse208625 <- subset(gse208625,
                    subset = doublet_class == "singlet" &
                      seurat_clusters != "3")

gse208625 <- gse208625 |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt")) |>
  RunPCA(npcs = 50, verbose = FALSE) |>
  FindNeighbors(dims = 1:n_pcs) |>
  FindClusters(resolution = 0.1) |>
  RunUMAP(dims = 1:n_pcs, verbose = FALSE)

saveRDS(gse208625, "results/objects/GSE208625_clean.rds")


#Celltype annotation

canonical <- list(
  Pluripotency  = c("POU5F1", "NANOG", "LIN28A", "EPCAM"),
  NSC           = c("PAX6", "SOX1", "SOX2", "NES", "HES5", "VIM"),
  Neuron        = c("DCX", "STMN2", "NEFL", "PHOX2A", "MAP2")
)
canonical <- lapply(canonical, function(g) g[g %in% rownames(gse208625)])

ggsave("results/figures/annotation/canonical_dotplot.png",
       DotPlot(gse208625, features = canonical,
               cols = c("lightgrey", "darkblue")) + RotatedAxis(),
       width = 14, height = 5, dpi = 300, bg = "white")

cell_type_map <- c("0" = "iPSC", "1" = "NSC", "2" = "Neuron", "3" = "iPSC")
gse208625$cell_type <- unname(cell_type_map[as.character(Idents(gse208625))])
Idents(gse208625) <- "cell_type"

ggsave("results/figures/annotation/umap_annotated.png",
       DimPlot(gse208625, reduction = "umap", group.by = "cell_type",
               label = TRUE, label.size = 5, pt.size = 0.3),
       width = 10, height = 8, dpi = 300, bg = "white")

saveRDS(gse208625, "results/objects/GSE208625_annotated.rds")


#DEGs

deg_nsc <- FindMarkers(gse208625, ident.1 = "NSC", ident.2 = "iPSC",
                       test.use = "wilcox",
                       min.pct = 0.1, logfc.threshold = 0.25) |>
  rownames_to_column("gene") |>
  arrange(desc(avg_log2FC))
write.csv(deg_nsc, "results/tables/DEG_NSC_vs_iPSC.csv", row.names = FALSE)

gse208625$lineage <- ifelse(gse208625$cell_type == "iPSC", "iPSC", "Neural")
Idents(gse208625) <- "lineage"
deg_neural <- FindMarkers(gse208625, ident.1 = "Neural", ident.2 = "iPSC",
                          test.use = "wilcox",
                          min.pct = 0.1, logfc.threshold = 0.25) |>
  rownames_to_column("gene") |>
  arrange(desc(avg_log2FC))
write.csv(deg_neural, "results/tables/DEG_Neural_vs_iPSC.csv", row.names = FALSE)


#Volcano plot

vd <- deg_nsc |>
  mutate(direction = case_when(
    p_val_adj < 0.05 & avg_log2FC >=  1 ~ "Up in NSC",
    p_val_adj < 0.05 & avg_log2FC <= -1 ~ "Up in iPSC",
    TRUE                                ~ "NS"),
    neg_log10_padj = -log10(p_val_adj + 1e-300))

label_genes <- vd |>
  filter(direction != "NS") |>
  group_by(direction) |>
  slice_max(abs(avg_log2FC), n = 10) |>
  pull(gene)

ggsave("results/figures/DEG/volcano_NSC_vs_iPSC.png",
       ggplot(vd, aes(avg_log2FC, neg_log10_padj, color = direction)) +
         geom_point(alpha = 0.5, size = 1) +
         scale_color_manual(values = c("Up in NSC"  = "#E64B35",
                                       "Up in iPSC" = "#3182bd",
                                       "NS"         = "grey80")) +
         geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
         geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
         geom_text_repel(data = filter(vd, gene %in% label_genes),
                         aes(label = gene), size = 3, max.overlaps = 30,
                         show.legend = FALSE) +
         theme_bw() +
         labs(title = "NSC vs iPSC",
              x = "log2 fold change (NSC / iPSC)",
              y = "-log10 adjusted p-value", color = NULL),
       width = 11, height = 8, dpi = 300, bg = "white")


