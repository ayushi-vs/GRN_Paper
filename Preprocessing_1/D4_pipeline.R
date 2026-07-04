#D4(PRJEB38269)
if (Sys.info()["sysname"] == "Darwin") try(mem.maxVSize(vsize = 64000), silent = TRUE)
options(Seurat.object.assay.version = "v3")

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(ggplot2); library(patchwork)
  library(dplyr); library(Matrix); library(methods)
})

setwd("~/Desktop/GRN_Paper_Deposit")
options(future.globals.maxSize = 8 * 1024^3)
set.seed(42)

IN_PATH   <- "results/objects/PRJEB38269_raw.rds"
OBJ_DIR   <- "results/objects"
FIG_DIR   <- "results/figures/PRJEB38269"
TABLE_DIR <- "results/tables"
N_VARFEATS <- 2000
N_PCS      <- 30
RESOLUTION <- 0.5


for (d in c(file.path(FIG_DIR, "QC"),
            file.path(FIG_DIR, "dimred"),
            file.path(FIG_DIR, "annotation"),
            file.path(FIG_DIR, "DEG"),
            TABLE_DIR, OBJ_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, path, w = 8, h = 6)
  ggsave(path, plot = p, width = w, height = h, dpi = 300, bg = "white")


#convert to Seurat

extract_counts_and_meta <- function(x) {
  if (inherits(x, "Seurat")) {
    return(list(counts = GetAssayData(x, assay = DefaultAssay(x), layer = "counts"),
                meta   = as.data.frame(x@meta.data)))
  }
  avail <- SummarizedExperiment::assayNames(x)
  counts_name <- intersect(c("counts", "X", "raw"), avail)[1]
  if (is.na(counts_name)) counts_name <- avail[1]
  cts <- SummarizedExperiment::assay(x, counts_name)
  if (!inherits(cts, "dgCMatrix")) {
    cts2 <- try(as(cts, "dgCMatrix"), silent = TRUE)
    if (inherits(cts2, "try-error") || !inherits(cts2, "dgCMatrix"))
      cts <- as(as(as.matrix(cts), "CsparseMatrix"), "dgCMatrix")
    else
      cts <- cts2
  }
  if (is.null(rownames(cts))) rownames(cts) <- rownames(x)
  if (is.null(colnames(cts))) colnames(cts) <- colnames(x)
  
  md <- as.data.frame(SummarizedExperiment::colData(x))
  for (j in seq_len(ncol(md)))
    if (isS4(md[[j]])) md[[j]] <- as.vector(md[[j]])
  
  list(counts = cts, meta = md)
}

extracted <- extract_counts_and_meta(readRDS(IN_PATH))
obj <- CreateSeuratObject(counts    = extracted$counts,
                          meta.data = extracted$meta,
                          project   = "PRJEB38269",
                          min.cells = 3, min.features = 200)
rm(extracted)

if (isS4(obj@assays) && !is.list(obj@assays))
  obj@assays <- if ("listData" %in% slotNames(obj@assays))
    obj@assays@listData else S4Vectors::as.list(obj@assays)

DefaultAssay(obj) <- "RNA"


#detect donor/stage columns

md_cols <- colnames(obj@meta.data)
DONOR_COLUMN <- intersect(
  c("donor_id","donor","pool_id","pool","line","cell_line",
    "sample_id","sample","experiment","Pool","Donor"), md_cols)[1]
STAGE_COLUMN <- intersect(
  c("day","Day","stage","timepoint","time","differentiation_day",
    "time_point","cell_type"), md_cols)[1]
HARMONY_VAR  <- DONOR_COLUMN

if (!"percent.mt" %in% md_cols)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")


#QC

qc_feats <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
group_var <- if (!is.na(STAGE_COLUMN)) STAGE_COLUMN
              {if (!is.na(DONOR_COLUMN)) DONOR_COLUMN
              else "orig.ident"
}
save_fig(VlnPlot(obj, features = qc_feats, group.by = group_var,
                 pt.size = 0, ncol = 3) &
           theme(axis.text.x = element_text(angle = 45, hjust = 1)),
         file.path(FIG_DIR, "QC/violin_pre_filter.png"), w = 14, h = 5)

nf_lo <- quantile(obj$nFeature_RNA, 0.02, na.rm = TRUE)
nf_hi <- max(quantile(obj$nFeature_RNA, 0.99, na.rm = TRUE), 200)
mt_hi <- quantile(obj$percent.mt,      0.95, na.rm = TRUE)

obj <- subset(obj, subset = nFeature_RNA >= nf_lo &
                nFeature_RNA <= nf_hi &
                percent.mt   <= mt_hi)
saveRDS(obj, file.path(OBJ_DIR, "PRJEB38269_filtered.rds"))

save_fig(VlnPlot(obj, features = qc_feats, group.by = group_var,
                 pt.size = 0, ncol = 3) &
           theme(axis.text.x = element_text(angle = 45, hjust = 1)),
         file.path(FIG_DIR, "QC/violin_post_filter.png"), w = 14, h = 5)


#Normalisation and pca

obj <- obj |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = N_VARFEATS, verbose = FALSE) |>
  ScaleData(features = NULL,
            vars.to.regress = c("nCount_RNA", "percent.mt"),
            verbose = FALSE) |>
  RunPCA(npcs = N_PCS, verbose = FALSE)

save_fig(LabelPoints(VariableFeaturePlot(obj),
                     points = head(VariableFeatures(obj), 10),
                     repel = TRUE),
         file.path(FIG_DIR, "dimred/variable_features.png"))
save_fig(ElbowPlot(obj, ndims = N_PCS),
         file.path(FIG_DIR, "dimred/elbow.png"))


#Preharmony UMAP

obj <- obj |>
  FindNeighbors(reduction = "pca", dims = 1:N_PCS, verbose = FALSE) |>
  FindClusters(resolution = RESOLUTION, verbose = FALSE) |>
  RunUMAP(reduction = "pca", dims = 1:N_PCS,
          reduction.name = "umap_prepca", verbose = FALSE)

colour_var <- if (!is.na(HARMONY_VAR)) HARMONY_VAR
                {if (!is.na(STAGE_COLUMN)) STAGE_COLUMN
                  else "seurat_clusters"
}
save_fig(DimPlot(obj, reduction = "umap_prepca", group.by = colour_var,
                 raster = TRUE, raster.dpi = c(1000, 1000)) + NoLegend(),
         file.path(FIG_DIR, "dimred/preintegration.png"), w = 8, h = 7)


#Harmony integration on donor_id

if (!is.na(HARMONY_VAR) &&
    length(unique(obj@meta.data[[HARMONY_VAR]])) > 1) {
  obj <- RunHarmony(obj, group.by.vars = HARMONY_VAR, reduction.use = "pca", 
                    dims.use = 1:N_PCS,
                    plot_convergence = FALSE,
                    project.dim = FALSE)
  obj <- obj |>
    RunUMAP(reduction = "harmony", dims = 1:N_PCS,
            reduction.name = "umap", verbose = FALSE) |>
    FindNeighbors(reduction = "harmony", dims = 1:N_PCS, verbose = FALSE) |>
    FindClusters(resolution = RESOLUTION, verbose = FALSE)
  
  save_fig(DimPlot(obj, reduction = "umap", group.by = colour_var,
                   raster = TRUE, raster.dpi = c(1000, 1000)) + NoLegend(),
           file.path(FIG_DIR, "dimred/postharmony.png"), w = 8, h = 7)
} else {
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:N_PCS,
                 reduction.name = "umap", verbose = FALSE)
}

saveRDS(obj, file.path(OBJ_DIR, "PRJEB38269_clustered.rds"))


#Canonical marker dot plot

markers <- intersect(
  c("TH","NR4A2","LMX1A","FOXA2","PITX3","EN1","CALB1","DDC",
    "SOX2","VIM","NES","HES1","ASCL1","NEUROG2"),
  rownames(obj))

save_fig(DotPlot(obj, features = markers, group.by = "seurat_clusters") +
           RotatedAxis(),
         file.path(FIG_DIR, "dimred/trajectory_markers.png"),
         w = 10, h = 6)


#Celltype annotation 

label_map <- c(
  FPP = "Progenitor", P_FPP = "Progenitor",
  DA  = "DA Neuron",
  Sert = "Other", P_Sert = "Other", Astro = "Other",
  Epen1 = "Other", Epen2 = "Other", NB = "Other",
  U_Neur1 = "Other", U_Neur2 = "Other", U_Neur3 = "Other"
)
obj$cell_type <- factor(
  unname(label_map[as.character(obj$celltype)]),
  levels = c("Progenitor", "Other", "DA Neuron")
)

save_fig(DimPlot(obj, reduction = "umap", group.by = "celltype",
                 label = TRUE, repel = TRUE,
                 raster = TRUE, raster.dpi = c(1000, 1000)),
         file.path(FIG_DIR, "annotation/annotated_umap.png"),
         w = 10, h = 7)

saveRDS(obj, file.path(OBJ_DIR, "PRJEB38269_annotated.rds"))


#Differential expression: DA Neuron vs Progenitor

Idents(obj) <- obj$cell_type
de <- FindMarkers(obj, ident.1 = "DA Neuron", ident.2 = "Progenitor",
                  test.use = "wilcox", min.pct = 0.1,
                  logfc.threshold = 0, verbose = FALSE)
de$gene   <- rownames(de)
de$signif <- ifelse(de$p_val_adj < 0.05 & abs(de$avg_log2FC) >= 1,
                    ifelse(de$avg_log2FC > 0, "Up", "Down"), "ns")

write.csv(de,
          file.path(TABLE_DIR, "PRJEB38269_DEG_Neuron_vs_Progenitor.csv"),
          row.names = FALSE)


#Volcano plot

top_up   <- de |> filter(signif == "Up")   |> arrange(desc(avg_log2FC)) |> head(8)
top_down <- de |> filter(signif == "Down") |> arrange(avg_log2FC)       |> head(8)

save_fig(
  ggplot(de, aes(avg_log2FC, -log10(p_val_adj + 1e-300), colour = signif)) +
    geom_point(alpha = 0.5, size = 0.7) +
    scale_colour_manual(values = c(Up = "#E31A1C", Down = "#1F78B4",
                                   ns = "grey70")) +
    geom_vline(xintercept = c(-1, 1), lty = 2, colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), lty = 2, colour = "grey40") +
    ggrepel::geom_text_repel(data = rbind(top_up, top_down),
                             aes(label = gene), size = 3,
                             max.overlaps = 30, box.padding = 0.3) +
    labs(title = "D4: DA Neuron vs Progenitor",
         x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"(adj. p)"),
         colour = "") +
    theme_minimal(base_size = 11),
  file.path(FIG_DIR, "DEG/volcano_neuron_vs_progenitor.png"))

