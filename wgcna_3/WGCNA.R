#WGCNA 

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(harmony); library(dplyr)
  library(tibble); library(tidyr); library(hdWGCNA); library(WGCNA); library(patchwork)
})
WGCNA::disableWGCNAThreads()

setwd("~/Desktop/GRN_Paper_Deposit")
mem.maxVSize(64000)
options(future.globals.maxSize = 4 * 1024^3)
set.seed(12345)

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/wgcna", recursive = TRUE, showWarnings = FALSE)

paths <- list(
  D1 = "results/objects/GSE208625_annotated.rds",
  D2 = "results/objects/GSE185275_annotated.rds",
  D3 = "results/objects/GSE138121_annotated.rds",
  D4 = "results/objects/PRJEB38269_annotated.rds"
)

for (nm in names(paths)) {
  o    <- readRDS(paths[[nm]])
  md   <- o@meta.data
  ctcol <- if ("cell_type" %in% colnames(md)) "cell_type" else "celltype"
  ct    <- as.character(md[[ctcol]])
  cnts  <- GetAssayData(o, assay = "RNA", layer = "counts")
  rm(o); gc(verbose = FALSE)
  s <- CreateSeuratObject(counts = cnts, project = nm)
  s$cell_type <- ct
  s$dataset   <- nm
  rm(cnts)
  saveRDS(s, sprintf("results/objects/_slim_%s.rds", nm))
  rm(s)
}


#Merge 

slim   <- lapply(names(paths),
                 function(nm) readRDS(sprintf("results/objects/_slim_%s.rds", nm)))
names(slim) <- names(paths)
merged <- merge(slim[[1]], y = slim[-1],
                add.cell.ids = names(paths),
                project = "iPSC_neural_4ds")
rm(slim)
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
saveRDS(merged, "results/objects/merged_4ds_raw.rds")

merged$analysis_group <- "neural"

merged <- merged |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 3000, verbose = FALSE) |>
  ScaleData(features = VariableFeatures(merged), verbose = FALSE) |>
  RunPCA(features = VariableFeatures(merged),
         npcs = 30, verbose = FALSE)

merged <- RunHarmony(merged,
                     group.by.vars = "dataset",
                     reduction.use = "pca",
                     dims.use = 1:30,
                     reduction.save = "harmony",
                     verbose = FALSE)

saveRDS(merged, "results/objects/merged_4ds_harmony.rds")


#hdWGCNA and construction of metacells

merged[["RNA"]]$scale.data <- NULL
gc(verbose = FALSE)

merged <- SetupForWGCNA(
  merged,
  gene_select = "fraction",
  fraction    = 0.05,
  group.by    = "dataset",
  wgcna_name  = "neural_consensus"
)

merged <- MetacellsByGroups(
  seurat_obj       = merged,
  group.by         = c("analysis_group", "dataset"),
  ident.group      = "analysis_group",
  reduction        = "harmony",
  k                = 25,
  max_shared       = 12,
  min_cells        = 50,
  target_metacells = 250
)
merged <- NormalizeMetacells(merged)
saveRDS(merged, "results/objects/merged_4ds_metacells.rds")

merged <- SetDatExpr(
  merged,
  group.by      = "analysis_group",
  group_name    = "neural",
  use_metacells = TRUE,
  assay         = "RNA",
  slot          = "data"
)
merged <- SetMultiExpr(
  merged,
  group_name     = "neural",
  group.by       = "analysis_group",
  multi.group.by = "dataset",
  multi_groups   = c("D1", "D2", "D3", "D4")
)

merged <- TestSoftPowersConsensus(
  merged,
  setDatExpr  = FALSE,
  networkType = "signed"
)
saveRDS(merged, "results/objects/merged_4ds_softpower.rds")

p_soft <- PlotSoftPowers(merged)
collect_ggplots <- function(x, acc = list()) {
  if (inherits(x, "ggplot")) return(c(acc, list(x)))
  if (is.list(x)) for (el in x) acc <- collect_ggplots(el, acc)
  acc
}
plots <- collect_ggplots(p_soft)
ggplot2::ggsave("results/figures/wgcna_softpower_consensus.png",
                wrap_plots(plots, ncol = 4) +
                  plot_annotation(title = "Soft thresholding: D1-D4"),
                width = 18, height = 14, dpi = 300, bg = "white")


#Construct network

merged <- ConstructNetwork(
  merged,
  soft_power    = c(10, 6, 14, 9),
  consensus     = TRUE,
  tom_name      = "neural_consensus",
  overwrite_tom = TRUE
)
saveRDS(merged, "results/objects/merged_4ds_network.rds")


#Compute eigengenes and kME 

mods <- GetModules(merged)
mc   <- GetMetacellObject(merged)
expr <- t(as.matrix(GetAssayData(mc, assay = "RNA", layer = "data")))
expr <- expr[, intersect(mods$gene_name, colnames(expr)), drop = FALSE]

gene2module <- setNames(as.character(mods$module), mods$gene_name)
gene2module <- gene2module[colnames(expr)]

me_res <- WGCNA::moduleEigengenes(expr, colors = gene2module, excludeGrey = TRUE)
MEs    <- me_res$eigengenes
kME    <- WGCNA::cor(expr, MEs, use = "pairwise.complete.obs")

for (m in colnames(MEs)) {
  col <- paste0("kME_", sub("^ME", "", m))
  mods[[col]] <- kME[mods$gene_name, m]
}

saveRDS(MEs,  "results/objects/MEs_metacells.rds")
saveRDS(kME,  "results/objects/kME_matrix.rds")
saveRDS(mods, "results/objects/modules_with_kME.rds")
saveRDS(expr, "results/objects/metacell_expr.rds")


#Cell type assignmennt

ct_lookup <- setNames(merged@meta.data$cell_type, rownames(merged@meta.data))

mc_md <- mc@meta.data |> rownames_to_column("mc_id")
parse_mc <- mc_md |>
  mutate(cells = strsplit(cells_merged, ",", fixed = TRUE)) |>
  select(mc_id, dataset, cells) |>
  unnest(cells) |>
  mutate(cell_type = ct_lookup[cells])

mc_ct <- parse_mc |>
  filter(!is.na(cell_type)) |>
  group_by(mc_id) |>
  summarise(
    cell_type = names(sort(table(cell_type), decreasing = TRUE))[1],
    purity    = max(table(cell_type)) / n(),
    .groups   = "drop"
  )
mc@meta.data$cell_type <-
  setNames(mc_ct$cell_type, mc_ct$mc_id)[rownames(mc@meta.data)]

stage_map <- c(
  iPSC              = "1_iPSC",
  Early_Progenitor  = "2_Progenitor",  Neural_Progenitor  = "2_Progenitor",
  MN_Progenitor     = "2_Progenitor",  Progenitor         = "2_Progenitor",
  Cycling_Progenitor= "2_Progenitor",  FPP                = "2_Progenitor",
  P_FPP             = "2_Progenitor",  NB                 = "2_Progenitor",
  d18_Progenitor    = "2_Progenitor",
  Transitional      = "3_Transitional", Proneural         = "3_Transitional",
  NSC               = "4_NSC",
  Neuron            = "5_Neuron",       Motor_Neuron      = "5_Neuron",
  DA_Neuron         = "5_Neuron",       DA                = "5_Neuron",
  Sert              = "5_Neuron",       P_Sert            = "5_Neuron",
  U_Neur1           = "5_Neuron",       U_Neur2           = "5_Neuron",
  Epen1             = "6_Other",        Astro             = "6_Other"
)
mc@meta.data$stage <- stage_map[mc@meta.data$cell_type]

saveRDS(mc, "results/objects/metacell_obj_with_celltype.rds")


#coorelation between module and trait

md     <- mc@meta.data
stages <- sort(unique(na.omit(md$stage)))

trait <- sapply(stages, function(s) as.integer(md$stage == s))
colnames(trait) <- stages
trait <- cbind(trait, diff_score = as.integer(sub("^([0-9]+).*", "\\1", md$stage)))
for (d in unique(md$dataset))
  trait <- cbind(trait, setNames(as.integer(md$dataset == d), paste0("ds_", d)))
rownames(trait) <- rownames(md)
trait <- trait[rownames(MEs), , drop = FALSE]

mt_cor  <- WGCNA::cor(MEs, trait, use = "pairwise.complete.obs")
mt_p    <- WGCNA::corPvalueStudent(mt_cor, nSamples = nrow(MEs))
mt_p_bh <- matrix(p.adjust(mt_p, method = "BH"),
                  nrow = nrow(mt_p), dimnames = dimnames(mt_p))

ds_cols <- grep("^ds_", colnames(mt_cor), value = TRUE)

ranked <- data.frame(
  module      = sub("^ME", "", rownames(mt_cor)),
  cor_diff    = mt_cor[, "diff_score"],
  p_diff_BH   = mt_p_bh[, "diff_score"],
  cor_iPSC    = mt_cor[, "1_iPSC"],
  cor_Neuron  = mt_cor[, "5_Neuron"],
  max_ds_cor  = apply(abs(mt_cor[, ds_cols, drop = FALSE]), 1, max),
  dominant_ds = apply(mt_cor[, ds_cols, drop = FALSE], 1,
                      function(v) ds_cols[which.max(abs(v))])
) |>
  mutate(dataset_dominated = max_ds_cor > 0.6) |>
  arrange(desc(abs(cor_diff)))

saveRDS(mt_cor,  "results/objects/module_trait_cor.rds")
saveRDS(mt_p_bh, "results/objects/module_trait_padj.rds")
saveRDS(ranked,  "results/objects/module_trait_ranked.rds")

TFs      <- read.csv("data/raw/Lambert2018_TFs.csv", stringsAsFactors = FALSE)
TF_set   <- unique(TFs$HGNC.symbol[TFs$Is.TF. == "Yes"])

# Top hub TF
top_hubs <- lapply(unique(mods$module), function(m) {
  gene_col <- paste0("kME_", m)
  if (!gene_col %in% names(mods)) return(NA_character_)
  sub_mod <- mods[mods$module == m & mods$gene_name %in% TF_set, ]
  if (nrow(sub_mod) == 0) return(NA_character_)
  sub_mod$gene_name[which.max(abs(sub_mod[[gene_col]]))]
})
names(top_hubs) <- unique(mods$module)

labels_df <- data.frame(
  module = unique(mods$module),
  top_hub = unlist(top_hubs),
  stringsAsFactors = FALSE
) |>
  left_join(ranked |> select(module, cor_diff), by = "module") |>
  mutate(direction = case_when(
    cor_diff >=  0.5 ~ "trajectory-up",
    cor_diff <= -0.5 ~ "trajectory-down",
    TRUE             ~ "trajectory-neutral"))

labels_df$biology <- NA_character_
labels_df$biology[labels_df$module == "turquoise"] <- "neural maturation"
labels_df$biology[labels_df$module == "green"]     <- "pluripotency"
labels_df$biology[labels_df$module == "pink"]      <- "cell cycle"
labels_df$biology[labels_df$module == "magenta"]   <- "glial / off-trajectory"

labels_df$label <- ifelse(
  !is.na(labels_df$biology),
  paste0(labels_df$top_hub, " module (", labels_df$biology, ")"),
  ifelse(labels_df$direction != "trajectory-neutral",
         paste0(labels_df$top_hub, " module (", labels_df$direction, ")"),
         paste0(labels_df$top_hub, " module"))
)

saveRDS(labels_df, "results/objects/module_labels.rds")

core4 <- readRDS("results/objects/consensus_core_4of4.rds")

focus_mods <- unique(mods$module[mods$gene_name %in% core4])
saveRDS(focus_mods, "results/objects/focus_modules.rds")

hub_tf_list <- lapply(focus_mods, function(m) {
  gene_col <- paste0("kME_", m)
  sub_mod  <- mods[mods$module == m & mods$gene_name %in% TF_set, ]
  if (!gene_col %in% names(sub_mod) || nrow(sub_mod) == 0) return(character(0))
  sub_mod$gene_name[order(abs(sub_mod[[gene_col]]), decreasing = TRUE)][1:min(20, nrow(sub_mod))]
})
names(hub_tf_list) <- focus_mods
saveRDS(hub_tf_list, "results/objects/hub_TFs_by_module.rds")

core_hubs <- intersect(unlist(hub_tf_list), core4)
saveRDS(core_hubs, "results/objects/INTERSECTION_consensus_hub_TFs.rds")

core_module_map <- mods |>
  filter(gene_name %in% core4) |>
  select(TF = gene_name, module) |>
  arrange(module, TF)
saveRDS(core_module_map, "results/objects/consensus_TFs_module_map.rds")
write.csv(core_module_map,
          "results/tables/INTERSECTION_4of4_WGCNA_hub_TFs.csv",
          row.names = FALSE)
