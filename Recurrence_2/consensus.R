#Cross dataset consensus meta analysis

USE_PRECOMPUTED <- TRUE

suppressPackageStartupMessages({
  library(tidyverse); library(UpSetR)
})

setwd("~/Desktop/GRN_Paper_Deposit")

dir.create("results/figures/consensus", recursive = TRUE, showWarnings = FALSE)

lambert_path <- "data/raw/Lambert2018_TFs.csv"
if (!file.exists(lambert_path)) {
  dir.create(dirname(lambert_path), showWarnings = FALSE, recursive = TRUE)
  download.file(
    "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv",
    lambert_path, mode = "wb", quiet = TRUE)
}
tf_table <- read.csv(lambert_path, stringsAsFactors = FALSE)
TF_set   <- unique(tf_table$HGNC.symbol[tf_table$Is.TF. == "Yes"])


if (USE_PRECOMPUTED) {
  d1_tf <- read.csv("data/intermediates/TF_candidates_NSC_vs_iPSC.csv",
                    stringsAsFactors = FALSE)
  tf_lists <- list(
    D1 = unique(d1_tf$gene[d1_tf$avg_log2FC >= 1]),
    D2 = unique(readRDS("data/intermediates/GSE185275_neural_TFs.rds")),
    D3 = unique(readRDS("data/intermediates/GSE138121_neural_TFs.rds")),
    D4 = unique(readRDS("data/intermediates/PRJEB38269_neural_TFs.rds"))
  )
} else {
  # Recompute from DEG CSVs produced by the minimal per-dataset pipelines
  deg_files <- list(
    D1 = "results/tables/DEG_NSC_vs_iPSC.csv",
    D2 = "results/tables/GSE185275_DEG_Neuron_vs_Progenitor.csv",
    D3 = "results/tables/GSE138121_DEG_d18_vs_d00.csv",
    D4 = "results/tables/PRJEB38269_DEG_Neuron_vs_Progenitor.csv"
  )
  
  extract_tfs <- function(path, direction = "up",
                          lfc_cut = 1, padj_cut = 0.05) {
    d <- read.csv(path, stringsAsFactors = FALSE)
    if ("avg_log2FC" %in% names(d)) names(d)[names(d) == "avg_log2FC"] <- "log2FC"
    if ("avg_logFC"  %in% names(d)) names(d)[names(d) == "avg_logFC"]  <- "log2FC"
    if ("p_val_adj"  %in% names(d)) names(d)[names(d) == "p_val_adj"]  <- "padj"
    if (!"gene" %in% names(d) && "X" %in% names(d))
      names(d)[names(d) == "X"] <- "gene"
    d <- d[!is.na(d$log2FC) & !is.na(d$padj) & d$padj < padj_cut, ]
    if (direction == "up")   d <- d[d$log2FC >=  lfc_cut, ]
    if (direction == "down") d <- d[d$log2FC <= -lfc_cut, ]
    unique(intersect(d$gene, TF_set))
  }
  
  tf_lists <- lapply(deg_files, extract_tfs, direction = "up")
}


universe <- sort(unique(unlist(tf_lists)))
rec <- tibble(
  TF = universe,
  D1 = universe %in% tf_lists$D1,
  D2 = universe %in% tf_lists$D2,
  D3 = universe %in% tf_lists$D3,
  D4 = universe %in% tf_lists$D4
) |>
  mutate(n_datasets = D1 + D2 + D3 + D4,
         tier = case_when(
           n_datasets == 4 ~ "Tier1_core_4of4",
           n_datasets == 3 ~ "Tier2_3of4",
           n_datasets == 2 ~ "Tier3_2of4",
           TRUE            ~ "Tier4_dataset_specific"
         )) |>
  arrange(desc(n_datasets), TF)

write.csv(rec, "results/tables/CONSENSUS_4way_neural_TFs.csv", row.names = FALSE)
saveRDS(rec$TF[rec$n_datasets == 4],
        "results/objects/consensus_core_4of4.rds")

png("results/figures/consensus/4way_upset.png",
    width = 2400, height = 1400, res = 220)
suppressWarnings(print(upset(
  fromList(setNames(tf_lists,
                    c("D1 iPSC-NSC", "D2 prog-neuron",
                      "D3 iPSC-MN",   "D4 prog-DA"))),
  order.by = "freq", nsets = 4,
  text.scale = 1.3, point.size = 3,
  mainbar.y.label = "Shared neural TFs",
  sets.x.label = "TFs per dataset")))
dev.off()

deg_files <- list(
  D1 = "results/tables/DEG_NSC_vs_iPSC.csv",
  D2 = "results/tables/GSE185275_DEG_Neuron_vs_Progenitor.csv",
  D3 = "results/tables/GSE138121_DEG_d18_vs_d00.csv",
  D4 = "results/tables/PRJEB38269_DEG_Neuron_vs_Progenitor.csv"
)

get_log2fc <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  if ("avg_log2FC" %in% names(d)) names(d)[names(d) == "avg_log2FC"] <- "log2FC"
  if ("avg_logFC"  %in% names(d)) names(d)[names(d) == "avg_logFC"]  <- "log2FC"
  if ("p_val_adj"  %in% names(d)) names(d)[names(d) == "p_val_adj"]  <- "padj"
  if (!"gene" %in% names(d) && "X" %in% names(d))
    names(d)[names(d) == "X"] <- "gene"
  d[!is.na(d$log2FC), c("gene", "log2FC", "padj")]
}

de_all <- lapply(deg_files, get_log2fc)

tf_universe_signed <- sort(unique(unlist(lapply(de_all, function(d) {
  d <- d[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FC) >= 1, ]
  intersect(d$gene, TF_set)
}))))

signed_mat <- sapply(de_all, function(d) {
  v <- setNames(rep(NA_real_, length(tf_universe_signed)), tf_universe_signed)
  d <- d[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FC) >= 1, ]
  d <- d[d$gene %in% tf_universe_signed, ]
  v[d$gene] <- d$log2FC
  v
})

n_up   <- rowSums(signed_mat > 0, na.rm = TRUE)
n_down <- rowSums(signed_mat < 0, na.rm = TRUE)
n_de   <- rowSums(!is.na(signed_mat))

class_v <- case_when(
  n_de == 1                             ~ "Single_dataset",
  n_up  == n_de & n_down == 0           ~ "Concordant_up",
  n_down == n_de & n_up == 0            ~ "Concordant_down",
  TRUE                                  ~ "Discordant"
)

signed <- tibble(gene = tf_universe_signed,
                 class = class_v,
                 n_up = n_up, n_down = n_down, n_de = n_de) |>
  bind_cols(as_tibble(signed_mat)) |>
  arrange(desc(n_de), class, gene)

saveRDS(signed, "results/objects/consensus_signed.rds")
write.csv(signed, "results/tables/consensus_signed_FULL.csv", row.names = FALSE)

