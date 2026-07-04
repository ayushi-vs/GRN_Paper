#RTN

suppressPackageStartupMessages({
  library(RTN); library(snow); library(WGCNA); library(dplyr)
})
WGCNA::disableWGCNAThreads()

setwd("~/Desktop/GRN_Paper_Deposit")
set.seed(12345)

TFs    <- read.csv("data/raw/Lambert2018_TFs.csv", stringsAsFactors = FALSE)
TF_set <- unique(TFs$HGNC.symbol[TFs$Is.TF. == "Yes"])

expr <- readRDS("results/objects/metacell_expr.rds")
expr <- t(expr)                                        

mc <- readRDS("results/objects/metacell_obj_with_celltype.rds")

rowAnnotation <- data.frame(
  PROBEID = rownames(expr),
  SYMBOL  = rownames(expr),
  row.names = rownames(expr),
  stringsAsFactors = FALSE
)

colAnnotation <- data.frame(
  Sample    = rownames(mc@meta.data),
  dataset   = mc@meta.data$dataset,
  cell_type = mc@meta.data$cell_type,
  stage     = mc@meta.data$stage,
  row.names = rownames(mc@meta.data),
  stringsAsFactors = FALSE
)
colAnnotation <- colAnnotation[colnames(expr), , drop = FALSE]

tfs_present <- intersect(TF_set, rownames(expr))


rtni <- tni.constructor(
  expData            = as.matrix(expr),
  regulatoryElements = tfs_present,
  rowAnnotation      = rowAnnotation,
  colAnnotation      = colAnnotation,
  cvfilter           = FALSE
)

options(cluster = snow::makeCluster(4, "SOCK"))
rtni <- tni.permutation(rtni, nPermutations = 1000,
                        pValueCutoff = 1e-5, verbose = TRUE)
saveRDS(rtni, "results/objects/RTN_after_permutation.rds")

rtni <- tni.bootstrap(rtni, nBootstraps = 100, consensus = 95, verbose = TRUE)
snow::stopCluster(getOption("cluster"))
saveRDS(rtni, "results/objects/RTN_after_bootstrap.rds")

rtni <- tni.dpi.filter(rtni, eps = 0)
saveRDS(rtni, "results/objects/RTN_final.rds")


signed <- readRDS("results/objects/consensus_signed.rds")

ds_mat   <- as.matrix(signed[, c("D1", "D2", "D3", "D4")])
mean_lfc <- rowMeans(ds_mat, na.rm = TRUE)
mean_lfc[is.nan(mean_lfc)] <- 0
names(mean_lfc) <- signed$gene

net_genes <- tni.get(rtni, what = "tnet")
all_genes <- rownames(net_genes)

phenotype <- setNames(rep(0.0, length(all_genes)), all_genes)
shared    <- intersect(names(mean_lfc), all_genes)
phenotype[shared] <- mean_lfc[shared]

hits <- c(
  signed$gene[signed$class == "Concordant_up"   & signed$n_de >= 3],
  signed$gene[signed$class == "Concordant_down" & signed$n_de >= 3]
)
hits <- intersect(hits, all_genes)

phenoIDs <- data.frame(PROBEID = all_genes, SYMBOL = all_genes)

rtna <- tni2tna.preprocess(
  object    = rtni,
  phenotype = phenotype,
  hits      = hits,
  phenoIDs  = phenoIDs
)

rtna <- tna.mra(rtna, pValueCutoff = 0.05, pAdjustMethod = "BH")
mra  <- tna.get(rtna, what = "mra")

rtna  <- tna.gsea2(rtna, nPermutations = 1000,
                   pValueCutoff = 0.05, pAdjustMethod = "BH")
gsea2 <- tna.get(rtna, what = "gsea2")

saveRDS(rtna,  "results/objects/RTN_tna.rds")
saveRDS(mra,   "results/objects/RTN_MRA.rds")
saveRDS(gsea2, "results/objects/RTN_GSEA2.rds")

