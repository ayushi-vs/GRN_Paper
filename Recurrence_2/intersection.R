#Crossmethod intersection

suppressPackageStartupMessages({
  library(RTN); library(igraph); library(dplyr)
})

setwd("~/Desktop/GRN_Paper_Deposit")

rtni      <- readRDS("results/objects/RTN_final.rds")
mra       <- readRDS("results/objects/RTN_MRA.rds")
signed    <- readRDS("results/objects/consensus_signed.rds")
regs_mode <- tni.get(rtni, what = "regulons.and.mode")

mr_dir <- lapply(mra$Regulon, function(tf) {
  targets    <- names(regs_mode[[tf]])
  signed_sub <- signed[signed$gene %in% targets, ]
  mlfc <- if (nrow(signed_sub) > 0)
    mean(rowMeans(as.matrix(signed_sub[, c("D1","D2","D3","D4")]),
                  na.rm = TRUE), na.rm = TRUE)
  else NA_real_
  
  data.frame(
    TF                   = tf,
    regulon_size         = length(targets),
    n_targets_in_signed  = nrow(signed_sub),
    n_conc_up            = sum(signed_sub$class == "Concordant_up"),
    n_conc_down          = sum(signed_sub$class == "Concordant_down"),
    mean_target_lfc      = round(mlfc, 3),
    direction            = ifelse(is.na(mlfc), "undetermined",
                                  ifelse(mlfc > 0, "drives_up (neural)",
                                         "drives_down (pluripotency)")),
    adj_pvalue           = mra$Adjusted.Pvalue[mra$Regulon == tf],
    stringsAsFactors = FALSE
  )
})
mr_dir_df <- do.call(rbind, mr_dir)
saveRDS(mr_dir_df, "results/objects/RTN_MRA_directional.rds")


#Network centrality via igraph

edges <- do.call(rbind, lapply(names(regs_mode), function(tf) {
  r <- regs_mode[[tf]]
  if (length(r) == 0) return(NULL)
  data.frame(from = tf, to = names(r),
             weight = abs(unname(r)),
             stringsAsFactors = FALSE)
}))

g       <- igraph::graph_from_data_frame(edges, directed = TRUE)
tfs     <- names(regs_mode)
centr   <- data.frame(
  TF           = tfs,
  out_degree   = igraph::degree(g, mode = "out")[tfs],
  betweenness  = igraph::betweenness(g, directed = TRUE, normalized = TRUE)[tfs],
  regulon_size = sapply(regs_mode, length)[tfs],
  stringsAsFactors = FALSE
)
centr <- centr[order(-centr$out_degree), ]
saveRDS(centr, "results/objects/RTN_centrality.rds")

consensus_4of4 <- readRDS("results/objects/consensus_core_4of4.rds")
core_hubs      <- readRDS("results/objects/INTERSECTION_consensus_hub_TFs.rds")
mra_sig        <- mra$Regulon
top_centrality <- head(centr$TF, 50)

triple <- data.frame(
  TF                    = consensus_4of4,
  recurrence_4of4       = TRUE,
  WGCNA_hub_NCAM1       = consensus_4of4 %in% core_hubs$gene_name,
  RTN_master_regulator  = consensus_4of4 %in% mra_sig,
  RTN_top50_centrality  = consensus_4of4 %in% top_centrality,
  stringsAsFactors = FALSE
)
triple$kME <- ifelse(triple$WGCNA_hub_NCAM1,
                     core_hubs$kME[match(triple$TF, core_hubs$gene_name)],
                     NA)
triple$regulon_size <- sapply(triple$TF, function(t)
  if (t %in% names(regs_mode)) length(regs_mode[[t]]) else 0)
triple$all_three_methods <- triple$recurrence_4of4 &
  triple$WGCNA_hub_NCAM1 &
  triple$RTN_master_regulator

saveRDS(triple, "results/objects/TRIPLE_INTERSECTION.rds")
