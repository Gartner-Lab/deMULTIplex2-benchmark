

# Read in all summary results

library(deMULTIplex2)
savePath <- "./data-raw/benchmark_res/summary_res/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)


test_text<- "benchmark_summary_realandsim_"
show_ratios <- c(1, 0.5, 0.1, 0.05, 0.01)
bc_path <- list(
    "sim1" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim1_n.bc_5_bc_all.rds",
    "sim1_0.5" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim1_n.bc_5__read_ds_0.5bc_all.rds",
    "sim1_0.1" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim1_n.bc_5__read_ds_0.1bc_all.rds",
    "sim1_0.05" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim1_n.bc_5__read_ds_0.05bc_all.rds",
    "sim1_0.01" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim1_n.bc_5__read_ds_0.01bc_all.rds",

    "sim2" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim2_n.bc_10_bc_all.rds",
    "sim2_0.5" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim2_n.bc_10__read_ds_0.5bc_all.rds",
    "sim2_0.1" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim2_n.bc_10__read_ds_0.1bc_all.rds",
    "sim2_0.05" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim2_n.bc_10__read_ds_0.05bc_all.rds",
    "sim2_0.01" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim2_n.bc_10__read_ds_0.01bc_all.rds",

    "sim3" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim3_n.bc_10_bc_all.rds",
    "sim3_0.5" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim3_n.bc_10__read_ds_0.5bc_all.rds",
    "sim3_0.1" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim3_n.bc_10__read_ds_0.1bc_all.rds",
    "sim3_0.05" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim3_n.bc_10__read_ds_0.05bc_all.rds",
    "sim3_0.01" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim3_n.bc_10__read_ds_0.01bc_all.rds",

    "sim4" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim4_n.bc_30_bc_all.rds",
    "sim4_0.5" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim4_n.bc_30__read_ds_0.5bc_all.rds",
    "sim4_0.1" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim4_n.bc_30__read_ds_0.1bc_all.rds",
    "sim4_0.05" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim4_n.bc_30__read_ds_0.05bc_all.rds",
    "sim4_0.01" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim4_n.bc_30__read_ds_0.01bc_all.rds",

    "sim5" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim5_n.bc_30_bc_all.rds",
    "sim5_0.5" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim5_n.bc_30__read_ds_0.5bc_all.rds",
    "sim5_0.1" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim5_n.bc_30__read_ds_0.1bc_all.rds",
    "sim5_0.05" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim5_n.bc_30__read_ds_0.05bc_all.rds",
    "sim5_0.01" ="./data-raw/benchmark_res/simulation_model/rds/simulation_model_update_sim5_n.bc_30__read_ds_0.01bc_all.rds",

    "McGinnis_8donor_MULTI" = "./data-raw/benchmark_res/McGinnis_8donor_MULTI/rds/McGinnis_8donor_MULTI_benchmarking_bc_all.rds",
    "McGinnis_8donor_MULTI_0.5" = "./data-raw/benchmark_res/ds_McGinnis_8donor_MULTI/rds/ds_McGinnis_8donor_MULTI_0.5bc_all.rds",
    "McGinnis_8donor_MULTI_0.1" = "./data-raw/benchmark_res/ds_McGinnis_8donor_MULTI/rds/ds_McGinnis_8donor_MULTI_0.1bc_all.rds",
    "McGinnis_8donor_MULTI_0.05" = "./data-raw/benchmark_res/ds_McGinnis_8donor_MULTI/rds/ds_McGinnis_8donor_MULTI_0.05bc_all.rds",
    "McGinnis_8donor_MULTI_0.01" = "./data-raw/benchmark_res/ds_McGinnis_8donor_MULTI/rds/ds_McGinnis_8donor_MULTI_0.01bc_all.rds",
    "McGinnis_8donor_SCMK" = "./data-raw/benchmark_res/McGinnis_8donor_SCMK/rds/McGinnis_8donor_SCMK_benchmarking_bc_all.rds",
    "McGinnis_8donor_SCMK_0.5" = "./data-raw/benchmark_res/ds_McGinnis_8donor_SCMK/rds/ds_McGinnis_8donor_SCMK_0.5bc_all.rds",
    "McGinnis_8donor_SCMK_0.1" = "./data-raw/benchmark_res/ds_McGinnis_8donor_SCMK/rds/ds_McGinnis_8donor_SCMK_0.1bc_all.rds",
    "McGinnis_8donor_SCMK_0.05" = "./data-raw/benchmark_res/ds_McGinnis_8donor_SCMK/rds/ds_McGinnis_8donor_SCMK_0.05bc_all.rds",
    "McGinnis_8donor_SCMK_0.01" = "./data-raw/benchmark_res/ds_McGinnis_8donor_SCMK/rds/ds_McGinnis_8donor_SCMK_0.01bc_all.rds",
    "Stoeckius_4cellline" = "./data-raw/benchmark_res/Stoeckius_4cellline/rds/benchmark_Stoeckius_4cellline_bc_all.rds",
    "Stoeckius_8donor" = "./data-raw/benchmark_res/Stoeckius_8donor/rds/benchmark_Stoeckius_8donor_bc_all.rds",
    "Gaublomme_8donor" = "./data-raw/benchmark_res/Gaublomme/human_st/rds/benchmark_Gaublomme_bc_all.rds",
    "lung_cellline" = "./data-raw/benchmark_res/lung_cellline/rds/benchmark_lungcellline_bc_cbn.rds",
    "bal_batch1" = "./data-raw/benchmark_res/BAL/rds/benchmark_bal_b1_bc_cbn.rds",
    "bal_batch2" = "./data-raw/benchmark_res/BAL/rds/benchmark_bal_b2_bc_cbn.rds",
    "bal_batch3" = "./data-raw/benchmark_res/BAL/rds/benchmark_bal_b3_bc_cbn.rds",
    "Winkler_PDX" = "./data-raw/benchmark_res/Winkler_PDX/rds/Winkler_PDX_benchmarking_bc_cbn.rds"
)

bc_res <- lapply(bc_path, function(x) {
    readRDS(x)
})

show_methods <- names(bc_res$McGinnis_8donor_SCMK)

singlet_stats <- lapply(bc_res, function(x) {
    sapply(show_methods, function(cur_method) {
        y = x[[cur_method]]
        if(is.null(y)) return(c(NA,NA,NA)) else y$singlet_avg_stats
    })
})

singlet_dfs <- lapply(1:length(singlet_stats), function(i){
    x <- singlet_stats[[i]]
    x <- reshape2::melt(x)
    colnames(x) <- c("stats", "method", "value")
    x$dataset <- names(singlet_stats)[i]
    return(x)
})

singlet_df_cbn <- do.call(rbind, singlet_dfs)

singlet_df_cbn$dataset <- factor(singlet_df_cbn$dataset, levels= names(bc_res))
plot_df <- singlet_df_cbn[singlet_df_cbn$stats == "f_score", ]
g1 <- ggplot(plot_df, aes(x = dataset, y = value, fill = method)) +
    geom_bar(stat="identity",position ="dodge")

ggsave(paste0(plotPath, test_text, "_barv1.pdf"), g1, width = 10, height =2)


# Heatmap
library(pheatmap)
singlet_mtx <- sapply(singlet_stats, function(x) {
    x["f_score",][show_methods]
})

breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "singlet_summary_mtx", ".pdf"), width = 15, height = 5)
pheatmap(singlet_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()

# Rank based heatmap
library(matrixStats)

#singlet_ranks <- t(colRanks(-round(singlet_mtx, 3), ties.method = "dense"))
#singlet_ranks <- t(colRanks(-round(singlet_mtx, 3), ties.method = "max"))
colnames(singlet_ranks) <- colnames(singlet_mtx)
rownames(singlet_ranks) <- rownames(singlet_mtx)
breaksList <- seq(0,10,by=1)
pdf(paste0(plotPath, test_text, "singlet_rank_mtx", ".pdf"), width = 15, height = 3.7)
pheatmap(singlet_ranks, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.0f", number_color ='black', fontsize = 8, fontsize_number=8, color = rev(get_gradient_color('RdBu', length(breaksList))), breaks = breaksList)
dev.off()


# Only show simulation with increasing noise
show_sims <- c("sim1", "sim2", "sim3", "sim4", "sim5")
singlet_ranks <- singlet_ranks[,colnames(singlet_ranks) %in% show_sims | !grepl("sim", colnames(singlet_ranks))]
breaksList <- seq(0,10,by=1)
pdf(paste0(plotPath, test_text, "singlet_rank_mtx_selected", ".pdf"), width = 7, height = 3.7)
pheatmap(singlet_ranks, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.0f", number_color ='black', fontsize = 8, fontsize_number=8, color = rev(get_gradient_color('RdBu', length(breaksList))), breaks = breaksList)
dev.off()


