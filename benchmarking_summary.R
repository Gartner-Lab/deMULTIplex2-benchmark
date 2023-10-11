

# Read in all summary results

library(deMULTIplex2)
savePath <- "./data-raw/benchmark_res/summary_res/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)


test_text<- "benchmark_summary_"

bc_path <- list(
    "McGinnis_8donor_MULTI" = "./data-raw/benchmark_res/McGinnis_8donor_MULTI/rds/McGinnis_8donor_MULTI_benchmarking_bc_all.rds",
    "McGinnis_8donor_SCMK" = "./data-raw/benchmark_res/McGinnis_8donor_SCMK/rds/McGinnis_8donor_SCMK_benchmarking_bc_all.rds",
    "Stoeckius_4cellline" = "./data-raw/benchmark_res/Stoeckius_4cellline/rds/benchmark_Stoeckius_4cellline_bc_all.rds",
    "Stoeckius_8donor" = "./data-raw/benchmark_res/Stoeckius_8donor/rds/benchmark_Stoeckius_8donor_bc_all.rds",
    "Gaublomme_8donor" = "./data-raw/benchmark_res/Gaublomme/human_st/rds/benchmark_Gaublomme_bc_all.rds",
    "lung_cellline" = "./data-raw/benchmark_res/lung_cellline/rds/benchmark_lungcellline_bc_cbn.rds",
    "bal_batch1" = "./data-raw/benchmark_res/BAL/rds/benchmark_bal_b1_bc_cbn.rds",
    "bal_batch2" = "./data-raw/benchmark_res/BAL/rds/benchmark_bal_b2_bc_cbn.rds",
    "bal_batch3" = "./data-raw/benchmark_res/BAL/rds/benchmark_bal_b3_bc_cbn.rds"
)

bc_res <- lapply(bc_path, function(x) {
    readRDS(x)
})

singlet_stats <- lapply(bc_res, function(x) {
    sapply(x, function(y) {
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
show_methods <- names(bc_res$McGinnis_8donor_MULTI)
singlet_mtx <- sapply(singlet_stats, function(x) {
    x["f_score",][show_methods]
})

breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "singlet_summary_mtx", ".pdf"), width = 4, height = 3.5)
pheatmap(singlet_mtx, cluster_rows = F, cluster_cols = T, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()



doublet_stats <- lapply(bc_res, function(x) {
    sapply(x, function(y) {
        if(is.null(y)) return(c(NA,NA,NA)) else y$doublet_avg_stats
    })
})

doublet_mtx <- sapply(doublet_stats, function(x) {
    x["recall",][show_methods]
})
doublet_mtx <- doublet_mtx[rownames(doublet_mtx)!= "GMM_Demux",]
breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "doublet_summary_mtx", ".pdf"), width = 4, height = 3.5)
pheatmap(doublet_mtx, cluster_rows = F, cluster_cols = T, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()



doublet_mtx <- sapply(doublet_stats, function(x) {
    x["doublet_called_singlet",][show_methods]
})
doublet_mtx <- doublet_mtx[rownames(doublet_mtx)!= "GMM_Demux",]
breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "doublet_called_singlet_summary_mtx", ".pdf"), width = 4, height = 3.5)
pheatmap(doublet_mtx, cluster_rows = F, cluster_cols = T, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()



# Demultiplex2 coefficients

dm2_coefs <- lapply(1:length(bc_res), function(i) {
    x = bc_res[[i]]
    df <- as.data.frame(do.call(rbind, x$deMULTIplex2$res$coefs))
    df$dataset <- names(bc_res)[i]
    return(df)
})

dm2_coef_tbl <- do.call(rbind, dm2_coefs)
colnames(dm2_coef_tbl) <- c("b0", "b1", "theta", "dataset")
dm2_coef_tbl_plot <- reshape2::melt(dm2_coef_tbl, id.var = "dataset")
g1 <- ggplot(dm2_coef_tbl_plot, aes(dataset, value)) +
    geom_boxplot(lwd=.2, fill = "dodger blue", outlier.shape = NA) +
    #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1),dotsize = .9,stroke = 0.2)+
    facet_wrap(~variable, scales = 'free')+
    guides(fill = 'none', color = 'none') +
    theme_bw() +
    ylab(NULL)+
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          #strip.text = element_blank(),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, 'bc_fitcoefts_boxplot.pdf'), g1, width =5, height = 2.5)



