




library(deMULTIplex2)

rm(list=ls())



savePath <- "./data-raw/benchmark_res/simulation_model/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

test_text <- "robustness_sim1_"
sim_res <- readRDS(paste0(rdsPath, paste0("simulation_model_update_", "sim", 1, "_", "n.bc_", 5, "_"), "sim_res.rds"))

# test_text <- "robustness_sim2_"
# sim_res <- readRDS(paste0(rdsPath, paste0("simulation_model_update_", "sim", 2, "_", "n.bc_", 10, "_"), "sim_res.rds"))

# test_text <- "robustness_sim3_"
# sim_res <- readRDS(paste0(rdsPath, paste0("simulation_model_update_", "sim", 3, "_", "n.bc_", 10, "_"), "sim_res.rds"))

# test_text <- "robustness_sim4_"
# sim_res <- readRDS(paste0(rdsPath, paste0("simulation_model_update_", "sim", 4, "_", "n.bc_", 30, "_"), "sim_res.rds"))

test_text <- "robustness_sim5_"
sim_res <- readRDS(paste0(rdsPath, paste0("simulation_model_update_", "sim", 5, "_", "n.bc_", 30, "_"), "sim_res.rds"))


tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
cell_comp <- sapply(strsplit(rownames(tag_mtx), "[|]"), function(x) {
    if(length(x) >= 2) {
        x = unlist(strsplit(x[[2]], "_"))
        if(x[1] == x[3]) x[1] else NA
    } else NA
})
true_label[!is.na(cell_comp)] <- cell_comp[!is.na(cell_comp)]
names(true_label) =rownames(tag_mtx)
tag_mapping <- data.frame(tag = colnames(tag_mtx), true_label = colnames(tag_mtx))



# Benchmark on cosine cut
cos.cuts <- seq(0.1,.9,.1)
cos_res_list <- list()
for(cos.cut in cos.cuts) {
    print(cos.cut)
    cos_res_list[[paste0("cut_", cos.cut)]] <- benchmark_demultiplex2(tag_mtx, true_label,
                                        tag_mapping = tag_mapping,
                                        true.multiplet = "doublet",
                                        plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                        seed = 1,
                                        init.cos.cut = cos.cut ,
                                        converge.threshold = 1e-3,
                                        prob.cut = 0.5,
                                        max.cell.fit = 5000,
                                        max.iter = 10,
                                        min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                        max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                        residual.type = "rqr",
                                        plot.umap = "none",
                                        plot.diagnostics = F)

}
saveRDS(cos_res_list, paste0(rdsPath, test_text, "cos_cut.rds"))

cos_res_list <- readRDS(paste0(rdsPath, test_text, "cos_cut.rds"))
summary_df <- sapply(cos_res_list, function(x) {
    x$singlet_avg_stats
})
plot_df <- as.data.frame(reshape2::melt(t(summary_df[c("f_score"), ,drop=F])))
colnames(plot_df) <- c("init.cos","category", "value")
plot_df$init.cos <- as.numeric(gsub("cut_", "", names(cos_res_list)))
g1<-ggplot(plot_df, aes(init.cos, value)) +
    geom_point(size=1, color = "#005A9C")+
    geom_line(color = "#005A9C") +
    guides(fill = 'none') +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(0,1, by = .2))+
    xlab("Cosine cutoff") +
    ylab("F-score") +
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, 'bc_cos.pdf'), g1, width = 1.9, height = 1.5)


# Benchmark on max cell and iter
max.cell.cuts.list <- 10^seq(1,4,.5)
max.iter.list <- seq(2,10,2)
res_list <- list()
for(max.cell.cut in max.cell.cuts.list) {
    for(max.iter in max.iter.list) {
        cur_name <- paste0(max.cell.cut, "_", max.iter)
        print(cur_name)
        res_list[[cur_name]] <- benchmark_demultiplex2(tag_mtx, true_label,
                                                                          tag_mapping = tag_mapping,
                                                                          true.multiplet = "doublet",
                                                                          plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                                                          seed = 1,
                                                                          init.cos.cut = .5 ,
                                                                          converge.threshold = 1e-3,
                                                                          prob.cut = 0.5,
                                                                          max.cell.fit = max.cell.cut,
                                                                          max.iter = max.iter,
                                                                          min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                                                          max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                                                          residual.type = "rqr",
                                                                          plot.umap = "none",
                                                                          plot.diagnostics = F)
        res_list[[]]
    }
}
saveRDS(res_list, paste0(rdsPath, test_text, "maxcell_maxiter_cut.rds"))

res_list <- readRDS(paste0(rdsPath, test_text, "maxcell_maxiter_cut.rds"))
summary_df <- sapply(res_list, function(x) {
    x$singlet_avg_stats
})
plot_df <- as.data.frame(reshape2::melt(t(summary_df[c("f_score"), ,drop=F])))
colnames(plot_df) <- c("condition","category", "value")
plot_df$max.cell.fit <- factor(log10(as.numeric(sapply(strsplit(as.character(plot_df$condition), "_"), function(x)x[1]))))
plot_df$max.iter <- as.numeric(sapply(strsplit(as.character(plot_df$condition), "_"), function(x)x[2]))
g1<-ggplot(plot_df, aes(max.iter, value)) +
    geom_point(size=1, aes(color = max.cell.fit))+
    geom_line(aes(color = max.cell.fit)) +
    scale_color_manual(values = get_gradient_color("BlueGreenRed", cnum = length(levels(plot_df$max.cell.fit)))) +
    guides(fill = 'none') +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme_bw() +
    scale_y_continuous(limits = c(0, 1))+
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, 'maxcell_maxiter_cut.pdf'), g1, width = 2.4, height = 1.5)


# Heatmap
library(pheatmap)
plot_df$max.cell.fit <- factor(plot_df$max.cell.fit, levels = rev(seq(1,4,.5)))
acc_mtx <- reshape2::dcast(plot_df[,c("max.cell.fit", "max.iter", "value")], max.iter~max.cell.fit)
rownames(acc_mtx) <- as.character(acc_mtx$max.iter)
acc_mtx$max.iter <- NULL
breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "singlet_summary_mtx", ".pdf"), width = 2.2, height = 1.5)
pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', angle_col = 0, fontsize = 7, fontsize_number=7, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()

















# McGinnis 8-donor MULTI
savePath <- "./data-raw/benchmark_res/McGinnis_8donor_MULTI/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

test_text <- "robustness_McGinnis_8donor_MULTI_"

rna_mtx <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "rna_mtx.rds"))
tag_mtx <- readRDS(paste0(rdsPath, "McGinnis_8donor_MULTI_benchmarking_", "tag_mtx.rds"))
true_label <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "true_label.rds"))
tag_mapping <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "tag_mapping.rds"))


# Benchmark on cosine cut
cos.cuts <- seq(0.1,.9,.1)
cos_res_list <- list()
for(cos.cut in cos.cuts) {
    print(cos.cut)
    cos_res_list[[paste0("cut_", cos.cut)]] <- benchmark_demultiplex2(tag_mtx, true_label,
                                                                      tag_mapping = tag_mapping,
                                                                      true.multiplet = "doublet",
                                                                      plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                                                      seed = 1,
                                                                      init.cos.cut = cos.cut ,
                                                                      converge.threshold = 1e-3,
                                                                      prob.cut = 0.5,
                                                                      max.cell.fit = 5000,
                                                                      max.iter = 30,
                                                                      min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                                                      max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                                                      residual.type = "rqr",
                                                                      plot.umap = "none",
                                                                      plot.diagnostics = F)

}
saveRDS(cos_res_list, paste0(rdsPath, test_text, "cos_cut.rds"))

cos_res_list <- readRDS(paste0(rdsPath, test_text, "cos_cut.rds"))
summary_df <- sapply(cos_res_list, function(x) {
    x$singlet_avg_stats
})
plot_df <- as.data.frame(reshape2::melt(t(summary_df[c("f_score"), ,drop=F])))
colnames(plot_df) <- c("init.cos","category", "value")
plot_df$init.cos <- as.numeric(gsub("cut_", "", names(cos_res_list)))
g1<-ggplot(plot_df, aes(init.cos, value)) +
    geom_point(size=1, color = "#005A9C")+
    geom_line(color = "#005A9C") +
    guides(fill = 'none') +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(0,1, by = .2))+
    xlab("Cosine cutoff") +
    ylab("F-score") +
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, 'bc_cos.pdf'), g1, width = 1.9, height = 1.5)



# Benchmark on max cell and iter
max.cell.cuts.list <- 10^seq(1,4,.5)
max.iter.list <- seq(2,10,2)
res_list <- list()
for(max.cell.cut in max.cell.cuts.list) {
    for(max.iter in max.iter.list) {
        cur_name <- paste0(max.cell.cut, "_", max.iter)
        print(cur_name)
        res_list[[cur_name]] <- benchmark_demultiplex2(tag_mtx, true_label,
                                                       tag_mapping = tag_mapping,
                                                       true.multiplet = "doublet",
                                                       plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                                       seed = 1,
                                                       init.cos.cut = .5 ,
                                                       converge.threshold = 1e-3,
                                                       prob.cut = 0.5,
                                                       max.cell.fit = max.cell.cut,
                                                       max.iter = max.iter,
                                                       min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                                       max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                                       residual.type = "rqr",
                                                       plot.umap = "none",
                                                       plot.diagnostics = F)
        res_list[[]]
    }
}
saveRDS(res_list, paste0(rdsPath, test_text, "maxcell_maxiter_cut.rds"))

res_list <- readRDS(paste0(rdsPath, test_text, "maxcell_maxiter_cut.rds"))
summary_df <- sapply(res_list, function(x) {
    x$singlet_avg_stats
})
plot_df <- as.data.frame(reshape2::melt(t(summary_df[c("f_score"), ,drop=F])))
colnames(plot_df) <- c("condition","category", "value")
plot_df$max.cell.fit <- factor(log10(as.numeric(sapply(strsplit(as.character(plot_df$condition), "_"), function(x)x[1]))))
plot_df$max.iter <- as.numeric(sapply(strsplit(as.character(plot_df$condition), "_"), function(x)x[2]))
g1<-ggplot(plot_df, aes(max.iter, value)) +
    geom_point(size=1, aes(color = max.cell.fit))+
    geom_line(aes(color = max.cell.fit)) +
    scale_color_manual(values = get_gradient_color("BlueGreenRed", cnum = length(levels(plot_df$max.cell.fit)))) +
    guides(fill = 'none') +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme_bw() +
    scale_y_continuous(limits = c(0, 1))+
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, 'maxcell_maxiter_cut.pdf'), g1, width = 2.4, height = 1.5)



# Heatmap
library(pheatmap)
plot_df$max.cell.fit <- factor(plot_df$max.cell.fit, levels = rev(seq(1,4,.5)))
acc_mtx <- reshape2::dcast(plot_df[,c("max.cell.fit", "max.iter", "value")], max.iter~max.cell.fit)
rownames(acc_mtx) <- as.character(acc_mtx$max.iter)
acc_mtx$max.iter <- NULL
breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "singlet_summary_mtx", ".pdf"), width = 2.3, height = 1.5)
pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', angle_col = 0, fontsize = 7, fontsize_number=7, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()







# McGinnis 8-donor SCMK
savePath <- "./data-raw/benchmark_res/McGinnis_8donor_SCMK/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

test_text <- "robustness_McGinnis_8donor_SCMK_"

rna_mtx <- readRDS(paste0(rdsPath,  "McGinnis_8donor_SCMK_benchmarking_", "rna_mtx.rds"))
tag_mtx <- readRDS(paste0(rdsPath, "McGinnis_8donor_SCMK_benchmarking_", "tag_mtx.rds"))
true_label <- readRDS(paste0(rdsPath,  "McGinnis_8donor_SCMK_benchmarking_", "true_label.rds"))
tag_mapping <- readRDS(paste0(rdsPath,  "McGinnis_8donor_SCMK_benchmarking_", "tag_mapping.rds"))


# Benchmark on cosine cut
cos.cuts <- seq(0.1,.9,.1)
cos_res_list <- list()
for(cos.cut in cos.cuts) {
    print(cos.cut)
    cos_res_list[[paste0("cut_", cos.cut)]] <- benchmark_demultiplex2(tag_mtx, true_label,
                                                                      tag_mapping = tag_mapping,
                                                                      true.multiplet = "doublet",
                                                                      plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                                                      seed = 1,
                                                                      init.cos.cut = cos.cut ,
                                                                      converge.threshold = 1e-3,
                                                                      prob.cut = 0.5,
                                                                      max.cell.fit = 5000,
                                                                      max.iter = 30,
                                                                      min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                                                      max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                                                      residual.type = "rqr",
                                                                      plot.umap = "none",
                                                                      plot.diagnostics = F)

}
saveRDS(cos_res_list, paste0(rdsPath, test_text, "cos_cut.rds"))

cos_res_list <- readRDS(paste0(rdsPath, test_text, "cos_cut.rds"))
summary_df <- sapply(cos_res_list, function(x) {
    x$singlet_avg_stats
})
plot_df <- as.data.frame(reshape2::melt(t(summary_df[c("f_score"), ,drop=F])))
colnames(plot_df) <- c("init.cos","category", "value")
plot_df$init.cos <- as.numeric(gsub("cut_", "", names(cos_res_list)))
g1<-ggplot(plot_df, aes(init.cos, value)) +
    geom_point(size=1, color = "#005A9C")+
    geom_line(color = "#005A9C") +
    guides(fill = 'none') +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(0,1, by = .2))+
    xlab("Cosine cutoff") +
    ylab("F-score") +
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, 'bc_cos.pdf'), g1, width = 1.9, height = 1.5)



# Benchmark on max cell and iter
max.cell.cuts.list <- 10^seq(1,4,.5)
max.iter.list <- seq(2,10,2)
res_list <- list()
for(max.cell.cut in max.cell.cuts.list) {
    for(max.iter in max.iter.list) {
        cur_name <- paste0(max.cell.cut, "_", max.iter)
        print(cur_name)
        res_list[[cur_name]] <- benchmark_demultiplex2(tag_mtx, true_label,
                                                       tag_mapping = tag_mapping,
                                                       true.multiplet = "doublet",
                                                       plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                                       seed = 1,
                                                       init.cos.cut = .5 ,
                                                       converge.threshold = 1e-3,
                                                       prob.cut = 0.5,
                                                       max.cell.fit = max.cell.cut,
                                                       max.iter = max.iter,
                                                       min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                                       max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                                       residual.type = "rqr",
                                                       plot.umap = "none",
                                                       plot.diagnostics = F)
        res_list[[]]
    }
}
saveRDS(res_list, paste0(rdsPath, test_text, "maxcell_maxiter_cut.rds"))

res_list <- readRDS(paste0(rdsPath, test_text, "maxcell_maxiter_cut.rds"))
summary_df <- sapply(res_list, function(x) {
    x$singlet_avg_stats
})
plot_df <- as.data.frame(reshape2::melt(t(summary_df[c("f_score"), ,drop=F])))
colnames(plot_df) <- c("condition","category", "value")
plot_df$max.cell.fit <- factor(log10(as.numeric(sapply(strsplit(as.character(plot_df$condition), "_"), function(x)x[1]))))
plot_df$max.iter <- as.numeric(sapply(strsplit(as.character(plot_df$condition), "_"), function(x)x[2]))
g1<-ggplot(plot_df, aes(max.iter, value)) +
    geom_point(size=1, aes(color = max.cell.fit))+
    geom_line(aes(color = max.cell.fit)) +
    scale_color_manual(values = get_gradient_color("BlueGreenRed", cnum = length(levels(plot_df$max.cell.fit)))) +
    guides(fill = 'none') +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme_bw() +
    scale_y_continuous(limits = c(0, 1))+
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, 'maxcell_maxiter_cut.pdf'), g1, width = 2.4, height = 1.5)



# Heatmap
library(pheatmap)
plot_df$max.cell.fit <- factor(plot_df$max.cell.fit, levels = rev(seq(1,4,.5)))
acc_mtx <- reshape2::dcast(plot_df[,c("max.cell.fit", "max.iter", "value")], max.iter~max.cell.fit)
rownames(acc_mtx) <- as.character(acc_mtx$max.iter)
acc_mtx$max.iter <- NULL
breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "singlet_summary_mtx", ".pdf"), width = 2.3, height = 1.5)
pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', angle_col = 0, fontsize = 7, fontsize_number=7, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()












# Revision, check the robustness of fit of glm with different number of downsampled cells
test_text <- "cellds_McGinnis_8donor_MULTI_"
savePath <- "./data-raw/benchmark_res/McGinnis_8donor_MULTI/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

rna_mtx <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "rna_mtx.rds"))
tag_mtx <- readRDS(paste0(rdsPath, "McGinnis_8donor_MULTI_benchmarking_", "tag_mtx.rds"))
true_label <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "true_label.rds"))
tag_mapping <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "tag_mapping.rds"))



max.cell.cuts.list <- 10^seq(1,4,.5)
res_list <- list()

bc <- "Bar12"
bc_donor <- tag_mapping$true_label[match(bc, tag_mapping$tag)]
df <- data.frame(bc.umi = tag_mtx[, bc], tt.umi = rowSums(tag_mtx), mem = (true_label == bc_donor)*1)
df <- df[df$tt.umi!=0,]
nsample = 100
library(MASS)
set.seed(1)
coef_list <- list()
post0_list <- list()
pi0 <- sum(df$mem == 0) / length(df$mem)
pi1 <- sum(df$mem == 1) / length(df$mem)
for(max.cell.fit in max.cell.cuts.list) {
    print(max.cell.fit)
    neg_num <- sum(df$mem ==0)
    pos_num <- sum(df$mem ==1)
    fit_coefs <- list()
    post0_cells <- list()
    post1_cells <- list()
    for(i in 1:nsample){
        df_fit0 = df[sample(which(df$mem ==0), min(neg_num, max.cell.fit)), ]
        fit0 = tryCatch(glm.nb("bc.umi~log(tt.umi)", df_fit0, link = log), error = function(e) return(NULL))

        df_fit1 = df[sample(which(df$mem ==1), min(pos_num, max.cell.fit)), ]
        fit1 = tryCatch(glm.nb("(tt.umi - bc.umi)~log(tt.umi)", df_fit1, link = log), error = function(e) return(NULL))

        if(is.null(fit0) || is.null(fit1)) {
            fit_coefs[[i]] <- c(NA,NA,NA,NA)
        } else {
            pred0 <- predict(fit0, df, type = "response", se.fit=FALSE)
            prob0 <- dnbinom(df$bc.umi, mu=pred0, size=fit0$theta)

            pred1 <- predict(fit1, df, type = "response", se.fit=FALSE)
            prob1 <- dnbinom(df$tt.umi - df$bc.umi, mu=pred1, size=fit1$theta)

            neg_c <- df$bc.umi < pred0
            prob0[neg_c] = 1 # Seems work best based on some super noisy test data
            pos_c <- (df$tt.umi - df$bc.umi) < pred1
            prob1[pos_c] = dnbinom(ceiling(pred1[pos_c]), mu=pred1[pos_c], size=fit1$theta)

            comp0 <- pi0 * prob0
            comp1 <- pi1 * prob1
            comp.sum <- comp0 + comp1

            post0 <- comp0/comp.sum
            post1 <- comp1/comp.sum

            fit_coefs[[i]] <- c(
                fit0_b0 = fit0$coefficients[1],
                fit0_b1 = fit0$coefficients[2],
                fit1_b0 = fit1$coefficients[1],
                fit1_b1 = fit1$coefficients[2]
            )
            post0_cells[[as.character(i)]] <- post0
            post1_cells[[as.character(i)]] <- post1
        }
    }
    fit_res <- as.data.frame(do.call(rbind, fit_coefs))
    colnames(fit_res) <- c("fit0_b0", "fit0_b1", "fit1_b0", "fit1_b1")
    fit_res$cell_num <- round(log10(max.cell.fit), 1)
    coef_list[[paste0("ds_", round(log10(max.cell.fit), 1))]] <- fit_res

    post0_list[[paste0("ds_", round(log10(max.cell.fit), 1))]] <- as.data.frame(do.call(cbind, post0_cells))
    post1_list[[paste0("ds_", round(log10(max.cell.fit), 1))]] <- as.data.frame(do.call(cbind, post1_cells))
}

saveRDS(coef_list, paste0(rdsPath,test_text, bc,  "_coef_list.rds"))
saveRDS(post0_list, paste0(rdsPath,test_text, bc,  "_post0_list.rds"))
saveRDS(post1_list, paste0(rdsPath,test_text, bc,  "_post1_list.rds"))

coef_list <- readRDS(paste0(rdsPath,test_text, bc,  "_coef_list.rds"))
coef_cbn <- do.call(rbind, coef_list)

coef_cbn$cell_num <- factor(coef_cbn$cell_num)
g1 <- ggplot(coef_cbn, aes(x = cell_num, y = fit0_b0, fill = cell_num)) +
    geom_boxplot(lwd=.2, outlier.shape = NA) +
    scale_fill_manual(values = get_gradient_color("BlueGreenRed", cnum = length(levels(coef_cbn$cell_num)))) +
    theme_bw()+
    guides(fill = 'none', color = 'none') +
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
ggsave(paste0(plotPath, test_text, bc, "_fit0_", "b0_", "maxcell.pdf"), g1, width = 1.7, height = 1.5)


g1 <- ggplot(coef_cbn, aes(x = cell_num, y = fit0_b1, fill = cell_num)) +
    geom_boxplot(lwd=.2, outlier.shape = NA) +
    scale_fill_manual(values = get_gradient_color("BlueGreenRed", cnum = length(levels(coef_cbn$cell_num)))) +
    theme_bw()+
    guides(fill = 'none', color = 'none') +
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
ggsave(paste0(plotPath, test_text, bc, "_fit0_", "b1_", "maxcell.pdf"), g1, width = 1.7, height = 1.5)


g1 <- ggplot(coef_cbn, aes(x = cell_num, y = fit1_b0, fill = cell_num)) +
    geom_boxplot(lwd=.2, outlier.shape = NA) +
    scale_fill_manual(values = get_gradient_color("BlueGreenRed", cnum = length(levels(coef_cbn$cell_num)))) +
    theme_bw()+
    guides(fill = 'none', color = 'none') +
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
ggsave(paste0(plotPath, test_text, bc, "_fit1_", "b0_", "maxcell.pdf"), g1, width = 1.7, height = 1.5)


g1 <- ggplot(coef_cbn, aes(x = cell_num, y = fit1_b1, fill = cell_num)) +
    geom_boxplot(lwd=.2, outlier.shape = NA) +
    scale_fill_manual(values = get_gradient_color("BlueGreenRed", cnum = length(levels(coef_cbn$cell_num)))) +
    theme_bw()+
    guides(fill = 'none', color = 'none') +
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
ggsave(paste0(plotPath, test_text, bc, "_fit1_", "b1_", "maxcell.pdf"), g1, width = 1.7, height = 1.5)




post0_list <- readRDS(paste0(rdsPath,test_text, bc,  "_post0_list.rds"))
post0_sd <- sapply(post0_list, function(x) apply(x, 1, sd))
colnames(post0_sd) <- gsub("ds_", "", colnames(post0_sd))
breaksList <- seq(0,.5,.05)
pdf(paste0(plotPath, test_text, bc, "post0_sd", ".pdf"), width = 4, height = 4)
#png(paste0(plotPath, test_text, bc, "post0_sd", ".png"), width = 4*300, height = 4*300, res = 300)
pheatmap(post0_sd, cluster_rows = T, cluster_cols = F, display_numbers =F, clustering_method = "ward.D2", border_color=NA,  angle_col = 0, fontsize = 8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()









