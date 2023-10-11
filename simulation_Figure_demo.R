

library(deMULTIplex2)
library(ggplot2)
library(ggrastr)
savePath <- "./data-raw/benchmark_res/simulation_b1/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

test_text <- "simulation_benchmarking_Fig1_"

sim_res <- simulateTags(n.cell = 1000,
                        n.bc = 2,
                        seed = 2023,
                        nb.theta = 10,
                        min.size.log = 2,
                        max.size.log = 7,
                        min.ambient.log = 3,
                        max.ambient.log = 4,
                        b0 = -3,
                        doublet.rate = 0,
                        separate.sim = T,
                        return.all = T)

tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
names(true_label) = rownames(tag_mtx)
tag_mapping <- data.frame(tag = colnames(tag_mtx), true_label = colnames(tag_mtx))

bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "Doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 1000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(sim_res, paste0(rdsPath, test_text, "sim_res.rds"))
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))

# Plot true label
sim_res <- readRDS(paste0(rdsPath, test_text, "sim_res.rds"))
bc_demux2 <- readRDS(paste0(rdsPath, test_text, "bc_demux2.rds"))
tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
names(true_label) = rownames(tag_mtx)


bc = 'bc2'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
plot_df$true_umi <- sim_res$cell.true.umi.mtx[,bc]
plot_df$cell_contam_umi <- sim_res$cell.contam.umi.mtx[,bc]
plot_df$ambient_contam_umi <- sim_res$bead.contam.umi.mtx[,bc]
map =c("log(tt.umi)", "log(bc.umi)", "true_label")


# g_crazyshow = ggplot(plot_df) +
#     geom_point_rast(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = 1) +
#     geom_line(aes_string(map[1], 'log(pred0)'),color = 'grey') +
#     geom_point_rast(aes_string(map[1], 'log(true_umi)'), color = "green", stroke = 0, size = 1) +
#     geom_point_rast(aes_string(map[1], 'log(cell_contam_umi)'), color = "blue", stroke = 0, size = 1) +
#     geom_point_rast(aes_string(map[1], 'log(ambient_contam_umi)'), color = "orange", stroke = 0, size = 1) +
#     geom_abline(slope = 1,intercept = 0) +
#     ggtitle(bc) +
#     theme_bw() +
#     labs(color = gsub("_","\n",map[3]))

use_color <- c(
    'bc1' = '#BE1E2D',
    'bc2' = '#00A651'
)
map =c("log(tt.umi)", "log(bc.umi)", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_true_label.pdf"), g1, width = 2.2, height =1.5)

map =c("log(tt.umi)", "log(bc.umi)", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    #guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, bc, "_post1.pdf"), g1, width = 2.2, height =1.5)


map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "true_label")
g2 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, bc, "_true_label_reverse.pdf"), g2, width = 2.2, height =1.5)

map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "true_label")
g2 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, bc, "_post1_reverse.pdf"), g2, width = 2.2, height =1.5)






# Simulate data with 0 ambient contaminaton
test_text <- "simulation_benchmarking_Fig1de_"
sim_res <- simulateTags(n.cell = 1000,
                        n.bc = 5,
                        seed = 2023,
                        nb.theta = 10,
                        min.size.log = 2,
                        max.size.log = 7,
                        min.ambient.log = 0,
                        max.ambient.log = 0,
                        b0 = -3,
                        doublet.rate = 0,
                        separate.sim = T,
                        return.all = T)

tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
names(true_label) = rownames(tag_mtx)
tag_mapping <- data.frame(tag = colnames(tag_mtx), true_label = colnames(tag_mtx))

bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "Doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 1000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(sim_res, paste0(rdsPath, test_text, "sim_res.rds"))
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))


sim_res <- readRDS(paste0(rdsPath, test_text, "sim_res.rds"))
bc_demux2 <- readRDS(paste0(rdsPath, test_text, "bc_demux2.rds"))
bc = 'bc1'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
plot_df$true_umi <- sim_res$cell.true.umi.mtx[,bc]
plot_df$cell_contam_umi <- sim_res$cell.contam.umi.mtx[,bc]
plot_df$ambient_contam_umi <- sim_res$bead.contam.umi.mtx[,bc]
map =c("log(tt.umi)", "log(bc.umi)", "true_label")

use_color <- get_factor_color(colnames(tag_mtx))
names(use_color) <- colnames(tag_mtx)
map =c("log(tt.umi)", "log(bc.umi)", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'grey') +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_true_label.pdf"), g1, width = 2.2, height =1.5)










# Simulate data with only ambient contaminaton (no cell attaching)
test_text <- "simulation_benchmarking_Fig1deonlycellcontam_"
sim_res <- simulateTags(n.cell = 1000,
                        n.bc = 5,
                        seed = 2023,
                        nb.theta = 10,
                        min.size.log = 2,
                        max.size.log = 7,
                        min.ambient.log = 3,
                        max.ambient.log = 4,
                        b0 = -3,
                        doublet.rate = 0,
                        separate.sim = T,
                        cell.contam = T,
                        ambient.contam = F,
                        return.all = T)

tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
names(true_label) = rownames(tag_mtx)
tag_mapping <- data.frame(tag = colnames(tag_mtx), true_label = colnames(tag_mtx))

bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "Doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 1000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(sim_res, paste0(rdsPath, test_text, "sim_res.rds"))
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))

bc = 'bc1'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
plot_df$true_umi <- sim_res$cell.true.umi.mtx[,bc]
plot_df$cell_contam_umi <- sim_res$cell.contam.umi.mtx[,bc]
plot_df$ambient_contam_umi <- sim_res$bead.contam.umi.mtx[,bc]
map =c("log(tt.umi)", "log(bc.umi)", "true_label")

use_color <- get_factor_color(colnames(tag_mtx))
names(use_color) <- colnames(tag_mtx)
map =c("log(tt.umi)", "log(bc.umi)", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'grey') +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(
        text=element_text(size=8),
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
ggsave(paste0(plotPath, test_text, "_true_label.pdf"), g1, width = 2.2, height =1.5)



test_text <- "simulation_benchmarking_Fig1deonlyambient_"
sim_res <- simulateTags(n.cell = 1000,
                        n.bc = 5,
                        seed = 2023,
                        nb.theta = 10,
                        min.size.log = 2,
                        max.size.log = 7,
                        min.ambient.log = 3,
                        max.ambient.log = 4,
                        b0 = -3,
                        doublet.rate = 0,
                        separate.sim = T,
                        cell.contam = F,
                        ambient.contam = T,
                        return.all = T)

tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
names(true_label) = rownames(tag_mtx)
tag_mapping <- data.frame(tag = colnames(tag_mtx), true_label = colnames(tag_mtx))

bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "Doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 1000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(sim_res, paste0(rdsPath, test_text, "sim_res.rds"))
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))

bc = 'bc1'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
plot_df$true_umi <- sim_res$cell.true.umi.mtx[,bc]
plot_df$cell_contam_umi <- sim_res$cell.contam.umi.mtx[,bc]
plot_df$ambient_contam_umi <- sim_res$bead.contam.umi.mtx[,bc]
map =c("log(tt.umi)", "log(bc.umi)", "true_label")

use_color <- get_factor_color(colnames(tag_mtx))
names(use_color) <- colnames(tag_mtx)
map =c("log(tt.umi)", "log(bc.umi)", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'grey') +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(
        text=element_text(size=8),
        legend.text=element_text(size=8),
        legend.title = element_blank(),
        axis.text = element_text(size=8),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_true_label.pdf"), g1, width = 2, height =1.3)





# Some ambient contamination

test_text <- "simulation_benchmarking_Fig1someambientsomecell_"
sim_res <- simulateTags(n.cell = 1000,
                        n.bc = 5,
                        seed = 2023,
                        nb.theta = 10,
                        min.size.log = 2,
                        max.size.log = 7,
                        min.ambient.log = 3,
                        max.ambient.log = 4,
                        b0 = -3,
                        doublet.rate = 0,
                        separate.sim = T,
                        cell.contam = T,
                        ambient.contam = F,
                        return.all = T)

tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
names(true_label) = rownames(tag_mtx)
tag_mapping <- data.frame(tag = colnames(tag_mtx), true_label = colnames(tag_mtx))

bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "Doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 1000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(sim_res, paste0(rdsPath, test_text, "sim_res.rds"))
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))

bc = 'bc1'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
plot_df$true_umi <- sim_res$cell.true.umi.mtx[,bc]
plot_df$cell_contam_umi <- sim_res$cell.contam.umi.mtx[,bc]
plot_df$ambient_contam_umi <- sim_res$bead.contam.umi.mtx[,bc]
map =c("log(tt.umi)", "log(bc.umi)", "true_label")

use_color <- get_factor_color(colnames(tag_mtx))
names(use_color) <- colnames(tag_mtx)
map =c("log(tt.umi)", "log(bc.umi)", "true_label")
map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "true_label")

map =c("log(tt.umi)", "log(bc.umi)", "cos.umi")
use_color = get_gradient_color("BlueGreenRed")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'grey') +
    ggtitle(bc) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    #scale_color_manual(values = use_color) +
    scale_color_gradientn(colors = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(
        text=element_text(size=8),
        legend.text=element_text(size=8),
        legend.title = element_blank(),
        axis.text = element_text(size=8),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_reg_cos.pdf"), g1, width = 2, height =2)















test_text <- "simulation_benchmarking_Fig2sigmoid_"
sim_res <- simulateTags(n.cell = 1000,
                        n.bc = 5,
                        seed = 2023,
                        nb.theta = 10,
                        min.size.log = log(10),
                        max.size.log = log(1000),
                        min.ambient.log = log(10),
                        max.ambient.log = log(100),
                        b0 = -4.5,
                        doublet.rate = 0.1,
                        separate.sim = T,
                        cell.contam = T,
                        ambient.contam = T,
                        return.all = T)

tag_mtx <- sim_res$final.umi.mtx
true_label = sapply(strsplit(rownames(tag_mtx), "_"), function(x){x[1]})
names(true_label) = rownames(tag_mtx)
tag_mapping <- data.frame(tag = colnames(tag_mtx), true_label = colnames(tag_mtx))

bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "Doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 1000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(sim_res, paste0(rdsPath, test_text, "sim_res.rds"))
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))

res <- bc_demux2$res
call_label <- res$assign_table$barcode_assign
call_label[res$assign_table$barcode_count == 0] = "Negative"
call_label[res$assign_table$barcode_count > 1] = "Multiplet"
names(call_label) <- rownames(res$assign_table)
true.multiplet = 'doublet'
call.multiplet = "Multiplet"
doublet_called_singlet_rate <- sum(true_label == true.multiplet & call_label %in% tag_mapping$tag, na.rm=T) / sum(true_label == true.multiplet, na.rm=T)
doublet_called_negative_rate <- sum(true_label == true.multiplet & !call_label %in% c(tag_mapping$tag, call.multiplet), na.rm=T) / sum(true_label == true.multiplet, na.rm=T)
doublet_called_doublet_rate <- sum(true_label == true.multiplet & call_label == call.multiplet, na.rm=T) / sum(true_label == true.multiplet, na.rm=T)
doublet_stats <- c(recall = doublet_called_doublet_rate , doublet_called_singlet = doublet_called_singlet_rate, doublet_called_negative = doublet_called_negative_rate)


bc = 'bc3'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
plot_df$is_positive <- ifelse(plot_df$true_label==bc, 'positive', ifelse(plot_df$true_label== 'doublet', ifelse(grepl(bc, rownames(plot_df)), 'positive(d)', 'negative'), 'negative'))

use_color <- c('positive' = '#e41a1c', 'negative' = '#4d4d4d', 'positive(d)' = '#fdae61')
map = c("log(bc.umi)", "cos.umi", "is_positive")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    #scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(text=element_text(size=8),
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

g1 <- ggMarginal(g1, groupColour = TRUE, groupFill = TRUE)

ggsave(paste0(plotPath, test_text, "_cosine_is_pos.pdf"), g1, width = 3, height =2)


map = c("log(bc.umi)", "cos.umi", "cos.umi")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(text=element_text(size=8),
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

g1 <- ggMarginal(g1, type = "histogram")

ggsave(paste0(plotPath, test_text, "_cosine_cos.pdf"), g1, width = 2.7, height =2)


use_color <- c('positive' = '#e41a1c', 'negative' = '#4d4d4d', 'positive(d)' = '#fdae61')
map =c("log(tt.umi)", "log(bc.umi)", "is_positive")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    geom_abline(slope = 1,intercept = 0) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(text=element_text(size=8),
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
g1 <- ggMarginal(g1, groupColour = TRUE, groupFill = TRUE, margins = 'y')

ggsave(paste0(plotPath, test_text, "_count_is_pos.pdf"), g1, width = 2.5, height =1.7)



use_color <- c('positive' = '#e41a1c', 'negative' = '#4d4d4d', 'positive(d)' = '#fdae61')
map =c("log(tt.umi)", "rqr", "is_positive")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .5) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(text=element_text(size=8),
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
g1 <- ggMarginal(g1, groupColour = TRUE, groupFill = TRUE, margins = 'y')

ggsave(paste0(plotPath, test_text, "_rqr_is_pos.pdf"), g1, width = 2.5, height =1.7)


p1 = ggplot(plot_df[plot_df$post1 < .5,], aes_string(sample='rqr'))+
    stat_qq(size = .5, stroke = 0) +
    stat_qq_line(alpha = .3) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(text=element_text(size=8),
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
ggsave(paste0(plotPath, test_text, "_rqr_qq_predicted.pdf"), p1, width = 1, height =.75)

p1 = ggplot(plot_df[plot_df$is_positive == 'negative',], aes_string(sample='rqr'))+
    stat_qq(size = .5, stroke = 0) +
    stat_qq_line(alpha = .3) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(text=element_text(size=8),
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
ggsave(paste0(plotPath, test_text, "_rqr_qq_trueneg.pdf"), p1, width = 1, height =.75)



