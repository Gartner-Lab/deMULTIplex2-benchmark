

library(deMULTIplex2)
library(ggrastr)
library(ggplot2)
savePath <- "./data-raw/benchmark_res/simulation_model/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)


library(truncnorm)
library(ggExtra)
library(ggrastr)
test_text <- "simulation_model_update_"


# Sim 1
set.seed(2023)
n.bc = 5

cell.staining.meanlog = rnorm(n.bc, mean = 7, sd = .1)
cell.staining.sdlog = .7

n.cell.min = 300
n.cell.max = 500
n.cell.logmean = log(mean(c(n.cell.min, n.cell.max)))
n.cell.logsd = .1
ncell <- round(rlnorm(n.bc,n.cell.logmean,n.cell.logsd))
ncell[ncell < n.cell.min] = n.cell.min
ncell[ncell > n.cell.max] = n.cell.max
hist(ncell)

nb.theta = 2
b0 = -5
ambient.meanlog = .1
doublet.rate = 0.02
dropout.lambda = 2

cur_text <- paste0(test_text, "sim1_", "n.bc_", n.bc, "_")



# Sim 2
set.seed(2023)
n.bc = 10

cell.staining.meanlog = rnorm(n.bc, mean = 5, sd = .1)
cell.staining.sdlog = 1

n.cell.min = 50
n.cell.max = 2000
n.cell.logmean = log(mean(c(n.cell.min, n.cell.max)))
n.cell.logsd = .8
ncell <- round(rlnorm(n.bc,n.cell.logmean,n.cell.logsd))
ncell[ncell < n.cell.min] = n.cell.min
ncell[ncell > n.cell.max] = n.cell.max
hist(ncell)

nb.theta = 7
b0 = -5
ambient.meanlog = .1
doublet.rate = 0.02
dropout.lambda = 2

cur_text <- paste0(test_text, "sim2_", "n.bc_", n.bc, "_")



# Sim 3
set.seed(2023)
n.bc = 10

cell.staining.meanlog = rnorm(n.bc, mean = 5, sd = .1)
cell.staining.sdlog = 1

n.cell.min = 50
n.cell.max = 2000
n.cell.logmean = log(mean(c(n.cell.min, n.cell.max)))
n.cell.logsd = .8
ncell <- round(rlnorm(n.bc,n.cell.logmean,n.cell.logsd))
ncell[ncell < n.cell.min] = n.cell.min
ncell[ncell > n.cell.max] = n.cell.max
hist(ncell)

nb.theta = 7
b0 = -4
ambient.meanlog = 1
doublet.rate = 0.1
dropout.lambda = 2

cur_text <- paste0(test_text, "sim3_", "n.bc_", n.bc, "_")


# Sim 4
set.seed(2023)
n.bc = 30

cell.staining.meanlog = rnorm(n.bc, mean = 5, sd = .1)
cell.staining.sdlog = 1

n.cell.min = 50
n.cell.max = 3000
n.cell.logmean = log(mean(c(n.cell.min, n.cell.max))/2)
n.cell.logsd = 1
ncell <- round(rlnorm(n.bc,n.cell.logmean,n.cell.logsd))
ncell[ncell < n.cell.min] = n.cell.min
ncell[ncell > n.cell.max] = n.cell.max
hist(ncell)

nb.theta = 7
b0 = -4
ambient.meanlog = 1
doublet.rate = 0.1
dropout.lambda = 2

cur_text <- paste0(test_text, "sim4_", "n.bc_", n.bc, "_")



# Sim 5
set.seed(2023)
n.bc = 30

cell.staining.meanlog = rnorm(n.bc, mean = 5, sd = .5)
cell.staining.sdlog = 1.5

n.cell.min = 50
n.cell.max = 3000
n.cell.logmean = log(mean(c(n.cell.min, n.cell.max))/2)
n.cell.logsd = 1
ncell <- round(rlnorm(n.bc,n.cell.logmean,n.cell.logsd))
ncell[ncell < n.cell.min] = n.cell.min
ncell[ncell > n.cell.max] = n.cell.max
hist(ncell)

nb.theta = 10
b0 = -4
ambient.meanlog = 3
doublet.rate = 0.2
dropout.lambda = .5

cur_text <- paste0(test_text, "sim5_", "n.bc_", n.bc, "_")




# ambient.meanlog = cell.staining.meanlog - rtruncnorm(n.bc, a = ambient.min.diff, b = ambient.max.diff, mean = ambient.meanlog.mean, sd = ambient.meanlog.sd)
# ambient.meanlog[ambient.meanlog < 0] = 0
sim_res <- simulateTags(n.cell = ncell,
                        n.bc = n.bc,
                        seed = 2023,
                        nb.theta = nb.theta,
                        cell.staining.meanlog = cell.staining.meanlog,
                        cell.staining.sdlog = cell.staining.sdlog,
                        ambient.meanlog = ambient.meanlog,
                        b0 = b0,
                        doublet.rate = doublet.rate,
                        dropout.lambda = dropout.lambda,
                        cell.contam = T,
                        ambient.contam = T,
                        return.all = T)
saveRDS(sim_res, paste0(rdsPath, cur_text, "sim_res.rds"))
print(dim(sim_res$final.umi.mtx))




#sim_res <- readRDS(paste0(rdsPath, cur_text, "sim_res.rds"))
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




set.seed(1)
umap_raw <- compute_umap(tag_mtx, use_dim = ncol(tag_mtx), n_component=2, n_neighbors = 30)
umap_raw <- as.data.frame(umap_raw)
umap_raw$true_label <- true_label
unq_bcs = colnames(tag_mtx)
use_color = get_factor_color(unq_bcs, "Set1")
names(use_color) = unq_bcs
use_color <- c(use_color, 'doublet' = 'black')
g1<- ggplot(umap_raw, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = "true_label"), stroke = 0, size =.1, raster.dpi=600) +
    scale_color_manual(values = use_color, na.value='lightgrey') +
    theme_bw() +
    # guides(colour = guide_legend(override.aes = list(size=3),
    #                              keywidth=0.1,
    #                              keyheight=0.12,
    #                              default.unit="inch"))+
    guides(color =F)+
    theme(text=element_text(size=8),
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
ggsave(paste0(plotPath, cur_text, "_umap_raw.pdf"), g1, width = 1.7, height =1.5)
ggsave(paste0(plotPath, cur_text, "_umap_raw.pdf"), g1, width = 2.5, height =1.6)

umap_raw$tt.umi <- rowSums(tag_mtx)
umap_raw$tt.umi[umap_raw$tt.umi > quantile(umap_raw$tt.umi,.975)] = quantile(umap_raw$tt.umi,.975)
g1<- ggplot(umap_raw, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = "log(tt.umi)"), stroke = 0, size = .5, raster.dpi=600) +
    scale_color_gradientn(colors = get_gradient_color("BlueGreenRed")) +
    theme_bw() +
    # guides(colour = guide_legend(override.aes = list(size=3),
    #                              keywidth=0.1,
    #                              keyheight=0.12,
    #                              default.unit="inch"))+
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
ggsave(paste0(plotPath, test_text, "_umap_raw_tt.count.pdf"), g1, width = 2.2, height =1.6)

# show_bc = names(which(table(umap_raw[["true_label"]]) > 0))
# label_data <- umap_raw %>% group_by_at("true_label") %>% summarize_at(c("UMAP_1", "UMAP_2"), median)
# label_data <- label_data[label_data[["true_label"]] %in% show_bc,,drop=F]
# g1 <- g1 + geom_label(
#     aes_string(
#         x="UMAP_1",y="UMAP_2",
#         label = "true_label"
#     ),
#     color = "black",
#     size = 2,
#     data = label_data
# )

bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = paste0(cur_text, 'demultiplex2_'), width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 5000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(bc_demux2, paste0(rdsPath, cur_text, "bc_demux2.rds"))

bc_demux2 <- readRDS(paste0(rdsPath, cur_text, "bc_demux2.rds"))
umap_rqr <- as.data.frame(bc_demux2$res$umap)
umap_rqr$tag_assign <- bc_demux2$res$final_assign
unq_bcs = colnames(tag_mtx)
use_color = get_factor_color(unq_bcs, "Set1")
names(use_color) = unq_bcs
use_color <- c(use_color, 'multiplet' = 'black', 'negative' = 'grey')
g1<- ggplot(umap_rqr, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = "tag_assign"), stroke = 0, size = .2, raster.dpi=600) +
    scale_color_manual(values = use_color, na.value='lightgrey') +
    theme_bw() +
    # guides(colour = guide_legend(override.aes = list(size=3),
    #                              keywidth=0.1,
    #                              keyheight=0.12,
    #                              default.unit="inch"))+
    guides(color =F)+
    theme(text=element_text(size=8),
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
ggsave(paste0(plotPath, cur_text, "_umap_rqr.pdf"), g1, width = 1.7, height =1.5)
ggsave(paste0(plotPath, cur_text, "_umap_rqr.pdf"), g1, width = 2.5, height =1.6)






# Benchmarking

# 1. deMULTIplex
library(deMULTIplex)
library(ggplot2)
bc_demux1 <- benchmark_demultiplex1(as.matrix(tag_mtx), true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'deMULTIplex1_'))




# 3. GMM-Demux
# The code to run this software is a bit cumbersome, generate script from R
# Step 1, save the tag matrix to the csv file
bcs <- make.names(colnames(tag_mtx))
tag_mtx_write <- as.matrix(tag_mtx)
colnames(tag_mtx_write) <- bcs
write.csv(tag_mtx_write, paste0(rdsPath, "tag_mtx_gmm.csv"), quote=F)
resPath <- paste0(savePath, test_text, "gmmdemux_res/"); dir.create(resPath)
cmds <- paste0("GMM-demux -c ./rds/tag_mtx_gmm.csv ", paste0(bcs, collapse = ","), " -x ", bcs, " -o ", cur_text, "gmmdemux_res/", bcs)
fileConn<-file(paste0(savePath, "run_gmm.sh"))
writeLines(cmds, fileConn)
close(fileConn)

# Get results from assignment
bcs <- make.names(colnames(tag_mtx))
gmm_res <- lapply(bcs, function(bc) {
    print(bc)
    barcode.path <- paste0(savePath, test_text,"gmmdemux_res/", bc, "/barcodes.tsv.gz")
    tryCatch({
        cells <<- read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
        colnames(cells) <- "cell_barcode"
        cells$assignment = bc
    }, error = function(e) {
        cells <<- data.frame(matrix(ncol = 2, nrow = 0))
    })
    return(cells)
})
gmm_res <- do.call(rbind, gmm_res)
gmm_res_call <- gmm_res$assignment
names(gmm_res_call) <- gmm_res$cell_barcode
gmm_res_call <- gsub("[.]", "-", gmm_res_call) # Correction
bc_gmm <- confusion_stats(gmm_res_call, true_label,
                          tag_mapping,
                          call.multiplet = 'Multiplet',
                          true.multiplet = "doublet",
                          plot.path =plotPath, plot.name = paste0(cur_text, 'gmmdemux_'),  width = 3.5, height = 2.5)


# 4. HTODemux
bc_seu <- benchmark_HTODemux(tag_mtx, rna_mtx = t(tag_mtx), true_label = true_label,
                             tag_mapping = tag_mapping,
                             true.multiplet = "doublet",
                             plot.path = plotPath, plot.name = paste0(cur_text, 'HTODemux_'), width = 3.5, height = 2.5)


# 5. hashedDrops
require(DropletUtils)
demux_res <- hashedDrops(t(tag_mtx))
demux_res$assign <- colnames(tag_mtx)[demux_res$Best]
demux_res$assign[!demux_res$Confident] = 'ambiguous'
demux_res$assign[demux_res$Doublet] = 'doublet'
hashed_call <- demux_res$assign
names(hashed_call) <- rownames(demux_res)
bc_hashed = confusion_stats(hashed_call, true_label,
                            tag_mapping = tag_mapping,
                            call.multiplet = 'doublet',
                            true.multiplet = "doublet",
                            plot.path =plotPath, plot.name = paste0(cur_text, 'hashedDrops_'), width = 3.5, height = 2.5)

# 6. demuxmix
require(demuxmix)
# bc_mixfull <- benchmark_demuxmix_full(tag_mtx, rna_mtx, true_label,
#                                       tag_mapping = tag_mapping,
#                                       true.multiplet = "doublet",
#                                       plot.path = plotPath, plot.name = paste0(cur_text, 'demuxmixfull_'), width = 3.5, height = 2.5)
bc_mixnaive <- benchmark_demuxmix_naive(tag_mtx, true_label = true_label,
                                        tag_mapping = tag_mapping,
                                        true.multiplet = "doublet",
                                        plot.path = plotPath, plot.name = paste0(cur_text, 'demuxmixnaive_'), width = 3.5, height = 2.5)


# 7. BFF raw and BFF cluster
library(cellhashR)
cellhashR_calls <- GenerateCellHashingCalls(barcodeMatrix = t(tag_mtx),c("bff_raw", "bff_cluster"),
                                            doTSNE = FALSE,
                                            doUMAP = FALSE, doHeatmap = FALSE, bff_raw.min_average_reads=0, bff_cluster.min_average_reads=0)
bff_raw_call <- cellhashR_calls$bff_raw; names(bff_raw_call) <- cellhashR_calls$cellbarcode
bff_cluster_call <-cellhashR_calls$bff_cluster; names(bff_cluster_call) <- cellhashR_calls$cellbarcode

bc_bffraw <- confusion_stats(bff_raw_call, true_label,
                             tag_mapping,
                             call.multiplet = 'Doublet',
                             true.multiplet = "doublet",
                             plot.path =plotPath, plot.name = paste0(cur_text, 'bffraw_'),  width = 3.5, height = 2.5)

bc_bffcluster <- confusion_stats(bff_cluster_call, true_label,
                                 tag_mapping,
                                 call.multiplet = 'Doublet',
                                 true.multiplet = "doublet",
                                 plot.path =plotPath, plot.name = paste0(cur_text, 'bffcluster_'),  width = 3.5, height = 2.5)


nsn <- function(x){
    if(exists(x[1])) {
        if(length(x) ==1) {
            return(get(x))
        } else {
            return(get(x[2], get(x[1])))
        }
    } else {
        return(NULL)
    }
}

# Combine all results
bc_all <- list(
    deMULTIplex2 = nsn("bc_demux2"),
    deMULTIplex = nsn("bc_demux1"),
    #demuxEM = nsn("bc_demuxEM"),
    GMM_Demux = nsn("bc_gmm"),
    HTODemux = nsn("bc_seu"),
    hashedDrops = nsn("bc_hashed"),
    #demuxmix = nsn("bc_mixfull"),
    demuxmix_naive = nsn("bc_mixnaive"),
    bff_raw = nsn("bc_bffraw"),
    bff_cluster = nsn("bc_bffcluster")
)
saveRDS(bc_all, paste0(rdsPath, cur_text, "bc_all.rds"))


# Plot summary plots
summary_df <- sapply(bc_all, function(x) {
    if(!is.null(x)) {
        x$singlet_avg_stats
    } else {
        return(c(0,0,0))
    }
})
write.csv(round(t(summary_df),3), paste0(plotPath, cur_text,"summary_df.csv"))
plot_df <- reshape2::melt(summary_df)
colnames(plot_df) <- c("stats", "method", "value")
g1<-ggplot(plot_df, aes(method, value)) +
    geom_col() +
    facet_wrap(~stats, scales = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(plotPath, cur_text, 'bc_all_bar.pdf'), g1, width = 8, height = 8)



bc_all <- readRDS(paste0(rdsPath, cur_text, "bc_all.rds"))
use_color <- gg_color_hue(length(bc_all[!names(bc_all) == "demuxmix"]))
names(use_color) <- names(bc_all)[!names(bc_all) == "demuxmix"]

bc_all[which(sapply(bc_all, is.null))] <- NULL
# Plot stats of each barcode
bc_stats <- lapply(1:length(bc_all), function(i) {
    x = bc_all[[i]]
    df <- as.data.frame(do.call(rbind, x$tag_stats))
    df$tag_id <- rownames(df)
    df$method =names(bc_all)[i]
    return(df)
})
bc_stats <- do.call(rbind, bc_stats)

plot_df <- bc_stats[,c('precision', 'recall', 'method')]
colnames(plot_df) <- c('Precision', 'Recall', 'Method')
plot_df <- reshape2::melt(plot_df, id.var = 'Method')
plot_df$Method <- factor(plot_df$Method, levels = names(bc_all))

g1 <- ggplot(plot_df, aes(x = Method, y = value, fill = Method)) +
    geom_boxplot(lwd=.2, outlier.shape = NA, alpha = .3) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1),dotsize = .5,stroke = 0.5)+
    stat_summary(fun.y=mean, geom="point", shape=23, size=1.5) +
    facet_wrap(~variable, nrow=2, strip.position = 'right')+
    scale_fill_manual(values = use_color) +
    ylim(c(0,1))+
    #coord_flip() +
    guides(fill = 'none', color = 'none') +
    theme_bw() +
    xlab(NULL) +
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

ggsave(paste0(plotPath, cur_text, 'bc_precision_recall_boxplot.pdf'), g1, width =2, height = 2.2)



# use_color <- gg_color_hue(length(names(bc_all)))
# names(use_color) <- names(bc_all)
summary_df <- sapply(bc_all, function(x) {
    x$singlet_avg_stats
})
plot_df <- as.data.frame(t(summary_df['f_score', ,drop=F]))
colnames(plot_df) <- c("F_score")
plot_df$Method <- factor(rownames(plot_df), levels = rev(names(bc_all)))
g1<-ggplot(plot_df, aes(Method, F_score)) +
    geom_bar(stat='identity',fill = '#1e90ff') +
    coord_flip() +
    guides(fill = 'none', color = 'none') +
    theme_bw() +
    xlab(NULL) +
    scale_color_manual(values = rev(use_color)) +
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
ggsave(paste0(plotPath, cur_text, 'bc_f_score_bar.pdf'), g1, width = 2, height = 1.8)


# Doublet assignment
doublet_df <- sapply(bc_all, function(x) {
    x$doublet_avg_stats
})
rownames(doublet_df) <- c("doublet", "singlet", "negative")
plot_df <- reshape2::melt(doublet_df)
colnames(plot_df) <- c("type", "method", "value")

use_color <- c(
    "negative" = "#878787",  "singlet"="#e31a1c", "doublet" = "#377eb8"
)
plot_df$type <- factor(plot_df$type, levels = names(use_color))
g1 <- ggplot(plot_df, aes(method, value, group = type)) +
    geom_bar(stat='identity',aes(fill = type))+
    ylab('Fraction') +
    xlab(NULL) +
    scale_fill_manual(values = use_color) +
    #guides(fill = guide_legend(override.aes = list(width = .1)))+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,1)) +
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          axis.text = element_text(size=8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    monocle:::monocle_theme_opts()
ggsave(paste0(plotPath, cur_text, 'bc_doublet_assign.pdf'), g1, width = 2.5, height = 1.8)










### Summarize across all simulations
sim_conditions <- paste0(test_text, "sim", 1:5, "_", "n.bc_", c(5,10,10,30,30), "_")

bc_res <- lapply(sim_conditions, function(x) {
    readRDS(paste0(rdsPath, x, "bc_all.rds"))
})
names(bc_res) <- paste0("sim_", 1:5)

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
singlet_df_cbn <- singlet_df_cbn[singlet_df_cbn$method!="demuxmix",]
singlet_df_cbn$dataset <- factor(singlet_df_cbn$dataset, levels= names(bc_res))
plot_df <- singlet_df_cbn[singlet_df_cbn$stats == "f_score", ]
g1 <- ggplot(plot_df, aes(x = dataset, y = value, fill = method)) +
    geom_bar(stat="identity",position ="dodge")

ggsave(paste0(plotPath, test_text, "_barv1.pdf"), g1, width = 10, height =2)


# Heatmap
library(pheatmap)
show_methods <- names(bc_res[[5]])
singlet_mtx <- sapply(singlet_stats, function(x) {
    x["f_score",][show_methods]
})

breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "singlet_summary_mtx", ".pdf"), width = 3, height = 1.8)
pheatmap(singlet_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
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
pdf(paste0(plotPath, test_text, "doublet_summary_mtx", ".pdf"), width = 3, height = 1.8)
pheatmap(doublet_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()



doublet_mtx <- sapply(doublet_stats, function(x) {
    x["doublet_called_singlet",][show_methods]
})
doublet_mtx <- doublet_mtx[rownames(doublet_mtx)!= "GMM_Demux",]
breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "doublet_called_singlet_summary_mtx", ".pdf"), width = 3, height = 1.8)
pheatmap(doublet_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
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











# Barcode case study

cur_text <- paste0(test_text, "sim", 1, "_", "n.bc_", 5, "_")
bc = 'bc1'

cur_text <- paste0(test_text, "sim", 2, "_", "n.bc_", 10, "_")
bc = 'bc9'

cur_text <- paste0(test_text, "sim", 3, "_", "n.bc_", 10, "_")
bc = 'bc4'

cur_text <- paste0(test_text, "sim", 4, "_", "n.bc_", 30, "_")
bc = 'bc7'

cur_text <- paste0(test_text, "sim", 5, "_", "n.bc_", 30, "_")
bc = 'bc2'



sim_res <- readRDS(paste0(rdsPath, cur_text, "sim_res.rds"))
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

bc_all <- readRDS(paste0(rdsPath, cur_text, "bc_all.rds"))
bc_demux2 <- bc_all$deMULTIplex2



plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
map =c("log(tt.umi)", "log(bc.umi)", "true_label")

plot_df$call_label <- paste0("call_", bc_demux2$res$assign_table$barcode_assign == bc)
plot_df$true_label_bin <- ifelse(grepl(paste0(bc, "_"), rownames(plot_df)), bc, "Other")

use_color <- c(
    '#8D4C6A',
    '#bdbdbd'
)
names(use_color) <- c(bc, "Other")


map =c("log(tt.umi)", "log(bc.umi)", "true_label_bin")
g1 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .2, raster.dpi=600) +
    geom_abline(slope = 1,intercept = 0, alpha = .2) +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_manual(values = use_color) +
    scale_x_continuous(limits = c(log(quantile(plot_df$tt.umi, .005)), log(quantile(plot_df$tt.umi, .995))))+
    # guides(colour = guide_legend(override.aes = list(size=3),
    #                              keywidth=0.1,
    #                              keyheight=0.12,
    #                              default.unit="inch"))+
    guides(color = F)+
    theme(text=element_text(size=7),
          legend.text=element_text(size=7),
          legend.title = element_blank(),
          axis.text = element_text(size=7),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, cur_text, bc, "_true_label.pdf"), g1, width = 1.4, height =1.1)


map =c("log(tt.umi)", "log(bc.umi)", "post1")
g1 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .2, raster.dpi=600) +
    geom_abline(slope = 1,intercept = 0, alpha = .2) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'lightsteelblue') +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    scale_x_continuous(limits = c(log(quantile(plot_df$tt.umi, .005)), log(quantile(plot_df$tt.umi, .995))))+
    guides(colour = F)+
    theme(text=element_text(size=7),
          legend.position = "right",
          legend.text=element_text(size=7),
          legend.title = element_blank(),
          axis.text = element_text(size=7),
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
#g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, cur_text, bc, "_post1.pdf"), g1, width = 1.35, height =.9)






# Plot for Figure S1

cur_text <- paste0(test_text, "sim", 1, "_", "n.bc_", 5, "_")
bc = 'bc1'

cur_text <- paste0(test_text, "sim", 5, "_", "n.bc_", 30, "_")
bc = 'bc3'



sim_res <- readRDS(paste0(rdsPath, cur_text, "sim_res.rds"))
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

bc_all <- readRDS(paste0(rdsPath, cur_text, "bc_all.rds"))
bc_demux2 <- bc_all$deMULTIplex2


plot_df <- bc_demux2$res$df_list[[bc]]


map =c("log(tt.umi)", "log(bc.umi)", "post1")
g1 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3, raster.dpi=600) +
    geom_abline(slope = 1,intercept = 0,  alpha = .2) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'lightsteelblue') +
    scale_x_continuous(limits = c(2.5, 10), breaks = seq(2.5,10,2.5))+
    scale_y_continuous(limits = c(0, 8), breaks = seq(0,8,2))+
    xlab("log(total)") +
    ylab("log(tag)") +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    guides(colour = guide_colorbar(barwidth = 1,
                                   barheight = .1,
                                   default.unit="inch"))+
    theme(text=element_text(size=8),
          legend.position = "top",
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, cur_text, bc, "_post1.pdf"), g1, width = 1.4, height =1.6)



map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "post1")
g2 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3, raster.dpi=600) +
    geom_abline(slope = 1,intercept = 0,  alpha = .2) +
    geom_line(aes_string(map[1], 'log(pred1)'),color = 'lightsteelblue') +
    scale_x_continuous(limits = c(2.5, 10), breaks = seq(2.5,10,2.5))+
    scale_y_continuous(limits = c(2.5, 10), breaks = seq(2.5,10,2.5))+
    xlab("log(total)") +
    ylab("log(total-tag)") +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    guides(colour = guide_colorbar(barwidth = 1,
                                   barheight = .1,
                                   default.unit="inch"))+
    theme(text=element_text(size=8),
          legend.position = "top",
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, cur_text, bc, "_post1_reverse.pdf"), g2, width = 1.5, height =1.6)



max.rqr = max(plot_df$rqr[is.finite(plot_df$rqr)]) + 1 # Best to cut inf?
plot_df$rqr[plot_df$rqr > max.rqr] = max.rqr
map =c("log(tt.umi)", "rqr", "true_label")
g1 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .3, raster.dpi=600) +
    #geom_abline(slope = 1,intercept = 0) +
    #geom_abline(slope = 0,intercept = 0, color = 'grey') +
    scale_x_continuous(limits = c(2.5, 10), breaks = seq(2.5,10,2.5))+
    xlab("log(total)") +
    ylab("RQR") +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    guides(colour = guide_colorbar(barwidth = 1,
                                   barheight = .1,
                                   default.unit="inch"))+
    theme(text=element_text(size=8),
          legend.position = "top",
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, cur_text, bc, "_rqr.pdf"), g1, width = 1.5, height =1.8)




p1 = ggplot(plot_df[plot_df$post1 < .5,], aes_string(sample='rqr'))+
    stat_qq(size = .5, stroke = 0) +
    stat_qq_line(alpha = .3) +
    xlab("Normal quantile") +
    ylab("RQR") +
    theme_bw() +
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, cur_text, bc, "_rqr_qq_predicted.pdf"), p1, width = 1.2, height =1)


max.rqr_p = max(plot_df$rqr_p[is.finite(plot_df$rqr_p)]) + 1 # Best to cut inf?
plot_df$rqr_p[plot_df$rqr_p > max.rqr_p] = max.rqr_p
map =c("log(tt.umi)", "rqr_p", "post1")
g1 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3, raster.dpi=600) +
    #geom_abline(slope = 1,intercept = 0) +
    #geom_abline(slope = 0,intercept = 0, color = 'grey') +
    scale_x_continuous(limits = c(2.5, 10), breaks = seq(2.5,10,2.5))+
    xlab("log(total)") +
    ylab("RQR") +
    theme_bw() +
    labs(color = gsub("_","\n",map[3]))+
    scale_color_gradientn(colors = get_gradient_color('BlueGreenRed')) +
    guides(colour = guide_colorbar(barwidth = 1,
                                   barheight = .1,
                                   default.unit="inch"))+
    theme(text=element_text(size=8),
          legend.position = "top",
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, cur_text, bc, "_rqrp.pdf"), g1, width = 1.5, height =1.8)


p1 = ggplot(plot_df[plot_df$post1 > .5,], aes_string(sample='rqr_p'))+
    stat_qq(size = .5, stroke = 0) +
    stat_qq_line(alpha = .3) +
    scale_x_continuous(limits = c(-3, 3))+
    xlab("Normal quantile") +
    ylab("RQR") +
    theme_bw() +
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, cur_text, bc, "_rqrp_qq_predicted.pdf"), p1, width = 1.2, height =1)


