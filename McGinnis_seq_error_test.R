

library(deMULTIplex2)
library(profvis)

test_text <- "McGinnis_seq_error_"
savePath <- "./data-raw/benchmark_res/McGinnis_seq_error/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/"); dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

library(deMULTIplex2)
fastq_dir <- "./benchmark_res/McGinnis_seq_error/BAR/"
read_table <- readTags(dir = fastq_dir,
                       name = "ACAGTG",
                       barcode.type = "MULTIseq",
                       assay = "RNA")

rna_mtx <- readRDS("benchmark_res/McGinnis_8donor_MULTI/rds/McGinnis_8donor_MULTI_benchmarking_rna_mtx.rds")

data(multiseq_oligos) # Current MULTI-seq oligo sequence provided with the package
tag.ref <- multiseq_oligos



tag_mtx <- alignTags(read_table, tag.ref, filter.cells = colnames(rna_mtx), max.dist = 1)
#cell_ids <- Matrix::rowSums(tag_mtx) > 100 # Filter for cells (it's better to use your own cell filter criteria based on RNA count)
tag_used <- Matrix::colSums(tag_mtx) > 1e4 # Filter for used tags
tag_mtx <- tag_mtx[,tag_used]

# True label processed from McGinnis_8donor_benchmarking_MULTI.R
true_label <- readRDS("benchmark_res/McGinnis_8donor_MULTI/rds/McGinnis_8donor_MULTI_benchmarking_true_label.rds")
#tag_mapping <- readRDS("benchmark_res/McGinnis_8donor_MULTI/rds/McGinnis_8donor_MULTI_benchmarking_tag_mapping.rds")
res <- demultiplexTags(tag_mtx,
                       plot.path = plotPath,
                       plot.name = "preassign_",
                       plot.umap = "none",
                       plot.diagnostics = F)
shared_cells <- intersect(names(res$final_assign), names(true_label))
donor_true <- true_label[shared_cells]
tag_assign <- res$final_assign[shared_cells]
conf_mtx <- as.data.frame.matrix(table(donor_true, tag_assign))
tag_mapping <- apply(conf_mtx[1:8,1:8], 1, function(x) colnames(conf_mtx)[1:8][which.max(x)])
tag_mapping <- data.frame(tag =tag_mapping, true_label = names(tag_mapping))

saveRDS(tag_mapping, paste0(rdsPath, test_text, "tag_mapping.rds"))

# tag_mapping needs to be remade based on one round of assignment because the name is different e.g., "Bar11" in the published file but "A11" in the original reference

tag_mapping <- readRDS(paste0(rdsPath, "McGinnis_seq_error_", "tag_mapping.rds"))






#### Run benchmarking for each distance cutoff #####

show_dist <- c(0,1,2,3)

cur_max_dist = show_dist[1]
#cur_max_dist = show_dist[2]
#cur_max_dist = show_dist[3]
#cur_max_dist = show_dist[4]

cur_text <- paste0(test_text, "hamming_dist_", cur_max_dist, "_")
tag_mtx <- alignTags(read_table, tag.ref, filter.cells = colnames(rna_mtx), max.dist = cur_max_dist)
tag_used <- Matrix::colSums(tag_mtx) > 1e4
tag_mtx <- as.matrix(tag_mtx[,tag_used])
saveRDS(tag_mtx, paste0(rdsPath, cur_text, "tag_mtx.rds"))

Matrix::colSums(tag_mtx)
print(dim(tag_mtx))

# Benchmark the performance of deMULTIplex2
bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = cur_text, width = 3.5, height = 2.5,
                                    seed = 1,
                                    init.cos.cut = .5,
                                    converge.threshold = 1e-3,
                                    prob.cut = 0.5,
                                    max.cell.fit = 5000,
                                    max.iter = 30,
                                    min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "residual",
                                    plot.diagnostics = T)
saveRDS(bc_demux2, paste0(rdsPath, cur_text, "bc_demux2.rds"))

# 1. deMULTIplex
bc_demux1 <- benchmark_demultiplex1(as.matrix(tag_mtx), true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = paste0(cur_text, 'deMULTIplex1_'))
saveRDS(bc_demux1, paste0(rdsPath, cur_text, "bc_demux1.rds"))


# 2. demuxEM, failed to run

# 3. GMM-Demux
# The code to run this software is a bit cumbersome, generate script from R
# Step 1, save the tag matrix to the csv file
bcs <- make.names(colnames(tag_mtx))
tag_mtx_write <- tag_mtx
colnames(tag_mtx_write) <- bcs
write.csv(tag_mtx_write, paste0(rdsPath, "tag_mtx_gmm.csv"), quote=F)
resPath <- paste0(savePath, cur_text, "gmmdemux_res/"); dir.create(resPath)
cmds <- paste0("GMM-demux -c ./rds/tag_mtx_gmm.csv ", paste0(bcs, collapse = ","), " -x ", bcs, " -o ", cur_text, "gmmdemux_res/", bcs)
fileConn<-file(paste0(savePath, "run_gmm.sh"))
writeLines(cmds, fileConn)
close(fileConn)

# Get results from assignment
bcs <- make.names(colnames(tag_mtx))
gmm_res <- lapply(bcs, function(bc) {
    print(bc)
    barcode.path <- paste0(savePath, cur_text,"gmmdemux_res/", bc, "/barcodes.tsv.gz")
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
bc_seu <- benchmark_HTODemux(tag_mtx, rna_mtx = rna_mtx, true_label = true_label,
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
bc_mixfull <- benchmark_demuxmix_full(tag_mtx, rna_mtx, true_label,
                                      tag_mapping = tag_mapping,
                                      true.multiplet = "doublet",
                                      plot.path = plotPath, plot.name = paste0(cur_text, 'demuxmixfull_'), width = 3.5, height = 2.5)
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





# Combine all results
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
    demuxEM = nsn("bc_demuxEM"),
    GMM_Demux = nsn("bc_gmm"),
    HTODemux = nsn("bc_seu"),
    hashedDrops = nsn("bc_hashed"),
    demuxmix = nsn("bc_mixfull"),
    demuxmix_naive = nsn("bc_mixnaive"),
    bff_raw = nsn("bc_bffraw"),
    bff_cluster = nsn("bc_bffcluster")
)
saveRDS(bc_all, paste0(rdsPath, cur_text, "bc_all.rds"))

demux2_call <- bc_demux2$res$assign_table$barcode_assign
demux2_call <- bc_demux2$res$assign_table$barcode_assign
demux2_call[bc_demux2$res$assign_table$barcode_count == 0] = "Negative"
demux2_call[bc_demux2$res$assign_table$barcode_count > 1] = "Multiplet"
names(demux2_call) <- rownames(bc_demux2$res$assign_table)
call_labels <- list(
    deMULTIplex2 = nsn("demux2_call"),
    deMULTIplex = nsn(c("bc_demux1", "res")),
    demuxEM = nsn("demuxem_call"),
    GMM_Demux = nsn("gmm_res_call"),
    HTODemux = nsn(c("bc_seu", "res")),
    hashedDrops = nsn("hashed_call"),
    demuxmix = nsn(c("bc_mixfull", "call")),
    demuxmix_naive = nsn(c("bc_mixnaive","call")),
    bff_raw = nsn("bff_raw_call"),
    bff_cluster = nsn("bff_cluster_call")
)
saveRDS(call_labels, paste0(rdsPath, cur_text, "call_labels.rds"))


# Plot summary plots
summary_df <- sapply(bc_all, function(x) {
    if(!is.null(x)) {
        x$singlet_avg_stats
    } else {
        return(c(0,0,0))
    }
})
write.csv(round(t(summary_df),3), paste0(plotPath, test_text,"summary_df.csv"))
plot_df <- reshape2::melt(summary_df)
colnames(plot_df) <- c("stats", "method", "value")
g1<-ggplot(plot_df, aes(method, value)) +
    geom_col() +
    facet_wrap(~stats, scales = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(plotPath, cur_text, 'bc_all_bar.pdf'), g1, width = 8, height = 8)
























bc_res <- lapply(show_dist, function(cur_dist) {
    print(cur_dist)
    readRDS(paste0(rdsPath, "McGinnis_seq_error_hamming_dist_",cur_dist, "_bc_all.rds"))
})
names(bc_res) <- as.character(show_dist)

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
show_methods <- names(bc_res[[1]])
show_methods <- show_methods[show_methods!= "demuxEM"]
singlet_df_cbn <- singlet_df_cbn[singlet_df_cbn$method %in% show_methods,]
singlet_df_cbn$method <- factor(as.character(singlet_df_cbn$method), levels = show_methods)
singlet_df_cbn$dataset <- factor(singlet_df_cbn$dataset, levels= names(bc_res))
plot_df <- singlet_df_cbn[singlet_df_cbn$stats == "f_score", ]
use_color <- get_gradient_color("BlueGreenRed",5)
names(use_color) <- names(bc_res)
g1 <- ggplot(plot_df, aes(x = method, y = value, fill = dataset)) +
    geom_bar(stat="identity",position ="dodge") +
    scale_fill_manual(values = use_color) +
    ylab("F score")+
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
          legend.text=element_text(size=8),
          axis.text = element_text(size=8),
          axis.text.x = element_text(angle = 15, vjust = 1, hjust=1),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))

ggsave(paste0(plotPath, test_text, "_barv1.pdf"), g1, width = 5, height =2)

# Line plot
plot_df$value[is.na(plot_df$value)] = 0
g1 <- ggplot(plot_df, aes(x = dataset, y = value, color = method, group = method)) +
    geom_line() +
    geom_point()+
    ylab("F score")+
    theme_bw()+
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.15,
                                 default.unit="inch"))+
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=8),
          axis.text.x = element_text(angle = 15, vjust = 1, hjust=1),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_linev1.pdf"), g1, width = 4, height =1.8)



# Heatmap
library(pheatmap)
singlet_mtx <- sapply(singlet_stats, function(x) {
    x["f_score",][show_methods]
})

breaksList <- seq(0,1,by=.1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "singlet_summary_mtx", ".pdf"), width = 2.8, height = 1.8)
pheatmap(singlet_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', angle_col = 0, fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()









# Analyze string distance and see how that relates to the pattern
library(stringdist)
tag_used <- tag_mapping$tag
ref_dist_mtx <-as.matrix(stringdistmatrix(tag.ref[tag_used], method = "hamming"))
colnames(ref_dist_mtx) <- names(tag.ref[tag_used])
rownames(ref_dist_mtx) <- names(tag.ref[tag_used])

pdf(paste0(plotPath, test_text, "multiseq_ref_dist_mtx", ".pdf"), width = 3, height = 2.6)
pheatmap(ref_dist_mtx, cluster_rows = T, cluster_cols = T, display_numbers =T, number_format ="%.0f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('viridis'), clustering_method = "ward.D2")
dev.off()

# Case study


# Individual barcode stats
true_label <- readRDS("./data-raw/benchmark_res/McGinnis_8donor_MULTI/rds/McGinnis_8donor_MULTI_benchmarking_true_label.rds")
cur_text <- paste0(test_text, "hamming_dist_", 3, "_")
tag_mtx <- readRDS(paste0(rdsPath, cur_text, "tag_mtx.rds"))
bc_all <- readRDS(paste0(rdsPath, "McGinnis_seq_error_hamming_dist_",3, "_bc_all.rds"))


bc_demux2 <- bc_all$deMULTIplex2
bc = 'A10'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label[rownames(plot_df)]
map =c("log(tt.umi)", "log(bc.umi)", "true_label")
td <- tag_mapping$true_label[tag_mapping$tag == bc]
plot_df$true_label_bin <- ifelse(plot_df$true_label == td,td, "Other")
plot_df$true_label_well <- tag_mapping$tag[match(plot_df$true_label, tag_mapping$true_label)]

plot_df$true_label_well <- factor(plot_df$true_label_well, levels = gtools::mixedsort(names(table(plot_df$true_label_well))))
map =c("log(tt.umi)", "log(bc.umi)", "true_label_well")
g1 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3, raster.dpi = 600) +
    geom_abline(slope = 1,intercept = 0, alpha =.2) +
    theme_bw() +
    xlab("log(total)") +
    ylab("log(tag)") +
    labs(color = gsub("_","\n",map[3]))+
    #scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
#g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, test_text, bc, "_true_label.pdf"), g1, width = 2.3, height =1.5)

map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "true_label_well")
g1 <- ggplot(plot_df) +
    geom_point_rast(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3, raster.dpi = 600) +
    geom_abline(slope = 1,intercept = 0, alpha =.2) +
    theme_bw() +
    xlab("log(total)") +
    ylab("log(total-tag)") +
    labs(color = gsub("_","\n",map[3]))+
    #scale_color_manual(values = use_color) +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
#g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, test_text, bc, "_true_label_reverse.pdf"), g1, width = 2.2, height =1.5)



map =c("log(tt.umi)", "log(bc.umi)", "post1")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3) +
    geom_abline(slope = 1,intercept = 0, alpha = .2) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'lightsteelblue') +
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
ggsave(paste0(plotPath, test_text, bc, "_post1.pdf"), g1, width = 1.8, height =1.8)



map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "post1")
g2 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3) +
    geom_abline(slope = 1,intercept = 0, alpha = .2) +
    geom_line(aes_string(map[1], 'log(pred1)'),color = 'lightsteelblue') +
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
ggsave(paste0(plotPath, test_text, bc, "_post1_reverse.pdf"), g2, width = 1.67, height =1.82)



max.rqr = max(plot_df$rqr[is.finite(plot_df$rqr)]) + 1 # Best to cut inf?
plot_df$rqr[plot_df$rqr > max.rqr] = max.rqr
map =c("log(tt.umi)", "rqr", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .3) +
    #geom_abline(slope = 1,intercept = 0) +
    #geom_abline(slope = 0,intercept = 0, color = 'grey') +
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
ggsave(paste0(plotPath, test_text, bc, "_rqr.pdf"), g1, width = 1.5, height =1.8)




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
ggsave(paste0(plotPath, test_text, "_rqr_qq_predicted.pdf"), p1, width = 1.3, height =1)


max.rqr_p = max(plot_df$rqr_p[is.finite(plot_df$rqr_p)]) + 1 # Best to cut inf?
plot_df$rqr_p[plot_df$rqr_p > max.rqr_p] = max.rqr_p
map =c("log(tt.umi)", "rqr_p", "post1")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3) +
    #geom_abline(slope = 1,intercept = 0) +
    #geom_abline(slope = 0,intercept = 0, color = 'grey') +
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
ggsave(paste0(plotPath, test_text, bc, "_rqrp.pdf"), g1, width = 1.5, height =1.8)


p1 = ggplot(plot_df[plot_df$post1 > .5,], aes_string(sample='rqr_p'))+
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
ggsave(paste0(plotPath, test_text, "_rqrp_qq_predicted.pdf"), p1, width = 1.2, height =1)



