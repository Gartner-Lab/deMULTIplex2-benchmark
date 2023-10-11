downsample_mtx <- function(mtx, ratio = .1, seed = 1) {
    set.seed(seed)
    ds_mtx <- t(apply(mtx, 1, function(x) {
        n <- floor(sum(x) * ratio)
        ds_reads <- sort(sample(seq_len(sum(x)), n))
        read_breaks <- c(0, cumsum(x))
        hist(ds_reads, breaks = read_breaks, plot = FALSE)$count
    }))
    colnames(ds_mtx) <- colnames(mtx)
    return(ds_mtx)
}

rm(list=ls())


library(deMULTIplex2)
savePath <- "./data-raw/benchmark_res/ds_McGinnis_8donor_SCMK/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

test_text <- "ds_McGinnis_8donor_SCMK_"


rna_mtx <- readRDS(paste0("./data-raw/benchmark_res/McGinnis_8donor_SCMK/rds/",  "McGinnis_8donor_SCMK_benchmarking_", "rna_mtx.rds"))
tag_mtx <- readRDS(paste0("./data-raw/benchmark_res/McGinnis_8donor_SCMK/rds/", "McGinnis_8donor_SCMK_benchmarking_", "tag_mtx.rds"))
true_label <- readRDS(paste0("./data-raw/benchmark_res/McGinnis_8donor_SCMK/rds/",  "McGinnis_8donor_SCMK_benchmarking_", "true_label.rds"))
tag_mapping <- readRDS(paste0("./data-raw/benchmark_res/McGinnis_8donor_SCMK/rds/",  "McGinnis_8donor_SCMK_benchmarking_", "tag_mapping.rds"))


#cur_ratio = 0.5
#cur_ratio = 0.1
#cur_ratio = 0.05
cur_ratio = 0.01

test_text <- paste0("ds_McGinnis_8donor_SCMK_",cur_ratio)

tag_mtx <- downsample_mtx(tag_mtx, ratio = cur_ratio)
head(tag_mtx)
bc_demux2 <- benchmark_demultiplex2(tag_mtx, true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'demultiplex2_'), width = 3.5, height = 2.5,
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
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))




# # No rna matrix simulated, just sample some random rna matrix
# hto12 <- readRDS(paste0(rdsPath, "hto12_processed.rds"))
# rna_mtx <- hto12@assays$RNA@counts
# rna_mtx <- rna_mtx[, sample(1:ncol(rna_mtx), nrow(tag_mtx), replace = T)]
# colnames(rna_mtx) <- rownames(tag_mtx)
# saveRDS(rna_mtx, paste0(rdsPath, test_text,"rna_mtx.rds"))


# 1. deMULTIplex
library(deMULTIplex)
library(ggplot2)
bc_demux1 <- benchmark_demultiplex1(as.matrix(tag_mtx), true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'deMULTIplex1_'))


# 2. demuxEM, failed
# Save data to required format

library(DropletUtils)
library(dplyr)
tag_mtx_write <- tag_mtx[rowSums(tag_mtx) > 0,]
rownames(tag_mtx_write) <- gsub("-1", "", rownames(tag_mtx_write))
tag_mtx_write <-as.data.frame(t(tag_mtx)) %>% add_rownames('tag')
write.csv(tag_mtx_write, paste0(rdsPath, test_text,"tag_mtx.csv"), quote = F, row.names = F)
rna_mtx_write<- rna_mtx[,rowSums(tag_mtx) > 0]
colnames(rna_mtx_write) <- gsub("-1", "", colnames(rna_mtx_write))
write10xCounts(paste0(rdsPath, test_text,"rna_mtx.h5"), rna_mtx_write, version='3')

# Run code on docker
paste0("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 100 --generate-diagnostic-plots ", test_text, "tag_mtx.csv ", test_text, "rna_mtx.h5 ", test_text, "demuxEM_res")

library(anndata)
adata <- read_h5ad(paste0(rdsPath, test_text, "demuxEM_res_demux.h5ad"))
demuxem_res <- adata$obs
demuxem_call <- as.character(demuxem_res$assignment)
demuxem_call[demuxem_res$demux_type!="singlet"] = as.character(demuxem_res$demux_type[demuxem_res$demux_type!="singlet"])
use_levels = sort(unique(demuxem_call))
use_levels <- c(use_levels[!use_levels %in% c("doublet", "unknown")], c("doublet", "unknown"))
demuxem_call <- factor(demuxem_call, levels = use_levels)
names(demuxem_call) <- paste0(rownames(demuxem_res))

bc_demuxEM = confusion_stats(demuxem_call, true_label,
                             tag_mapping,
                             call.multiplet = 'doublet',
                             true.multiplet = "doublet",
                             plot.path =plotPath, plot.name = paste0(test_text, 'demuxEM_'), width = 3.5, height = 2.5)


# 3. GMM-Demux
# The code to run this software is a bit cumbersome, generate script from R
# Step 1, save the tag matrix to the csv file
bcs <- make.names(colnames(tag_mtx))
tag_mtx_write <- as.matrix(tag_mtx)
colnames(tag_mtx_write) <- bcs
write.csv(tag_mtx_write, paste0(rdsPath, "tag_mtx_gmm.csv"), quote=F)
resPath <- paste0(savePath, test_text, "gmmdemux_res/"); dir.create(resPath)
cmds <- paste0("GMM-demux -c ./rds/tag_mtx_gmm.csv ", paste0(bcs, collapse = ","), " -x ", bcs, " -o ", test_text, "gmmdemux_res/", bcs)
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
                          plot.path =plotPath, plot.name = paste0(test_text, 'gmmdemux_'),  width = 3.5, height = 2.5)


# 4. HTODemux
bc_seu <- benchmark_HTODemux(tag_mtx, rna_mtx = t(tag_mtx), true_label = true_label,
                             tag_mapping = tag_mapping,
                             true.multiplet = "doublet",
                             plot.path = plotPath, plot.name = paste0(test_text, 'HTODemux_'), width = 3.5, height = 2.5)


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
                            plot.path =plotPath, plot.name = paste0(test_text, 'hashedDrops_'), width = 3.5, height = 2.5)

# 6. demuxmix
require(demuxmix)
bc_mixfull <- benchmark_demuxmix_full(tag_mtx, rna_mtx, true_label,
                                      tag_mapping = tag_mapping,
                                      true.multiplet = "doublet",
                                      plot.path = plotPath, plot.name = paste0(test_text, 'demuxmixfull_'), width = 3.5, height = 2.5)
bc_mixnaive <- benchmark_demuxmix_naive(tag_mtx, true_label = true_label,
                                        tag_mapping = tag_mapping,
                                        true.multiplet = "doublet",
                                        plot.path = plotPath, plot.name = paste0(test_text, 'demuxmixnaive_'), width = 3.5, height = 2.5)


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
                             plot.path =plotPath, plot.name = paste0(test_text, 'bffraw_'),  width = 3.5, height = 2.5)

bc_bffcluster <- confusion_stats(bff_cluster_call, true_label,
                                 tag_mapping,
                                 call.multiplet = 'Doublet',
                                 true.multiplet = "doublet",
                                 plot.path =plotPath, plot.name = paste0(test_text, 'bffcluster_'),  width = 3.5, height = 2.5)


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
    demuxmix = nsn("bc_mixfull"),
    demuxmix_naive = nsn("bc_mixnaive"),
    bff_raw = nsn("bc_bffraw"),
    bff_cluster = nsn("bc_bffcluster")
)
saveRDS(bc_all, paste0(rdsPath, test_text, "bc_all.rds"))


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
ggsave(paste0(plotPath, test_text, 'bc_all_bar.pdf'), g1, width = 8, height = 8)






show_ratios <- c(1, 0.5, 0.1, 0.05, 0.01)
test_text <- "ds_McGinnis_8donor_SCMK_"

bc_res <- lapply(show_ratios, function(cur_ratio) {
    if(cur_ratio == 1) {
        res = readRDS(paste0("./data-raw/benchmark_res/McGinnis_8donor_SCMK/rds/", "McGinnis_8donor_SCMK_benchmarking_bc_all.rds"))
        res$demuxEM <- NULL
    } else {
        res = readRDS(paste0(rdsPath, paste0("ds_McGinnis_8donor_SCMK_",cur_ratio), "bc_all.rds"))
    }
    return(res)
})
names(bc_res) <- as.character(show_ratios)

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
pheatmap(singlet_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.2f", number_color ='black', fontsize = 8, fontsize_number=8, color = get_gradient_color('RdBu', length(breaksList)), breaks = breaksList)
dev.off()




