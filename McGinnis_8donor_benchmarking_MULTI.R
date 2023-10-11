
library(deMULTIplex2)
savePath <- "./data-raw/benchmark_res/McGinnis_8donor_MULTI/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

test_text <- "McGinnis_8donor_MULTI_benchmarking_"


# Preprocess data
library(Seurat)
rna_mtx <- Read10X(paste0(savePath, 'data/AH_RNA'))
colnames(rna_mtx) <- gsub("-1", "", colnames(rna_mtx))

tag_mtx <- read.csv(paste0(savePath, 'data/AH_MULTI/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv'), row.names = 1)
tag_mtx <- tag_mtx[,!colnames(tag_mtx) %in% c("nUMI", "nUMI_total")]
tag_orig_assign <- read.csv(paste0(savePath, 'data/AH_MULTI/GSM4904942_8donor_PBMC_AH_MULTI_barcodes.csv'), row.names = 1)

# True label
souporcell <- read.table(paste0(savePath,"/souporcell/clusters.tsv"), sep="\t", header=T)
souporcell$barcode <- unlist(strsplit(as.character(souporcell$barcode), split="-1"))
rownames(souporcell) <- souporcell$barcode

sum(colnames(rna_mtx) %in% rownames(souporcell))
sum(rownames(tag_mtx) %in% rownames(souporcell))

# Use cells which have all three info
use_cells <- Reduce(intersect, list(colnames(rna_mtx), rownames(tag_mtx), rownames(souporcell)))
rna_mtx = rna_mtx[,use_cells]
tag_mtx = tag_mtx[use_cells,]
souporcell = souporcell[use_cells,]
tag_orig_assign <- tag_orig_assign[use_cells,,drop=F]

true_label = paste0("D", as.numeric(souporcell$assignment)+1)
true_label[souporcell$status != 'singlet'] = souporcell$status[souporcell$status != 'singlet']
names(true_label) = rownames(souporcell)
tag_mapping <- as.data.frame.matrix(table(true_label, tag_orig_assign$sample_classification))
tag_mapping <- apply(tag_mapping, 1, function(x)colnames(tag_mapping)[which.max(x)])
tag_mapping <- data.frame(tag = tag_mapping, true_label = names(tag_mapping))
tag_mapping <- tag_mapping[!rownames(tag_mapping) %in% c("doublet", "unassigned"), ]


saveRDS(rna_mtx, paste0(rdsPath, test_text, "rna_mtx.rds"))
saveRDS(tag_mtx, paste0(rdsPath, test_text, "tag_mtx.rds"))
saveRDS(true_label, paste0(rdsPath, test_text, "true_label.rds"))
saveRDS(tag_mapping, paste0(rdsPath, test_text, "tag_mapping.rds"))

# Create Seurat object of RNA
seu_rna <- CreateSeuratObject(counts = rna_mtx, min.features = 300)
# Normalize data
seu_rna <- NormalizeData(seu_rna)
seu_rna <- FindVariableFeatures(seu_rna, selection.method = "vst", nfeatures = 2000)
seu_rna <- ScaleData(seu_rna)
seu_rna <- RunPCA(seu_rna, features = VariableFeatures(object = seu_rna))
seu_rna <- FindNeighbors(seu_rna, dims = 1:10)
seu_rna <- FindClusters(seu_rna, resolution = 0.5)
seu_rna <- RunUMAP(seu_rna, dims = 1:10)
DimPlot(seu_rna, reduction = "umap", label = TRUE)
FeaturePlot(seu_rna, features = c("CD68" , "CD3D", "MS4A1", "NKG7"))
seu_rna$cell_type <- 'Myeloid cells'
seu_rna$cell_type[seu_rna$seurat_clusters %in% c(1)] = 'NK cells'
seu_rna$cell_type[seu_rna$seurat_clusters %in% c(2,7,5,3,4,9)] = 'T cells'
seu_rna$cell_type[seu_rna$seurat_clusters %in% c(13,8)] = 'B cells'
DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'cell_type')
saveRDS(seu_rna, paste0(rdsPath, test_text, "seu_rna.rds"))



set.seed(1)
umap_raw <- compute_umap(tag_mtx, use_dim = ncol(tag_mtx), n_component=2, n_neighbors = 30)
umap_raw <- as.data.frame(umap_raw)
umap_raw$true_label <- true_label[rownames(umap_raw)]
use_color = get_factor_color(tag_mapping$true_label, "Set1")
names(use_color) = tag_mapping$true_label
#umap_raw <- umap_raw[umap_raw$true_label!='doublet',] # Skip doublet
use_color <- c(use_color, 'doublet' = 'black', 'unassigned' = 'grey')
g1<- ggplot(umap_raw, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = "true_label"), stroke = 0, size = .5, raster.dpi=600) +
    scale_color_manual(values = use_color, na.value='lightgrey') +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    #guides(color = 'none') +
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

ggsave(paste0(plotPath, test_text, "_umap_tagraw.pdf"), g1, width = 2.6, height =1.7)


rna_mtx<- readRDS(paste0(rdsPath, test_text, "rna_mtx.rds"))
tag_mtx <- readRDS(paste0(rdsPath, test_text, "tag_mtx.rds"))
true_label <- readRDS(paste0(rdsPath, test_text, "true_label.rds"))
tag_mapping <- readRDS(paste0(rdsPath, test_text, "tag_mapping.rds"))

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
                                    max.quantile.fit = 0.9, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                    residual.type = "rqr",
                                    plot.umap = "none",
                                    plot.diagnostics = F)
saveRDS(bc_demux2, paste0(rdsPath, test_text, "bc_demux2.rds"))



bc_demux2 <- readRDS(paste0(rdsPath, test_text, "bc_demux2.rds"))
umap_rqr <- as.data.frame(bc_demux2$res$umap)
umap_rqr$tag_assign <- bc_demux2$res$assign_table$barcode_assign
umap_rqr$droplet_type <- bc_demux2$res$assign_table$droplet_type
umap_rqr$tag_plot <- tag_mapping$true_label[match(umap_rqr$tag_assign, tag_mapping$tag)]
umap_rqr$tag_plot[umap_rqr$droplet_type == 'multiplet'] = 'multiplet'
umap_rqr$tag_plot[umap_rqr$droplet_type == 'negative'] = 'negative'
use_color = get_factor_color(tag_mapping$true_label, "Set1")
names(use_color) = tag_mapping$true_label
use_color <- c(use_color, 'multiplet' = 'black', 'negative' = 'grey')
#use_color <- c(use_color, 'negative' = 'grey')
g1<- ggplot(umap_rqr, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(data = umap_rqr[umap_rqr$tag_plot == 'negative', ], aes_string(color = "tag_plot"), stroke = 0, size = .5, raster.dpi=600) +
    geom_point_rast(data = umap_rqr[umap_rqr$tag_plot != 'negative', ], aes_string(color = "tag_plot"), stroke = 0, size = .5, raster.dpi=600) +
    scale_color_manual(values = use_color, na.value='lightgrey') +
    theme_bw() +
    # guides(colour = guide_legend(override.aes = list(size=3),
    #                              keywidth=0.1,
    #                              keyheight=0.12,
    #                              default.unit="inch"))+
    guides(color = 'none') +
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
ggsave(paste0(plotPath, test_text, "_umap_rqr.pdf"), g1, width = 1.7, height =1.5)



acc_mtx<- bc_demux2$acc_mtx
rownames(acc_mtx)[which(rownames(acc_mtx) %in% tag_mapping$tag)] <- tag_mapping$true_label[match(rownames(acc_mtx)[which(rownames(acc_mtx) %in% tag_mapping$tag)], tag_mapping$tag)]
acc_mtx <- acc_mtx[!rownames(acc_mtx) %in% c("NA"), !colnames(acc_mtx) %in% c("unassigned")]
acc_mtx <- acc_mtx[order(rownames(acc_mtx)), ]
breaksList <- seq(0,300,by=1e1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "accmtx", ".pdf"), width = 3, height = 2)
pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.0f", number_color ='white', fontsize = 8, fontsize_number=8, color = get_gradient_color('viridis', length(breaksList)), breaks = breaksList)
dev.off()



# 1. deMULTIplex
bc_demux1 <- benchmark_demultiplex1(as.matrix(tag_mtx), true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'deMULTIplex1_'))


# 2. demuxEM, failed
# Save data to required format
library(DropletUtils)
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
#adata <- read_h5ad(paste0(rdsPath, test_text, "demuxEM_res_demux.h5ad"))
adata <- read_h5ad(paste0(rdsPath, "experiment1_human_st_demux.h5ad"))
demuxem_res <- adata$obs
demuxem_call <- as.character(demuxem_res$assignment)
demuxem_call[demuxem_res$demux_type!="singlet"] = as.character(demuxem_res$demux_type[demuxem_res$demux_type!="singlet"])
use_levels = sort(unique(demuxem_call))
use_levels <- c(use_levels[!use_levels %in% c("doublet", "unknown")], c("doublet", "unknown"))
demuxem_call <- factor(demuxem_call, levels = use_levels)
names(demuxem_call) <- paste0(rownames(demuxem_res), "-1")

bc_demuxEM = confusion_stats(demuxem_call, true_label,
                             tag_mapping,
                             call.multiplet = 'doublet',
                             true.multiplet = "doublet",
                             plot.path =plotPath, plot.name = paste0(test_text, 'demuxEM_'), width = 3.5, height = 2.5)


# 3. GMM-Demux
# The code to run this software is a bit cumbersome, generate script from R
# Step 1, save the tag matrix to the csv file
bcs <- make.names(colnames(tag_mtx))
tag_mtx_write <- tag_mtx
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
                                            doUMAP = FALSE, doHeatmap = FALSE)
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
saveRDS(bc_all, paste0(rdsPath, test_text, "bc_all.rds"))

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
saveRDS(call_labels, paste0(rdsPath, test_text, "call_labels.rds"))



# Plot summary plots
summary_df <- sapply(bc_all, function(x) {
    x$singlet_avg_stats
})
write.csv(round(t(summary_df),3), paste0(plotPath, test_text,"summary_df.csv"))
plot_df <- reshape2::melt(summary_df)
colnames(plot_df) <- c("stats", "method", "value")
g1<-ggplot(plot_df, aes(method, value)) +
    geom_col() +
    facet_wrap(~stats, scales = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(plotPath, test_text, 'bc_all_bar.pdf'), g1, width = 8, height = 8)






bc_all <- readRDS(paste0(rdsPath, test_text, "bc_all.rds"))

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
use_color <- gg_color_hue(length(bc_all))
names(use_color) <- names(bc_all)
g1 <- ggplot(plot_df, aes(x = Method, y = value, fill = Method)) +
    geom_boxplot(lwd=.2, outlier.shape = NA) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1),dotsize = .9,stroke = 0.2)+
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
ggsave(paste0(plotPath, test_text, 'bc_precision_recall_boxplot.pdf'), g1, width =2, height = 2)



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
ggsave(paste0(plotPath, test_text, 'bc_f_score_bar.pdf'), g1, width = 2, height = 1.8)


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
ggsave(paste0(plotPath, test_text, 'bc_doublet_assign.pdf'), g1, width = 2.5, height = 1.8)







# Individual barcode stats
rna_mtx <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "rna_mtx.rds"))
tag_mtx <- readRDS(paste0(rdsPath, "McGinnis_8donor_MULTI_benchmarking_", "tag_mtx.rds"))
true_label <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "true_label.rds"))
tag_mapping <- readRDS(paste0(rdsPath,  "McGinnis_8donor_MULTI_benchmarking_", "tag_mapping.rds"))
bc_all <- readRDS(paste0(rdsPath, "McGinnis_8donor_MULTI_benchmarking_",  "bc_all.rds"))


bc_demux2 <- readRDS(paste0(rdsPath, test_text, "bc_demux2.rds"))
bc = 'Bar11'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label[rownames(plot_df)]
map =c("log(tt.umi)", "log(bc.umi)", "true_label")
td <- tag_mapping$true_label[tag_mapping$tag == bc]
plot_df$true_label_bin <- ifelse(plot_df$true_label == td,td, "Other")

use_color <- c(
    '#e41a1c',
    '#bdbdbd'
)
names(use_color) <- c(td, "Other")

map = c("log(bc.umi)", "cos.umi", "true_label_bin")
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
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))

g1 <- ggMarginal(g1,type = "histogram")

ggsave(paste0(plotPath, test_text, bc, "_cosine_is_pos.pdf"), g1, width = 2.7, height =2)



map =c("log(tt.umi)", "log(bc.umi)", "true_label_bin")
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
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, test_text, bc, "_true_label.pdf"), g1, width = 2.5, height =1.8)



map =c("log(tt.umi)", "log(bc.umi)", "post1")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'grey') +
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
ggsave(paste0(plotPath, test_text, bc, "_post1.pdf"), g1, width = 1.4, height =1.6)



map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "post1")
g2 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = map[3]), stroke = 0, size = .3) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred1)'),color = 'grey') +
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
ggsave(paste0(plotPath, test_text, bc, "_post1_reverse.pdf"), g2, width = 1.4, height =1.6)



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






# ROC analysis
bc_demux2 <- readRDS(paste0(rdsPath, test_text, "bc_demux2.rds"))

true_label_tag <- tag_mapping$tag[match(true_label, tag_mapping$true_label)]
true_label_tag[is.na(true_label_tag)] = "n/d"
names(true_label_tag) <- names(true_label)

dmmFull <- demuxmix(as.matrix(t(tag_mtx)), rna = Matrix::colSums(rna_mtx > 0))
full_post <- t(dmmFull@posteriorProb)

dmmNaive <- demuxmix(as.matrix(t(tag_mtx)), model = "naive")
naive_post <- t(dmmNaive@posteriorProb)

full_roc <- prob_to_mean_roc(full_post, true_label_tag)
full_roc$method <- "demuxmix_full"
print(compute_auc(full_roc$mean_fpr,full_roc$mean_tpr))

naive_roc <- prob_to_mean_roc(naive_post, true_label_tag)
naive_roc$method <- "demuxmix_naive"
print(compute_auc(naive_roc$mean_fpr,naive_roc$mean_tpr))

demux2_roc <- prob_to_mean_roc(bc_demux2$res$prob_mtx, true_label_tag[rownames(bc_demux2$res$prob_mtx)])
demux2_roc$method <- "deMULTIplex2"
print(compute_auc(demux2_roc$mean_fpr,demux2_roc$mean_tpr))

roc_bind <- rbind(demux2_roc, full_roc, naive_roc)
saveRDS(roc_bind, paste0(rdsPath, test_text, "roc_bind.rds"))



method_color <- c(
    "deMULTIplex2" = "#F8766D",
    "demuxmix_full" = "#00B6EB",
    "demuxmix_naive" = "#9590FF"
)
g1 <- ggplot(roc_bind, aes(mean_fpr, mean_tpr, color = method)) +
    geom_line(size = .5)+
    scale_color_manual(values = method_color) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black") +
    guides(colour = guide_legend(override.aes = list(size=3),
                                 keywidth=0.1,
                                 keyheight=0.12,
                                 default.unit="inch"))+
    theme_bw() +
    theme(text=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_blank(),
          axis.text = element_text(size=8),
          legend.margin=margin(0,0,0,0),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_roc.pdf"), g1, width = 2.5, height =1)






