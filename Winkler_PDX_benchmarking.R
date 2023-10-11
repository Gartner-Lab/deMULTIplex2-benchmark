





library(deMULTIplex2)
savePath <- "./data-raw/benchmark_res/Winkler_PDX/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)

test_text <- "Winkler_PDX_benchmarking_"

# Preprocess data
library(Seurat)
rna_mtx <- fread("./data-raw/benchmark_res/Winkler_PDX/data/GSE211145_MULTI_raw_counts.csv")
# Convert to sparse matrix and save
cell_barcode <- rna_mtx[,V1]
rna_mtx <- rna_mtx[,V1:=NULL]
gene_name = colnames(rna_mtx)
rna_mtx <- Matrix(t(as.matrix(rna_mtx)), sparse = T)
rownames(rna_mtx) <- gene_name
colnames(rna_mtx) <- cell_barcode
saveRDS(rna_mtx, paste0(rdsPath, test_text, "rna_mtx_global.rds"))
cell_meta <- read.csv("./data-raw/benchmark_res/Winkler_PDX/data/GSE211145_MULTI_metadata.csv", row.names = 1)
tag_meta <- read.csv("./data-raw/benchmark_res/Winkler_PDX/data/PDX_MULTI-SEQ_Metadata.csv", row.names = 1)
tag_mapping_global <- tag_meta[,c("MULTI", "Tumor")]
colnames(tag_mapping_global) <- c("tag", "true_label")

tag_mtx_pdx1 <- read.csv("./data-raw/benchmark_res/Winkler_PDX/data/GSE211145_MULTI_barcode_count_PDX1.csv", row.names = 1)
tag_mtx_pdx2 <- read.csv("./data-raw/benchmark_res/Winkler_PDX/data/GSE211145_MULTI_barcode_count_PDX2.csv", row.names = 1)
tag_mtx_pdx2 <- tag_mtx_pdx2[,!grepl('nUMI',colnames(tag_mtx_pdx2))]
tag_mtx_pdx3 <- read.csv("./data-raw/benchmark_res/Winkler_PDX/data/GSE211145_MULTI_barcode_count_PDX3.csv", row.names = 1)
tag_mtx_pdx3 <- tag_mtx_pdx3[,!grepl('nUMI',colnames(tag_mtx_pdx3))]
use_bc <- names(which(colSums(tag_mtx_pdx3) > 1e4))
tag_mtx_pdx3 <- tag_mtx_pdx3[,use_bc]
sum(rownames(tag_mtx_pdx1) %in% cell_barcode)
sum(rownames(tag_mtx_pdx2) %in% cell_barcode)
sum(rownames(tag_mtx_pdx3) %in% cell_barcode)

# Try initial deMULTIplex2
res_pdx1 <- demutiplexTags(tag_mtx_pdx1,
                      plot.diagnostics = T,
                      plot.path = plotPath,
                      plot.name = paste0(test_text, 'pdx1'))
res_pdx2 <- demutiplexTags(tag_mtx_pdx2,
                           plot.diagnostics = T,
                           plot.path = plotPath,
                           plot.name = paste0(test_text, 'pdx2'))
res_pdx3 <- demutiplexTags(tag_mtx_pdx3,
                           plot.diagnostics = T,
                           plot.path = plotPath,
                           plot.name = paste0(test_text, 'pdx3'))



# Create Seurat object of RNA
seu_rna <- CreateSeuratObject(counts = rna_mtx, min.features = 0)
seu_rna <- NormalizeData(seu_rna)
seu_rna <- FindVariableFeatures(seu_rna, selection.method = "vst", nfeatures = 2000)
seu_rna <- ScaleData(seu_rna)
seu_rna <- RunPCA(seu_rna, features = VariableFeatures(object = seu_rna))
seu_rna <- FindNeighbors(seu_rna, dims = 1:10)
seu_rna <- FindClusters(seu_rna, resolution = 0.5)
seu_rna <- RunUMAP(seu_rna, dims = 1:10)
seu_rna@meta.data <- cbind(seu_rna@meta.data, cell_meta[colnames(seu_rna),])
DimPlot(seu_rna, reduction = "umap", label = TRUE)

DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'MULTI_barcode')
DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'MULTI_tumor_ID')
DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'Tumor_ID')
DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'run_id')
DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'sequencing_batch')
DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'leiden')
DimPlot(seu_rna, reduction = "umap", label = TRUE, group.by = 'louvain')

saveRDS(seu_rna, paste0(rdsPath, test_text, "seu_rna.rds"))

table(seu_rna$Tumor_ID, seu_rna$MULTI_barcode)



seu_rna <- readRDS(paste0(rdsPath,"Winkler_PDX_benchmarking_", "seu_rna.rds"))

# Benchmark
# test_text <- paste0("Winkler_PDX_benchmarking_", "pdx1_")
# tag_mtx = tag_mtx_pdx1

# test_text <- paste0("Winkler_PDX_benchmarking_", "pdx2_")
# tag_mtx = tag_mtx_pdx2

test_text <- paste0("Winkler_PDX_benchmarking_", "pdx3_")
tag_mtx = tag_mtx_pdx3


cur_seu <- seu_rna[,match(rownames(tag_mtx), colnames(seu_rna))]
true_label <- cur_seu$Tumor_ID
names(true_label) = rownames(tag_mtx)
rna_mtx <- cur_seu@assays$RNA@counts
tag_mapping = tag_mapping_global[tag_mapping_global$tag %in% colnames(tag_mtx),]



set.seed(1)
umap_raw <- compute_umap(tag_mtx, use_dim = ncol(tag_mtx), n_component=2, n_neighbors = 30)
umap_raw <- as.data.frame(umap_raw)
umap_raw$true_label <- true_label[rownames(umap_raw)]
#use_color = get_factor_color(tag_mapping$true_label, "Set1")
#names(use_color) = tag_mapping$true_label
use_color = get_factor_color(unique(tag_mapping_global$true_label), "Set1")
names(use_color) <- unique(tag_mapping_global$true_label)
#umap_raw <- umap_raw[umap_raw$true_label!='doublet',] # Skip doublet
use_color <- c(use_color, 'doublet' = 'black', 'unassigned' = 'grey')
g1<- ggplot(umap_raw, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = "true_label"), stroke = 0, size = .3, raster.dpi=600) +
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

ggsave(paste0(plotPath, test_text, "_umap_tagraw.pdf"), g1, width = 2.1, height =1.8)


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
rownames(acc_mtx)[which(rownames(acc_mtx) %in% tag_mapping$tag)] <- paste0(tag_mapping$true_label, "-", tag_mapping$tag)[match(rownames(acc_mtx)[which(rownames(acc_mtx) %in% tag_mapping$tag)], tag_mapping$tag)]
#acc_mtx <- acc_mtx[!rownames(acc_mtx) %in% c("NA"), !colnames(acc_mtx) %in% c("unassigned")]
acc_mtx <- acc_mtx[order(rownames(acc_mtx)), order(colnames(acc_mtx))]
breaksList <- seq(0,1000,by=1e1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "accmtx", ".pdf"), width = 3, height = 2)
pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.0f", number_color ='white', fontsize = 8, fontsize_number=8, color = get_gradient_color('viridis', length(breaksList)), breaks = breaksList)
dev.off()



bc_demux2 <- readRDS(paste0(rdsPath, test_text, "bc_demux2.rds"))
demux2_res <- bc_demux2$res
demux2_res$final_assign <- demux2_res$assign_table$barcode_assign
demux2_res$final_assign[is.na(demux2_res$final_assign)] <- demux2_res$assign_table$droplet_type[is.na(demux2_res$final_assign)]
names(demux2_res$final_assign) <- rownames(demux2_res$assign_table)
demux2_knn_res <-knn_refinement(demux2_res, k = 3)
bc_demux2_3nn <- confusion_stats(demux2_knn_res$final_assign, true_label,
                                 tag_mapping,
                                 call.multiplet = 'multiplet',
                                 true.multiplet = "doublet",
                                 plot.path =plotPath, plot.name = paste0(test_text, 'demux2_3nn_'),  width = 3.5, height = 2.5)


# 1. deMULTIplex
bc_demux1 <- benchmark_demultiplex1(as.matrix(tag_mtx), true_label,
                                    tag_mapping = tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = plotPath, plot.name = paste0(test_text, 'deMULTIplex1_'))


# 2. demuxEM
# Save data to required format
library(DropletUtils)
tag_mtx_write <- tag_mtx[rowSums(tag_mtx) > 0,]
tag_mtx_write <-as.data.frame(t(tag_mtx)) %>% add_rownames('tag')
write.csv(tag_mtx_write, paste0(rdsPath, test_text,"tag_mtx.csv"), quote = F, row.names = F)
rna_mtx_write<- rna_mtx[,rowSums(tag_mtx) > 0]
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
bc_gmm <- confusion_stats(gmm_res_call, true_label,
                          tag_mapping,
                          call.multiplet = 'Multiplet',
                          true.multiplet = "doublet",
                          plot.path =plotPath, plot.name = paste0(test_text, 'gmmdemux_'),  width = 3.5, height = 2.5)


# 4. HTODemux
bc_seu <- benchmark_HTODemux(tag_mtx, rna_mtx = rna_mtx, true_label = true_label,
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



bc_all <- readRDS(paste0(rdsPath, test_text, "bc_all.rds"))

# Plot stats of each barcode
bc_stats <- lapply(1:length(bc_all), function(i) {
    x = bc_all[[i]]
    if(!is.null(x)) {
        df <- as.data.frame(do.call(rbind, x$tag_stats))
        df$tag_id <- rownames(df)
        df$method =names(bc_all)[i]
    } else {
        df <- NULL
    }
    return(df)
})
bc_stats <- do.call(rbind, bc_stats[sapply(bc_stats, function(x)!is.null(x))])

plot_df <- bc_stats[,c('precision', 'recall', 'method')]
colnames(plot_df) <- c('Precision', 'Recall', 'Method')
plot_df <- reshape2::melt(plot_df, id.var = 'Method')
plot_df$Method <- factor(plot_df$Method, levels = names(bc_all))
use_color <- gg_color_hue(length(bc_all))
names(use_color) <- names(bc_all)
g1 <- ggplot(plot_df, aes(x = Method, y = value, fill = Method)) +
    geom_boxplot(lwd=.2, outlier.shape = NA) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1),dotsize = 1.3,stroke = 0.5)+
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
    if(!is.null(x)) {
        x$singlet_avg_stats
    } else {
        return(c(0,0,0))
    }
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








# GLOBAL analysis on all results

test_text <- "Winkler_PDX_benchmarking_"
call_pdx1 <- readRDS(paste0(rdsPath, test_text, "pdx1_", "call_labels.rds"))
call_pdx2 <- readRDS(paste0(rdsPath, test_text, "pdx2_", "call_labels.rds"))
call_pdx3 <- readRDS(paste0(rdsPath, test_text, "pdx3_", "call_labels.rds"))

seu_rna <- readRDS(paste0(rdsPath,"Winkler_PDX_benchmarking_", "seu_rna.rds"))

tag_meta <- read.csv("./data-raw/benchmark_res/Winkler_PDX/data/PDX_MULTI-SEQ_Metadata.csv", row.names = 1)
tag_mapping_global <- tag_meta[,c("MULTI", "Tumor")]
colnames(tag_mapping_global) <- c("tag", "true_label")
tag_mapping_global$true_label[grepl("naïve", tag_mapping_global$true_label)] = "naïve"
# For each method get the single assignment vector
bc_methods = c("paper",names(call_pdx1))
assign_res_bc <-sapply(bc_methods, function(x) {
    if(x == "paper") {
        y = seu_rna$MULTI_barcode
    } else {
        y <- c(call_pdx1[[x]], call_pdx2[[x]], call_pdx3[[x]])
        if(is.null(y)) y = NA
        y = as.character(y[colnames(seu_rna)])
    }
     y[y %in% c('Multiplet', 'multiplet', 'Doublet', 'doublet')] = 'Doublet'
     y[y %in% c('Negative', 'negative', 'unknown', 'uncertain', 'ambiguous') | is.na(y)] = 'Unclassified'
     #y[which(y %in% tag_mapping_global$tag)] = tag_mapping_global$true_label[match(y[which(y %in% tag_mapping_global$tag)], tag_mapping_global$tag)]
    return(y)
})
assign_res_all <- apply(assign_res_bc, 2, function(y) {
    y[which(y %in% tag_mapping_global$tag)] = tag_mapping_global$true_label[match(y[which(y %in% tag_mapping_global$tag)], tag_mapping_global$tag)]
    return(y)
})
saveRDS(as.data.frame(assign_res_bc), paste0(rdsPath, test_text, "assign_res_bc.rds"))
saveRDS(assign_res_all, paste0(rdsPath, test_text, "assign_res_all.rds"))


assign_res_all <- readRDS(paste0(rdsPath, test_text, "assign_res_all.rds"))
umap_rna <- as.data.frame(cbind(seu_rna@reductions$umap@cell.embeddings, seu_rna@meta.data, assign_res_all))

use_color <- get_factor_color(unique(tag_mapping_global$true_label), "Set1")
names(use_color) <- unique(tag_mapping_global$true_label)
use_color <- c(use_color, 'Doublet' = '#636363', 'Unclassified' = 'black')

show_method = "paper"
show_method = "Tumor_ID"
#show_method = "deMULTIplex2"
#show_method = "deMULTIplex"
g1<- ggplot(umap_rna, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = show_method), stroke = 0, size = .4, raster.dpi=600) +
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

ggsave(paste0(plotPath, test_text, show_method, "_umap.pdf"), g1, width = 2.8, height =1.9)

umap_rna$deMULTIplex_bin <- umap_rna$deMULTIplex; umap_rna$deMULTIplex_bin[!umap_rna$deMULTIplex_bin %in% c("Unclassified", "Doublet")] = "Classified"
umap_rna$deMULTIplex2_bin <- umap_rna$deMULTIplex2; umap_rna$deMULTIplex2_bin[!umap_rna$deMULTIplex2_bin %in% c("Unclassified", "Doublet")] = "Classified"

#show_method = "deMULTIplex2_bin"
#show_method = "deMULTIplex_bin"
use_color <-  c('Doublet' = '#636363', 'Unclassified' = 'black',  "Classified" = '#7bccc4')
g1<- ggplot(umap_rna, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = show_method), stroke = 0, size = .4, raster.dpi=600) +
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
ggsave(paste0(plotPath, test_text, show_method, "_umap.pdf"), g1, width = 2.8, height =1.9)


umap_rna$sequencing_batch <- paste0("batch_", umap_rna$sequencing_batch)
plot_col<- "sequencing_batch"
g1<- ggplot(umap_rna, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = plot_col), stroke = 0, size = .2, raster.dpi=600) +
    #scale_color_manual(values = use_color, na.value='lightgrey') +
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

ggsave(paste0(plotPath, test_text, plot_col, "_umap.pdf"), g1, width = 2.6, height =1.9)




# For each of the tumor show recovered proportion
assign_res_all <- readRDS(paste0(rdsPath, test_text, "assign_res_all.rds"))
#assign_res_bc <- readRDS(paste0(rdsPath, test_text, "assign_res_bc.rds"))
true_label <- seu_rna$Tumor_ID

call_df<- as.data.frame(apply(assign_res_all,2,function(y){
    y[y == true_label] = 'Correct'
    y[y != true_label & !y %in% c("Correct", "Doublet", "Unclassified")] = 'Incorrect'
    return(y)
}))

call_df_bin <- as.data.frame(call_df=="Correct")

call_df$Tumor_ID <- true_label
call_df_bin$Tumor_ID <- true_label
library(data.table)
call_df_bin <- as.data.table(call_df_bin)

call_df_bin[, lapply(.SD, mean, na.rm=TRUE), by=Tumor_ID ]

plot_df <- reshape2::melt(call_df, id.var = "Tumor_ID")
colnames(plot_df) <- c("cell_type", "method", "result")

use_color <-  c('Unclassified' = '#737373', "Doublet" = '#4daf4a',  "Incorrect" = '#e41a1c', "Correct" = '#377eb8')
plot_df$result <- factor(plot_df$result, levels = names(use_color))
plot_df <- plot_df[!plot_df$method %in% c("paper", "demuxmix", "demuxmix_naive"),]
g1 <- ggplot(plot_df, aes(x = cell_type, fill = result)) +
    geom_bar(position="fill") +
    facet_wrap(~method, nrow = 2)+
    ylab('Fraction') +
    xlab(NULL) +
    scale_fill_manual(values = use_color) +
    #guides(fill = guide_legend(override.aes = list(width = .1)))+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,1)) +
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
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
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_tumor_assign_stats_show.pdf"), g1, width = 5, height =2.3)


# For each of the tumor show recovered proportion, version 2
assign_res_all <- readRDS(paste0(rdsPath, test_text, "assign_res_all.rds"))
#assign_res_bc <- readRDS(paste0(rdsPath, test_text, "assign_res_bc.rds"))
true_label <- seu_rna$Tumor_ID

call_df<- as.data.frame(apply(assign_res_all,2,function(y){
    y[y == true_label] = 'Correct'
    y[y != 'Correct'] = 'Un/Mis-classified'
    return(y)
}))

call_df$Tumor_ID <- true_label
sapply(call_df, function(x) {
    mean(x == "Correct")
})

plot_df <- reshape2::melt(call_df, id.var = "Tumor_ID")
colnames(plot_df) <- c("cell_type", "method", "result")

use_color <-  c('Un/Mis-classified' = '#999999', "Correct" = '#005A9C')
plot_df$result <- factor(plot_df$result, levels = names(use_color))
plot_df <- plot_df[!plot_df$method %in% c("paper", "demuxmix", "demuxmix_naive"),]
g1 <- ggplot(plot_df, aes(x = cell_type, fill = result)) +
    geom_bar(position="fill") +
    facet_wrap(~method, nrow = 2)+
    ylab('Fraction') +
    xlab(NULL) +
    scale_fill_manual(values = use_color) +
    #guides(fill = guide_legend(override.aes = list(width = .1)))+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,1)) +
    theme(text=element_text(size=8),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
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
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
ggsave(paste0(plotPath, test_text, "_tumor_assign_stats_show_v2.pdf"), g1, width = 6, height =2.3)





assign_res_all <- readRDS(paste0(rdsPath, test_text, "assign_res_all.rds"))
assign_res_bc <- readRDS(paste0(rdsPath, test_text, "assign_res_bc.rds"))
assign_res_cbn <- sapply(colnames(assign_res_bc), function(x){
    x <- ifelse(!assign_res_bc[,x] %in% c("Unclassified", "Doublet"), paste0(assign_res_bc[,x], "_", assign_res_all[,x]), assign_res_bc[,x])
    return(x)
})
rownames(assign_res_cbn) <- rownames(assign_res_bc)
umap_plot <- as.data.frame(cbind(umap_raw, seu_rna@meta.data[rownames(umap_raw),], assign_res_cbn[rownames(umap_raw), ]))
#show_method = "paper"
#show_method = "deMULTIplex2"
#show_method = "deMULTIplex"
tag_mapping = tag_mapping_global[tag_mapping_global$tag %in% colnames(tag_mtx_pdx3),]
unq_labels = gtools::mixedsort(unique(paste0(tag_mapping$tag, "_", tag_mapping$true_label)))
use_color <- get_factor_color(unq_labels, "Set1")
names(use_color) <- unq_labels
use_color <- c(use_color, 'Doublet' = '#bdbdbd', 'Unclassified' = '#bdbdbd')

umap_plot[[show_method]] <- factor(umap_plot[[show_method]], levels = names(use_color))
g1<- ggplot(umap_plot, aes_string("UMAP_1", "UMAP_2")) +
    geom_point_rast(aes_string(color = show_method), stroke = 0, size = .2, raster.dpi=600) +
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
          legend.position="top",
          axis.text = element_text(size=8),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))

ggsave(paste0(plotPath, test_text, show_method, "_umapzoom_bc.pdf"), g1, width = 2.6, height =1.7)



#
# seu_J53 = seu_J53[,seu_J53$Tumor_ID == "J53353"]
# seu_J53 <- FindVariableFeatures(seu_J53, selection.method = "vst", nfeatures = 2000)
# seu_J53 <- ScaleData(seu_J53)
# seu_J53 <- RunPCA(seu_J53, features = VariableFeatures(object = seu_J53))
# seu_J53 <- FindNeighbors(seu_J53, dims = 1:10)
# seu_J53 <- FindClusters(seu_J53, resolution = 0.5)
# seu_J53 <- RunUMAP(seu_J53, dims = 1:10)
# DimPlot(seu_J53, reduction = "umap", label = TRUE)
#
# umap_J53 <- as.data.frame(cbind(seu_J53@reductions$umap@cell.embeddings, seu_J53@meta.data, assign_res_all[colnames(seu_J53),]))
#
#
# show_method = "paper"
# #show_method = "deMULTIplex2"
# #show_method = "deMULTIplex"
# g1<- ggplot(umap_J53, aes_string("UMAP_1", "UMAP_2")) +
#     geom_point_rast(aes_string(color = show_method), stroke = 0, size = .2, raster.dpi=600) +
#     scale_color_manual(values = use_color, na.value='lightgrey') +
#     theme_bw() +
#     guides(colour = guide_legend(override.aes = list(size=3),
#                                  keywidth=0.1,
#                                  keyheight=0.12,
#                                  default.unit="inch"))+
#     #guides(color = 'none') +
#     theme(text=element_text(size=8),
#           legend.text=element_text(size=8),
#           legend.title = element_blank(),
#           axis.text = element_text(size=8),
#           axis.ticks.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.y=element_blank(),
#           axis.text.y=element_blank(),
#           legend.margin=margin(0,0,0,0),
#           #legend.box.margin=margin(-10,-10,-10,-10),
#           plot.margin = unit(c(0,0.4,0,0.1), "cm"))
#
# ggsave(paste0(plotPath, test_text, show_method, "_J53umapredo.pdf"), g1, width = 3, height =2)



# Summary stats
assign_res_all <- readRDS(paste0(rdsPath, test_text, "assign_res_all.rds"))
seu_rna <- readRDS(paste0(rdsPath,"Winkler_PDX_benchmarking_", "seu_rna.rds"))
true_label_cbn <- seu_rna$Tumor_ID
names(true_label_cbn) <- colnames(seu_rna)
true_label_cbn <-true_label_cbn[rownames(assign_res_all)]
tag_mapping_corrected <- tag_mapping_global
tag_mapping_corrected$tag = tag_mapping_corrected$true_label
tag_mapping_corrected <- tag_mapping_corrected[!duplicated(tag_mapping_corrected),]
bc_cbn <- lapply(colnames(assign_res_all), function(x) {
    call_label <- assign_res_all[,x]
    confusion_stats(call_label, true_label_cbn,
                    tag_mapping_corrected,
                    call.multiplet = 'Doublet',
                    true.multiplet = "Doublet",
                    plot.path =plotPath, plot.name = paste0(test_text, x, '_cbn_'), width = 3.5, height = 2.5)
})
names(bc_cbn) <- colnames(assign_res_all)

saveRDS(bc_cbn, paste0(rdsPath, test_text, "bc_cbn.rds"))

summary_df <- sapply(bc_cbn, function(x) {
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









# Barcode case study

test_text <- paste0("Winkler_PDX_benchmarking_", "pdx3_")
tag_mtx = tag_mtx_pdx3


cur_seu <- seu_rna[,match(rownames(tag_mtx), colnames(seu_rna))]
true_label <- cur_seu$Tumor_ID
names(true_label) = rownames(tag_mtx)
rna_mtx <- cur_seu@assays$RNA@counts
tag_mapping = tag_mapping_global[tag_mapping_global$tag %in% colnames(tag_mtx),]
bc_demux2 <- readRDS(paste0(rdsPath, test_text, "bc_demux2.rds"))



bc = 'Bar2'
plot_df <- bc_demux2$res$df_list[[bc]]
plot_df$true_label <- true_label
map =c("log(tt.umi)", "log(bc.umi)", "true_label")

plot_df$call_label <- paste0("call_", bc_demux2$res$assign_table$barcode_assign == bc)
plot_df$true_label_bin <- ifelse(plot_df$true_label == "HCI011", "HCI011", "Other")

use_color <- c(
    'HCI011' = '#8D4C6A',
    'Other' = '#bdbdbd'
)


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
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
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
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, test_text, bc, "_true_label.pdf"), g1, width = 2.5, height =1.8)







map =c("log(tt.umi)", "log(bc.umi)", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .3) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred0)'),color = 'grey') +
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
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
#g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, test_text, bc, "_post1.pdf"), g1, width = 1.6, height =1.9)



map =c("log(tt.umi)", "log(tt.umi-bc.umi)", "true_label")
g2 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .3) +
    geom_abline(slope = 1,intercept = 0) +
    geom_line(aes_string(map[1], 'log(pred1)'),color = 'grey') +
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
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
#g2 <- ggMarginal(g2, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, test_text, bc, "_post1_reverse.pdf"), g2, width = 1.6, height =1.9)



max.rqr = max(plot_df$rqr[is.finite(plot_df$rqr)]) + 1 # Best to cut inf?
plot_df$rqr[plot_df$rqr > max.rqr] = max.rqr
map =c("log(tt.umi)", "rqr", "true_label")
g1 <- ggplot(plot_df) +
    geom_point(aes_string(map[1], map[2], color = 'post1'), stroke = 0, size = .3) +
    geom_abline(slope = 1,intercept = 0) +
    geom_abline(slope = 0,intercept = 0, color = 'grey') +
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
          # axis.ticks.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.y=element_blank(),
          # axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm"))
g1 <- ggMarginal(g1, margins = 'y', type = "histogram")
ggsave(paste0(plotPath, test_text, bc, "_rqr.pdf"), g1, width = 1.6, height =2)



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
ggsave(paste0(plotPath, test_text, "_rqr_qq_predicted.pdf"), p1, width = 1.2, height =1)


