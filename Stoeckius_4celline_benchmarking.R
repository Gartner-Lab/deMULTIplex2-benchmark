

# The processed hto object was obtained following the Seurat tutorial: https://satijalab.org/seurat/articles/hashing_vignette


library(deMULTIplex2)
library(Matrix)
savePath <- "./data-raw/benchmark_res/Stoeckius_4cellline/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
dataPath <- paste0(savePath, "data/"); dir.create(plotPath)
test_text <- "benchmark_Stoeckius_4cellline_"

# Data preprocessing to obtain tag_mtx, rna_mtx, true_label for benchmarking
hto12 <- readRDS(paste0(rdsPath, "hto12_processed.rds"))

tag_mtx <- t(hto12@assays$HTO@counts)
rna_mtx <- hto12@assays$RNA@counts

true_label = hto12$assign_rna
true_label[true_label == "Doublet_rna"] = "doublet"
names(true_label) <- colnames(hto12)

tag_mapping <- data.frame(
    tag = colnames(tag_mtx)
)
tag_mapping$true_label <- sapply(strsplit(tag_mapping$tag, "-"), function(x)x[1])



# Some stats plot
plot_df<- data.frame(tt.rna.umi = hto12$nCount_RNA, tt.rna.gene = hto12$nFeature_RNA, tt.tag.umi = hto12$nCount_HTO, cell.type = hto12$assign_rna)
plot_df <- plot_df[!plot_df$cell.type == "Doublet_rna",]
g1 = ggplot(plot_df, aes(log(tt.rna.gene), log(tt.tag.umi),color= cell.type)) +
    geom_point_rast(aes(color= cell.type), size = .5, stroke = 0, alpha = .4) +
    geom_smooth(method = "lm") +
    ylim(c(3,7))+
    theme(text=element_text(family = "Helvetica", size=8),
          legend.text=element_text(size=8),
          axis.text = element_text(size=9),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(0,0.4,0,0.1), "cm")) +
    theme_bw()
g1 <- ggMarginal(g1, groupColour = TRUE, groupFill = TRUE)

ggsave(paste0(plotPath, test_text, "tag_vs_gene.pdf"), g1, width = 3.5, height = 2.5)





#Methods to be compared against:
#deMULTIplex1,
#demuxEM, https://www.nature.com/articles/s41467-019-10756-2#Sec8
#GMM-Demux, https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02084-2
#HTODemux, Seurat
#hashedDrops, https://rdrr.io/github/MarioniLab/DropletUtils/man/hashedDrops.html
#demuxmix, https://www.biorxiv.org/content/10.1101/2023.01.27.525961v1.full.pdf

# Methods for SNP-based classification
#1. demuxlet
#2. vireo
#3. Souporcell, https://www.nature.com/articles/s41592-020-0820-1


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






acc_mtx<- bc_demux2$acc_mtx
rownames(acc_mtx)[which(rownames(acc_mtx) %in% tag_mapping$tag)] <- paste0(tag_mapping$true_label, "-", tag_mapping$tag)[match(rownames(acc_mtx)[which(rownames(acc_mtx) %in% tag_mapping$tag)], tag_mapping$tag)]
#acc_mtx <- acc_mtx[!rownames(acc_mtx) %in% c("NA"), !colnames(acc_mtx) %in% c("unassigned")]
acc_mtx <- acc_mtx[order(rownames(acc_mtx)), order(colnames(acc_mtx))]
breaksList <- seq(0,1000,by=1e1)
graphics.off() # Shut down open plotting devices if any left over
pdf(paste0(plotPath, test_text, "accmtx", ".pdf"), width = 3, height = 2)
pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, border_color=NA, number_format ="%.0f", number_color ='white', fontsize = 8, fontsize_number=8, color = get_gradient_color('viridis', length(breaksList)), breaks = breaksList)
dev.off()



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
gmm_res_call <- gsub("[.]", "-", gmm_res_call) # Correction
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










