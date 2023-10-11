

# Files can be obtained from the docker files at https://hub.docker.com/r/regevlab/demuxem.

library(deMULTIplex2)

savePath <- "./data-raw/benchmark_res/Gaublomme/human_st/"; dir.create(savePath)
rdsPath <- paste0(savePath, "rds/");dir.create(rdsPath)
plotPath <- paste0(savePath, "plots/"); dir.create(plotPath)
dataPath <- paste0(savePath, "data/"); dir.create(plotPath)

test_text <- "benchmark_Gaublomme_human_st_"


# Process barcode data
adt_count = data.table::fread('./data-raw/benchmark_res/Gaublomme/docker_files/software/inputs/experiment1_human_st_ADT.csv')
antibody = adt_count[,Antibody]
# Convert to sparse matrix
tag_mtx<- t(Matrix(as.matrix(adt_count[,Antibody:=NULL]), sparse=T))
colnames(tag_mtx) <- antibody
rownames(tag_mtx) <- paste0(rownames(tag_mtx), "-1")

saveRDS(tag_mtx, paste0(rdsPath, "tag_mtx.rds"))

# Process count data
rna_mtx <-Seurat::Read10X('./data-raw/benchmark_res/Gaublomme/docker_files/software/inputs/experiment1_human_st_raw_GRCh38_premrna/')
table(sapply(strsplit(colnames(rna_mtx), "-"), function(x)x[2]))
# Simple filter on cells
total_umi_rna <- colSums(rna_mtx)
min.rna.umi <- 100
rna_mtx <- rna_mtx[,total_umi_rna >= min.rna.umi]
saveRDS(rna_mtx, paste0(rdsPath, "rna_mtx.rds"))

cells <- intersect(rownames(tag_mtx), colnames(rna_mtx))
tag_mtx <- tag_mtx[cells,]

min.tt.umi = 0
# Filter out cells which have too few tag counts which could be beads that can mess up with fitting
tag_mtx <- tag_mtx[rowSums(tag_mtx) > min.tt.umi, ]
print(dim(tag_mtx))
rna_mtx <- rna_mtx[, cells]
print(dim(rna_mtx))


init.cos.cut = .5
converge.threshold = 1e-3
prob.cut = 0.5
max.cell.fit = 3000
max.iter = 30

res <- demutiplexTags(tag_mtx,
                      init.cos.cut = init.cos.cut,
                      converge.threshold = converge.threshold,
                      max.iter = max.iter,
                      prob.cut = prob.cut,
                      min.cell.fit = 10,
                      max.cell.fit = max.cell.fit,
                      residual.type = c("rqr", 'pearson'), # ONLY use RQR for future
                      plot.umap = c("residual", "umi"),
                      plot.diagnostics = T,
                      plot.path = plotPath,
                      plot.name = paste0(test_text),
                      umap.nn = 30L,
                      seed = 1,
                      point.size = 1,
                      label.size = 3,
                      min.bc.show = 50)
saveRDS(res, paste0(rdsPath, test_text, "cosc_", init.cos.cut, "probc_",  prob.cut, "maxc_",max.cell.fit, "maxi_",max.iter, "_res_list.rds"))


# Benchmarking with Donor data, following and reproducing the python code
demuxlet_res <- fread('./data-raw/benchmark_res/Gaublomme/docker_files/software/inputs/experiment1_demuxlet.best')

demuxlet_res[['demux_type']] = 'doublet'
idx_amb <- grepl("AMB", demuxlet_res$BEST)
demuxlet_res[['demux_type']][idx_amb] <-"unknown"
idx_sng = (demuxlet_res[['PRB.DBL']] <= 0.99) & (!idx_amb)
demuxlet_res[['demux_type']][idx_sng] = 'singlet'
demuxlet_res[['assignment']] = ''
demuxlet_res[['assignment']][idx_sng] = demuxlet_res[['SNG.1ST']][idx_sng]
idx_dbt = demuxlet_res[['demux_type']]== 'doublet'
demuxlet_res[['assignment']][idx_dbt] = paste0(demuxlet_res[['SNG.1ST']][idx_dbt], ",", demuxlet_res[['SNG.2ND']][idx_dbt])

demuxlet_res$assignment_simplified <- demuxlet_res$assignment
demuxlet_res$assignment_simplified[demuxlet_res$demux_type == "doublet"] <- "doublet"
demuxlet_res$assignment_simplified[demuxlet_res$demux_type == "unknown"] <- "unknown"
use_levels = sort(unique(demuxlet_res$assignment_simplified))
use_levels <- c(use_levels[!use_levels %in% c("doublet", "unknown")], c("doublet", "unknown"))
demuxlet_res$assignment_simplified <- factor(demuxlet_res$assignment_simplified, levels = use_levels)


saveRDS(demuxlet_res, paste0(rdsPath, "demuxlet_res.rds"))


shared_cells <- intersect(demuxlet_res[,BARCODE], rownames(res$assign_table))

subset_col = c('demux_type', 'assignment')
combine_res <- cbind(res$umap[shared_cells, ], res$assign_table[shared_cells, ], demuxlet_res[demuxlet_res$BARCODE %in% shared_cells, ..subset_col])

combine_res$assignment_singlet <- combine_res$assignment
combine_res$assignment_singlet[combine_res$demux_type == "doublet"] <- "doublet"
combine_res$assignment_singlet[combine_res$demux_type == "unknown"] <- "unknown"
use_levels = sort(unique(combine_res$assignment_singlet))
use_levels <- c(use_levels[!use_levels %in% c("doublet", "unknown")], c("doublet", "unknown"))
combine_res$assignment_singlet <- factor(combine_res$assignment_singlet, levels = use_levels)

combine_res$barcode_assign[is.na(combine_res$barcode_assign)] <- "assign_NA"
acc_mtx <- as.data.frame.matrix(table(combine_res$barcode_assign, combine_res$assignment_singlet))

breaksList <- seq(0,1e3,by=1e2)
pdf(paste0(plotPath, test_text, "_demux2_assign_donor_bc_accmtx", ".pdf"), width = 3.5, height = 2.5)
pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, number_format ="%.0f", fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList)
dev.off()


use_row <- rownames(acc_mtx)[!rownames(acc_mtx) %in% c("assign_NA")]
use_col <- colnames(acc_mtx)[!colnames(acc_mtx) %in% c("doublet", "unknown")]
sub_mtx <- as.matrix(acc_mtx[use_row,use_col])
pc <- sum(diag(sub_mtx )) / sum(acc_mtx[use_row,use_col])
rc = sum(diag(sub_mtx )) / sum(acc_mtx[use_row,])
st <- sum(diag(sub_mtx )) / sum(acc_mtx)


