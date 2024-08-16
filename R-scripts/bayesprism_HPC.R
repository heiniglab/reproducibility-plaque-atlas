suppressWarnings(library(BayesPrism))
#suppressWarnings(library(TSIclient))
suppressWarnings(library(dplyr))
suppressWarnings(library(SummarizedExperiment))
library(SingleCellExperiment)
#suppressPackageStartupMessages(library(rtracklayer)) # Parse GTF
suppressPackageStartupMessages(library(ggplot2))
library(gplots) # heatmap.2
library("openxlsx")
library(zellkonverter)


setwd("/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/bulk")
output_dir = "cell_type_deconvolution/lars/sign/onlylevel2V3"; dir.create(file.path(output_dir), showWarnings = FALSE)


###########################################################
# Extract data of interest
###########################################################

print("Reading in scdata...")
# Rows: genes (ENSGID remove versioning with stringr::str_replace(colnames(bulk.count_matrix), "\\.[0-9]+", ""))
# Column: cells
#sc.count_matrix = readH5AD('/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Altas_level2-noprototypes-allgenes-names-mse3-reload-Bayes-uncorrected.h5ad')
#sc.count_matrix = readH5AD('/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Altas_level2-noprototypes-allgenes-names-mse3-reload-Bayes-uncorrected-V2.h5ad') # with new EC subtypes
#sc.count_matrix = readH5AD('/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Altas_level2-noprototypes-allgenes-names-mse3-reload-Bayes-uncorrected-V2-swap.h5ad') # with swapped EC subtypes
#sc.count_matrix = readH5AD('/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Altas_level2-noprototypes-allgenes-names-mse3-reload-Bayes-uncorrected-V2-swap-carotid.h5ad') # with swapped EC subtypes only carotid
#sc.count_matrix = readH5AD('/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Altas_level2-noprototypes-allgenes-names-mse3-reload-Bayes-uncorrected-V3.h5ad') # V2 but with new SMC subtypes
sc.count_matrix = readH5AD('/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Altas_level2-noprototypes-allgenes-names-mse3-reload-Bayes-uncorrected-V3-miller-ECswap.h5ad') # v2-swap but with new SMC subtypes from C. Miller and modSMC from jessica


original_colData <- colData(sc.count_matrix)
original_rowData <- rowData(sc.count_matrix)

counts <- assays(sc.count_matrix)$X # Extract the transposed count matrix
transposed_counts <- as.matrix(t(counts))

# Rows: genes (ENSGID, remove versioning with stringr::str_replace(colnames(bulk.count_matrix), "\\.[0-9]+", ""))
# Column: samples
print("Reading in bulk..")

#Lars
bulk.count_matrix <- readRDS("/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/bulk/Maegdefessel/gene_summarized_experiment.RDS")
rownames(bulk.count_matrix) <- stringr::str_replace(rownames(bulk.count_matrix), "\\.[0-9]+", "")
bulk.count_matrix <- as.matrix(t(assays(bulk.count_matrix)$counts))



# Vector of cell types, same length as columns in sc.count_matrix
cell_types.level1 <- as.character(original_colData$cell_type_level1)
cell_types.level2 <- as.character(original_colData$cell_type_level2)


print("Filter genes...")
# Filter genes
sc.count_matrix.filtered = cleanup.genes(input = transposed_counts,
                                         input.type = "count.matrix",
                                         species = "hs", 
                                         gene.group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                         exp.cells = 5)


# Filter for protein coding genes
print("Filter for protein coding genes")
sc.count_matrix.filtered.pc = select.gene.type(sc.count_matrix.filtered, gene.type = "protein_coding")


print("Filter for signature genes")
# Optionally, in cases where cell types are defined in a way that some of them
# show very similar transcription or severe batch effects exist, e.g., reference
# and mixture are from very different assays (ribo-depleted RNA-seq vs poly-A
# tail RNA-seq or PRO-seq (nascent RNA-seq) vs RNA-seq (steady state RNA)),
# selecting signature genes can be beneficial. This is because the selection of
# signature genes can enrich for genes informative for deconvolution while
# reducing the impact of noise caused by technical batch effects.
diff.exp.stat = get.exp.stat(sc.count_matrix.filtered.pc[,colSums(sc.count_matrix.filtered.pc > 0) > 3], # filter genes to reduce memory use
                             cell.type.labels = cell_types.level2,
                             cell.state.labels = cell_types.level2, # change this line to include cell states
                             pseudo.count = 0.1, # a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                             cell.count.cutoff = 50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                             n.cores = 1 # number of threads
)

sc.count_matrix.filtered.pc.sig = select.marker(sc.dat = sc.count_matrix.filtered.pc,
                                                stat = diff.exp.stat,
                                                pval.max = 0.01,
                                                lfc.min = 0.1)


###########################################################
# Prism
###########################################################
print("Construct Prism")
myPrism = new.prism(
  reference = sc.count_matrix.filtered.pc.sig, 
  mixture = bulk.count_matrix,
  input.type = "count.matrix", 
  cell.type.labels = cell_types.level2,
  cell.state.labels = cell_types.level2, # change this line to include cell states
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
  key = NULL
)

print("Run Prism")
bp.res = run.prism(prism = myPrism, n.cores=50)
saveRDS(bp.res, file.path(output_dir, "bayesprism_output.rds"))



###########################################################
# Extract info
###########################################################

print("Get ct fraction and CV")
# Extract deconvolved cell type fraction (rowSums = 1)
theta = get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")
saveRDS(theta, file.path(output_dir, "theta.rds"))
# theta = readRDS(file.path(output_dir, "theta.rds"))

# Extract coefficient of variation (CV) of cell type fraction
theta.cv = bp.res@posterior.theta_f@theta.cv
saveRDS(theta.cv, file.path(output_dir, "theta.cv.rds"))