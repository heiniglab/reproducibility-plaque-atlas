library(scPower)
library(SingleCellExperiment)
library(reshape2) 
library(ggplot2)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)
library(zellkonverter) 

#environment location: /Users/korbinian.traeuble/Library/Caches/org.R-project.R/R/basilisk/1.14.3/0
#environment location: /Users/korbinian.traeuble/Library/Caches/org.R-project.R/R/basilisk/1.14.3/zellkonverter/1.12.1/zellkonverterAnnDataEnv-0.10.2


# ------------------------------------------------Grid calcs Power-----------------------------------------------------------------------------------------


gamma.linear.fit.new <- readRDS("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/big-atlas/scPower/corrected2/all_gamma_linear_fits.rds")
disp.fun.general.new <- readRDS("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/big-atlas/scPower/corrected2/all_disp_fun_general.rds")


# pilot data
read.umis <- data.frame(
  transcriptome.mapped.reads = c(42918.738 ,1134244.8 ,957056.115 ,753901.0 ,151165.035 ,96032.88900000001 ,144224.952 ,123504.68 ,348684.804 ,171858.988 ,396504.79 ,66962.89600000001 ,703458.368 ,446701.79199999996 ,139964.30800000002 ,144603.433 ,129465.252 ,61029.166),
  mean.umi = c(1703 ,1833 ,1893 ,1561 ,1789 ,1685 ,1460 ,1096 ,2287 ,1468 ,5552 ,2940 ,3553 ,5247 ,1282 ,1360 ,322 ,283)
)

read.umi.fit.new<-umi.read.relation(read.umis)



cell_types = c("T cell", "Macrophage", "EC", "Smooth Muscle Cell", "Fibromyocyte", "NK cell", "Fibroblast", "Monocyte", "Plasma cell", "Mast cell", "B cell", "Dendritic cell", "Neutrophil")
frequencies = c(0.245075, 0.151300, 0.162333, 0.163909, 0.059102, 0.054374, 0.053586, 0.041765, 0.019701, 0.015760, 0.010244, 0.014184, 0.008668)
ct_vector <- setNames(frequencies, cell_types)
#genesNr <- c(100, 1000, 5000)
FoldChanges <- c(1.1,1.5,2.5)

# Read in ranks with cell_type and pathway info
ranks_all = read.csv("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/big-atlas/scPower/corrected2/all_gene_ranks6.csv")

ranks_all <- ranks_all[!(ranks_all$pathway %in% c("Neovascularization and angiogenesis", "CCL21")), ]
#ranks_all$pathway[ranks_all$pathway == "Toll like receptor pathway"] <- "TLR signaling pathway"

# Identify the unique pathways in the ranks_all data
pathways <- unique(ranks_all$pathway)

# Prepare results dataframe
resultsDF <- data.frame(
  Celltype = character(0),
  `Mean Fold Change` = numeric(0),
  power = numeric(0),
  pathway = character(0)
)

for (ctName in names(ct_vector)) {
  # Frequency for current cell type
  ctFreq <- ct_vector[ctName]
  
  # Loop through each pathway scenario
  for (path in pathways) {
    # Filter ranks for current cell type and pathway
    gene_set <- ranks_all %>% 
      filter(cell_type == ctName, pathway == path)
    
    if (nrow(gene_set) == 0) {
      next  # If no genes for this combination, skip
    }
    
    # Extract the ranks from the filtered gene set
    ranks <- gene_set$rank
    
    # For each fold change scenario
    for (FCmean in FoldChanges) {
      # Simulate fold changes for all genes in the scenario
      numGenes <- length(ranks)
      foldChange <- effectSize.DE.simulation(mean = FCmean, sd = 0.5, numGenes = numGenes)
      
      simulated.de.genes <- data.frame(ranks = ranks, FoldChange = foldChange, name = "Simulated")
      
      power <- power.general.withDoublets(
        nSamples = 22,
        nCells = 200,
        readDepth = 334015,
        ct.freq = ctFreq,
        nGenes = 29000,
        type = "de",
        ref.study = simulated.de.genes,
        ref.study.name = "Simulated",
        samplesPerLane = 1,
        read.umi.fit = read.umi.fit.new,
        gamma.mixed.fits = gamma.linear.fit.new,
        ct = ctName,
        disp.fun.param = disp.fun.general.new,
        mappingEfficiency = 0.43,
        min.UMI.counts = 3,
        perc.indiv.expr = 0.5,
        sign.threshold = 0.05,
        MTmethod = "FDR"
      )
      
      resultpower <- power$power
      
      new_row <- data.frame(
        Celltype = ctName,
        `Mean Fold Change` = FCmean,
        power = resultpower,
        pathway = path
      )
      
      resultsDF <- rbind(resultsDF, new_row)
    }
  }
}

# Add Frequency column
resultsDF$Frequency <- ct_vector[as.character(resultsDF$Celltype)]

# Aggregate by cell type to get Frequency for plotting
aggDF <- resultsDF %>%
  group_by(Celltype) %>%
  summarise(Frequency = dplyr::first(Frequency))

rotate_x <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
theme_set(theme_bw())

ct_order <- order(ct_vector)
resultsDF$Celltype <- factor(resultsDF$Celltype, levels = names(ct_vector)[rev(ct_order)])

mat <- ggplot(data = resultsDF, aes(x = Celltype, y = power)) + 
  geom_point(aes(col = factor(`Mean.Fold.Change`))) + 
  geom_col(data = aggDF, aes(x = Celltype, y = Frequency), fill = "gray", 
           show.legend = FALSE, alpha = 0.5) +
  facet_grid(. ~ pathway) +
  scale_y_continuous(limits = c(0, 1),sec.axis = sec_axis(~., name = "Frequency")) +
  rotate_x +
  labs(col = "Fold Change")

print(mat)

ggsave(
  filename = "/Users/korbinian.traeuble/PhD-local/projects/main_Roche/manuscript/revision1/corrected2/matthiasplots_RemodelAndIFN-yAndCCL19.pdf",
  plot = mat,
  width = 10,
  height = 6,
  device = pdf
)

write.csv(resultsDF, "resultsDF_newIFN.csv", row.names = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------
# Power to detect rare cell types
power.detect.celltype(nCells=500,min.num.cells = 10, cell.type.frac=0.05,nSamples=20)
number.cells.detect.celltype(prob.cut=0.95,min.num.cells = 10, cell.type.frac=0.05,nSamples=20)