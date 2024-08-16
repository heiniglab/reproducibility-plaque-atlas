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


gamma.linear.fit.new <- readRDS("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/big-atlas/scPower/all_gamma_linear_fits.rds")
disp.fun.general.new <- readRDS("/Users/korbinian.traeuble/PhD-local/projects/main_Roche/data/big-atlas/scPower/all_disp_fun_general.rds")


# pilot data
read.umis <- data.frame(
  transcriptome.mapped.reads = c(42918.738 ,1134244.8 ,957056.115 ,753901.0 ,151165.035 ,96032.88900000001 ,144224.952 ,123504.68 ,348684.804 ,171858.988 ,396504.79 ,66962.89600000001 ,703458.368 ,446701.79199999996 ,139964.30800000002 ,144603.433 ,129465.252 ,61029.166),
  mean.umi = c(1703 ,1833 ,1893 ,1561 ,1789 ,1685 ,1460 ,1096 ,2287 ,1468 ,5552 ,2940 ,3553 ,5247 ,1282 ,1360 ,322 ,283)
)

read.umi.fit.new<-umi.read.relation(read.umis)


cell_types = c("T cell", "Macrophage", "EC", "SMC", "Fibromyocyte", "NK", "Fibroblast", "Monocyte", "Plasma cell", "Mast cell", "B cell", "DC")
frequencies = c(0.253521,0.169797,0.161189,0.121283,0.113459,0.047731,0.039906,0.035994,0.019562,0.017997,0.010172,0.009390) 
ct_vector <- setNames(frequencies, cell_types)
genesNr <- c(10, 20, 50)
FoldChanges <- c(1.1,1.5,2)

# Define an empty dataframe with specified column names
resultsDF <- data.frame(
  Celltype = character(0),
  `Number of Genes` = integer(0),
  `Mean Fold Change` = integer(0),
  `power` = numeric(0)
)


for (ctName in names(ct_vector)){
  print(ctName)
  # ctName already iter
  ctFreq <- ct_vector[ctName]
  for (numGenes in genesNr){
    for (FCmean in FoldChanges){
      #Uniform distributed gene ranks for 100 genes in the interval from 1-2,000
      ranks<-uniform.ranks.interval(start=1,end=numGenes,numGenes=20) 
      #Simulation of fold changes 
      foldChange<-effectSize.DE.simulation(mean=FCmean,sd=0.5,numGenes=20)
      simulated.de.genes<-data.frame(ranks=ranks, FoldChange=foldChange, name="Simulated")
      
      power<-power.general.withDoublets(nSamples=22,
                                        nCells=200,
                                        readDepth=334015, 
                                        ct.freq=ctFreq,
                                        nGenes = 29000,
                                        type="de",
                                        ref.study=simulated.de.genes, 
                                        ref.study.name="Simulated", 
                                        samplesPerLane=1,
                                        read.umi.fit = read.umi.fit.new,
                                        gamma.mixed.fits = gamma.linear.fit.new, 
                                        ct=ctName, 
                                        disp.fun.param=disp.fun.general.new, 
                                        mappingEfficiency = 0.43,
                                        min.UMI.counts = 3,
                                        perc.indiv.expr = 0.5,
                                        sign.threshold = 0.05,
                                        MTmethod="FDR")
      resultpower <- power$power
      
      #results <- c(ctName, numGenes, FCmean, resultpower)
      new_row <- data.frame(ctName, numGenes, FCmean, resultpower)
      colnames(new_row) <- names(resultsDF)
      resultsDF <- rbind(resultsDF, new_row)
    }
  }
}


# Ensure that each cell type in resultsDF has the correct frequency value from ct_vector
resultsDF$Frequency <- ct_vector[as.character(resultsDF$Celltype)]



aggDF <- resultsDF %>%
  group_by(Celltype) %>%
  summarise(Frequency = dplyr::first(Frequency))




rotate_x <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
theme_set(theme_bw())

ct_order <- order(ct_vector)
resultsDF$Celltype <- factor(resultsDF$Celltype, levels=names(ct_vector)[rev(ct_order)])


mat <- ggplot(data=resultsDF, aes(x=Celltype, y=power)) + 
  geom_point(aes(col=factor(Mean.Fold.Change))) + 
  geom_col(data = aggDF, aes(x = Celltype, y = Frequency), fill = "gray", show.legend = FALSE, alpha = 0.5) +
  facet_grid(.~Number.of.Genes, labeller = labeller(Number.of.Genes = function(x) paste("Expression rank <", x))) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Frequency")) +
  rotate_x + 
  labs(col="Fold Change")
print(mat)

ggsave(filename = "/Users/korbinian.traeuble/PhD-local/projects/main_Roche/figures/scPower/big-atlas/matthiasplots_generanks.pdf", 
       plot = mat, 
       width = 7, 
       height = 4, 
       device = pdf)

write.csv(resultsDF, "resultsDF_new.csv", row.names = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------
# Power to detect rare cell types
power.detect.celltype(nCells=500,min.num.cells = 10, cell.type.frac=0.05,nSamples=20)
number.cells.detect.celltype(prob.cut=0.95,min.num.cells = 10, cell.type.frac=0.05,nSamples=20)