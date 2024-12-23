library(scPower)
library(SingleCellExperiment)
library(Matrix)
library(reshape2) 
library(ggplot2)
library(zellkonverter)
library(scuttle)

sce <- readH5AD("/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Big-atlas-uncorrected_counts-nootherlayers-corrected2-ensembleNO-noCITE-withlevel2.h5ad")


cell_types = c("T cell", "Macrophage", "NK cell", "Monocyte", "Smooth Muscle Cell","B cell","EC", "Fibroblast", "Fibromyocyte", "Mast cell", "Dendritic cell", "Plasma cell", "Neutrophil") #

gene_lists <- list(
"Neovascularization and angiogenesis" = c("VEGFA", "FGD6", "TGFB1", "ANKS1A", 
                                            "CCM2", "BCAS3", "ZFPM2", "ARHGEF12", 
                                            "DAB2IP", "ARHGEF26", "TSPAN14", 
                                            "SMAD7", "SMAD3", "TCF21"),
"Vascular remodelling" = c("PECAM1", "HTRA1", "COL6A3", "P116", "FGD6", 
                              "MIA3", "SH3PXD2A", "TNS1", "ZC3HC1", "SERPINH1", 
                              "SEMA5A", "9p21", "COL4A2", "FURIN", "CALCRL", 
                              "ITGB5", "BMP1", "REST-NOA1", "MMP9", "RPL17", 
                              "SVEP1", "MYH11", "NEDD4", "ADORAZA", "ITGA1", 
                              "WWP2", "ADAMTS7", "FN1", "PDGFD", "BCAS3", 
                              "LOX", "SWAP70", "FLT1", "COL4A1", "COL4A2", 
                              "KSR2"),


"INF-y" = c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2",
  "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5",
  "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH",
  "CIITA", "CMKLR1", "CMPK2", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60",
  "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "CMTR1",
  "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA",
  "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44",
  "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA",
  "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4",
  "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1",
  "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "TMT1B", "MT2A",
  "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5",
  "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14",
  "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "HELZ2",
  "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2",
  "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31",
  "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28",
  "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4",
  "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3",
  "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26",
  "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1",
  "XCL1", "ZBP1", "ZNFX1"),


"CCL19" = c("IL23A","LTBR","PEG10","CCND1","CDKN1A","MYC","TNF","FOS","VEGFA","IRF2BPL","ID1","SMAD3","PIM1","NFKBIA","TP53","UBC","DDIT4",
            "BCL3","JUNB","GADD45B","CDKN1B","IL6","RARA","PLEC","IER3","SERPINE1","ICAM1","GATA4","PTMA","NR1D1","HIST1H4H","JUN","TRAF4",
            "SMAD7","BCL2L1","HES1","ZFP36","BRCA1","EPHA2","BHLHE40","ID3","BCL2","DUSP1","ATF3","BCL6","HSP90AB1","EGR1","IL1B","PTGS2","LMNA"),

"CCL21" = c("MAPK1","CDK1","CCNB1","CCNA2","CCND1","CDKN1A","VEGFA","MYC","IRF2BPL","FOS","ID1","SMAD3","TNF","UBC","PIM1","TP53","BCL3",
            "PLEC","SMAD7","HES1","SERPINE1","GADD45B","EPHA2","DDIT4","BCL2L1","TRAF4","BCL6","RARA","JUNB","HSP90AB1","JUN","PTMA","IER3",
            "GATA4","NFKBIA","CDKN1B","ID3","BRCA1","NR1D1","ATF3","HIST1H4H","BHLHE40","BCL2","CEBPB","IL6","EDN1","FN1","TGIF1","ICAM1","ID2")
)

# Initialize an empty list to store results
all_results <- list()

for (celltype in cell_types) {
  counts_matrix <- assay(sce, "X")
  counts_matrix <- as(counts_matrix, "matrix")
  
  print(paste0("Processing cell type: ", celltype))
  
  # Subset to specific cell type
  cells_indices <- sce$cell_type_level1 == celltype
  example.matrix <- counts_matrix[, cells_indices]
  
  # Remove genes with all zeros
  example.matrix <- example.matrix[rowSums(example.matrix) > 0, ]
  
  # Normalize by count per cell
  example.matrix <- t(t(example.matrix) / colSums(example.matrix))
  
  for (pathway in names(gene_lists)) {
    sign.genes <- gene_lists[[pathway]]
    
    print(paste0("Calculating gene ranks for pathway: ", pathway))
    
    # Calculate gene ranks
    gene.ranks <- gene.rank.calculation(example.matrix, sign.genes)
    
    # Convert gene ranks to a data frame and add cell type and pathway columns
    gene.ranks <- as.data.frame(gene.ranks)
    gene.ranks$gene_symbol <- rownames(gene.ranks)
    gene.ranks$cell_type <- celltype
    gene.ranks$pathway <- pathway
    
    # Append to results list
    all_results <- append(all_results, list(gene.ranks))
  }
}

# Combine all results into a single data frame
combined_results <- do.call(rbind, all_results)

# Save combined results to a single CSV file
write.csv(combined_results, "output/scPower/corrected2/all_gene_ranks.csv", row.names = FALSE)

print("Gene rank calculations completed and saved.")