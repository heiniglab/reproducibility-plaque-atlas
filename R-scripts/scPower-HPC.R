library(scPower)
library(SingleCellExperiment)
library(Matrix)
library(reshape2) 
library(ggplot2)
library(zellkonverter)
library(scuttle)

sce <- readH5AD("/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/Altas_level2-noprototypes-allgenes-names-mse3-reload-scPower-uncorrected.h5ad")

cell_types = c("T cell", "Macrophage", "NK", "Monocyte", "SMC","B cell","EC", "Fibroblast", "Fibromyocyte", "Mast cell", "DC", "Plasma cell") #
  

all_gamma_linear_fits <- data.frame()
all_disp_fun_general <- data.frame()

# Function to subsample a matrix to a desired total count
subsample_matrix <- function(mat, target_total_count) {
  # Get the current total count of the matrix
  current_total_count <- sum(mat)
  
  # Calculate the probability of each entry being selected
  prob <- mat / current_total_count
  
  # Perform multinomial subsampling
  subsampled_counts <- rmultinom(1, target_total_count, prob)
  
  # Reshape the subsampled counts back into the original matrix dimensions
  subsampled_matrix <- matrix(subsampled_counts, nrow(mat), ncol(mat))
  rownames(subsampled_matrix) <- rownames(mat)
  colnames(subsampled_matrix) <- colnames(mat)
  
  return(subsampled_matrix)
}

sparseToMatrix <- function(sparseMatrix) {
  # Get indices and values from sparse matrix
  i <- sparseMatrix@i + 1
  j <- rep(seq_len(ncol(sparseMatrix)), diff(sparseMatrix@p))
  v <- sparseMatrix@x
  
  # Create a zero matrix of the correct dimensions
  denseMatrix <- matrix(0, nrow = sparseMatrix@Dim[1], ncol = sparseMatrix@Dim[2])
  
  # Populate the non-zero entries
  denseMatrix[cbind(i, j)] <- v
  
  # Assign row and column names
  row.names(denseMatrix) <- sparseMatrix@Dimnames[[1]]
  colnames(denseMatrix) <- sparseMatrix@Dimnames[[2]]
  
  return(denseMatrix)
}

subsampleIntoList <- function(counts) {
  tmp <- vector("list", 4)
  
  tmp[[1]] <- counts

  # Define the proportions
  proportions <- c(0.75, 0.5, 0.25)
  
  # Loop over the proportions and fill the list
  for(i in seq_along(proportions)){
    subsample <- downsampleMatrix(counts, prop = proportions[i], bycol = TRUE)
    subsample <- sparseToMatrix(subsample)
    
    tmp[[i + 1]] <- subsample
  }

  # Name the list elements
  names(tmp) <- c("Complete", "Subsampled_75", "Subsampled_50", "Subsampled_25")

  print("Subsampling process done successfully.")
  return(tmp)
}




for(celltype in cell_types){ 

counts_matrix <- assay(sce, "X")
counts_matrix <- as(counts_matrix, "matrix")  

print(paste0("Calculating gamma fits for ", celltype))
  
# only one cell_type.. FOR CLUSTER ONLY
cells_indices <- sce$cell_type_level1 == celltype
counts_matrix <- counts_matrix[, cells_indices]


#subsampling

if (celltype == "T cell"){
all_matrices <- subsampleIntoList(counts_matrix)
}
else {
# Determine the target total counts for 75%, 50%, and 25% subsampling
target_counts <- sum(counts_matrix) * c(0.75, 0.5, 0.25)

# Generate subsampled matrices
subsampled_75 <- subsample_matrix(counts_matrix, target_counts[1])
subsampled_50 <- subsample_matrix(counts_matrix, target_counts[2])
subsampled_25 <- subsample_matrix(counts_matrix, target_counts[3])

all_matrices <- list(Complete = counts_matrix, Subsampled_75 = subsampled_75, Subsampled_50 = subsampled_50, Subsampled_25 = subsampled_25)
}




# ----------------------------------------------------------------------------------------------------------------------------------------------------------

norm.mean.values<-NULL
disp.param<-NULL
for(name in names(all_matrices)){
  temp<-nbinom.estimation(all_matrices[[name]], sizeFactorMethod = "poscounts")
  #Save the normalized mean values
  norm.mean.values.temp<-temp[[1]] 
  norm.mean.values.temp$matrix<-name 
  norm.mean.values<-rbind(norm.mean.values,norm.mean.values.temp)
  #Save the parameter of the mean-dispersion function
  disp.param.temp<-temp[[3]] 
  disp.param.temp$matrix<-name 
  disp.param<-rbind(disp.param,disp.param.temp)
}
#First rows of the data frame with normalized mean values
#head(norm.mean.values)
#print("Dispersion Params")
#print(disp.param)


gamma.fits<-NULL
for(name in names(all_matrices)){
  #Number of cells per cell type as censoring point
  censoredPoint<- 1 / ncol(all_matrices[[name]])
  norm.mean.values.temp<-norm.mean.values[norm.mean.values$matrix==name,] 
  gamma.fit.temp<-mixed.gamma.estimation(norm.mean.values.temp$mean,num.genes.kept = 28000, censoredPoint = censoredPoint)
  gamma.fit.temp$matrix<-name
  gamma.fits<-rbind(gamma.fits,gamma.fit.temp) }

print(gamma.fits)

directory_path <- paste0("output/scPower/big-atlas/", celltype, "/")
if (!dir.exists(directory_path)) {
  dir.create(path = directory_path, recursive = TRUE)
}

png(filename = paste0("output/scPower/big-atlas/", celltype, "/gammafit.png"), width = 800, height = 600)

g<-visualize.gamma.fits(norm.mean.values$mean[norm.mean.values$matrix=="Complete"], gamma.fits[gamma.fits$matrix=="Complete",],
                        nGenes=28000)
print(g)

dev.off()



umi.values<-NULL
for(name in names(all_matrices)){
  mean.umi<-meanUMI.calculation(all_matrices[[name]])
  umi.values<-rbind(umi.values,data.frame(mean.umi,matrix=name)) 
  }

#print(umi.values)

gamma.fits<-merge(gamma.fits,umi.values,by="matrix")
#Convert the gamma fits from the shape-rate parametrization to the mean-sd parametrization
gamma.fits<-convert.gamma.parameters(gamma.fits)

#Fit relationship between gamma parameters and UMI values
gamma.linear.fit.new<-umi.gamma.relation(gamma.fits) 

disp.fun.general.new<-dispersion.function.estimation(disp.param) 

gamma.linear.fit.new$ct<-celltype 
disp.fun.general.new$ct<-celltype

# Append new data frames to the aggregated data frames
all_gamma_linear_fits <- rbind(all_gamma_linear_fits, gamma.linear.fit.new)
all_disp_fun_general <- rbind(all_disp_fun_general, disp.fun.general.new)


}

print("Done")
print("Gamma fits")
print(all_gamma_linear_fits)
print("Dispersions")
print(all_disp_fun_general)


saveRDS(all_gamma_linear_fits, file = "output/scPower/big-atlas/all_gamma_linear_fits.rds")
saveRDS(all_disp_fun_general, file = "output/scPower/big-atlas/all_disp_fun_general.rds")
