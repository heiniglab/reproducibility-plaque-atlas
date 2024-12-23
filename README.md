## Reproducibility code for "Integrated single-cell atlas of human atherosclerotic plaques"

This is the repository that includes the code used for the scRNA-seq human atherosclerotic plaque atlas (preprint: https://www.biorxiv.org/content/10.1101/2024.09.11.612431v1)

## How to run the code

We provide environment yml files for the python and R packages, required to run the notebooks/scripts. 

The notebooks can be used in the given order to reproduce the atlas from scratch and the calculations for the plots and results in the preprint. 

To execute the first notebook, 1_preprocessing.ipynb, the corresponding raw datasets need to be downloaded from GEO. The corresponding accession numbers are: 

| Dataset           | Plaque Location | # Cells | # Samples | Accession Number |
|-------------------|-----------------|---------|-----------|------------------|
| Pan et al.      | carotid         | 8,867   | 3         | GSE155512        |
| Alsaigh et al.  | carotid         | 34,716  | 3         | GSE159677        |
| Jaiswal et al.  | carotid         | 5,567   | 2         | GSE179159        |
| Dib et al.       | carotid         | 25,194  | 6         | GSE210152        |
| Fernandez et al. | carotid         | 9,984   | 7         | GSE224273        |
| Slysz et al.     | carotid         | 27,540  | 4         | GSE234077        |
| Pauli et al.    | carotid         | 1,278   | 8         | GSE247238        |
| Bashore et al.  | carotid         | 77,112  | 18        | GSE253904        |
| Wirka et al.     | coronary        | 11,756  | 8         | GSE131778        |
| Emoto et al.     | coronary        | 1,624   | 2         | GSE184073        |
| Chowdhury et al. | coronary        | 22,671  | 12        | GSE196943        |
| Slysz et al.    | femoral         | 35,438  | 9         | GSE234077        |
| **Total**         |                 | **261,747** | **80** |                  |

Some of the scripts have a "HPC" suffix, which means that this script requires substantial amount of memory and is hence recommended to be ran on a cluster.

## Usefull links

preprint: https://www.biorxiv.org/content/10.1101/2024.09.11.612431v1 \
\
The resulting atlas can also be downloaded on CELLxGENE: https://cellxgene.cziscience.com/collections/db70986c-7d91-49fe-a399-a4730be394ac. (is currently being updated to the new atlas version) \
\
Repo for automatic cell type annotation: https://github.com/kotr98/plaque-atlas-mapping

