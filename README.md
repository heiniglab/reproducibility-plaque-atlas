## Reproducibility code for "Integrated single-cell atlas of human atherosclerotic plaques"

This is the repository that includes the code used for the scRNA-seq human atherosclerotic plaque atlas (https://www.nature.com/articles/s41467-025-63202-x). 
Because not all cells in the revisions were executed in these notebooks, some output figures need to be recomputed.

## Usefull links

Paper: https://www.nature.com/articles/s41467-025-63202-x \
\
The resulting atlas can also be downloaded on CELLxGENE: https://cellxgene.cziscience.com/collections/db70986c-7d91-49fe-a399-a4730be394ac. \
\
Repo for automatic cell type annotation: https://github.com/heiniglab/plaque-atlas-mapping

## Changelog
### 2026-01-29
- **Location of counts in CELLxGENE:** Due to several email requests, we want to clarify that the gene counts of the CELLxGENE object are at adata.raw.X.

### 2025-10-27
- **Duplicate cells identified:** During a recent CELLxGENE data audit, duplicated cells were found in the atlas caused by overlapping barcodes between certain samples.  
  - `donor_id:9_Chowdhury` (230 cells) is fully contained within `donor_id:2_Wirka`.  
  - `donor_id:8_Chowdhury` (295 cells) is fully contained within the combination of `donor_id:3_Wirka` and `donor_id:4_Wirka`.
  - Of note, the Chowdhury dataset was sorted for T-cells.  
- **Action:** The CELLxGENE submission was updated to remove `donor_id:9_Chowdhury` and `donor_id:8_Chowdhury`. 

## How to run the code

We provide environment yml files for the python and R packages, required to run the notebooks/scripts. 

The notebooks can be used in the given order to reproduce the atlas from scratch and the calculations for the plots and results in the preprint. 

To execute the first notebook, 1_preprocessing.ipynb, the corresponding raw datasets need to be downloaded from GEO. The corresponding accession numbers are: 

| Dataset           | Plaque Location | # Cells | # Samples | Accession Number |
|-------------------|-----------------|---------|-----------|------------------|
| Pan et al.      | carotid         | 8,850   | 3         | GSE155512        |
| Alsaigh et al.  | carotid         | 34,626  | 3         | GSE159677        |
| Jaiswal et al.  | carotid         | 5,554   | 2         | GSE179159        |
| Dib et al.       | carotid         | 25,147  | 6         | GSE210152        |
| Fernandez et al. | carotid         | 9,935   | 7         | GSE224273        |
| Slysz et al.     | carotid         | 27,245  | 3         | GSE234077        |
| Pauli et al.    | carotid         | 1,269   | 8         | GSE247238        |
| Bashore et al.  | carotid         | 75,582  | 18        | GSE253904        |
| Wirka et al.     | coronary        | 11,750  | 8         | GSE131778        |
| Emoto et al.     | coronary        | 1,621   | 2         | GSE184073        |
| Chowdhury et al. | coronary        | 22,543  | 12        | GSE196943        |
| Slysz et al.    | femoral         | 35,371  | 7         | GSE234077        |
| **Total**         |                 | **259,493** | **79** |                  |

Some of the scripts have a "HPC" suffix, which means that this script requires substantial amount of memory and is hence recommended to be ran on a cluster.


## Citation

If you use this resource, please cite our paper:

Traeuble, K., Munz, M., Pauli, J. et al. Integrated single-cell atlas of human atherosclerotic plaques. Nat Commun 16, 8255 (2025). https://doi.org/10.1038/s41467-025-63202-x

### BibTeX
```bibtex
@article{Traeuble2025,
  title = {Integrated single-cell atlas of human atherosclerotic plaques},
  volume = {16},
  ISSN = {2041-1723},
  url = {http://dx.doi.org/10.1038/s41467-025-63202-x},
  DOI = {10.1038/s41467-025-63202-x},
  number = {1},
  journal = {Nature Communications},
  publisher = {Springer Science and Business Media LLC},
  author = {Traeuble,  Korbinian and Munz,  Matthias and Pauli,  Jessica and Sachs,  Nadja and Vafadarnejad,  Eshan and Carrillo-Roa,  Tania and Maegdefessel,  Lars and Kastner,  Peter and Heinig,  Matthias},
  year = {2025},
  month = sep 
}