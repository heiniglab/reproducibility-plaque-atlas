import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker
import jax

print("Jax GPU:", jax.devices() )

# ## Evaluation

adata = sc.read_h5ad("/home/icb/korbinian.traeuble/projects/Roche/hpc-data-transfer/Roche_main/data/Plaque-atlas/scarches-noroche-mse.h5ad")

print("Starting evaluations...")

embedding_obsm_keys_list = ['PCA', 'scVI', 'scANVI', "scGen",'scPoli','Harmony', 'LIGER']


bm = Benchmarker(
	adata,
	batch_key="sample",
	label_key="cell_type_level1",
	embedding_obsm_keys=embedding_obsm_keys_list,
	n_jobs=7,
	)
bm.benchmark()


print("Finished benchmarker..")
print("Writing data..")
bm.plot_results_table(save_dir="output_small_atlas/scarches/mse")
df = bm.get_results(min_max_scale=False)
df.to_csv("output_small_atlas/scarches/mse/results_scIB.csv")
df1 = df.transpose()
df1.to_csv("output_small_atlas/scarches/mse/results_transposed_sciB.csv")
print("Finished without errors")



