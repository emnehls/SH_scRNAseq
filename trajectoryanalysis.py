## A script to maybe run the trajectory analysis on naive.integrated.

import scvelo as scv

adata = scv.read("wt_vel.h5ad")
print(adata)


scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
#scv.tl.velocity(adata)
#scv.tl.velocity_graph(adata)
#scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters")

