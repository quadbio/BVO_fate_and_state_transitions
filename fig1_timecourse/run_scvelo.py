import scvelo as scv
from scipy.io import mmwrite

adata = scv.read('../data/timecourse_velocity.loom', cache=True) # all cells

scv.pp.filter_genes(adata, min_shared_counts=10)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=3000)
scv.pp.log1p(adata)

scv.pp.moments(adata, n_neighbors=30, use_rep="css_embedding") # all cells
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.tl.rank_velocity_genes(adata, match_with='clusters', resolution=.4)

adata.write("../data/timecourse_scvelo.h5ad")

tm = scv.tl.transition_matrix(adata)
mmwrite("../data/timecourse_scvelo_transition_mat.mtx", tm)

