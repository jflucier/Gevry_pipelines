import os
import scanpy as sc
import matplotlib.pyplot as plt


work_dir = '/home/def-gevrynic/programs/scenicplus/pbmc_tutorial'

#scRNA-seq preprocessing using Scanpy
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor='white')

if not os.path.exists(os.path.join(work_dir, 'scRNA')):
    os.makedirs(os.path.join(work_dir, 'scRNA'))

adata = sc.read_10x_h5(os.path.join(work_dir, 'data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5'))
adata.var_names_make_unique()

#Basic quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
#Optionally, predict and filter out doublets using Scrublet.
sc.external.pp.scrublet(adata) #estimates doublets
adata = adata[adata.obs['predicted_doublet'] == False] #do the actual filtering
adata

#Filter based on mitochondrial counts and total counts.
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

mito_filter = 25
n_counts_filter = 4300
fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax = axs[0], show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax = axs[1], show = False)
#draw horizontal red lines indicating thresholds.
axs[0].hlines(y = mito_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
axs[1].hlines(y = n_counts_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
fig.tight_layout()
#plt.show()
plt.savefig(work_dir + '/fig1.png')

adata = adata[adata.obs.n_genes_by_counts < n_counts_filter, :]
adata = adata[adata.obs.pct_counts_mt < mito_filter, :]
adata

# Data normalization
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

#Cell type annotation
adata_ref = sc.datasets.pbmc3k_processed() #use the preprocessed data from the Scanpy tutorial as reference
var_names = adata_ref.var_names.intersection(adata.var_names) #use genes which are present in both assays
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]
sc.pp.pca(adata_ref) #calculate PCA embedding
sc.pp.neighbors(adata_ref) #calculate neighborhood graph
sc.tl.umap(adata_ref) #calculate umap embedding
sc.tl.ingest(adata, adata_ref, obs='louvain') #run label transfer
adata.obs.rename({'louvain': 'ingest_celltype_label'}, inplace = True, axis = 1)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig(work_dir + '/fig2.png')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)
sc.pl.umap(adata, color = 'ingest_celltype_label')
plt.savefig(work_dir + '/fig3.png')

sc.tl.leiden(adata, resolution = 0.8, key_added = 'leiden_res_0.8')
sc.pl.umap(adata, color = 'leiden_res_0.8')
plt.savefig(work_dir + '/fig4.png')

tmp_df = adata.obs.groupby(['leiden_res_0.8', 'ingest_celltype_label']).size().unstack(fill_value=0)
tmp_df = (tmp_df / tmp_df.sum(0)).fillna(0)
leiden_to_annotation = tmp_df.idxmax(1).to_dict()
leiden_to_annotation

leiden_to_annotation['7'] = 'B cells 1'
leiden_to_annotation['12'] = 'B cells 2'
leiden_to_annotation = {cluster: leiden_to_annotation[cluster].replace(' ', '_') for cluster in leiden_to_annotation.keys()}
leiden_to_annotation

adata.obs['celltype'] = [leiden_to_annotation[cluster_id] for cluster_id in adata.obs['leiden_res_0.8']]
#del(leiden_to_annotation)
#del(tmp_df)
sc.pl.umap(adata, color = 'celltype')
#fig, axs = plt.subplots(ncols = 2, figsize = (8,8))
plt.savefig(work_dir + '/fig5.png', bbox_inches="tight")

#Save results
adata.write(os.path.join(work_dir, 'scRNA/adata.h5ad'), compression='gzip')




