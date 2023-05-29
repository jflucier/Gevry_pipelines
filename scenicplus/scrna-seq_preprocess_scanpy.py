import os
import argparse
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor='white')

def run_scrna_precprocess(work_dir, in_data, mtx_dir, mito_filter, n_counts_filter, min_genes, min_cells, min_mean,
                          max_mean,
                          min_disp, max_value, n_neighbors, leiden_res, cpu, mem):

    # scRNA-seq preprocessing using Scanpy
    print(f"initialisation of variable and output folder {work_dir}/scRNA")
    if not os.path.exists(os.path.join(work_dir, 'scRNA')):
        os.makedirs(os.path.join(work_dir, 'scRNA'))

    print(f"Setting execution to {cpu} threads and {mem}G")
    sc.settings.n_jobs = cpu
    sc.settings.max_memory = mem
    outdir = f"{work_dir}/scRNA"

    print(f"reading input 10x h5 file")
    adata = sc.read_10x_h5(in_data)
    # adata.var_names_make_unique()

    # Basic quality control
    print(f"filtering input by min_genes {min_genes} & min_cells {min_cells}")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    # Optionally, predict and filter out doublets using Scrublet.
    # estimates doublets
    print(f"Predicting and filtering out doublets using Scrublet.")
    sc.external.pp.scrublet(adata)
    # do the actual filtering
    adata = adata[adata.obs['predicted_doublet'] == False]
    # adata

    # Filter based on mitochondrial counts and total counts.
    # annotate the group of mitochondrial genes as 'mt'
    print(f"Annotating the group of mitochondrial genes as 'mt'")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    print(f"Generating figure qc_counts_filter_mt{mito_filter}_cnt{n_counts_filter}.png")
    fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axs[0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axs[1], show=False)
    # draw horizontal red lines indicating thresholds.
    axs[0].hlines(y=mito_filter, xmin=0, xmax=max(adata.obs['total_counts']), color='red', ls='dashed')
    axs[1].hlines(y=n_counts_filter, xmin=0, xmax=max(adata.obs['total_counts']), color='red', ls='dashed')
    fig.tight_layout()
    plt.savefig(f"{outdir}/qc_counts_filter_mt{mito_filter}_cnt{n_counts_filter}.png")

    print(f"Filtering all count data by total counts {n_counts_filter} and mitochondrial counts {mito_filter}")
    adata = adata[adata.obs.n_genes_by_counts < n_counts_filter, :]
    adata = adata[adata.obs.pct_counts_mt < mito_filter, :]

    # OPTIONNAL??? Data normalization for visualisation purpose
    print(f"Backup filtered data.")
    adata.raw = adata
    print(f"Performing data normalisation")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=max_value)

    # Cell type annotation
    # use the preprocessed data from the Scanpy tutorial as reference
    # adata_ref = sc.datasets.pbmc3k_processed()
    print(f"Perform cell type annotation using mtx dir {mtx_dir}")
    adata_ref = sc.read_10x_mtx(
        # the directory with the `.mtx` file. Tuto data taken from
        # https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_raw_gene_bc_matrices.tar.gz
        mtx_dir,
        # use gene symbols for the variable names (variables-axis index)
        var_names='gene_symbols',
        # write a cache file for faster subsequent reading
        cache=True,
    )
    adata_ref.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
    # adata.write('write/pbmc3k_raw.h5ad', compression='gzip')

    # use genes which are present in both assays
    print(f"Crossing count data with annotation data using gene symbols")
    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names]
    adata = adata[:, var_names]
    # calculate PCA embedding
    print(f"Calculate PCA")
    sc.pp.pca(adata_ref)
    # calculate neighborhood graph
    print(f"Calculate neighborhood graph")
    sc.pp.neighbors(adata_ref)
    # calculate umap embedding
    print(f"Calculate umap embedding")
    sc.tl.umap(adata_ref)
    # run label transfer
    print(f"Run label transfer")
    sc.tl.ingest(adata, adata_ref, obs='louvain')
    adata.obs.rename({'louvain': 'ingest_celltype_label'}, inplace=True, axis=1)

    # Let’s visualize the labels
    print(f"Generate cell type annotation plot cell_type_annotation_variance_ratio.png")
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    plt.savefig(f"{outdir}/cell_type_annotation_variance_ratio.png")

    print(f"Generate cell subgroups plot by label cell_subgroups.label.png")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=10)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='ingest_celltype_label')
    plt.savefig(f"{outdir}/cell_subgroups.label.png")

    # Cluster cells into subgroups using the Leiden algorithm
    print(f"Clustering cells into subgroups using the Leiden algorithm")
    sc.tl.leiden(adata, resolution=leiden_res, key_added=f'leiden_res_{leiden_res}')
    sc.pl.umap(adata, color=f'leiden_res_{leiden_res}')
    plt.savefig(f'{outdir}/cell_subgroups.leiden{leiden_res}.png')

    print(f"TEST: reannotate some groups")
    tmp_df = adata.obs.groupby([f'leiden_res_{leiden_res}', 'ingest_celltype_label']).size().unstack(fill_value=0)
    tmp_df = (tmp_df / tmp_df.sum(0)).fillna(0)
    leiden_to_annotation = tmp_df.idxmax(1).to_dict()
    # leiden_to_annotation

    leiden_to_annotation['7'] = 'B cells 1'
    leiden_to_annotation['12'] = 'B cells 2'
    leiden_to_annotation = {cluster: leiden_to_annotation[cluster].replace(' ', '_') for cluster in
                            leiden_to_annotation.keys()}
    # leiden_to_annotation
    print(f"Generate cell subgroups plot by cell type cell_subgroups.celltype.png")
    adata.obs['celltype'] = [leiden_to_annotation[cluster_id] for cluster_id in adata.obs['leiden_res_0.8']]
    sc.pl.umap(adata, color='celltype')
    plt.savefig(f"{outdir}/cell_subgroups.celltype.png", bbox_inches="tight")

    # Save results
    print(f"Saving data object in {outdir}/adata.h5ad")
    adata.write(os.path.join(outdir, 'adata.h5ad'), compression='gzip')


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-w", "--workdir", help="your working directory", required=True)
    argParser.add_argument("-i", "--input", help="your h5 input file", required=True)
    argParser.add_argument("-mtx", "--mtx_dir", help="directory with the mtx file from 10x genomic", required=True)
    argParser.add_argument("--cpu", nargs='?',
                           help="Number of cpu to use", const=24,
                           type=int, default=24)
    argParser.add_argument("--mem", nargs='?',
                           help="Max memery usage in Gigabyte", const=30,
                           type=int, default=30)

    # qc
    argParser.add_argument("--qc_min_genes", nargs='?',
                           help="Only keep cells with at least <<min_genes>> genes expressed", const=200,
                           type=int, default=200)
    argParser.add_argument("--qc_min_cells", nargs='?',
                           help="Only keep genes which are expressed in at least <<min_cells>> cells", const=3,
                           type=int, default=3)

    # filter
    argParser.add_argument("--filter_mito", nargs='?', help="Filter based on mitochondrial counts", const=25, type=int,
                           default=25)
    argParser.add_argument("--filter_n_counts", nargs='?', help="Filter based on total counts", const=4300, type=int,
                           default=4300)

    # normalisation
    argParser.add_argument("--norm_min_mean", nargs='?', help="Normalisation min mean value", const=0.0125, type=float,
                           default=0.0125)
    argParser.add_argument("--norm_max_mean", nargs='?', help="Normalisation max mean value", const=3, type=float,
                           default=3)
    argParser.add_argument("--norm_min_disp", nargs='?', help="Normalisation min dispersion value", const=0.5,
                           type=float,
                           default=0.5)
    argParser.add_argument("--norm_max_value", nargs='?', help="Clip (truncate) to this value after scaling", const=10,
                           type=float,
                           default=10)

    # Cell type annot
    argParser.add_argument("--ct_annot_n_neighbors", nargs='?',
                           help="The size of local neighborhood (in terms of number of "
                                "neighboring data points) used for manifold approximation. "
                                "Larger values result in more global views of the manifold, "
                                "while smaller values result in more local data being "
                                "preserved. ",
                           const=10, type=int, default=10)
    argParser.add_argument("--ct_annot_leiden_res", nargs='?', help="parameter value controlling the coarseness of the "
                                                                    "clustering. Higher values lead to more clusters.",
                           const=0.8, type=float, default=0.8)

    args = argParser.parse_args()
    # print("args=%s\n" % args)
    # print(vars(args))
    for k, v in vars(args).items():
        print(f"Input {k}: {v}")

    # for p in args:
    #     print(f"{p}")
    # print("Input workdir:%s" % args.workdir)
    # print("Input data:%s" % args.input)

    run_scrna_precprocess(
        args.workdir,
        args.input,
        args.mtx_dir,
        args.filter_mito,
        args.filter_n_counts,
        args.qc_min_genes,
        args.qc_min_cells,
        args.norm_min_mean,
        args.norm_max_mean,
        args.norm_min_disp,
        args.norm_max_value,
        args.ct_annot_n_neighbors,
        args.ct_annot_leiden_res,
        args.cpu,
        args.mem
    )
