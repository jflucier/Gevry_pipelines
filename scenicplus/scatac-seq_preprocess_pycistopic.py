import os
import pycisTopic
import scanpy as sc
import pyranges as pr
import requests
import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
import pickle
from pycisTopic.pseudobulk_peak_calling import peak_calling
from pycisTopic.iterative_peak_calling import *
import pybiomart as pbm
from pycisTopic.qc import *
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
import matplotlib.pyplot as plt


def run_scattac_precprocess(work_dir, tmp_dir, atac_frag_file, scrna_data):
    print(f"initialisation of variable and output folder {work_dir}/scATAC")
    if not os.path.exists(os.path.join(work_dir, 'scATAC')):
        os.makedirs(os.path.join(work_dir, 'scATAC'))

    outdir = f"{work_dir}/scATAC"

    # atac_frag_file = os.path.join(work_dir, 'data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz')
    fragments_dict = {'10x_pbmc': atac_frag_file}

    # Generate pseudobulk ATAC-seq profiles, call peaks and generate a consensus peak set
    adata = sc.read_h5ad(scrna_data)
    cell_data = adata.obs
    cell_data['sample_id'] = '10x_pbmc'
    cell_data['celltype'] = cell_data['celltype'].astype(
        str)  # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
    del (adata)

    # Get chromosome sizes (for hg38 here)
    target_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
    chromsizes = pd.read_csv(target_url, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0] * chromsizes.shape[0]
    chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
    # Exceptionally in this case, to agree with CellRangerARC annotations
    chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in
                                range(len(chromsizes['Chromosome']))]
    chromsizes['Chromosome'] = [
        chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else
        chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
    chromsizes = pr.PyRanges(chromsizes)

    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data,
        variable='celltype',
        # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
        sample_id_col='sample_id',
        chromsizes=chromsizes,
        bed_path=os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/'),
        # specify where pseudobulk_bed_files should be stored
        bigwig_path=os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bw_files/'),
        # specify where pseudobulk_bw_files should be stored
        path_to_fragments=fragments_dict,  # location of fragment fiels
        n_cpu=8,  # specify the number of cores to use, we use ray for multi processing
        normalize_bigwig=True,
        remove_duplicates=True,
        _temp_dir=os.path.join(tmp_dir, 'ray_spill'),
        split_pattern='-'
    )

    pickle.dump(bed_paths,
                open(os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'wb'))
    pickle.dump(bw_paths,
                open(os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'wb'))

    bed_paths = pickle.load(
        open(os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'rb'))
    bw_paths = pickle.load(
        open(os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'rb'))
    macs_path = 'macs2'
    # Run peak calling
    narrow_peaks_dict = peak_calling(macs_path,
                                     bed_paths,
                                     os.path.join(work_dir, 'scATAC/consensus_peak_calling/MACS/'),
                                     genome_size='hs',
                                     n_cpu=8,
                                     input_format='BEDPE',
                                     shift=73,
                                     ext_size=146,
                                     keep_dup='all',
                                     q_value=0.05,
                                     _temp_dir=os.path.join(tmp_dir, 'ray_spill'))

    pickle.dump(narrow_peaks_dict,
                open(os.path.join(work_dir, 'scATAC/consensus_peak_calling/MACS/narrow_peaks_dict.pkl'), 'wb'))

    # Other param
    peak_half_width = 250
    path_to_blacklist = work_dir + '/hg38-blacklist.v2.bed'
    # Get consensus peaks
    consensus_peaks = get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=chromsizes,
                                          path_to_blacklist=path_to_blacklist)

    consensus_peaks.to_bed(
        path=os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
        keep=True,
        compression='infer',
        chain=False)

    ###Quality control

    dataset = pbm.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name',
                                      'transcript_biotype'])
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype=str)
    filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot = annot[annot.Transcript_type == 'protein_coding']
    path_to_regions = {'10x_pbmc': os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed')}

    metadata_bc, profile_data_dict = compute_qc_stats(
        fragments_dict=fragments_dict,
        tss_annotation=annot,
        stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'frip'],
        label_list=None,
        path_to_regions=path_to_regions,
        n_cpu=1,
        valid_bc=None,
        n_frag=100,
        n_bc=None,
        tss_flank_window=1000,
        tss_window=50,
        tss_minimum_signal_window=100,
        tss_rolling_window=10,
        remove_duplicates=True,
        _temp_dir=os.path.join(tmp_dir + 'ray_spill'))

    if not os.path.exists(os.path.join(work_dir, 'scATAC/quality_control')):
        os.makedirs(os.path.join(work_dir, 'scATAC/quality_control'))

    pickle.dump(metadata_bc,
                open(os.path.join(work_dir, 'scATAC/quality_control/metadata_bc.pkl'), 'wb'))

    pickle.dump(profile_data_dict,
                open(os.path.join(work_dir, 'scATAC/quality_control/profile_data_dict.pkl'), 'wb'))

    # Filter cell barcodes.
    # [min,  #max]
    QC_filters = {
        'Log_unique_nr_frag': [3.3, None],
        'FRIP': [0.45, None],
        'TSS_enrichment': [5, None],
        'Dupl_rate': [None, None]

    }

    # Return figure to plot together with other metrics, and cells passing filters. Figure will be saved as pdf.
    FRIP_NR_FRAG_fig, FRIP_NR_FRAG_filter = plot_barcode_metrics(metadata_bc['10x_pbmc'],
                                                                 var_x='Log_unique_nr_frag',
                                                                 var_y='FRIP',
                                                                 min_x=QC_filters['Log_unique_nr_frag'][0],
                                                                 max_x=QC_filters['Log_unique_nr_frag'][1],
                                                                 min_y=QC_filters['FRIP'][0],
                                                                 max_y=QC_filters['FRIP'][1],
                                                                 return_cells=True,
                                                                 return_fig=True,
                                                                 plot=False)
    # Return figure to plot together with other metrics, and cells passing filters
    TSS_NR_FRAG_fig, TSS_NR_FRAG_filter = plot_barcode_metrics(metadata_bc['10x_pbmc'],
                                                               var_x='Log_unique_nr_frag',
                                                               var_y='TSS_enrichment',
                                                               min_x=QC_filters['Log_unique_nr_frag'][0],
                                                               max_x=QC_filters['Log_unique_nr_frag'][1],
                                                               min_y=QC_filters['TSS_enrichment'][0],
                                                               max_y=QC_filters['TSS_enrichment'][1],
                                                               return_cells=True,
                                                               return_fig=True,
                                                               plot=False)
    # Return figure to plot together with other metrics, but not returning cells (no filter applied for the duplication rate  per barcode)
    DR_NR_FRAG_fig = plot_barcode_metrics(metadata_bc['10x_pbmc'],
                                          var_x='Log_unique_nr_frag',
                                          var_y='Dupl_rate',
                                          min_x=QC_filters['Log_unique_nr_frag'][0],
                                          max_x=QC_filters['Log_unique_nr_frag'][1],
                                          min_y=QC_filters['Dupl_rate'][0],
                                          max_y=QC_filters['Dupl_rate'][1],
                                          return_cells=False,
                                          return_fig=True,
                                          plot=False,
                                          plot_as_hexbin=True)

    # Plot barcode stats in one figure
    fig = plt.figure(figsize=(10, 10))
    plt.subplot(1, 3, 1)
    img = fig2img(FRIP_NR_FRAG_fig)
    plt.imshow(img)
    plt.axis('off')
    plt.subplot(1, 3, 2)
    img = fig2img(TSS_NR_FRAG_fig)
    plt.imshow(img)
    plt.axis('off')
    plt.subplot(1, 3, 3)
    img = fig2img(DR_NR_FRAG_fig)
    plt.imshow(img)
    plt.axis('off')
    plt.savefig(work_dir + '/fig6.png')

    bc_passing_filters = {'10x_pbmc': []}
    bc_passing_filters['10x_pbmc'] = list((set(FRIP_NR_FRAG_filter) & set(TSS_NR_FRAG_filter)))
    pickle.dump(bc_passing_filters,
                open(os.path.join(work_dir, 'scATAC/quality_control/bc_passing_filters.pkl'), 'wb'))
    print(f"{len(bc_passing_filters['10x_pbmc'])} barcodes passed QC stats")

    ### Creating a cisTopic object and topic modeling
    adata = sc.read_h5ad(os.path.join(work_dir, 'scRNA/adata.h5ad'))
    scRNA_bc = adata.obs_names
    cell_data = adata.obs
    cell_data['sample_id'] = '10x_pbmc'
    cell_data['celltype'] = cell_data['celltype'].astype(
        str)  # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
    del (adata)

    fragments_dict = {'10x_pbmc': os.path.join(work_dir, 'data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz')}
    path_to_regions = {'10x_pbmc': os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed')}
    path_to_blacklist = work_dir + '/hg38-blacklist.v2.bed'
    metadata_bc = pickle.load(open(os.path.join(work_dir, 'scATAC/quality_control/metadata_bc.pkl'), 'rb'))
    bc_passing_filters = pickle.load(
        open(os.path.join(work_dir, 'scATAC/quality_control/bc_passing_filters.pkl'), 'rb'))
    print(
        f"{len(list(set(bc_passing_filters['10x_pbmc']) & set(scRNA_bc)))} cell barcodes pass both scATAC-seq and scRNA-seq based filtering")

    key = '10x_pbmc'
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments=fragments_dict[key],
        path_to_regions=path_to_regions[key],
        path_to_blacklist=path_to_blacklist,
        metrics=metadata_bc[key],
        valid_bc=list(set(bc_passing_filters[key]) & set(scRNA_bc)),
        n_cpu=10,
        project=key,
        split_pattern='-')
    cistopic_obj.add_cell_data(cell_data, split_pattern='-')
    print(cistopic_obj)

    pickle.dump(cistopic_obj,
                open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

    # Run topic modeling.
    cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
    models = run_cgs_models(cistopic_obj,
                            n_topics=[2, 4, 10, 16, 32, 48],
                            n_cpu=25,
                            n_iter=500,
                            random_state=555,
                            alpha=50,
                            alpha_by_topic=True,
                            eta=0.1,
                            eta_by_topic=False,
                            save_path=None,
                            _temp_dir=os.path.join(tmp_dir + '/ray_spill'))

    if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
        os.makedirs(os.path.join(work_dir, 'scATAC/models'))

    pickle.dump(models,
                open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'wb'))

    models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'rb'))
    cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
    model = evaluate_models(models,
                            select_model=16,
                            return_model=True,
                            metrics=['Arun_2010', 'Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                            plot_metrics=False)
    plt.savefig(work_dir + '/fig7.png', bbox_inches="tight")

    cistopic_obj.add_LDA_model(model)
    pickle.dump(cistopic_obj,
                open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

    ### Visualization

    # We can use the cell-topic probabilities to generate dimensionality reductions.
    run_umap(cistopic_obj, target='cell', scale=True)
    plot_metadata(cistopic_obj, reduction_name='UMAP', variables=['celltype'])
    plt.savefig(work_dir + '/fig8.png', bbox_inches="tight")

    # We can also plot the cell-topic probabilities on the UMAP, to visualize their cell type specifiticy.
    plot_topic(cistopic_obj, reduction_name='UMAP', num_columns=4)
    plt.savefig(work_dir + '/fig9.png', bbox_inches="tight")

    ### Inferring candidate enhancer regions

    # First we will binarize the topics using the otsu method and by taking the top 3k regions per topic.
    region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
    region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop=3000)

    # Next we will calculate DARs per cell type
    imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None,
                                           scale_factor=10 ** 6)
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10 ** 4)
    variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot=False)
    markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='celltype', var_features=variable_regions,
                                      split_pattern='-')

    if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
        os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))

    pickle.dump(region_bin_topics_otsu,
                open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
    pickle.dump(region_bin_topics_top3k,
                open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
    pickle.dump(markers_dict, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-w", "--workdir", help="your working directory", required=True)
    argParser.add_argument("--frag_file", help="Path to ATAC fragments file", required=True)
    argParser.add_argument("--tmp", nargs='?',
                           help="Temp directory, defaults to /tmp", const="/tmp", default="/tmp")

    argParser.add_argument("--scrna", nargs='?',
                           help=f"Scanpy scRNA data, defaults to <<workdir>>/scRNA/adata.h5ad", const="", default="")

    args = argParser.parse_args()
    # print("args=%s\n" % args)
    # print(vars(args))
    for k, v in vars(args).items():
        print(f"Input {k}: {v}")

    if args.scrna == "":
        args.scrna = args.workdir + "/scRNA/adata.h5ad"

    run_scattac_precprocess(
        args.workdir,
        args.tmp,
        args.frag_file,
        args.scrna
    )
