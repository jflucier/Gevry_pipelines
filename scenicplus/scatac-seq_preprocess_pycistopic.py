import argparse
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

target_urls = {
    'hs': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
    'mm': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes',
    'dm': 'http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes'
}

pbm_datasets = {
    'hs': 'hsapiens_gene_ensembl',
    'mm': 'mmusculus_gene_ensembl',
    'dm': 'dmelanogaster_gene_ensembl'
}


def get_chromsizes(sp):
    target_url = target_urls[sp]
    chromsizes = pd.read_csv(target_url, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0] * chromsizes.shape[0]
    chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
    # Exceptionally in this case, to agree with CellRangerARC annotations
    chromsizes['Chromosome'] = [
        chromsizes['Chromosome'][x].replace('v', '.')
        for x in range(len(chromsizes['Chromosome']))
    ]

    chromsizes['Chromosome'] = [
        chromsizes['Chromosome'][x].split('_')[1]
        if len(chromsizes['Chromosome'][x].split('_')) > 1
        else chromsizes['Chromosome'][x]
        for x in range(len(chromsizes['Chromosome']))
    ]
    return pr.PyRanges(chromsizes)


def run_scattac_preprocess(work_dir, tmp_dir, atac_frag_file, scrna_data, scrna_sample_id, cpu, shift, ext_size, qval,
                            peak_half_width, path_to_blacklist, sp, overwrite):
    print(f"initialisation of variable and output folder {work_dir}/scATAC")
    if not os.path.exists(os.path.join(work_dir, 'scATAC')):
        os.makedirs(os.path.join(work_dir, 'scATAC'))

    outdir = f"{work_dir}/scATAC"

    # atac_frag_file = os.path.join(work_dir, 'data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz')
    print(f"reading ATAC fragments tsv")
    fragments_dict = {'10x_pbmc': atac_frag_file}

    # Generate pseudobulk ATAC-seq profiles, call peaks and generate a consensus peak set
    print(f"reading scRNA data ")
    adata = sc.read_h5ad(scrna_data)
    scRNA_bc = adata.obs_names
    cell_data = adata.obs
    cell_data['sample_id'] = scrna_sample_id
    # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
    cell_data['celltype'] = cell_data['celltype'].astype(str)
    # del (adata)

    # Get chromosome sizes (for hg38 here)
    print(f"Get chromsize from UCSC")
    chromsizes = get_chromsizes(sp)

    if os.path.isfile(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl')) \
            and not overwrite:
        print(f"Reusing previously computed pseudobulk ATAC-seq profiles")
        bed_paths = pickle.load(
            open(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'rb'))
        bw_paths = pickle.load(
            open(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'rb'))
    else:
        print(f"Generate pseudobulk ATAC-seq profiles")
        if not os.path.exists(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files')):
            os.makedirs(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files'))
        if not os.path.exists(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bw_files')):
            os.makedirs(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bw_files'))
        if not os.path.exists(os.path.join(tmp_dir, 'ray_spill')):
            os.makedirs(os.path.join(tmp_dir, 'ray_spill'))

        bw_paths, bed_paths = export_pseudobulk(
            input_data=cell_data,
            variable='celltype',
            # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
            sample_id_col='sample_id',
            chromsizes=chromsizes,
            bed_path=os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files/'),
            # specify where pseudobulk_bed_files should be stored
            bigwig_path=os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bw_files/'),
            # specify where pseudobulk_bw_files should be stored
            path_to_fragments=fragments_dict,  # location of fragment files
            n_cpu=cpu,  # specify the number of cores to use, we use ray for multi processing
            normalize_bigwig=True,
            remove_duplicates=True,
            _temp_dir=os.path.join(tmp_dir, 'ray_spill'),
            split_pattern='-'
        )

        print(f"Dump pseudobulk ATAC-seq profiles bed and bigwig filepaths into pickle")
        pickle.dump(
            bed_paths,
            open(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'wb')
        )
        pickle.dump(
            bw_paths,
            open(os.path.join(outdir, 'consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'wb')
        )

    if os.path.isfile(os.path.join(outdir, 'consensus_peak_calling/MACS/narrow_peaks_dict.pkl')) \
            and not overwrite:
        print(f"Reusing macs2 peak call")
        narrow_peaks_dict = pickle.load(
            open(os.path.join(outdir, 'consensus_peak_calling/MACS/narrow_peaks_dict.pkl'), 'rb'))
    else:
        print(f"Call peaks using macs2")
        macs_path = 'macs2'
        # Run peak calling
        narrow_peaks_dict = peak_calling(
            macs_path,
            bed_paths,
            os.path.join(outdir, 'consensus_peak_calling/MACS/'),
            genome_size=sp,
            n_cpu=cpu,
            input_format='BEDPE',
            shift=shift,
            ext_size=ext_size,
            keep_dup='all',
            q_value=qval,
            _temp_dir=os.path.join(tmp_dir, 'ray_spill')
        )

        print(f"Dump peaks dict into pickle")
        pickle.dump(
            narrow_peaks_dict,
            open(os.path.join(outdir, 'consensus_peak_calling/MACS/narrow_peaks_dict.pkl'), 'wb')
        )

    if not os.path.isfile(os.path.join(outdir, 'consensus_peak_calling/consensus_regions.bed')) \
            and not overwrite:
        # Get consensus peaks
        print(f"Get consensus peaks")
        consensus_peaks = get_consensus_peaks(
            narrow_peaks_dict,
            peak_half_width,
            chromsizes=chromsizes,
            path_to_blacklist=path_to_blacklist
        )

        print(f"Output consensus peaks to BED")
        consensus_peaks.to_bed(
            path=os.path.join(outdir, 'consensus_peak_calling/consensus_regions.bed'),
            keep=True,
            compression='infer',
            chain=False
        )

    if os.path.isfile(os.path.join(outdir, 'quality_control/metadata_bc.pkl')) \
            and not overwrite:
        print(f"Reusing QC")
        metadata_bc = pickle.load(open(os.path.join(outdir, 'quality_control/metadata_bc.pkl'), 'rb'))
        profile_data_dict = pickle.load(open(os.path.join(outdir, 'quality_control/profile_data_dict.pkl'), 'rb'))
        path_to_regions = {scrna_sample_id: os.path.join(outdir, 'consensus_peak_calling/consensus_regions.bed')}
    else:
        ###Quality control
        print(f"QC using ensembl {pbm_datasets[sp]}")
        dataset = pbm.Dataset(name=pbm_datasets[sp], host='http://www.ensembl.org')
        annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name',
                                          'transcript_biotype'])
        annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype=str)
        filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
        annot = annot[~filter]
        annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
        annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
        annot = annot[annot.Transcript_type == 'protein_coding']
        path_to_regions = {scrna_sample_id: os.path.join(outdir, 'consensus_peak_calling/consensus_regions.bed')}

        print(f"Compute QC stats")
        metadata_bc, profile_data_dict = compute_qc_stats(
            fragments_dict=fragments_dict,
            tss_annotation=annot,
            stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'frip'],
            label_list=None,
            path_to_regions=path_to_regions,
            n_cpu=cpu,
            valid_bc=None,
            n_frag=100,
            n_bc=None,
            tss_flank_window=1000,
            tss_window=50,
            tss_minimum_signal_window=100,
            tss_rolling_window=10,
            remove_duplicates=True,
            _temp_dir=os.path.join(tmp_dir + 'ray_spill')
        )

        if not os.path.exists(os.path.join(outdir, 'quality_control')):
            os.makedirs(os.path.join(outdir, 'quality_control'))

        print(f"Dump QC into pickle")
        pickle.dump(
            metadata_bc,
            open(os.path.join(outdir, 'quality_control/metadata_bc.pkl'), 'wb')
        )

        pickle.dump(
            profile_data_dict,
            open(os.path.join(outdir, 'quality_control/profile_data_dict.pkl'), 'wb')
        )

    # Filter cell barcodes.
    if os.path.isfile(os.path.join(outdir, 'quality_control/bc_passing_filters.pkl')) \
            and not overwrite:
        print(f"Reusing barcode passing filters")
        bc_passing_filters = pickle.load(open(os.path.join(outdir, 'quality_control/bc_passing_filters.pkl'), 'rb'))
    else:
        # Plot barcode stats in one figure
        print(f"Filter cell barcodes and export figure")
        QC_filters = {
            'Log_unique_nr_frag': [3.3, None],
            'FRIP': [0.45, None],
            'TSS_enrichment': [5, None],
            'Dupl_rate': [None, None]

        }

        fig = plt.figure(figsize=(10, 10))
        c = 0
        tmp_passing_filters = []
        img_lst = []
        for t in QC_filters:
            if t == 'Log_unique_nr_frag':
                continue

            print(f"exporting figure panel {t}")
            if t == 'Dupl_rate':
                t_fig = plot_barcode_metrics(
                    metadata_bc[scrna_sample_id],
                    var_x='Log_unique_nr_frag',
                    var_y='Dupl_rate',
                    min_x=QC_filters['Log_unique_nr_frag'][0],
                    max_x=QC_filters['Log_unique_nr_frag'][1],
                    min_y=QC_filters['Dupl_rate'][0],
                    max_y=QC_filters['Dupl_rate'][1],
                    return_cells=False,
                    return_fig=True,
                    plot=False,
                    plot_as_hexbin=True
                )
            else:
                # Return figure to plot together with other metrics, and cells passing filters. Figure will be saved as pdf.
                t_fig, t_filter = plot_barcode_metrics(
                    metadata_bc[scrna_sample_id],
                    var_x='Log_unique_nr_frag',
                    var_y=t,
                    min_x=QC_filters['Log_unique_nr_frag'][0],
                    max_x=QC_filters['Log_unique_nr_frag'][1],
                    min_y=QC_filters[t][0],
                    max_y=QC_filters[t][1],
                    return_cells=True,
                    return_fig=True,
                    plot=False
                )

                if len(tmp_passing_filters) == 0:
                    tmp_passing_filters = set(t_filter)
                else:
                    tmp_passing_filters = list((tmp_passing_filters & set(t_filter)))

            # c += 1
            # plt.subplot(1, 3, c)
            # img = fig2img(t_fig)
            img_lst.append(t_fig)

        fig = plt.figure(figsize=(15, 10))
        plt.subplot(1, 3, 1)
        img = fig2img(img_lst[0])
        plt.imshow(img)
        plt.axis('off')
        plt.subplot(1, 3, 2)
        img = fig2img(img_lst[1])
        plt.imshow(img)
        plt.axis('off')
        plt.subplot(1, 3, 3)
        img = fig2img(img_lst[2])
        plt.imshow(img)
        plt.axis('off')
        plt.savefig(outdir + '/fig6.png', dpi=500)

        print(f"{len(tmp_passing_filters)} cell barcodes pass scATAC-seq filtering")
        bc_passing_filters = {scrna_sample_id: []}
        bc_passing_filters[scrna_sample_id] = tmp_passing_filters
        print(f"Outputting barcodes passed QC stats to pickle")
        pickle.dump(
            bc_passing_filters,
            open(os.path.join(outdir, 'quality_control/bc_passing_filters.pkl'), 'wb')
        )

    print(
        f"{len(list(set(bc_passing_filters['10x_pbmc']) & set(scRNA_bc)))} cell barcodes pass both scATAC-seq and scRNA-seq based filtering")

    if os.path.isfile(os.path.join(outdir, 'cistopic_obj.tmp.pkl')) \
            and not overwrite:
        print(f"Reusing cistopic object created from fragments")
        cistopic_obj = pickle.load(open(os.path.join(outdir, 'cistopic_obj.pkl'), 'rb'))
    else:
        print(f"Creating cistopic object from fragments")
        # key = '10x_pbmc'
        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=fragments_dict[scrna_sample_id],
            path_to_regions=path_to_regions[scrna_sample_id],
            path_to_blacklist=path_to_blacklist,
            metrics=metadata_bc[scrna_sample_id],
            valid_bc=list(set(bc_passing_filters[scrna_sample_id]) & set(scRNA_bc)),
            n_cpu=cpu,
            project=scrna_sample_id,
            split_pattern='-'
        )
        cistopic_obj.add_cell_data(cell_data, split_pattern='-')
        print(cistopic_obj)

        print(f"Dumping cistopic object to pickle")
        pickle.dump(
            cistopic_obj,
            open(os.path.join(outdir, 'cistopic_obj.tmp.pkl'), 'wb')
        )

    # Run topic modeling.
    if os.path.isfile(os.path.join(outdir, 'models/10x_pbmc_models_500_iter_LDA.pkl')) \
            and not overwrite:
        print(f"Reusing models")
        models = pickle.load(open(os.path.join(outdir, 'models/10x_pbmc_models_500_iter_LDA.pkl'), 'rb'))
    else:
        print(f"Run topic modeling.")
        if not os.path.exists(os.path.join(outdir, 'models')):
            os.makedirs(os.path.join(outdir, 'models'))

        models = run_cgs_models(
            cistopic_obj,
            n_topics=[2, 4, 10, 16, 32, 48],
            n_cpu=cpu,
            n_iter=500,
            random_state=555,
            alpha=50,
            alpha_by_topic=True,
            eta=0.1,
            eta_by_topic=False,
            save_path=os.path.join(outdir, 'models'),
            _temp_dir=os.path.join(tmp_dir + '/ray_spill')
        )

        print(f"Dumping models to pickle")
        pickle.dump(
            models,
            open(os.path.join(outdir, 'models/10x_pbmc_models_500_iter_LDA.pkl'), 'wb')
        )

    if os.path.isfile(os.path.join(outdir, 'cistopic_obj.pkl')) \
            and not overwrite:
        print(f"Reusing models")
        cistopic_obj = pickle.load(open(os.path.join(outdir, 'cistopic_obj.pkl'), 'rb'))
    else:
        print(f"Evaluate models")
        model = evaluate_models(
            models,
            select_model=None,
            return_model=True,
            metrics=['Arun_2010', 'Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
            plot_metrics=False,
            save=outdir + '/model_selection.pdf'
        )
        plt.savefig(outdir + '/model_selection.png', bbox_inches="tight")

        print(f"Save models in cistopic_obj")
        cistopic_obj.add_LDA_model(model)

        print(f"Dumping cistopic_obj to pickle")
        pickle.dump(
            cistopic_obj,
            open(os.path.join(outdir, 'cistopic_obj.pkl'), 'wb')
        )

    ### Visualization

    # We can use the cell-topic probabilities to generate dimensionality reductions.
    print(f"Generate dimensionality reductions")
    run_umap(cistopic_obj, target='cell', scale=True)
    plot_metadata(
        cistopic_obj,
        reduction_name='UMAP',
        variables=['celltype'],
        save=outdir + '/dimensionality_reduction_label.pdf'
    )
    plt.savefig(outdir + '/dimensionality_reduction_label.png', bbox_inches="tight")

    # We can also plot the cell-topic probabilities on the UMAP, to visualize their cell type specifiticy.
    print(f"Plotting the cell-topic probabilities on the UMAP, to visualize their cell type specifiticy. Output: "
          f"dimensionality_reduction_topic_uncorrected.*")
    plot_topic(
        cistopic_obj,
        reduction_name='UMAP',
        num_columns=4,
        save=outdir + '/dimensionality_reduction_topic_uncorrected.pdf'
    )
    plt.savefig(outdir + '/dimensionality_reduction_topic_uncorrected.png', bbox_inches="tight")

    ### Inferring candidate enhancer regions

    # First we will binarize the topics using the otsu method and by taking the top 3k regions per topic.
    if os.path.isfile(os.path.join(outdir, 'candidate_enhancers/region_bin_topics_otsu.pkl')) \
            and not overwrite:
        print(f"Reusing models")
        region_bin_topics_otsu = pickle.load(
            open(os.path.join(outdir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
        region_bin_topics_top3k = pickle.load(
            open(os.path.join(outdir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
        markers_dict = pickle.load(open(os.path.join(outdir, 'candidate_enhancers/markers_dict.pkl'), 'rb'))
    else:
        print(f"Binarizing the topics using the otsu method")
        region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
        region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop=3000)

        # Next we will calculate DARs per cell type
        print(f"Calculate Differentially Accessible Regions (DAR) per cell type")
        imputed_acc_obj = impute_accessibility(
            cistopic_obj,
            selected_cells=None,
            selected_regions=None,
            scale_factor=10 ** 6
        )

        print(f"log-normalize the DAR imputed data.")
        normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10 ** 4)
        print(f"identify highly variable regions.")
        variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot=False)
        print(f"identify differentially accessible regions by celltype")
        markers_dict = find_diff_features(
            cistopic_obj,
            imputed_acc_obj,
            variable='celltype',
            var_features=variable_regions,
            split_pattern='-'
        )

        if not os.path.exists(os.path.join(outdir, 'candidate_enhancers')):
            os.makedirs(os.path.join(outdir, 'candidate_enhancers'))

        print(f"Dumping binarized topics to pickle")
        pickle.dump(
            region_bin_topics_otsu,
            open(os.path.join(outdir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb')
        )
        print(f"Dumping binarized topics top3k to pickle")
        pickle.dump(
            region_bin_topics_top3k,
            open(os.path.join(outdir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb')
        )
        print(f"Dumping DAR by celltype to pickle")
        pickle.dump(
            markers_dict,
            open(os.path.join(outdir, 'candidate_enhancers/markers_dict.pkl'), 'wb')
        )


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-w", "--workdir", help="your working directory", required=True)
    argParser.add_argument("--sample", help="The sample id", required=True)
    argParser.add_argument("--frag_file", help="Path to ATAC fragments file", required=True)
    argParser.add_argument("--tmp", nargs='?',
                           help="Temp directory. Defaults to /tmp", const="/tmp", default="/tmp")
    argParser.add_argument("--scrna", nargs='?',
                           help=f"Scanpy scRNA data. Defaults to <<workdir>>/scRNA/adata.h5ad", const="", default="")
    argParser.add_argument(
        "--cpu",
        nargs='?',
        help="Number of cpu to use",
        const=24,
        type=int,
        default=24
    )

    argParser.add_argument(
        "--shift",
        nargs='?',
        help="To set an arbitrary shift in bp. For finding enriched cutting sites(such as in ATAC - seq) a shift of 73bp is recommended.Default: 73",
        const=73,
        type=int,
        default=73
    )

    argParser.add_argument(
        "--ext_size",
        nargs='?',
        help="To extend reads in 5’->3’ direction to fix-sized fragment. For ATAC-seq data, a extension of 146 bp is recommended. Default: 146",
        const=146,
        type=int,
        default=146
    )

    argParser.add_argument(
        "--q_value",
        nargs='?',
        help="The q-value (minimum FDR) cutoff to call significant regions. Default: 0.05",
        const=0.05,
        type=float,
        default=0.05
    )

    argParser.add_argument(
        "--peak_half_width",
        nargs='?',
        help="Number of base pairs that each summit will be extended in each direction. Default: 250",
        const=250,
        type=int,
        default=250
    )

    argParser.add_argument(
        "--blacklist_regions",
        nargs='?',
        help=f"Path to bed file containing blacklist regions (Amemiya et al., 2019). Defaults to <<workdir>>/data/hg38-blacklist.v2.bed",
        const="",
        default=""
    )

    argParser.add_argument(
        "--overwrite",
        nargs='?',
        help=f"Recalculate all steps even if they completed sucessfully.",
        const=True,
        default=False
    )

    argParser.add_argument(
        "--specie",
        nargs='?',
        help="Species from from which genome size will be inputted to MACS2, options are: homo_sapiens, mus_musculus, drosophila_melanogaster",
        const="homo_sapiens",
        default="homo_sapiens"
    )

    args = argParser.parse_args()

    if args.specie == "homo_sapiens":
        args.specie = "hs"
    elif args.specie == "mus_musculus":
        args.specie = "mm"
    elif args.specie == "drosophila_melanogaster":
        args.specie = "dm"
    else:
        print(
            f"Unrecongnised specie provided: {args.specie}. Please provide one of the following: homo_sapiens, mus_musculus, drosophila_melanogaster.")
        argParser.print_help()
        exit(1)

    if args.scrna == "":
        args.scrna = args.workdir + "/scRNA/adata.h5ad"

    if args.blacklist_regions == "":
        args.blacklist_regions = args.workdir + "/data/hg38-blacklist.v2.bed"

    for k, v in vars(args).items():
        print(f"Input {k}: {v}")


    # sample_id = '10x_pbmc'
    # cpu = args.cpu

    # To set an arbitrary shift in bp. For finding enriched cutting sites(such as in ATAC - seq) a shift of 73bp is recommended.Default: 73.
    # shift = 73

    # To extend reads in 5’->3’ direction to fix-sized fragment. For ATAC-seq data, a extension of 146 bp is recommended. Default: 146.
    # ext_size = 146

    # The q-value (minimum FDR) cutoff to call significant regions. Default: 0.05.
    # q_value = 0.05

    # for consensus peaks step
    # Number of base pairs that each summit will be extended in each direction.
    # peak_half_width = 250
    # Path to bed file containing blacklist regions (Amemiya et al., 2019).
    # path_to_blacklist = args.workdir + '/data/hg38-blacklist.v2.bed'

    # overwrite = False

    run_scattac_preprocess(
        args.workdir,
        args.tmp,
        args.frag_file,
        args.scrna,
        args.sample,
        args.cpu,
        args.shift,
        args.ext_size,
        args.q_value,
        args.peak_half_width,
        args.blacklist_regions,
        args.specie,
        args.overwrite
    )
