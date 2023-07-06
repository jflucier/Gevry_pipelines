import argparse
import os
import pickle
# import matplotlib.pyplot as plt
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget
import dill
from IPython.core.display import HTML


def run_motif_enrichment(work_dir, tmp_dir, otsu, top3k, markers, scores_db, rankings_db, motif_annotation, motifs_version, species, cpu, overwrite):

    #fig = plt.figure(figsize=(10, 10))

    ### Cistarget databases

    ## custom db creation
    # You can choose to compute this database yourself by scoring the consensus peaks generated
    # in the scATAC-seq analysis using a set of motifs. The advantage of creating a sample
    # specific database is that you can potentially pick up more target regions, given that only
    # regions included/overlappig with regions in the cistarget database will be used
    # for the SCENIC+ analysis.
    # https://github.com/aertslab/create_cisTarget_databases

    # if os.path.isfile(os.path.join(work_dir, 'motifs/menr.pkl')) \
    #         and not overwrite:
    #     print(f"Reusing motif results")
    #     menr = pickle.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))
    # else:

    # use precomputed db
    print(f"loading region_bin_topics_otsu: {otsu}")
    region_bin_topics_otsu = pickle.load(open(otsu, 'rb'))
    print(f"loading region_bin_topics_top3k: {top3k}")
    region_bin_topics_top3k = pickle.load(open(top3k, 'rb'))
    print(f"loading markers_dict: {markers}")
    markers_dict = pickle.load(open(markers, 'rb'))

    # Convert to dictionary of pyranges objects.
    print(f"convert region to pyrange")
    region_sets = {'topics_otsu': {}, 'topics_top_3': {}, 'DARs': {}}
    for topic in region_bin_topics_otsu.keys():
        regions = region_bin_topics_otsu[topic].index[
            region_bin_topics_otsu[topic].index.str.startswith('chr')]  # only keep regions on known chromosomes
        region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

    for topic in region_bin_topics_top3k.keys():
        regions = region_bin_topics_top3k[topic].index[
            region_bin_topics_top3k[topic].index.str.startswith('chr')]  # only keep regions on known chromosomes
        region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

    for DAR in markers_dict.keys():
        regions = markers_dict[DAR].index[
            markers_dict[DAR].index.str.startswith('chr')]  # only keep regions on known chromosomes
        region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

    for key in region_sets.keys():
        print(f'{key}: {region_sets[key].keys()}')

    # Define rankings, score and motif annotation database.

    # rankings_db = os.path.join(work_dir, 'data/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
    # scores_db = os.path.join(work_dir, 'data/hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
    # # motif_annotation = os.path.join(work_dir, 'data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')
    # motif_annotation = os.path.join(work_dir, 'data/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')

    # run cistarget based and DEM based motif enrichment analysis with or without promoter regions.
    if not os.path.exists(os.path.join(work_dir, 'motifs')):
        os.makedirs(os.path.join(work_dir, 'motifs'))

    print(f"running pycistarget")
    run_pycistarget(
        region_sets=region_sets,
        species=species,
        save_path=os.path.join(work_dir, 'motifs'),
        ctx_db_path=rankings_db,
        dem_db_path=scores_db,
        path_to_motif_annotations=motif_annotation,
        run_without_promoters=True,
        n_cpu=cpu,
        _temp_dir=os.path.join(tmp_dir, 'ray_spill'),
        annotation_version=motifs_version
    )

    # explore some of the results. Below we show the motifs found for topic 8 (specific to B-cells) using DEM.
    # menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

    # print(f"Loading motif results")
    # menr = pickle.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

    # print(f"Output motif for each topics")
    # menr['DEM_topics_otsu_All'].DEM_results('Topic8')
    # with open(os.path.join(work_dir, 'motifs/fig10.html'), 'w') as f:
    #     f.write(menr['DEM_topics_otsu_All'].DEM_results('Topic8').data)


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-w", "--workdir", help="your working directory", required=True)
    argParser.add_argument("--tmp", nargs='?',
                           help="Temp directory (default /tmp)", const="/tmp", default="/tmp")

    argParser.add_argument(
        "--specie",
        nargs='?',
        help="Species from which genomic coordinates come from, options are: homo_sapiens, mus_musculus, drosophila_melanogaster (default homo_sapiens)",
        const="homo_sapiens",
        default="homo_sapiens"
    )

    argParser.add_argument("--otsu", nargs='?', help="Path to region bin topic otsu pickle (default <<workdir>>/scATAC/candidate_enhancers/region_bin_topics_otsu.pkl)"
                                                     "<<workdir>>/scATAC/candidate_enhancers/region_bin_topics_otsu.pkl",
                           const="", default="")

    argParser.add_argument("--top3k", nargs='?', help="Path to region bin topic top3k pickle (default <<workdir>>/scATAC/candidate_enhancers/region_bin_topics_top3k.pkl)"
                                                      "<<workdir>>/scATAC/candidate_enhancers/region_bin_topics_top3k.pkl",
                           const="", default="")

    argParser.add_argument("--markers", nargs='?', help="Path to marker dictionary pickle (default <<workdir>>/scATAC/candidate_enhancers/markers_dict.pkl)"
                                                        "<<workdir>>/scATAC/candidate_enhancers/markers_dict.pkl",
                           const="", default="")

    argParser.add_argument("--scores_db", nargs='?', help="Path to score feather file (default <<workdir>>/data/hg38_screen_v10_clust.regions_vs_motifs.scores.feather)",
                           const="", default="")

    argParser.add_argument("--rank_db", nargs='?', help="Path to ranking feather file (default <<workdir>>/data/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather)",
                           const="", default="")

    argParser.add_argument("--motifs", nargs='?', help="Path to motif annotation table (default <<workdir>>/data/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl)",
                           const="", default="")

    argParser.add_argument("--motifs_version", nargs='?', help="Motif annotation version (default v10nr_clust)",
                           const="v10nr_clust", default="v10nr_clust")

    argParser.add_argument("--cpu", nargs='?',
                           help="Number of cpu to use (default 24)", const=24,
                           type=int, default=24)
    argParser.add_argument(
        "--overwrite",
        nargs='?',
        help=f"Recalculate all steps even if they completed sucessfully.",
        const=True,
        default=False
    )

    args = argParser.parse_args()

    if args.otsu == "":
        args.otsu = args.workdir + "/scATAC/candidate_enhancers/region_bin_topics_otsu.pkl"

    if args.top3k == "":
        args.top3k = args.workdir + "/scATAC/candidate_enhancers/region_bin_topics_top3k.pkl"

    if args.markers == "":
        args.markers = args.workdir + "/scATAC/candidate_enhancers/markers_dict.pkl"

    if args.scores_db == "":
        args.scores_db = args.workdir + "/data/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"

    if args.rank_db == "":
        args.rank_db = args.workdir + "/data/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"

    if args.motifs == "":
        args.motifs = args.workdir + "/data/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl"

    for k, v in vars(args).items():
        print(f"Input {k}: {v}")
    
    # Species from which genomic coordinates come from, options are: homo_sapiens, mus_musculus,
    # drosophila_melanogaster and gallus_gallus.
    # sp = 'homo_sapiens'
    # cpu = args.cpu
    # overwrite = False

    run_motif_enrichment(
        args.workdir,
        args.tmp,
        args.otsu,
        args.top3k,
        args.markers,
        args.scores_db,
        args.rank_db,
        args.motifs,
        args.motifs_version,
        args.specie,
        args.cpu,
        args.overwrite
    )
