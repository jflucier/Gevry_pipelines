import argparse
import pickle

import dill
import scanpy as sc
import os
import sys
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus
import numpy as np
import pybiomart as pbm
import warnings

warnings.filterwarnings("ignore")
# import pandas
# import pyranges
# Set stderr to null to avoid strange messages from ray
_stderr = sys.stderr
null = open(os.devnull, 'wb')

ensembl_version_dict = {
    '105': 'http://www.ensembl.org',
    '104': 'http://may2021.archive.ensembl.org/',
    '103': 'http://feb2021.archive.ensembl.org/',
    '102': 'http://nov2020.archive.ensembl.org/',
    '101': 'http://aug2020.archive.ensembl.org/',
    '100': 'http://apr2020.archive.ensembl.org/',
    '99': 'http://jan2020.archive.ensembl.org/',
    '98': 'http://sep2019.archive.ensembl.org/',
    '97': 'http://jul2019.archive.ensembl.org/',
    '96': 'http://apr2019.archive.ensembl.org/',
    '95': 'http://jan2019.archive.ensembl.org/',
    '94': 'http://oct2018.archive.ensembl.org/',
    '93': 'http://jul2018.archive.ensembl.org/',
    '92': 'http://apr2018.archive.ensembl.org/'
}


def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species + '_gene_ensembl', host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name',
                                      'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov = len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov


def infer_enhancer_driven_gene(work_dir, scrna_path, cistopic_path, menr_path, tmp, sample_id, ensembl_specie, tf_file, cpu, overwrite):

    # loading prev analysis
    print(f"reading scrna data {scrna_path}")
    adata = sc.read_h5ad(scrna_path)
    print(f"reading cistopic data {cistopic_path}")
    cistopic_obj = pickle.load(open(cistopic_path, 'rb'))
    print(f"reading motifs data {menr_path}")
    menr = pickle.load(open(menr_path, 'rb'))

    if os.path.isfile(os.path.join(work_dir, 'scenicplus/scplus_obj.1.pkl')) \
            and not overwrite:
        print(f"Reusing scenicplus object")
        scplus_obj = pickle.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.1.pkl'), 'rb'))
    else:
        # Create the SCENIC+ object.
        print(f"Creating scenicplus object")
        scplus_obj = create_SCENICPLUS_object(
            GEX_anndata=adata.raw.to_adata(),
            cisTopic_obj=cistopic_obj,
            menr=menr,
            bc_transform_func=lambda x: f'{x}-{sample_id}'  # function to convert scATAC-seq barcodes to scRNA-seq ones
        )
        scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
        if not os.path.exists(os.path.join(work_dir, 'scenicplus')):
            os.makedirs(os.path.join(work_dir, 'scenicplus'))

        print(f"Dumping scenicplus object to pickle")
        pickle.dump(
            scplus_obj,
            open(os.path.join(work_dir, 'scenicplus/scplus_obj.1.pkl'), 'wb')
        )

    print(f"Finding biggest gene overlap in multiple ensembl releases")
    n_overlap = {}
    for version in ensembl_version_dict.keys():
        print(f'host: {version}')
        try:
            n_overlap[version] = test_ensembl_host(scplus_obj, ensembl_version_dict[version], ensembl_specie)
        except:
            print('Host not reachable')

    v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
    print(f"version: {v} has the largest overlap, will use {ensembl_version_dict[v]} as biomart host")

    # Lets take the rensembl release with most overlap
    biomart_host = ensembl_version_dict[v]

    # only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
    scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]
    scplus_obj.dr_cell['GEX_rep'] = scplus_obj.dr_cell['GEX_rep'].iloc[:, 0:2]

    if not os.path.exists(os.path.join(work_dir, 'bedToBigBed')):
        os.makedirs(os.path.join(work_dir, 'bedToBigBed'))

    # run the analysis
    try:
        print(f"Running scenicplus")
        run_scenicplus(
            scplus_obj=scplus_obj,
            variable=['GEX_celltype'],
            species=ensembl_specie,
            assembly='hg38',
            tf_file=tf_file,
            save_path=os.path.join(work_dir, 'scenicplus'),
            biomart_host=biomart_host,
            upstream=[1000, 150000],
            downstream=[1000, 150000],
            calculate_TF_eGRN_correlation=True,
            calculate_DEGs_DARs=True,
            export_to_loom_file=True,
            export_to_UCSC_file=True,
            path_bedToBigBed='/ucsc.v386',
            n_cpu=cpu,
            _temp_dir=tmp
        )
    except Exception as e:
        # in case of failure, still save the object
        pickle.dump(
            scplus_obj,
            open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb')
        )
        raise (e)

    print(f"Ouputting results to tsv in {work_dir}/scenicplus")
    counts_df = scplus_obj.to_df('EXP')
    counts_df.to_csv(os.path.join(work_dir, 'scenicplus/gene_expression_counts.tsv'), sep="\t")

    acc_df = scplus_obj.to_df('ACC')
    acc_df.to_csv(os.path.join(work_dir, 'scenicplus/chromatin_accessibility.tsv'), sep="\t")

    scplus_obj.metadata_cell.to_csv(os.path.join(work_dir, 'scenicplus/cells_metadata.tsv'), sep="\t")

    scplus_obj.metadata_regions.to_csv(os.path.join(work_dir, 'scenicplus/regions_metadata.tsv'), sep="\t")

    scplus_obj.metadata_genes.to_csv(os.path.join(work_dir, 'scenicplus/genes_metadata.tsv'), sep="\t")

    scplus_obj.uns['eRegulon_metadata'].to_csv(os.path.join(work_dir, 'scenicplus/eRegulon_metadata.tsv'), sep="\t")


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-w", "--workdir", help="your working directory", required=True)
    argParser.add_argument("-tf", "--tf_file", help="Path to file containing genes that are TFs", required=True)
    argParser.add_argument("--sample", help="The sample id", required=True)
    argParser.add_argument("--scrna", nargs='?',
                           help=f"Scanpy scRNA data. Defaults to <<workdir>>/scRNA/adata.h5ad", const="", default="")
    argParser.add_argument("--cistopic", nargs='?',
                           help=f"cistopic object data pickle data. Defaults to <<workdir>>/scATAC/cistopic_obj2.pkl", const="", default="")
    argParser.add_argument("--menr", nargs='?',
                           help=f"menr object data pickle data. Defaults to <<workdir>>/motifs/menr.pkl",
                           const="", default="")
    argParser.add_argument("--tmp", nargs='?',
                           help="Temp directory. Defaults to /tmp", const="/tmp", default="/tmp")

    argParser.add_argument("--cpu", nargs='?',
                           help="Number of cpu to use", const=24,
                           type=int, default=24)

    argParser.add_argument(
        "--specie",
        nargs='?',
        help="Species from which data comes from. options are: homo_sapiens, mus_musculus, drosophila_melanogaster.",
        const="homo_sapiens",
        default="homo_sapiens"
    )

    argParser.add_argument(
        "--overwrite",
        nargs='?',
        help=f"Recalculate all steps even if they completed sucessfully.",
        const=True,
        default=False
    )

    args = argParser.parse_args()

    if args.scrna == "":
        args.scrna = args.workdir + "/scRNA/adata.h5ad"
    if args.cistopic == "":
        args.cistopic = args.workdir + "/scATAC/cistopic_obj.pkl"
    if args.menr == "":
        args.menr = args.workdir + "/motifs/menr.pkl"

    if args.specie == "homo_sapiens":
        args.specie = "hsapiens"
        assembly = "hg38"
    elif args.specie == "mus_musculus":
        args.specie = "mmusculus"
        assembly = "mm10"
    elif args.specie == "drosophila_melanogaster":
        args.specie="dmelanogaster"
        assembly = "dm6"
    else:
        print(
            f"Unrecongnised specie provided: {args.specie}. Please provide one of the following: homo_sapiens, mus_musculus, drosophila_melanogaster.")
        argParser.print_help()
        exit(1)

    for k, v in vars(args).items():
        print(f"Input {k}: {v}")

    # cpu = args.cpu
    # overw= False
    # sample_id = '10x_pbmc'
    # Species from which data comes from. Possible values: 'hsapiens', 'mmusculus', 'dmelanogaster'
    # ensembl_sp = 'hsapiens'

    infer_enhancer_driven_gene(
        args.workdir,
        args.scrna,
        args.cistopic,
        args.menr,
        args.tmp,
        args.sample,
        args.specie,
        args.tf_file,
        args.cpu,
        args.overwrite
    )
