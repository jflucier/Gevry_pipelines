# ILL Humann custom pipeline User Manual

----

## Contents ##

* [Requirements](#requirements)
* [Installation](#installation)
* [How to run](#how-to-run)
    * [Run scenicplus pipelines](#Run-scenicplus-pipelines)
      * [Step1: scRNA-seq preprocessing using Scanpy](#Step1-scRNA-seq-preprocessing-using-Scanpy)
      * [Step2: scATAC-seq preprocessing using pycisTopic](#Step2-scATAC-seq-preprocessing-using-pycisTopic)
      * [Step3: Motif enrichment analysis using pycistarget](#Step3-Motif-enrichment-analysis-using-pycistarget)
      * [Step4: Inferring enhancer-driven Gene Regulatory Networks using SCENICplus](#Step4-Inferring-enhancer-driven-Gene-Regulatory-Networks-using-SCENICplus)

----

## Requirements ##

All pipelines are self contained. The only requirements needed is [Apptainer](https://apptainer.org/). The apptainer executable "singularity" should be available in your path.

**Note**: On interactive node run ``module load StdEnv/2020 apptainer/1.1.5 ``. You can include this command in your ~/.bashrc file

----

## Installation ##

**If you are running on ip34, the installation is already available. Skip to next section in documentation.**

To install Gevry_pipelines you need to:

* Create a clone of the repository:

    ``git clone https://github.com/jflucier/Gevry_pipelines.git ``

    Note: Creating a clone of the repository requires [Github](https://github.com/) to be installed.

* For convenience, set environment variable G_PIPELINES in your ~/.bashrc:

    ``export G_PIPELINES=/path/to/Gevry_pipelines ``

    Note: On ip34, Gevry pipelines path is /home/def-gevrynic/programs/Gevry_pipelines

* Go to $G_PIPELINES/containers and run these commands:
```
module load StdEnv/2020 apptainer/1.1.5
cd $G_PIPELINES/containers
singularity build --force --fakeroot scenicplus.sif scenicplus.def
```


----

## Run scenicplus pipelines ##

Scenicplus analysis is composed of these steps and was developped based on the following tutorial available [here](https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html):

1. [scrna-seq_preprocess_scanpy.py](scenicplus%2Fscrna-seq_preprocess_scanpy.py): scRNA-seq preprocessing using Scanpy
2. [scatac-seq_preprocess_pycistopic.py](scenicplus%2Fscatac-seq_preprocess_pycistopic.py): scATAC-seq preprocessing using pycisTopic
3. [motif_enrichment_pycistarget.py](scenicplus%2Fmotif_enrichment_pycistarget.py): Motif enrichment analysis using pycistarget
4. [infer_enhancer-driven_gene.py](scenicplus%2Finfer_enhancer-driven_gene.py): inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
5. [downstream_analysis.py](scenicplus%2Fdownstream_analysis.py): **Still in development**

For each of these scripts, you can access script help by providing -h option like below:

```
$ singularity exec \
-B /fast/:/fast/ \
-B /home:/home \
-e $G_PIPELINES/containers/scenicplus.sif \
/miniconda3/envs/scenicplus/bin/python3 -W ignore scrna-seq_preprocess_scanpy.py -h
INFO:    underlay of /usr/share/zoneinfo/America/New_York required more than 50 (51) bind mounts
usage: scrna-seq_preprocess_scanpy.py [-h] -w WORKDIR -i INPUT [--cpu [CPU]] [--mem [MEM]] [--qc_min_genes [QC_MIN_GENES]] [--qc_min_cells [QC_MIN_CELLS]] [--filter_mito [FILTER_MITO]]
                                      [--filter_n_counts [FILTER_N_COUNTS]] [--norm_min_mean [NORM_MIN_MEAN]] [--norm_max_mean [NORM_MAX_MEAN]] [--norm_min_disp [NORM_MIN_DISP]]
                                      [--norm_max_value [NORM_MAX_VALUE]] [--ct_annot_n_neighbors [CT_ANNOT_N_NEIGHBORS]] [--ct_annot_leiden_res [CT_ANNOT_LEIDEN_RES]]

optional arguments:
  -h, --help            show this help message and exit
  -w WORKDIR, --workdir WORKDIR
                        your working directory
  -i INPUT, --input INPUT
                        your h5 input file
  --cpu [CPU]           Number of cpu to use
  --mem [MEM]           Max memery usage in Gigabyte
  --qc_min_genes [QC_MIN_GENES]
                        Only keep cells with at least <<min_genes>> genes expressed
  --qc_min_cells [QC_MIN_CELLS]
                        Only keep genes which are expressed in at least <<min_cells>> cells
  --filter_mito [FILTER_MITO]
                        Filter based on mitochondrial counts
  --filter_n_counts [FILTER_N_COUNTS]
                        Filter based on total counts
  --norm_min_mean [NORM_MIN_MEAN]
                        Normalisation min mean value
  --norm_max_mean [NORM_MAX_MEAN]
                        Normalisation max mean value
  --norm_min_disp [NORM_MIN_DISP]
                        Normalisation min dispersion value
  --norm_max_value [NORM_MAX_VALUE]
                        Clip (truncate) to this value after scaling
  --ct_annot_n_neighbors [CT_ANNOT_N_NEIGHBORS]
                        The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller
                        values result in more local data being preserved.
  --ct_annot_leiden_res [CT_ANNOT_LEIDEN_RES]
                        parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.

```

Most default values are the same as used in the tutorial. Lets now look at each step in more details:

### Step1 scRNA-seq preprocessing using Scanpy ###

This step processes single cell rna-seq. An example call would be:

```

> singularity exec \
-B /fast/:/fast/ \
-B /home:/home \
-e $G_PIPELINES/containers/scenicplus.sif \
/miniconda3/envs/scenicplus/bin/python3 -W ignore scrna-seq_preprocess_scanpy.py \
-w $G_PIPELINES/scenicplus/test \
-i $G_PIPELINES/scenicplus/test/data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5

```

The following parameters are mandatory:
- -w WORKDIR, --workdir WORKDIR: your working directory
- -i INPUT, --input INPUT: your h5 input file

These options are optional. **Please make sure default values are ok prior to running**
- cpu: Number of cpu to use (default 24)
- mem: Max memory usage in Gigabyte (default 30)
- qc_min_genes: Only keep cells with at least <<min_genes>> genes expressed (default 200)
- qc_min_cells: Only keep genes which are expressed in at least <<min_cells>> cells (default 3)
- filter_mito: Filter based on mitochondrial counts (default 25)
- filter_n_counts: Filter based on total counts (default 4300)
- norm_min_mean: Normalisation min mean value (default 0.0125)
- norm_max_mean: Normalisation max mean value (default 3)
- norm_min_disp: Normalisation min dispersion value (default 0.5)
- norm_max_value: Clip (truncate) to this value after scaling (default 10)
- ct_annot_n_neighbors: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller  values result in more local data being preserved (default 10)
- ct_annot_leiden_res: Parameter value controlling the coarseness of the clustering. Higher values lead to more clusters (default 0.8)

**Notice** 
In its current form cell type annotation uses the preprocessed data from the Scanpy tutorial as reference 
(see [here](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k_processed.html) for more information.
Please contact to change this.

### Step2 scATAC-seq preprocessing using pycisTopic ###

This step processes single cell atac-seq. An example call would be:

```
singularity exec \
-B /fast/:/fast/ \
-B /home:/home \
-e $G_PIPELINES/containers/scenicplus.sif \
/miniconda3/envs/scenicplus/bin/python3 -W ignore scatac-seq_preprocess_pycistopic.py \
-w $G_PIPELINES/scenicplus/test \
--frag_file $G_PIPELINES/scenicplus/test/data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz \
--tmp /fast/tmp --sample '10x_pbmc'

```

The following parameters are mandatory:
- workdir: your working directory
- frag_file: Path to ATAC fragments file
- sample: The sample id

These options are optional. **Please make sure default values are ok prior to running**
- tmp: Temp directory (default /tmp)
- scrna: Scanpy scRNA data (default <<workdir>>/scRNA/adata.h5ad)
- cpu: Number of cpu to use (default 24)
- shift; To set an arbitrary shift in bp. For finding enriched cutting sites(such as in ATAC - seq) a shift of 73bp is recommended (default 73)
- ext_size: To extend reads in 5’->3’ direction to fix-sized fragment. For ATAC-seq data, a extension of 146 bp is recommended (default 146)
- q_value: The q-value (minimum FDR) cutoff to call significant regions (default 0.05)
- peak_half_width: Number of base pairs that each summit will be extended in each direction (default 250)
- blacklist_regions: Path to bed file containing blacklist regions (Amemiya et al., 2019) (default <<workdir>>/data/hg38-blacklist.v2.bed)
- overwrite: Recalculate all steps even if they completed sucessfully.
- specie: Species from from which genome size will be inputted to MACS2, options are: homo_sapiens, mus_musculus, drosophila_melanogaster (default homo_sapiens)

### Step3 Motif enrichment analysis using pycistarget ###

This step identifies motif enrichment. An example call would be:

```
singularity exec \
-B /fast/:/fast/ \
-B /home:/home \
-e $G_PIPELINES/containers/scenicplus.sif \
/miniconda3/envs/scenicplus/bin/python3 -W ignore motif_enrichment_pycistarget.py \
-w $G_PIPELINES/scenicplus/test \
--tmp /fast/tmp

```

The following parameters are mandatory:
- workdir: your working directory

These options are optional. **Please make sure default values are ok prior to running**
- tmp: Temp directory (default /tmp)
- specie: Species from which genomic coordinates come from, options are: homo_sapiens, mus_musculus, drosophila_melanogaster (default homo_sapiens)
- otsu: Path to region bin topic otsu pickle (default <<workdir>>/scATAC/candidate_enhancers/region_bin_topics_otsu.pkl)<<workdir>>/scATAC/candidate_enhancers/region_bin_topics_otsu.pkl
- top3k: Path to region bin topic top3k pickle (default <<workdir>>/scATAC/candidate_enhancers/region_bin_topics_top3k.pkl)<<workdir>>/scATAC/candidate_enhancers/region_bin_topics_top3k.pkl
- markers: Path to marker dictionary pickle (default <<workdir>>/scATAC/candidate_enhancers/markers_dict.pkl)<<workdir>>/scATAC/candidate_enhancers/markers_dict.pkl
- scores_db: Path to score feather file (default <<workdir>>/data/hg38_screen_v10_clust.regions_vs_motifs.scores.feather)
- rank_db: Path to ranking feather file (default <<workdir>>/data/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather)
- motifs: Path to motif annotation table (default <<workdir>>/data/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl)
- motifs_version: Motif annotation version (default v10nr_clust)
- cpu: Number of cpu to use (default 24)
- overwrite: Recalculate all steps even if they completed sucessfully.

### Step4 Inferring enhancer-driven Gene Regulatory Networks using SCENICplus ###

This step identifies eGRN. An example call would be:

```
singularity exec \
-B /fast/:/fast/ \
-B /home:/home \
-e $G_PIPELINES/containers/scenicplus.sif \
/miniconda3/envs/scenicplus/bin/python3 -W ignore infer_enhancer-driven_gene.py \
-w $G_PIPELINES/scenicplus/test \
-tf $G_PIPELINES/scenicplus/test/data/utoronto_human_tfs_v_1.01.txt \
--tmp /fast/tmp \
--cpu 20 --sample '10x_pbmc'
```

The following parameters are mandatory:
- workdir: your working directory
- tf_file: Path to file containing genes symbols that are TFs
- sample: The sample id

These options are optional. **Please make sure default values are ok prior to running**
- tmp: Temp directory (default /tmp)
- specie: Species from which genomic coordinates come from, options are: homo_sapiens, mus_musculus, drosophila_melanogaster (default homo_sapiens)
- otsu: Path to region bin topic otsu pickle (default <<workdir>>/scATAC/candidate_enhancers/region_bin_topics_otsu.pkl)<<workdir>>/scATAC/candidate_enhancers/region_bin_topics_otsu.pkl
- top3k: Path to region bin topic top3k pickle (default <<workdir>>/scATAC/candidate_enhancers/region_bin_topics_top3k.pkl)<<workdir>>/scATAC/candidate_enhancers/region_bin_topics_top3k.pkl
- markers: Path to marker dictionary pickle (default <<workdir>>/scATAC/candidate_enhancers/markers_dict.pkl)<<workdir>>/scATAC/candidate_enhancers/markers_dict.pkl
- scores_db: Path to score feather file (default <<workdir>>/data/hg38_screen_v10_clust.regions_vs_motifs.scores.feather)
- rank_db: Path to ranking feather file (default <<workdir>>/data/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather)
- motifs: Path to motif annotation table (default <<workdir>>/data/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl)
- motifs_version: Motif annotation version (default v10nr_clust)
- cpu: Number of cpu to use (default 24)
- overwrite: Recalculate all steps even if they completed sucessfully.
