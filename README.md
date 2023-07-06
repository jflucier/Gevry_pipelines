# ILL Humann custom pipeline User Manual

----

## Contents ##

* [Requirements](#requirements)
* [Installation](#installation)
* [How to run](#how-to-run)
    * [Run scenicplus pipelines](#Run-scenicplus-pipelines)
      * [Step1: scRNA-seq preprocessing using Scanpy](Step1:-scRNA-seq-preprocessing-using-Scanpy)

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

### Step1: scRNA-seq preprocessing using Scanpy ###

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
- --cpu: Number of cpu to use (default 24)
- --mem: Max memery usage in Gigabyte (default 30)
- --qc_min_genes: Only keep cells with at least <<min_genes>> genes expressed (default )

**Notice** that preprocess script generates sample tsv file (i.e. precocess/preprocessed_reads.sample.tsv) that should be used
for the taxonomic profile and the functionnal profile pipeline.

Finally, the preprocess script can be executed on a single sample.
Use -h option to view usage:

```

$ bash $ILL_PIPELINES/scripts/preprocess.kneaddata.sh -h

Usage: preprocess.kneaddata.sh -s sample_name -o /path/to/out [--db] [--trimmomatic_options "trim options"] [--bowtie2_options "bowtie2 options"]
Options:

	-s STR	sample name
	-o STR	path to output dir
	-tmp STR	path to temp dir (default output_dir/temp)
	-t	# of threads (default 8)
	-m	memory (default 40G)
	-fq1	path to fastq1
	-fq2	path to fastq2
	--db	kneaddata database path (default /net/nfs-ip34/fast/def-ilafores/host_genomes/GRCh38_index/grch38_1kgmaj)
	--trimmomatic_options	options to pass to trimmomatic (default ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:30 MINLEN:100)
	--bowtie2_options	options to pass to trimmomatic (default --very-sensitive-local)

  -h --help	Display help


```

### Run Sourmash taxonomic abundance per sample ###

For full list of options:

```
$ bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.sourmash.sh  -h

Usage: generateslurm_taxonomic_abundance.sourmash.sh --sample_tsv /path/to/tsv --out /path/to/out [--SM_db /path/to/sourmash/db] [--SM_db_prefix sourmash_db_prefix] [--kmer kmer_size]
Options:

   --sample_tsv STR     path to sample tsv (3 columns: sample name<tab>fastq1 path<tab>fastq2 path)
   --out STR    path to output dir
   --SM_db sourmash databases directory path (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/sourmash_db/)
   --SM_db_prefix  sourmash database prefix, allowing wildcards (default genbank-2022.03)
   --kmer  choice of k-mer, dependent on database choices (default 51, make sure to have them available)

Slurm options:
   --slurm_alloc STR    slurm allocation (default def-ilafores)
   --slurm_log STR      slurm log file output directory (default to output_dir/logs)
   --slurm_email "your@email.com"       Slurm email setting
   --slurm_walltime STR slurm requested walltime (default 24:00:00)
   --slurm_threads INT  slurm requested number of threads (default 12)
   --slurm_mem STR      slurm requested memory (default 62G)

   -h --help    Display help

```

**Notice** that preprocess script generates sample tsv file needed here (i.e. precocess/preprocessed_reads.sample.tsv).

The sourmash taxonomic abundance script can also be executed on a single sample.
Use -h option to view usage:

```

$ bash $ILL_PIPELINES/scripts/taxonomic_abundance.sourmash.sh -h

Usage: taxonomic_abundance.sourmash.sh -s sample_name -o /path/to/out [-t threads] -fq1 /path/to/fastq1 -fq2 /path/to/fastq2 [--SM_db /path/to/sourmash/db] [--SM_db_prefix sourmash_db_prefix] [--kmer kmer_size]
Options:

        -s STR  sample name
        -o STR  path to output dir
        -tmp STR        path to temp dir (default output_dir/temp)
        -t      # of threads (default 8)
        -fq1    path to fastq1
        -fq2    path to fastq2
        --SM_db sourmash databases directory path (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/sourmash_db/)
        --SM_db_prefix  sourmash database prefix, allowing wildcards (default genbank-2022.03)
        --kmer  choice of k-mer size, dependent on available databases (default 51, make sure database is available)

  -h --help     Display help



```

### Run Metaphlan taxonomic abundance ###

For full list of options:

```
$ bash $ILL_PIPELINES/generateslurm_taxonomic_abundance.metaphlan.sh -h

Usage: generateslurm_taxonomic_abundance.metaphlan.sh --sample_tsv /path/to/tsv --out /path/to/out [--db /path/to/metaphlan/db]
Options:

   --sample_tsv STR     path to sample tsv (3 columns: sample name<tab>fastq1 path<tab>fastq2 path)
   --out STR    path to output dir
   --db   metaphlan db path (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212)

Slurm options:
   --slurm_alloc STR    slurm allocation (default def-ilafores)
   --slurm_log STR      slurm log file output directory (default to output_dir/logs)
   --slurm_email "your@email.com"       Slurm email setting
   --slurm_walltime STR slurm requested walltime (default 24:00:00)
   --slurm_threads INT  slurm requested number of threads (default 12)
   --slurm_mem STR      slurm requested memory (default 25G)

   -h --help    Display help
```

**Notice** that preprocess script generates sample tsv file needed here (i.e. precocess/preprocessed_reads.sample.tsv).

The metaphlan taxonomic abundance script can also be executed on a single sample.
Use -h option to view usage:

```
$ bash $ILL_PIPELINES/scripts/taxonomic_abundance.metaphlan.sh -h

Usage: taxonomic_abundance.metaphlan.sh -s sample_name -o /path/to/out [-db /path/to/metaphlan/db] -fq1 /path/to/fastq1 -fq2 /path/to/fastq2 [-fq1_single /path/to/single1.fastq] [-fq2_single /path/to/single2.fastq]
Options:

        -s STR  sample name
        -o STR  path to output dir
        -tmp STR        path to temp dir (default output_dir/temp)
        -t      # of threads (default 8)
        -fq1    path to fastq1
        -fq1_single     path to fastq1 unpaired reads
        -fq2    path to fastq2
        -fq2_single     path to fastq2 unpaired reads
        -db     metaphlan db path (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/metaphlan4_db/mpa_vOct22_CHOCOPhlAnSGB_202212)

  -h --help     Display help

```

Once metaphlan as run on all samples, you can merge results table by running the following script:

```
$ bash $ILL_PIPELINES/scripts/taxonomic_abundance.metaphlan.all.sh -h

Usage: taxonomic_abundance.metaphlan.all.sh -profiles /path/to/metaphlan_out/*_profile.txt -o /path/to/out
Options:

        -profiles Path to metaphlan outputs (i.e. /path/to/metaphlan_out/*_profile.txt)
        -o STR  path to output dir
        -tmp STR        path to temp dir (default output_dir/temp)

  -h --help     Display help

```

### Run Kraken2 taxonomic profile per sample ###

For full list of options:

```
$ bash $ILL_PIPELINES/generateslurm_taxonomic_profile.sample.sh -h

Usage: generateslurm_taxonomic_profile.sample.sh --sample_tsv /path/to/tsv --out /path/to/out [--kraken_db "kraken database"]
Options:

	--sample_tsv STR	path to sample tsv (3 columns: sample name<tab>fastq1 path<tab>fastq2 path)
	--out STR	path to output dir
	--kraken_db	kraken2 database path (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/kraken2_dbs/k2_pluspfp_16gb_20210517)
	--bracken_readlen	bracken read length option (default 150)

Slurm options:
	--slurm_alloc STR	slurm allocation (default def-ilafores)
	--slurm_log STR	slurm log file output directory (default to output_dir/logs)
	--slurm_email "your@email.com"	Slurm email setting
	--slurm_walltime STR	slurm requested walltime (default 6:00:00)
	--slurm_threads INT	slurm requested number of threads (default 24)
	--slurm_mem STR	slurm requested memory (default 125)

  -h --help	Display help

```

**Notice** that preprocess script generates sample tsv file needed here (i.e. precocess/preprocessed_reads.sample.tsv).

The taxonomic profile script can also be executed on a single sample.
Use -h option to view usage:

```

$ bash $ILL_PIPELINES/scripts/taxonomic_profile.sample.sh -h

Usage: taxonomic_profile.sample.sh [--kraken_db /path/to/krakendb] [--bracken_readlen int] [--confidence float] [-t thread_nbr] [-m mem_in_G] -fq1 /path/fastq1 -fq2 /path/fastq2 -o /path/to/out
Options:

	-s STR	sample name
	-o STR	path to output dir
	-tmp STR	path to temp dir (default output_dir/temp)
	-t	# of threads (default 8)
	-m	memory (default 40G)
	-fq1	path to fastq1
	-fq2	path to fastq2
	--kraken_db	kraken2 database path (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/kraken2_dbs/k2_pluspfp_16gb_20210517)
	--bracken_readlen	bracken read length option (default 150)
    --confidence	kraken confidence level to reduce false-positive rate (default 0.05)
    
  -h --help	Display help

```

### Taxonomy profile on all samples for all taxonomic level ###

For full list of options:

```
$ bash $ILL_PIPELINES/generateslurm_taxonomic_profile.allsamples.sh -h

Usage: generateslurm_taxonomic_profile.allsamples.sh --kreports 'kraken_report_regex' --out /path/to/out --bowtie_index_name idx_nbame
Options:

        --kreports STR  base path regex to retrieve species level kraken reports (i.e.: /home/def-ilafores/programs/ILL_pipelines/taxonomic_profile/*/*_bracken/*_bracken_S.kreport).
        --out STR       path to output dir
        --bowtie_index_name  name of the bowtie index that will be generated
        --chocophlan_db path to the full chocoplan db (default: /net/nfs-ip34/fast/def-ilafores/humann_dbs/chocophlan)

Slurm options:
        --slurm_alloc STR       slurm allocation (default def-ilafores)
        --slurm_log STR slurm log file output directory (default to output_dir/logs)
        --slurm_email "your@email.com"  Slurm email setting
        --slurm_walltime STR    slurm requested walltime (default 24:00:00)
        --slurm_threads INT     slurm requested number of threads (default 48)
        --slurm_mem STR slurm requested memory (default 251G)

  -h --help     Display help


```

The kreports parameter is a regular expression that points to all kraken report generated at specie level.
The analysis begins and creates the buglist and then creates the bowtie index on the buglist. It finishes by 
generating the taxonomic table for each taxonomic level.

### Generate HUMAnN bugs list ###


For full list of options:

```
$ bash $ILL_PIPELINES/generateslurm_taxonomic_profile.allsamples.sh -h

Usage: generateslurm_taxonomic_profile.allsamples.sh --kreports 'kraken_report_regex' --out /path/to/out --bowtie_index_name idx_nbame
Options:

	--kreports STR	base path regex to retrieve species level kraken reports (i.e.: /home/jflucier/tmp/taxonomic_profile/*/*_bracken/*_bracken_S.kreport).
	--out STR	path to output dir
	--bowtie_index_name  name of the bowtie index that will be generated
	--chocophlan_db	path to the full chocoplan db (default: /net/nfs-ip34/fast/def-ilafores/humann_dbs/chocophlan)

Slurm options:
	--slurm_alloc STR	slurm allocation (default def-ilafores)
	--slurm_log STR	slurm log file output directory (default to output_dir/logs)
	--slurm_email "your@email.com"	Slurm email setting
	--slurm_walltime STR	slurm requested walltime (default 24:00:00)
	--slurm_threads INT	slurm requested number of threads (default 48)
	--slurm_mem STR	slurm requested memory (default 251G)

  -h --help	Display help

```

The kreports parameter is a regular expression that points to all species level kraken report files that will be used in analysis.

If you wish, the humann bug list generation can be directly runned locally on serve using similar code as below:

```
kreports="/net/nfs-ip34/home/def-ilafores//projet_PROVID19/taxKB_conf01_jfl/*/*_bracken/*_S.kreport"
out=test
tmp=test/temp
threads=24
bowtie_idx_name=my_bt_idx
choco_db=/net/nfs-ip34/fast/def-ilafores/humann_dbs/chocophlan

bash $ILL_PIPELINES//scripts/taxonomic_profile.allsamples.sh \
--kreports "$kreports" \
--out ${out} \
--tmp $tmp \
--threads ${threads} \
--bowtie_index_name $bowtie_idx_name \
--chocophlan_db $choco_db

```

### Run HUMAnN functionnal profile ###

Before running this pipeline, make sure [HUMAnN](https://huttenhower.sph.harvard.edu/humann/) environment is acessible.

For full list of options:

```
$ bash $ILL_PIPELINES/generateslurm_functionnal_profile.humann.sh -h

Usage: generateslurm_functionnal_profile.humann.sh --sample_tsv /path/to/tsv --out /path/to/out --nt_db "nt database path" [--search_mode "search mode"] [--prot_db "protein database path"]
Options:

  --sample_tsv STR      path to sample tsv (5 columns: sample name<tab>fastq1 path<tab>fastq2 path<tab>fastq1 single path<tab>fastq2 single path). Generated in preprocess step.
        --out STR       path to output dir
        --search_mode   Search mode. Possible values are: dual, nt, prot (default prot)
        --nt_db the nucleotide database to use (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/humann_dbs/chocophlan)
        --prot_db       the protein database to use (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/humann_dbs/uniref)
        --utility_map_db        the protein database to use (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/humann_dbs/utility_mapping)

Slurm options:
        --slurm_alloc STR       slurm allocation (default def-ilafores)
        --slurm_log STR slurm log file output directory (default to output_dir/logs)
        --slurm_email "your@email.com"  Slurm email setting
        --slurm_walltime STR    slurm requested walltime (default 24:00:00)
        --slurm_threads INT     slurm requested number of threads (default 24)
        --slurm_mem STR slurm requested memory (default 30G)

  -h --help     Display help



```

The sample_tsv that can be used was created in the preprocess step (i.e. precocess/preprocessed_reads.sample.tsv).

The functionnal profile script can also be executed on a single sample.
Use -h option to view usage:

```

$ bash $ILL_PIPELINES/scripts/functionnal_profile.humann.sh -h

Usage: functionnal_profile.humann.sh -s sample_name -o /path/to/out --nt_db "nt database path" [--search_mode "search mode"] [--prot_db "protein database path"]
Options:

        -s STR  sample name
        -o STR  path to output dir
        -tmp STR        path to temp dir (default output_dir/temp)
        -t      # of threads (default 8)
        -fq1    path to fastq1
        -fq1_single     path to fastq1 unpaired reads
        -fq2    path to fastq2
        -fq2_single     path to fastq2 unpaired reads
        --search_mode   Search mode. Possible values are: dual, nt, prot (default prot)
        --nt_db the nucleotide database to use (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/humann_dbs/chocophlan)
        --prot_db       the protein database to use (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/humann_dbs/uniref)
        --utility_map_db        the protein database to use (default /cvmfs/datahub.genap.ca/vhost34/def-ilafores/humann_dbs/utility_mapping)

  -h --help     Display help


```

### Run MetaWRAP assembly, binning and bin refinement ###

Before running this pipeline, make sure singularity and BBmap executables are in your path. On ip29, just do the following:

```
module load singularity mugqic/BBMap/38.90

```

For full list of options:

```

$ bash ${ILL_PIPELINES}/generateslurm_assembly_bin_refinement.metawrap.sh -h

Usage: generateslurm_denovo_assembly_bin_refinement.metawrap.sh --sample_tsv /path/to/tsv --out /path/to/out [--assembly] [--binning] [--refinement]
Options:

    --sample_tsv STR	path to sample tsv (3 columns: sample name<tab>fastq1 path<tab>fastq2 path)
	--out STR	path to output dir
	--assembly	perform assembly
	--binning	perform binning step
	--refinement	perform refinement step

Metawrap options:
	--metaspades	use metaspades for assembly (default: true)
	--megahit	use megahit for assembly (default: true)
	--metabat2	use metabat2 for binning (default: true)
	--maxbin2	use maxbin2 for binning (default: true)
	--concoct	use concoct for binning (default: true)
	--run-checkm	run checkm for binning (default: true)
	--refinement_min_compl INT	refinement bin minimum completion percent (default 50)
	--refinement_max_cont INT	refinement bin maximum contamination percent (default 10)

Slurm options:
	--slurm_alloc STR	slurm allocation (default def-ilafores)
	--slurm_log STR	slurm log file output directory (default to output_dir/logs)
	--slurm_email "your@email.com"	Slurm email setting
	--slurm_walltime STR	slurm requested walltime (default 24:00:00)
	--slurm_threads INT	slurm requested number of threads (default 48)
	--slurm_mem STR	slurm requested memory (default 251G)

  -h --help	Display help


```

Most default values should be ok in a cluster environment. Make sure you specify sample_tsv, ouput path and steps you wich to execute (assembly and/or binning and/or refinement). Obviously, before running binning, you must perform assembly step.

The sample_tsv that can be used was created in the preprocess step (i.e. precocess/preprocessed_reads.sample.tsv).

Here a re some example commands you can perform for this pipeline:
```
# on ip29, load singularity and bbmap in path
module load singularity mugqic/BBMap/38.90

# Run all steps with defualt parameters
$ bash ${ILL_PIPELINES}/generateslurm_assembly_bin_refinement.metawrap.sh \
--out path/to/out --sample_tsv /path/to/tsv

# Run only the assembly step
$ bash ${ILL_PIPELINES}/generateslurm_assembly_bin_refinement.metawrap.sh \
--out path/to/out --sample_tsv /path/to/tsv \
--assembly

# Run only the assembly step using only megahit assembler
$ bash ${ILL_PIPELINES}/generateslurm_assembly_bin_refinement.metawrap.sh \
--out path/to/out --sample_tsv /path/to/tsv \
--assembly --megahit

# Run the assembly step using only megahit assembler and binning step with
# concoct and maxbin binner software
$ bash ${ILL_PIPELINES}/generateslurm_assembly_bin_refinement.metawrap.sh \
--out path/to/out --sample_tsv /path/to/tsv \
--assembly --megahit \
--binning --maxbin2 --concoct

# Run assembly and binning with default paramters and refinement step using
# specific bin completion and contamination values
$ bash ${ILL_PIPELINES}/generateslurm_assembly_bin_refinement.metawrap.sh \
--out path/to/out --sample_tsv /path/to/tsv \
--assembly --binning \
--refinement --refinement_min_compl 90 --refinement_max_cont 5

```

Finally, the assembly, binning and refinement script can be executed on a single sample.
Use -h option to view usage:

```
# on ip29, load singularity and bbmap in path
module load singularity mugqic/BBMap/38.90


## assembly script usage:
$ bash /home/jflucier/localhost/projet/ILL_pipelines/scripts/assembly.metawrap.sh -h

Usage: assembly.metawrap.sh [-tmp /path/tmp] [-t threads] [-m memory] [--metaspades] [--megahit] -s sample_name -o /path/to/out -fq1 /path/to/fastq1 -fq2 /path/to/fastq2
Options:

	-s STR	sample name
	-o STR	path to output dir
	-tmp STR	path to temp dir (default output_dir/temp)
	-t	# of threads (default 8)
	-m	memory (default 40G)
	-fq1	path to fastq1
	-fq2	path to fastq2
	--metaspades	use metaspades for assembly (default: true)
	--megahit	use megahit for assembly (default: true)

  -h --help	Display help

## Binning script usage:
$ bash $ILL_PIPELINES/scripts/binning.metawrap.sh -h

Usage: binning.metawrap.sh [-tmp /path/tmp] [-t threads] [-m memory] [--metabat2] [--maxbin2] [--concoct] [--run-checkm] -s sample_name -o /path/to/out -a /path/to/assembly -fq1 /path/to/fastq1 -fq2 /path/to/fastq2
Options:

	-s STR	sample name
	-o STR	path to output dir
	-tmp STR	path to temp dir (default output_dir/temp)
	-t	# of threads (default 8)
	-m	memory (default 40G)
	-a	assembly fasta filepath
	-fq1	path to fastq1
	-fq2	path to fastq2
	--metabat2	use metabat2 for binning (default: true)
	--maxbin2	use maxbin2 for binning (default: true)
	--concoct	use concoct for binning (default: true)
	--run-checkm	run checkm on bins (default: true)

  -h --help	Display help


## Binning refinement script usage:
$ $ILL_PIPELINES/scripts/bin_refinement.metawrap.sh -h

Usage: bin_refinement.metawrap.sh [-tmp /path/tmp] [-t threads] [-m memory] [--metaspades] [--megahit] -s sample_name -o /path/to/out -fq1 /path/to/fastq1 -fq2 /path/to/fastq2
Options:

	-s STR	sample name
	-o STR	path to output dir
	-tmp STR	path to temp dir (default output_dir/temp)
	-t	# of threads (default 8)
	-m	memory (default 40G)
	--metabat2_bins	path to metabats bin direcotry
	--maxbin2_bins	path to maxbin2 bin direcotry
	--concoct_bins	path to concoct bin direcotry
	--refinement_min_compl INT	refinement bin minimum completion percent (default 50)
	--refinement_max_cont INT	refinement bin maximum contamination percent (default 10)

  -h --help	Display help

```
