import os
import pickle
import matplotlib.pyplot as plt
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget
import dill
from IPython.core.display import HTML

work_dir = '/home/def-gevrynic/programs/scenicplus/pbmc_tutorial'
tmp_dir = '/tmp'

fig=plt.figure(figsize=(10,10))

### Cistarget databases

## custom db creation
# You can choose to compute this database yourself by scoring the consensus peaks generated 
# in the scATAC-seq analysis using a set of motifs. The advantage of creating a sample 
# specific database is that you can potentially pick up more target regions, given that only 
# regions included/overlappig with regions in the cistarget database will be used 
# for the SCENIC+ analysis.
# https://github.com/aertslab/create_cisTarget_databases

##to TEST!!!##

# use precomputed db
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))

# Convert to dictionary of pyranges objects.
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

# Define rankings, score and motif annotation database.
rankings_db = os.path.join(work_dir, 'data/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(work_dir, 'data/hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
#motif_annotation = os.path.join(work_dir, 'data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')
motif_annotation = os.path.join(work_dir, 'data/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')

# run cistarget based and DEM based motif enrichment analysis with or without promoter regions.
if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))

run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 25,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust',
    )

# explore some of the results. Below we show the motifs found for topic 8 (specific to B-cells) using DEM.
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))
menr['DEM_topics_otsu_All'].DEM_results('Topic8')
with open(os.path.join(work_dir, 'fig10..html'), 'w') as f:
    f.write(menr['DEM_topics_otsu_All'].DEM_results('Topic8').data)

