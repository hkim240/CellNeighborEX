import CellNeighborEX.neighbors
import CellNeighborEX.categorization
import CellNeighborEX.DEanalysis
import CellNeighborEX.visualization
import squidpy as sq
import pandas as pd


#### (1) Import spatial transcriptomics (ST) data
# Image-based data: seqFISH
adata = sq.datasets.seqfish()  # 19416 cells x 351 genes


#### (2) Define immediate neighbors for cell contact & Categorize cells into heterotypic neighbors and homotypic neighbors
df =  CellNeighborEX.neighbors.create_dataframe(adata, coord_key='spatial', celltype_key='celltype_mapped_refined')

# The closest_distances function is optional. If it is needed to caculate closest distances, this function can be used.
#closest_distances = CellNeighborEX.neighbors.calculate_closest_distance(df)

matrix = CellNeighborEX.neighbors.detect_neighbors(adata, coord_key='spatial', type='generic', knn=None, radius_value=None, delaunay=True)
neiNum = CellNeighborEX.neighbors.get_neighbors(matrix)

# CellNeighborEX.neighbors.process_dataframe (save=True): df_processed is saved in the neighbor_info folder of the root directory
df_processed = CellNeighborEX.neighbors.process_dataframe(df, matrix, neiNum, save=True) 

# CellNeighborEX.categorization.generate_input_files: all categorzied input files are saved in the categorized_data folder of the root directory
CellNeighborEX.categorization.generate_input_files(data_type = "Image", df = df_processed, sample_size=30, min_sample_size=1)


#### (3) Perform neighbor-dependent gene expression analysis
# set paths of input data: categorized data files and expression data
# you can download 'cell_id.txt', 'gene_name.txt', and 'log_data.txt' from the github page.
path_lognormalized_data = '/Users/kimh15/expression_data/seqFISH_embryo/'
df_cell_id = pd.read_csv(path_lognormalized_data + "cell_id.txt", delimiter="\t", header=None)
df_gene_name = pd.read_csv(path_lognormalized_data + "gene_name.txt", delimiter="\t", header=None)
df_log_data = pd.read_csv(path_lognormalized_data + "log_data.txt", delimiter=",", header=0)

path_categorization = '/Users/kimh15/categorized_data/'  

# remove categorized data files with a specific keyword: this function can be used to exclude unwanted cell types.
specific_keyword = 'Low-quality' # In the seqFISH data, 'Low-quality' is a unwanted cell type.
CellNeighborEX.DEanalysis.delete_files_with_keyword(path_categorization, specific_keyword)

# set argument values for CellNeighborEX.DEanalysis.analyze_data()
data_type = "Image"  # Image: image-based ST data, NGS: NGS-based ST data
lrCutoff = 0.4 # log ratio 
pCutoff = 0.01 # p-value 
pCutoff2 = 0.05 # false discovery rate
direction = 'up' # up: up-reguated genes, down: down-regulated genes
normality_test = False # True: the statistical test for DE analysis is determined between parameteric and non-parametric. 
                       # False: when sample size (number of cells/spots) is larger than 30, the parameteric tests are used.
top_genes = 10 # For only top 10 DEGs, gene names are annotated in the volcano plot.

DEG_list = CellNeighborEX.DEanalysis.analyze_data(df_cell_id, df_gene_name, df_log_data, path_categorization, data_type, lrCutoff, pCutoff, pCutoff2, direction, normality_test, top_genes, save=True)


#### (4) Visualize the neighbor-dependent gene expression for spatial validation
# select a cell type and DEG for spatial visualization and load the data
path_selected = '/Users/kimh15/DE_results/Gut-tube+Neural-crest/'
column_names = ['barcode', 'logdata', 'zscore']
heterotypic = pd.read_csv(path_selected + "Gut-tube+Neural-crest_Pitx1.txt", delimiter=",", names = column_names)
homotypic1 = pd.read_csv(path_selected + "Gut-tube+Gut-tube_Pitx1.txt", delimiter=",", names = column_names)
homotypic2 = pd.read_csv(path_selected + "Neural-crest+Neural-crest_Pitx1.txt", delimiter=",", names = column_names)
heterotypic['type'] = 'Gut-tube+Neural-crest'
homotypic1['type'] = 'Gut-tube+Gut-tube'
homotypic2['type'] = 'Neural-crest'
df_exp = pd.concat([heterotypic, homotypic1, homotypic2])

# set parameter values
# seqFISH: beadsize_bg=5, beadsize_red=60, beadsize_blue=30, beadsize_black=30 
df_bg, df_red, df_blue, df_black = CellNeighborEX.visualization.set_parameters(df_processed, df_exp, beadsize_bg=5, edgecolor_bg=(0.85,0.85,0.85), beadcolor_bg=(0.85,0.85,0.85), beadsize_red=60, beadsize_blue=30, beadsize_black=30, type_red='Gut-tube+Neural-crest', type_blue='Gut-tube+Gut-tube', type_black='Neural-crest') 

# get spatial map
# seqFISH:(20,28)
CellNeighborEX.visualization.get_spatialPlot(df_bg, df_red, df_blue, df_black, label_red='Gut-tube+Neural-crest', label_blue='Gut-tube+Gut-tube', label_black='Neural-crest', label_gene='Pitx1', figsize=(20,28), save=True)