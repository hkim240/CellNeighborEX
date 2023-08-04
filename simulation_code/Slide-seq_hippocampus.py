import CellNeighborEX.categorization
import CellNeighborEX.DEanalysis
import CellNeighborEX.visualization
import pandas as pd


# To run CellNeighborEX, information on cell type annotation, spatial location, and expression values is required.
# For preparation of input data, please refer to simulation_data on the github page.
# The CellNeighborEX python scripts provide the description of parameters and returns for each function: 
# https://github.com/hkim240/CellContact/tree/main/CellNeighborEX

#### (1) Import spatial transcriptomics (ST) data
# Slide-seq data: mouse hippocampus
# you can download 'Slide-seq_hippocampus.csv' and 'hippocampus_abbrev.csv' in simulation_data/annotated_data on the github page.
# 'Slide-seq_hippocampus.csv' was produced based on the RCTD (deconvolution tool) result.
path = '/Users/kimh15/annotated_data/'
df_processed = pd.read_csv(path + 'Slide-seq_hippocampus.csv', header=0)

# add abbreviation of cell typs to df_processed.
df_abbrev = pd.read_csv(path + 'hippocampus_abbrev.csv', header=0)
for i in range(len(df_processed)):
    
    print(i)
    ct1 = df_processed['first_type'][i]
    ct2 = df_processed['second_type'][i]
    for j in range(len(df_abbrev)):
        
        if ct1 == df_abbrev['Cluster'][j]:
            df_processed['celltype1'][i] = df_abbrev['Abbrev'][j]
            
        if ct2 == df_abbrev['Cluster'][j]:
            df_processed['celltype2'][i] = df_abbrev['Abbrev'][j]


#### (2) Categorize cells into heterotypic neighbors and homotypic neighbors
# CellNeighborEX.categorization.generate_input_files: all categorzied input files are saved in the categorized_data folder of the root directory.
CellNeighborEX.categorization.generate_input_files(data_type = "NGS", df = df_processed, sample_size=30, min_sample_size=1)


#### (3) Perform neighbor-dependent gene expression analysis
# set paths of input data: categorized data files and expression data.
# you can download 'cell_id.txt', 'gene_name.txt', and 'log_data.txt' in simulation_data/expression_data/Slide-seq_liver_cancer on the github page.
# cell_id.txt: cell or spot barcodes
# gene_name.txt: genes of interest
# log_data.txt: log-normalized data
path_lognormalized_data = '/Users/kimh15/expression_data/Slide-seq_hippocampus/'
df_cell_id = pd.read_csv(path_lognormalized_data + "cell_id.txt", delimiter="\t", header=None) # The length of df_processed must be the same as the length of df_cell_id!
df_gene_name = pd.read_csv(path_lognormalized_data + "gene_name.txt", delimiter="\t", header=None)
df_log_data = pd.read_csv(path_lognormalized_data + "log_data.txt", delimiter="\t", header=None)

path_categorization = '/Users/kimh15/categorized_data/'  

# set argument values for CellNeighborEX.DEanalysis.analyze_data().
data_type = "NGS" # Image: image-based ST data, NGS: NGS-based ST data
lrCutoff = 0.4 # log ratio 
pCutoff = 0.01 # p-value 
pCutoff2 = 0.01 # false discovery rate
direction = 'up' # up: up-reguated genes, down: down-regulated genes
normality_test = False # True: depending on the result of the normality test, the statistical test is determined. If the data is normal, the parametric test is used. Otherwise, the non-parametric test is used.
                       # False: when sample size (number of cells/spots) is larger than 30, the parameteric test is used. Otherwise, the non-parametric test is used.
top_genes = 10 # Top 10 DEGs are annotated in the volcano plot.

DEG_list = CellNeighborEX.DEanalysis.analyze_data(df_cell_id, df_gene_name, df_log_data, path_categorization, data_type, lrCutoff, pCutoff, pCutoff2, direction, normality_test, top_genes, save=True)


#### (4) Visualize the neighbor-dependent gene expression for spatial validation
# select a cell type and DEG for spatial visualization and load the data.
# for example, Fabp7 is one of up-regulated genes identified from the heterotypic spots of EnT+A.
path_selected = '/Users/kimh15/DE_results/EnT+A/'
column_names = ['barcode', 'logdata', 'zscore']
heterotypic = pd.read_csv(path_selected + "EnT+A_Fabp7.txt", delimiter=",", names = column_names)
homotypic1 = pd.read_csv(path_selected + "EnT+EnT_Fabp7.txt", delimiter=",", names = column_names)
homotypic2 = pd.read_csv(path_selected + "A+A_Fabp7.txt", delimiter=",", names = column_names)
heterotypic['type'] = 'EnT+A'
homotypic1['type'] = 'EnT+EnT'
homotypic2['type'] = 'A+A'
df_exp = pd.concat([heterotypic, homotypic1, homotypic2])

# set parameter values.
# Slide-seq: beadsize_bg=10, beadsize_red=600, beadsize_blue=200, beadsize_black=200
df_bg, df_red, df_blue, df_black = CellNeighborEX.visualization.set_parameters(df_processed, df_exp, beadsize_bg=10, edgecolor_bg=(0.85,0.85,0.85), beadcolor_bg=(0.85,0.85,0.85), beadsize_red=600, beadsize_blue=200, beadsize_black=200, type_red='EnT+A', type_blue='EnT+EnT', type_black='A+A') 

# get spatial map.
#Slide-seq:(28,28)
# zorder_red, zorder_blue, and zorder_black are parameters that determine the drawing order in the spatial map.
# CellNeighborEX.visualization.get_spatialPlot (save=True): The spatial map is saved in the spatialMap folder of the root directory.
CellNeighborEX.visualization.get_spatialPlot(df_bg, df_red, df_blue, df_black, label_red='Endothelial_Tip+Astrocyte', label_blue='Endothelial_Tip', label_black='Astrocyte', label_gene='Fabp7', zorder_red=4.0, zorder_blue=3.0, zorder_black=2.0, figsize=(28,28), save=True)