# CellNeighborEX
<p align="justify">CellNeighborEX is a computational approach to identify genes up-regulated or down-regulated by immediate neighbors from spatial transcriptomics (ST) data. It works for both image-based and NGS-based ST data. For image-based ST data where exact cell locations are available, CellNeighborEX uses various algorithms including Delaunay triangulation and KNN to find immediate neighbors. For NGS-based ST data where exact cell locations are not available, CellNeighborEX leverages the mixture of transcriptomes in each spot. CellNeighborEX dissects cells or spots based on the cell types of the immediate neighbors. Carrying out differential expression analysis for the categorized cells or spots, CellNeighborEX detects neighbor-dependent genes. The expression of neighbor-dependent genes is validated in the spatial context.</p> 

The figure below shows the workflow of CellNeighborEX:

![Fig 1](https://user-images.githubusercontent.com/99720939/229945240-2c9a2ef9-2566-496f-9981-0823cd95b813.png)

# Running CellNeighborEX
CellNeighborEX was implemented in Python and Matlab. The scripts and input data files have been uploaded to the following folders:

(1) <code>Python/seqFISH_neighbor</code> 

- get the information of neighboring cells in the mouse embryo seqFISH data (output: seqFISH_neiInfo.csv)
- extract cells that have the same cell type or one different cell type as the cell type of immediate neighbors)  (output: seqFISH_embryo.csv)

(2) <code>Python/categorization</code>

- categorize cells or spots based on the cell types of immediate neighbors (output: index_.csv, matchComb_.csv, neiCombUnique_.csv, prop_.csv)
- For seqFISH data, seqFISH_embryo.csv obtained by (1) was used as input.
- For Slide-seq data, result files obtained from deconvolution tool RCTD were used as input.

(3) <code>Matlab/input_logdata</code>

- generate input log data (example of output: mouseLiver_Slideseq_RCTD_top2000_plus.mat)

(4) <code>Python/input_markers</code>

- generate input marker data (example of output: liver_matchMarker_plus_selected.csv)

(5) <code>Matlab/input_category</code>

- generate input categorization data from four pairs of csv files obtained by (2)(example of output: HepatocyteI+HepaticStellateCell.mat)

(6) <code>Matlab/DE_analysis</code>

- find neighbor-dependent genes from input files obtained by (3),(4) and (5) (output: DEG_lists, expression data, heat maps)

(7) <code>Python/visualization</code>

- get spatial mapping plots from expression data obtained by (6) (output: neighbor-dependent gene.png)

To reduce the file size, all the output files for each folder have been stored in Results.zip.

# Citation
Hyobin Kim, Cecilia Lövkvist, Patrick C.N. Martin, Junil Kim, Kyoung Jae Won, Detecting Cell Contact-dependent Gene Expression from Spatial Transcriptomics Data, bioRxiv, 2022.
