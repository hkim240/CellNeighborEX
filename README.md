# CellNeighborEX
<p align="justify">CellNeighborEX is a computational approach to identify genes up-regulated or down-regulated by immediate neighbors from spatial transcriptomics (ST) data. It works for both image-based and NGS-based ST data. For image-based ST data where exact cell locations are available, CellNeighborEX uses various algorithms including Delaunay triangulation and KNN to find immediate neighbors. For NGS-based ST data where exact cell locations are not available, CellNeighborEX leverages the mixture of transcriptomes in each spot. CellNeighborEX dissects cells or spots based on the cell types of the immediate neighbors. Carrying out differential expression analysis for the categorized cells or spots, CellNeighborEX detects neighbor-dependent genes. The expression of neighbor-dependent genes is validated in the spatial context.</p> 

The figure below shows the workflow of CellNeighborEX:

![Fig 1](https://user-images.githubusercontent.com/99720939/229945240-2c9a2ef9-2566-496f-9981-0823cd95b813.png)

# Install
Install CellNeighborEX via PyPI by running:
<code>pip install CellNeighborEX</code> 

# Citation
Hyobin Kim, Amit Kumar, Cecilia Lövkvist, António M. Palma, Patrick Martin, Junil Kim, Praveen Bhoopathi, Jose Trevino, Paul Fisher, Esha Madan, Rajan Gogna, and Kyoung Jae Won, CellNeighborEX: Deciphering Neighbor-Dependent Gene Expression from Spatial Transcriptomics Data, bioRxiv.