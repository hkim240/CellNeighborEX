# CellNeighborEX
<p align="justify">CellNeighborEX is a computational approach to identify genes up-regulated or down-regulated by immediate neighbors from spatial transcriptomics (ST) data at single cell or near cellular resolution. It works for both image-based and NGS-based ST data. For image-based ST data where exact cell locations are available, CellNeighborEX uses various algorithms including Delaunay triangulation and KNN to find immediate neighbors. For NGS-based ST data where exact cell locations are not available, CellNeighborEX leverages the mixture of transcriptomes in each spot. CellNeighborEX dissects cells or spots based on the cell types of the immediate neighbors. Carrying out differential expression analysis for the categorized cells or spots, CellNeighborEX detects neighbor-dependent genes. The expression of neighbor-dependent genes is validated in the spatial context.</p> 

The figure below shows the workflow of CellNeighborEX:

![Fig 1](https://user-images.githubusercontent.com/99720939/229945240-2c9a2ef9-2566-496f-9981-0823cd95b813.png)

# Installation
<p align="justify">CellNeighborEX requires Python version >=3.8, <3.11. We recommend using conda environment to avoid dependency conflicts. The dependencies are listed in requirements.txt.</p> 

Create conda environment “myenv”:
<code><p>conda create -n myenv python=3.10</p>
<p>conda activate myenv</p></code> 

Navigate into the directory where requirements.txt is located. Then, install dependencies:
<code>pip install -r requirements.txt</code> 

Install CellNeighborEX from PyPI:
<code>pip install CellNeighborEX</code> 

# Python API documentation and tutorials
Please see this [Read the Docs page](https://CellNeighborEX.readthedocs.io/en/latest/).