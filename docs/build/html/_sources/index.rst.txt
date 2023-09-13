|PyPI|

CellNeighborEX 
==========================================================

CellNeighborEX is a computational approach to identify genes up-regulated or down-regulated by immediate neighbors from spatial transcriptomics (ST) data at sinlge cell or near cellular resolution. It works for both image-based and NGS-based ST data. For image-based ST data where exact cell locations are available, CellNeighborEX uses various algorithms including Delaunay triangulation and KNN to find immediate neighbors. For NGS-based ST data where exact cell locations are not available, CellNeighborEX leverages the mixture of transcriptomes in each spot. CellNeighborEX dissects cells or spots based on the cell types of the immediate neighbors. Carrying out differential expression analysis for the categorized cells or spots, CellNeighborEX detects neighbor-dependent genes. The expression of neighbor-dependent genes is validated in the spatial context.

The figure below shows the workflow of CellNeighborEX:

.. image:: https://user-images.githubusercontent.com/99720939/229945240-2c9a2ef9-2566-496f-9981-0823cd95b813.png
    :alt: CellNeighborEX title figure
    :align: center
    :target: https://github.com/hkim240/CellNeighborEX

.. toctree::
   :maxdepth: 2
   :hidden:
   
   installation
   api
   tutorials

.. |PyPI| image:: https://img.shields.io/pypi/v/CellNeighborEX.svg
    :target: https://pypi.org/project/CellNeighborEX/
    :alt: PyPI  