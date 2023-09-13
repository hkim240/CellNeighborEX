.. _api:

API
===
Import CellNeighborEX::

    import CellNeighborEX

Immediate neighbor detection
~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: CellNeighborEX.neighbors
.. currentmodule:: CellNeighborEX

.. autosummary::
    :toctree: api

    neighbors.create_dataframe
    neighbors.detect_neighbors
    neighbors.calculate_closest_distance
    neighbors.get_neighbors
    neighbors.process_dataframe
    
Cell type categorization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: CellNeighborEX.categorization
.. currentmodule:: CellNeighborEX

.. autosummary::
    :toctree: api

    categorization.generate_input_files

Neighbor-dependent gene expression analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: CellNeighborEX.DEanalysis
.. currentmodule:: CellNeighborEX

.. autosummary::
    :toctree: api

    DEanalysis.two_sample_f_test
    DEanalysis.create_nullmodel
    DEanalysis.find_nullDEGs
    DEanalysis.find_contactDEGs
    DEanalysis.get_heatmap
    DEanalysis.get_volcano_plot
    DEanalysis.delete_files_with_keyword
    DEanalysis.analyze_data

Spatial visualization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: CellNeighborEX.visualization
.. currentmodule:: CellNeighborEX

.. autosummary::
    :toctree: api

    visualization.import_expdata
    visualization.set_parameters
    visualization.get_spatialPlot