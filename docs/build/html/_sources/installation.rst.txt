.. _installation:

Installation
============
CellNeighborEX requires Python version >=3.8, <3.11. 
We recommend using conda environment to avoid dependency conflicts. 
The dependencies are listed in requirements.txt, which can be downloaded from https://github.com/hkim240/CellNeighborEX. 
Cloning the CellNeighborEX repository from GitHub is also available. 

Create conda environment "myenv"::
    
    conda create -n myenv python=3.10
    conda activate myenv

Navigate into the directory where requirements.txt is located. Then, install dependencies:: 

    pip install -r requirements.txt

Install CellNeighborEX from PyPI::

    pip install CellNeighborEX