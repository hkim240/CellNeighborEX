import squidpy as sq
import pandas as pd
from collections import Counter
import anndata 
from scipy.sparse import csr_matrix
import numpy as np
import os
import math
from tqdm import tqdm
import matplotlib.pyplot as plt


def create_dataframe(adata:anndata, coord_key:str, celltype_key:str):
    """
    Create a DataFrame from the provided AnnData object.

    Parameters:
        adata (AnnData): AnnData object containing the dataset.
        celltype_key (str): Key to access the cell type information in `adata.obs`.
        coord_key (str): Key to access the spatial coordinates in `adata.obsm`.

    Returns:
        df (DataFrame): Constructed DataFrame with barcode, cell type, coordinates, and additional columns.
    """
    # Create a DataFrame for barcode
    df_barcode = pd.DataFrame(adata.obs.index, columns=["barcode"])

    # Create a DataFrame for cell type
    df_celltype = pd.DataFrame(adata.obs[celltype_key])
    df_celltype.index = df_barcode.index
    df_celltype = df_celltype.rename(columns={celltype_key: "celltype1"})

    # Extract x and y values from spatial coordinates
    xVal_list = [coord[0] for coord in adata.obsm[coord_key]]
    yVal_list = [coord[1] for coord in adata.obsm[coord_key]]
    df_coord = pd.DataFrame(list(zip(xVal_list, yVal_list)), columns=["x", "y"])

    # Concatenate the barcode, cell type, and coordinate DataFrames
    df = pd.concat([df_barcode, df_celltype, df_coord], axis=1)

    # Add additional columns
    df['first_type'] = -1
    df['second_type'] = -1
    df['celltype2'] = 'NA'
    df['prop1'] = 1.0
    df['prop2'] = 0.0

    return df


def detect_neighbors(adata:anndata, coord_key:str, type:str, knn:int, radius_value:float, delaunay:bool):
     
    """
    Detect the spatial neighbors using squidpy and retrieve the spatial connectivity matrix.

    Parameters:
        adata (AnnData): Loaded dataset.
        coord_key (str): Key to access the spatial coordinates in `adata.obsm`.
        type (str): Type of spatial coordinates, either 'generic' or 'grid'.
        knn (int): Number of nearest neighbors to consider for connectivity.
        radius_value (float): Radius to define the spatial connectivity.
        delaunay (bool): Flag indicating whether to use Delaunay triangulation for connectivity.

    Returns:
        matrix (sparse matrix): Spatial connectivity matrix.
    """
    # Detect neighbors via Delaunay triangulation, KNN, or radius
    if type == 'generic' and delaunay == True:
        sq.gr.spatial_neighbors(adata, spatial_key=coord_key, coord_type=type, delaunay=True)
    elif knn != 0 and radius_value is None and delaunay is False:
        sq.gr.spatial_neighbors(adata, spatial_key=coord_key, coord_type=type, n_neighs=knn, delaunay=False)
    elif type == 'generic' and radius_value != 0.0 and knn is None and delaunay is False:
        sq.gr.spatial_neighbors(adata, spatial_key=coord_key, coord_type=type, radius=radius_value, delaunay=False)
    else:
        print("Error: Please check the arguments of CellNeighborEX.neighbors.detect_neighbors.")

    matrix = adata.obsp["spatial_connectivities"]

    return matrix


def calculate_closest_distance(df:pd.DataFrame):

    # Calcuate minimum distance between cells to determine the radius value
    points = []
    for i in range(len(df)):
        x = df['x'][i]
        y = df['y'][i]
        points.append((x, y))
    
    closest_distances = []
    for i in tqdm(range(len(points)), desc="Calculating closest distances"):
        closest_distance = float('inf')
        for j in range(len(points)):
            if i != j:  
                distance = math.sqrt((points[j][0] - points[i][0])**2 + (points[j][1] - points[i][1])**2)
                closest_distance = min(closest_distance, distance)
        closest_distances.append(closest_distance)

    fig = plt.figure(figsize=(8, 6))
    sns.histplot(closest_distances, kde=True, color='steelblue')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title('Distribution of the closest distance')
    plt.show()   

    return closest_distances


def get_neighbors(matrix:csr_matrix):
    """
    Calculate the number of neighbors for each cell.

    Parameters:
        matrix (sparse matrix): Spatial connectivity matrix.

    Returns:
        neiNum (Counter): Counter object containing the count of neighbors for each cell.
    """
    # Calculate the number of neighbors for each cell
    neiNum = Counter(matrix.tocoo().row)

    return neiNum


def process_dataframe(df:pd.DataFrame, matrix:csr_matrix, neiNum:Counter, save:bool, root ='neighbor_info/'):
    """
    Processes the dataframe by adding additional columns based on the neighbor matrix and neighbor counts.

    Parameters:
        - df (DataFrame): The input dataframe containing the barcode, cell type, coordinates, and other information.
        - matrix (sparse matrix): The neighbor matrix representing the spatial connectivities.
        - neiNum (Counter): The neighbor counts for each cell index.
        - save (bool): Flag indicating whether to save the processed dataframe.
        - root (str): Root directory for saving the dataframe (default is 'neighbor_info/').

    Returns:
        - df (DataFrame): The processed dataframe with additional columns for neighbor information.
    """
    # Add 'neiNum' column with neighbor counts
    df['neiNum'] = [neiNum[idx] for idx in range(len(df))]

    # Create a list to store neighbor types
    neiType_list = []
    cumulativeIdx = 0
    for idx in range(len(df)):
        start_idx = cumulativeIdx
        end_idx = cumulativeIdx + df['neiNum'][idx]
        neiIdx_list = matrix.tocoo().col[start_idx:end_idx]
        celltype_list = df['celltype1'][neiIdx_list].tolist()
        neiType_list.append(celltype_list)
        cumulativeIdx += df['neiNum'][idx]

    # Assign neighbor types to 'neiType' column
    df['neiType'] = neiType_list

    # Calculate 'catNum' column
    df['catNum'] = [len(set(neiType)) - (df['celltype1'][idx] in neiType) for idx, neiType in enumerate(df['neiType'])]
    
    # Filter rows based on conditions
    df = df[(df['neiNum'] > 0) & ((df['catNum'] == 0) | (df['catNum'] == 1))]
    df.reset_index(drop=True, inplace=True)
    
    def get_neighboring_cell_type(neiType, celltype1):
        """
        Retrieve the neighboring cell type based on the neighbor types and current cell type.

        Parameters:
            neiType (str): Neighbor types.
            celltype1 (str): Current cell type.

        Returns:
            celltype2 (str): Neighboring cell type.
        """
        # Extract neighboring cell type from the list
        temp = str(neiType).split("'")       
        for temp_neiType in temp:
            if temp_neiType != celltype1 and temp_neiType != ', ' and temp_neiType != '[' and temp_neiType != ']':
                return temp_neiType
        return 'NA'

    # Assign values to 'celltype2' column
    for i in range(len(df)):
        if df['catNum'][i] == 0:
            df['celltype2'][i] = df['celltype1'][i]
        else:    
            df['celltype2'][i] = get_neighboring_cell_type(neiType=df['neiType'][i], celltype1=df['celltype1'][i])

    # Assign integer codes to 'first_type' and 'second_type' columns
    df['first_type'] = df['celltype1'].cat.codes
    df['celltype2'] = df['celltype2'].astype('category')
    df['second_type'] = df['celltype2'].cat.codes

    # Replace spaces and slashes in cell type columns
    df['celltype1'] = df['celltype1'].str.replace(' ', '-')
    df['celltype1'] = df['celltype1'].str.replace('/', '-')
    df['celltype2'] = df['celltype2'].str.replace(' ', '-')
    df['celltype2'] = df['celltype2'].str.replace('/', '-')

    # Drop multiple columns
    columns_to_drop = ['neiNum', 'neiType', 'catNum']
    df = df.drop(columns_to_drop, axis=1)

    # Arrange columns in a specific order
    desired_order = ['barcode', 'first_type', 'second_type', 'celltype1', 'celltype2', 'x', 'y', 'prop1', 'prop2']
    df = df[desired_order]

    # Check if saving is requested
    if save == True:
        if not os.path.exists(root):
            os.makedirs(root)  
        df.to_csv(root + 'df_processed.csv', index = None)

    return df