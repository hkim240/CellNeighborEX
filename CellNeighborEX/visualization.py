import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def import_expdata(path:str):
    """
    Import expression data from CSV files in the specified directory.

    Parameters:
        path (str): Path to the expression data directory.

    Returns:
        df (pd.DataFrame): DataFrame containing the merged expression data.
    """

    # Initialize an empty DataFrame to store the merged expression data
    df = pd.DataFrame()

    # Iterate over each file in the specified directory
    for filename in os.listdir(path):

        # Exclude the '.DS_Store' file commonly found on macOS
        if filename != '.DS_Store':
            print(filename)

            # Read the CSV file into a temporary DataFrame
            df_temp = pd.read_csv(path + '/' + filename, sep=',', header=0)

            # Create a new DataFrame for spatial data with selected columns
            df_spatial = pd.DataFrame({
                'barcode': df_temp['barcode'],
                'type': filename[:-4],  # Remove ".csv" extension from the filename
                'zscore': df_temp['zscore']
            })

            # Append the spatial DataFrame to the total DataFrame
            df = df.append(df_spatial, ignore_index=True)

    # Return the merged DataFrame containing all the expression data
    return df


def set_parameters(df_processed:pd.DataFrame, df_exp:pd.DataFrame, beadsize_bg:float, edgecolor_bg:tuple, beadcolor_bg:tuple, beadsize_red:float, beadsize_blue:float, beadsize_black:float, type_red:str , type_blue:str, type_black:str):
    """
    Set parameters for the processed DataFrame based on the expression data.

     Parameters:
        df_processed (pd.DataFrame): DataFrame to be processed.
        df_exp (pd.DataFrame): DataFrame containing the expression data.
        beadsize_bg (float): Bead size for the background.
        edgecolor_bg (tuple): Edge color for the background.
        beadcolor_bg (tuple): Bead color for the background.
        beadsize_red (float): Bead size for cells of type red.
        beadsize_blue (float): Bead size for cells of type blue.
        beadsize_black (float): Bead size for cells of type black.
        type_red (str): Type name for cells of type red.
        type_blue (str): Type name for cells of type blue.
        type_black (str): Type name for cells of type black.

    Returns:
        df_processed (pd.DataFrame): Processed DataFrame with updated parameters.
        df_red (pd.DataFrame): Subset of df_processed containing cells of type red.
        df_blue (pd.DataFrame): Subset of df_processed containing cells of type blue.
        df_black (pd.DataFrame): Subset of df_processed containing cells of type black.
    """
    
    df_processed['beadsize'] = beadsize_bg 
    df_processed['edgecolor'] = 'NA'
    df_processed['beadcolor'] = 'NA'
    for idx in range(len(df_processed)):
        df_processed['edgecolor'][idx] = edgecolor_bg  
        df_processed['beadcolor'][idx] = edgecolor_bg  

    df_processed['type'] = 'NA'
    df_processed['expression'] = (np.max(df_exp['zscore'])-np.min(df_exp['zscore']))/2  # The colors of cells or beads (background) with gray edges are set to white
    
    for idx in range(len(df_processed)):

        barcode = df_processed['barcode'][idx]
        matching_rows = df_exp[df_exp['barcode'] == barcode]

        if not matching_rows.empty:

            row = matching_rows.iloc[0]
            df_processed['expression'][idx] = row['zscore']

            if row['type'] == type_red:
                df_processed['edgecolor'][idx] = (1, 0, 0)  # Red
                df_processed['type'][idx] = type_red
                df_processed['beadsize'][idx] = beadsize_red  
            elif row['type'] == type_blue:
                df_processed['edgecolor'][idx] = (0, 0, 1)  # Blue
                df_processed['type'][idx] = type_blue
                df_processed['beadsize'][idx] = beadsize_blue  
            elif row['type'] == type_black:
                df_processed['edgecolor'][idx] = (0, 0, 0)  # Black
                df_processed['type'][idx] = type_black
                df_processed['beadsize'][idx] =  beadsize_black  

    df_red = df_processed[df_processed['type'] == type_red].reset_index()
    df_blue = df_processed[df_processed['type'] == type_blue].reset_index()
    df_black= df_processed[df_processed['type'] == type_black].reset_index()

    return df_processed, df_red, df_blue, df_black            


def get_spatialPlot(df_bg:pd.DataFrame, df_red:pd.DataFrame, df_blue:pd.DataFrame, df_black:pd.DataFrame, label_red:str, label_blue:str, label_black:str, label_gene:str, figsize:tuple, save:bool, root ='spatialMap/'):
    """
    Generate a spatial plot using the provided DataFrames and parameters.

     Parameters:
        df_bg (pd.DataFrame): DataFrame for the background.
        df_red (pd.DataFrame): DataFrame for cells of type red.
        df_blue (pd.DataFrame): DataFrame for cells of type blue.
        df_black (pd.DataFrame): DataFrame for cells of type black.
        label_red (str): Label for cells of type red in the legend.
        label_blue (str): Label for cells of type blue in the legend.
        label_black (str): Label for cells of type black in the legend.
        label_gene (str): Label for the gene expression in the colorbar title.
        figsize (tuple): Size of the figure (width, height).
        save (bool): Flag indicating whether to save the plot.
        root (str): Root directory for saving the plot.

    Returns:
        None
    """

    fig, ax = plt.subplots(figsize=figsize)  

    # Gray - background
    cbar_flag = ax.scatter(df_bg['x'], df_bg['y'], s=df_bg['beadsize'], c=df_bg['expression'].to_numpy(), cmap='bwr',
                           edgecolor=df_bg['edgecolor'], linewidths=2, zorder=2.0)

    # Red
    ax.scatter(df_red['x'], df_red['y'], s=df_red['beadsize'], c=df_red['expression'].to_numpy(),
               cmap='bwr', edgecolor=df_red['edgecolor'], linewidths=2, zorder=4.0)
    ax.scatter(df_red['x'], df_red['y'], s=df_red['beadsize'], c='white',
               edgecolor=df_red['edgecolor'], linewidths=2, label=label_red)

    # Blue
    ax.scatter(df_blue['x'], df_blue['y'], s=df_blue['beadsize'], c=df_blue['expression'].to_numpy(),
               cmap='bwr', edgecolor=df_blue['edgecolor'], linewidths=2, zorder=3.0)
    ax.scatter(df_blue['x'], df_blue['y'], s=df_blue['beadsize'], c='white',
               edgecolor=df_blue['edgecolor'], linewidths=2, label=label_blue)

    # Black
    ax.scatter(df_black['x'], df_black['y'], s=df_black['beadsize'], c=df_black['expression'].to_numpy(),
               cmap='bwr', edgecolor=df_black['edgecolor'], linewidths=2, zorder=5.0)
    ax.scatter(df_black['x'], df_black['y'], s=df_black['beadsize'], c='white',
               edgecolor=df_black['edgecolor'], linewidths=2, label=label_black)

    # Legend setting
    ax.legend(bbox_to_anchor=(0.98, 0.42), prop={'size': 30}, frameon=False, fontsize=30)
    legend = ax.get_legend()
    legend.set_title('Cell type', prop={'size': 30})
    legend._legend_box.align = "left"
    legend.legendHandles[0]._sizes = [80]
    legend.legendHandles[1]._sizes = [40]
    legend.legendHandles[2]._sizes = [40]

    # Spatial plot setting
    ax.axis('off')
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    cbar = plt.colorbar(cbar_flag, cax=fig.add_axes([0.91, 0.45, 0.03, 0.2]), ticks=[-0.5, 0, 1, 2, 3, 4, 5])
    cbar.ax.tick_params(labelsize=30)
    cbar.ax.set_title('%s \n expression' %label_gene, fontsize=30)

    # Check if saving is requested
    if save == True:
        if not os.path.exists(root):
            os.makedirs(root)  
        plt.savefig('spatialMap/%s.pdf' %label_gene, bbox_inches='tight')    
        
    #plt.show()