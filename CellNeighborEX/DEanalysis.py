import os
import pandas as pd
import numpy as np
from scipy.stats import ranksums, ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from adjustText import adjust_text
from scipy.stats import f
from scipy import stats
from tqdm import tqdm


def two_sample_f_test(data1, data2):
    
    """
    Calculates the F-statistic and p-value for a two-sample F-test.

    Parameters:
        data1 (array-like): The first sample data.
        data2 (array-like): The second sample data.

    Returns:
        tuple: A tuple containing the F-statistic and p-value.
            - f_stat (float): The F-statistic value.
            - p_value (float): The p-value associated with the F-statistic.
    """
    
    # Calculate the variances
    var1 = np.var(data1, ddof=1)
    var2 = np.var(data2, ddof=1)
    
    if var2 == 0:
        
        return None, np.nan
    
    # Calculate the F-statistic
    f_stat = var1 / var2
    
    # Calculate the degrees of freedom
    df1 = len(data1) - 1
    df2 = len(data2) - 1
    
    # Calculate the p-value
    p_value = 2 * min(f.cdf(f_stat, df1, df2), 1 - f.cdf(f_stat, df1, df2))
    
    return f_stat, p_value


def create_nullmodel(seedNumber, randSize, prop, matchComb, clusterSelect, log_data):
    
    """
    Generates artificial heterotypic spots and computes normalized properties.

    Parameters:
        seedNumber (int): Seed number for randomization.
        randSize (int): Size of randomization.
        prop (numpy.ndarray): Proportions of cell types in heterotypic beads.
        matchComb (numpy.ndarray): Array indicating the matching combinations.
        clusterSelect (numpy.ndarray): Array specifying cluster numbers.
        log_data (numpy.ndarray): Log data of cells.

    Returns:
        log_data_artificialHeteroSpots (numpy.ndarray): Artificial heterotypic spots.
        normalized_props (numpy.ndarray): Normalized properties of heterotypic beads.
    """
    
    # Set the seed for the random number generator to ensure reproducibility
    np.random.seed(seedNumber)
    
    # Calculate the number of elements in 'matchComb' that are equal to 'clusterSelect[0]'
    true_heteroSpot_size = np.sum(matchComb == clusterSelect[0])
    
    # Create an array of zeros with dimensions (number of rows in 'log_data', true_heteroSpot_size * randSize)
    log_data_artificialHeteroSpots = np.zeros((log_data.shape[0], true_heteroSpot_size * randSize))
    
    # Create an array of zeros with dimensions (true_heteroSpot_size, 2)
    normalized_props = np.zeros((true_heteroSpot_size, 2))
    
    # Iterate over the range of 'true_heteroSpot_size'
    for idx in range(true_heteroSpot_size):
        
        # Assign the value at index 1 of 'clusterSelect' to 'clusterIndex1'
        clusterIndex1 = clusterSelect[1]
        
        # Assign the value at index 2 of 'clusterSelect' to 'clusterIndex2'
        clusterIndex2 = clusterSelect[2]
        
        # Find the indices in 'matchComb' where the value is equal to 'clusterIndex1'
        cellIndex1 = np.where(matchComb == clusterIndex1)[0]
        
        # Find the indices in 'matchComb' where the value is equal to 'clusterIndex2'
        cellIndex2 = np.where(matchComb == clusterIndex2)[0]
        
        # Calculate the sum of the values at 'prop1[idx]' and 'prop2[idx]'
        sum_prop = prop['prop1'][idx] + prop['prop2'][idx]
        
        # Calculate the normalized proportion for 'prop1'
        prop1 = prop['prop1'][idx] / sum_prop
        
        # Calculate the normalized proportion for 'prop2'
        prop2 = prop['prop2'][idx] / sum_prop
        
        # Assign 'prop1' and 'prop2' to the corresponding positions in 'normalized_props'
        normalized_props[idx, 0] = prop1
        normalized_props[idx, 1] = prop2
        
        # Generate random indices for 'cellIndex1' with a size of 'randSize'
        randIndex1 = np.random.randint(0, len(cellIndex1), size=randSize)
        
        # Generate random indices for 'cellIndex2' with a size of 'randSize'
        randIndex2 = np.random.randint(0, len(cellIndex2), size=randSize)
        
        # Select elements from 'cellIndex1' based on the generated random indices
        selected1 = cellIndex1[randIndex1]
        
        # Convert the selected indices to strings
        stringIDX1 = [str(num) for num in selected1]
        
        # Select elements from 'cellIndex2' based on the generated random indices
        selected2 = cellIndex2[randIndex2]
        
        # Convert the selected indices to strings
        stringIDX2 = [str(num) for num in selected2]
        
        # Select columns from 'log_data' based on 'stringIDX1' and multiply them by 'prop1'
        df1 = prop1 * log_data.loc[:, stringIDX1]
        
        # Select columns from 'log_data' based on 'stringIDX2' and multiply them by 'prop2'
        df2 = prop2 * log_data.loc[:, stringIDX2]
        
        # Create new column names as strings from 0 to the number of columns in 'df1'
        new_column_names = [str(i) for i in range(len(df1.columns))]
        
        # Assign the new column names to the columns of 'df1' and 'df2'
        df1.columns = new_column_names
        df2.columns = new_column_names
        
        # Add 'df1' and 'df2' element-wise and assign the result to 'tmp_log_data'
        tmp_log_data = df1 + df2
        
        # Assign the values of 'tmp_log_data' to the corresponding columns in 'log_data_artificialHeteroSpots'
        log_data_artificialHeteroSpots[:, idx * randSize:(idx + 1) * randSize] = tmp_log_data

        
    return pd.DataFrame(log_data_artificialHeteroSpots), normalized_props


def find_nullDEGs(center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, log_data_artificialHeteroSpots, gene_name, pCutoff, pCutoff2, lrCutoff, direction, normality_test):
    
    """
    Find Null Differentially Expressed Genes (DEGs) based on the provided parameters.

    Parameters:
    - center_celltype: The center cell type.
    - clusterSelect: The selected cluster.
    - matchComb: The combination used for matching.
    - neiCombUnique: Unique neighboring combination.
    - log_data: The log-transformed data.
    - log_data_artificialHeteroSpots: Artificial heterogeneity spots in log-transformed data.
    - gene_name: Gene names.
    - pCutoff: P-value cutoff.
    - pCutoff2: Secondary P-value cutoff.
    - lrCutoff: Log ratio cutoff.
    - direction: Direction flag for DEG filtering.

    Returns:
    - null_DEGs: A list of gene names identified as Null DEGs.
    - logRatio_null: The log ratio values of Null DEGs.
    - pvalue_null: The p-values of Null DEGs.
    - fdr_null: The False Discovery Rate (FDR) values of Null DEGs.
    """

    cellSelect1 = matchComb == clusterSelect[0]  # Selects cells based on matching combinations.

    max_leng = 1  # Maximum length for arrays.
    pvalue_total = [None] * max_leng  # Placeholder array for p-values.
    fdr_total = [None] * max_leng  # Placeholder array for false discovery rates (FDR).
    logRatio_total = [None] * max_leng  # Placeholder array for log ratios.

    pvalue = np.ones(log_data.shape[0])  # Initializes p-values with ones.
    fdr = np.ones(log_data.shape[0])  # Initializes FDR with ones.
    logRatio = np.zeros(log_data.shape[0])  # Initializes log ratios with zeros.

    sample_size = 30
    alpha = 0.05
    for i in range(log_data.shape[0]):
        
        if (log_data_artificialHeteroSpots.shape[1] > 2) and (np.sum(cellSelect1) > 2):
            
            if normality_test == True:
                
                _, p_normal2 = stats.shapiro(log_data_artificialHeteroSpots.loc[i, :])
                _, p_normal1 = stats.shapiro(log_data.loc[i, cellSelect1])
                
                if (p_normal2 > alpha) and (p_normal1 > alpha):
                    
                    # Variance test: Calculates the test statistic and heteroscedasticity variance.
                    _, h_var = two_sample_f_test(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1])
                    
                    if np.isnan(h_var):
                        p = 1  # If heteroscedasticity variance is NaN, set p-value to 1.
                    elif h_var < alpha:
                        # Welch's t-test: unequal variance
                        _, p = ttest_ind(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1], equal_var=False)
                    else:
                        # Student's t-test: equal variance
                        _, p = ttest_ind(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1])
                        
                else: 
                    # Wilcoxon rank sum test
                    _, p = ranksums(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1])            
            else:
                
                if (log_data_artificialHeteroSpots.shape[1] > sample_size) and (np.sum(cellSelect1) > sample_size):
                    
                    # Variance test: Calculates the test statistic and heteroscedasticity variance.
                    _, h_var = two_sample_f_test(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1])
                    
                    if np.isnan(h_var):
                        p = 1  # If heteroscedasticity variance is NaN, set p-value to 1.
                    elif h_var < alpha:
                        # Welch's t-test: unequal variance
                        _, p = ttest_ind(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1], equal_var=False)
                    else:
                        # Student's t-test: equal variance
                        _, p = ttest_ind(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1])
                        
                else: 
                    # Wilcoxon rank sum test
                    _, p = ranksums(log_data_artificialHeteroSpots.loc[i, :], log_data.loc[i, cellSelect1])  
                    
        pvalue[i] = p # Assign calculated p-value to array.
        logRatio[i] = np.mean(log_data.loc[i, cellSelect1] + 1) - np.mean(log_data_artificialHeteroSpots.loc[i, :] + 1) # Calculate log ratio and assign it to array.
        
    _, fdr, _, _ = multipletests(pvalue, method='fdr_tsbh', alpha=0.15) # Adjust p-values using the Benjamini-Hochberg procedure.

    pvalue_total[0] = pvalue  # Store p-values in total array.
    fdr_total[0] = fdr  # Store FDR values in total array.
    logRatio_total[0] = logRatio  # Store log ratios in total array.

    max_leng = len(pvalue_total)  # Update maximum length based on total array size.
    DEGindex = np.zeros((gene_name.shape[0], max_leng))  # Initialize DEG index array.
    DEGindex = np.zeros((gene_name.shape[0], max_leng))  # Initialize DEG index array.
    
    if direction == 'up':
        # If direction is 1 (positive), check for significant upregulated DEGs.
        DEGindex[:, 0] = (fdr_total[0] < pCutoff2) & (logRatio_total[0] > lrCutoff)
    elif direction == 'down':
        # Otherwise (negative), check for significant downregulated DEGs.
        DEGindex[:, 0] = (fdr_total[0] < pCutoff2) & (logRatio_total[0] < -lrCutoff)
    else:
        DEGindex = np.zeros((gene_name.shape[0], max_leng)) 
        
    geneIndexTemp = np.where(DEGindex == 1)[0]  # Get indices of significant DEGs.
    sortIndex = np.argsort(logRatio_total[0][geneIndexTemp])[::-1]  # Sort indices based on log ratios.
    geneIndex_sorted = geneIndexTemp[sortIndex]  # Sort gene indices based on sorted indices.

    logRatio_null = logRatio_total[0]  # Get log ratios for null DEGs.
    logRatio_null = logRatio_null[geneIndex_sorted]  # Sort log ratios based on sorted gene indices.
    logRatio_null_total = logRatio_total[0]
    pvalue_null = pvalue_total[0]  # Get p-values for null DEGs.
    pvalue_null = pvalue_null[geneIndex_sorted]  # Sort p-values based on sorted gene indices.
    fdr_null = fdr_total[0]  # Get FDR values for null DEGs.
    fdr_null = fdr_null[geneIndex_sorted]  # Sort FDR values based on sorted gene indices.
    fdr_null_total = fdr_total[0]
    null_DEGs = list(gene_name[0][geneIndex_sorted])  # Get names of null DEGs.

    return null_DEGs, logRatio_null, pvalue_null, fdr_null, logRatio_null_total, fdr_null_total


def find_contactDEGs(data_type, center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, gene_name, pCutoff, pCutoff2, lrCutoff, null_DEGs, fdr_null, logRatio_null, direction, normality_test):
    """
    Find differentially expressed genes (DEGs) related to cell-cell contacts.
    
    Parameters:
    - data_type: Type of data. It has to be 'Image' or 'NGS'.
    - center_celltype: The center cell type used for comparison.
    - clusterSelect: List of cluster indices to select for comparison.
    - matchComb: Combination of cluster indices.
    - neiCombUnique: Unique combination of neighboring cluster indices.
    - log_data: Log-transformed data for gene expression.
    - gene_name: Names of genes.
    - pCutoff: Cutoff value for p-value.
    - lrCutoff: Cutoff value for log-ratio.
    - null_DEGs: DEGs identified from the null model.
    - fdr_null: False discovery rate (FDR) for null DEGs.
    - logRatio_null: Log-ratio values for null DEGs.
    - direction: Direction of differential expression. Can be 1 for up-regulated or -1 for down-regulated.
    
    Return:
    - cellContact_DEGs_IDX: Indices of cell-contact DEGs.
    - cellContact_DEGs: Names of cell-contact DEGs.
    - logRatio1_cellContact: Log-ratio values for cell-contact DEGs (comparison 1).
    - pvalue1_cellContact: p-values for cell-contact DEGs (comparison 1).
    - fdr1_cellContact: False discovery rate (FDR) for cell-contact DEGs (comparison 1).
    - logRatio2_cellContact: Log-ratio values for cell-contact DEGs (comparison 2) (only for data_type = 'NGS').
    - pvalue2_cellContact: p-values for cell-contact DEGs (comparison 2) (only for data_type = 'NGS').
    - fdr2_cellContact: False discovery rate (FDR) for cell-contact DEGs (comparison 2) (only for data_type = 'NGS').
    - logRatio_null_cellContact: Log-ratio values for null DEGs (cell-contact) (only for data_type = 'NGS').
    - fdr_null_cellContact: False discovery rate (FDR) for null DEGs (cell-contact) (only for data_type = 'NGS').
    - logRatio1_total: Log-ratio values for all genes.
    - pvalue1_total: p-values for all genes.
    - fdr1_total: False discovery rate (FDR) for all genes 
    """
    
    cellSelect1 = matchComb == clusterSelect[0]  # Selects cells that match clusterSelect[0]
    cellSelect1 = np.array(cellSelect1)
    cellSelect2 = matchComb == clusterSelect[1]  # Selects cells that match clusterSelect[1]
    cellSelect2 = np.array(cellSelect2)
    
    if data_type == 'NGS':
        cellSelect3 = matchComb == clusterSelect[2]  # Selects cells that match clusterSelect[2]
        cellSelect3 = np.array(cellSelect3)   
    
    max_leng = 1  # Maximum length for result arrays
    pvalue1_total = [None] * max_leng  # Array to store p-values for comparison 1
    fdr1_total = [None] * max_leng  # Array to store FDR values for comparison 1
    logRatio1_total = [None] * max_leng  # Array to store log ratios for comparison 1
    pvalue1 = np.ones(log_data.shape[0])  # Array to store p-values for comparison 1
    logRatio1 = np.zeros(log_data.shape[0])  # Array to store log ratios for comparison 1
    
    if data_type == 'NGS':
        pvalue2_total = [None] * max_leng  # Array to store p-values for comparison 2
        fdr2_total = [None] * max_leng  # Array to store FDR values for comparison 2
        logRatio2_total = [None] * max_leng  # Array to store log ratios for comparison 2
        pvalue2 = np.ones(log_data.shape[0])  # Array to store p-values for comparison 2
        logRatio2 = np.zeros(log_data.shape[0])  # Array to store log ratios for comparison 2

    sample_size = 30  # Sample size for statistical tests
    alpha = 0.05
    for i in range(log_data.shape[0]):
        
        if (np.sum(cellSelect2) > 2) and (np.sum(cellSelect1) > 2): 
            
            if normality_test == True:
                
                _, p_normal2 = stats.shapiro(log_data.loc[i, cellSelect2])
                _, p_normal1 = stats.shapiro(log_data.loc[i, cellSelect1])
                
                if (p_normal2 > alpha) and (p_normal1 > alpha):
                    
                    # Variance test
                    _, h_var = two_sample_f_test(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1])
                    
                    if np.isnan(h_var):
                        p = 1
                    elif h_var < alpha:
                        # Welch's t-test: unequal variance
                        _, p = ttest_ind(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1], equal_var=False)
                    else:
                        # Student's t-test: equal variance
                        _, p = ttest_ind(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1])
                
                else:
                    
                    # Wilcoxon rank sum test
                    _, p = ranksums(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1])
                    
            else:
                
                #Comparing heterotypic beads with the first cell type of homotypic beads
                if (np.sum(cellSelect2) > sample_size) and (np.sum(cellSelect1) > sample_size):
                    _, h_var = two_sample_f_test(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1])
                
                    if np.isnan(h_var):
                        p = 1
                    elif h_var < alpha:
                        # Welch's t-test: unequal variance
                        _, p = ttest_ind(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1], equal_var=False)
                    else:
                        # Student's t-test: equal variance
                        _, p = ttest_ind(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1])
                
                else:
                    # Wilcoxon rank sum test
                    _, p = ranksums(log_data.loc[i, cellSelect2], log_data.loc[i, cellSelect1])
            pvalue1[i] = p
            logRatio1[i] = np.mean(log_data.loc[i, cellSelect1] + 1) - np.mean(log_data.loc[i, cellSelect2] + 1)
        
        if data_type == 'NGS':
            
            if (np.sum(cellSelect3) > 2) and (np.sum(cellSelect1) > 2):
                
                if normality_test == True:
                    _, p_normal3 = stats.shapiro(log_data.loc[i, cellSelect3])
                    _, p_normal1 = stats.shapiro(log_data.loc[i, cellSelect1])
                    
                    if (p_normal3 > alpha) and (p_normal1 > alpha):
                        
                        _, h_var2 = two_sample_f_test(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1])
                        if np.isnan(h_var2):
                            p2 = 1
                        elif h_var2 < alpha:
                            _, p2 = ttest_ind(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1], equal_var=False)
                        else:
                            _, p2 = ttest_ind(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1])
                            
                    else:
                        _, p2 = ranksums(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1])
                
                else:
                    
                    #Comparing heterotypic beads with the second cell type of homotypic beads
                    
                    if (np.sum(cellSelect3) > sample_size) and (np.sum(cellSelect1) > sample_size):
                        _, h_var2 = two_sample_f_test(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1])
                        if np.isnan(h_var2):
                            p2 = 1
                        elif h_var2 < alpha:
                            _, p2 = ttest_ind(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1], equal_var=False)
                        else:
                            _, p2 = ttest_ind(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1])
                    else: 
                        _, p2 = ranksums(log_data.loc[i, cellSelect3], log_data.loc[i, cellSelect1])
                
                pvalue2[i] = p2
                logRatio2[i] = np.mean(log_data.loc[i, cellSelect1] + 1) - np.mean(log_data.loc[i, cellSelect3] + 1)
            
    _, fdr1, _, _ = multipletests(pvalue1, method='fdr_tsbh', alpha=0.15) # Adjust p-values using the Benjamini-Hochberg procedure.
    
    pvalue1_total[0] = pvalue1 # Store p-values in total array.
    fdr1_total[0] = fdr1 # Store FDR values in total array.
    logRatio1_total[0] = logRatio1 # Store log ratios in total array.
    
    if data_type == 'NGS':
        _, fdr2, _, _ = multipletests(pvalue2, method='fdr_tsbh', alpha=0.15)
        pvalue2_total[0] = pvalue2
        fdr2_total[0] = fdr2
        logRatio2_total[0] = logRatio2

    # Finding DEGs: comparison between heterotypic beads and homotypic beads
    max_leng = len(pvalue1_total) # Update maximum length based on total array size.
    DEGindex = np.zeros((gene_name.shape[0], max_leng)) # Initialize DEG index array.

    if data_type == 'Image':
        
        if direction == 'up':  # up-regulated genes
            DEGindex[:, 0] = (pvalue1_total[0] < pCutoff) & (logRatio1_total[0] > lrCutoff) & (fdr1_total[0] < pCutoff2)
        elif direction == 'down':  # down-regulated genes
            DEGindex[:, 0] = (pvalue1_total[0] < pCutoff) & (logRatio1_total[0] < -lrCutoff) & (fdr1_total[0] < pCutoff2)    
        else:
            DEGindex = np.zeros((gene_name.shape[0], max_leng))
            
        cellContact_DEGs_IDX = np.where(DEGindex == 1)[0]  # Get indices of significant DEGs.
        cellContact_DEGs = [] # Initialize DEG list.
        for i in range(len(cellContact_DEGs_IDX)):
            index = cellContact_DEGs_IDX[i]
            cellContact_DEGs.append(gene_name[0][index]) # Get the names of DEGs
    
        # Confirm the significance of DEGs: intersection between DEGx and DEGy
        # DEGy (temp_DEGs): DEGs identified from the comparison between heterotypic and homotypic beads
        # Final cell-contact DEGs
        logRatio1_cellContact = logRatio1_total[0][cellContact_DEGs_IDX]
        logRatio1_total = logRatio1_total[0]
        pvalue1_cellContact = pvalue1_total[0][cellContact_DEGs_IDX]
        pvalue1_total = pvalue1_total[0]
        fdr1_cellContact = fdr1_total[0][cellContact_DEGs_IDX]
        fdr1_total = fdr1_total[0]
        
        return cellContact_DEGs_IDX, cellContact_DEGs, logRatio1_cellContact, pvalue1_cellContact, fdr1_cellContact, logRatio1_total, pvalue1_total, fdr1_total
                 
    elif data_type == 'NGS':
        
        if direction == 'up':  # Up-regulated genes
            DEGindex[:, 0] = (
                (pvalue1_total[0] < pCutoff) & (logRatio1_total[0] > lrCutoff) & (pvalue2_total[0] < pCutoff) & (logRatio2_total[0] > lrCutoff)
            )
        elif direction == 'down':  # Down-regulated genes
            DEGindex[:, 0] = (
                (pvalue1_total[0] < pCutoff) & (logRatio1_total[0] < -lrCutoff) & (pvalue2_total[0] < pCutoff) & (logRatio2_total[0] < -lrCutoff)
            )
        else:
            DEGindex = np.zeros((gene_name.shape[0], max_leng))
            
        geneIndexTemp2 = np.where(DEGindex == 1)[0]  # Get indices of significant DEGs. 
        temp_DEGs = list(gene_name[0][geneIndexTemp2])  # Get names of DEGs.
    
        # Confirming the significance of DEGs: intersection between DEGx and DEGy
        # DEGx (null_DEGs): DEGs identified from the null model
        # DEGy (temp_DEGs): DEGs identified from the comparison between heterotypic and homotypic beads
    
        # Final cell-contact DEGs
        cellContact_DEGs = list(np.intersect1d(null_DEGs, temp_DEGs)) # Get names of DEGs.
        cellContact_DEGs_IDX = np.where(np.isin(gene_name, cellContact_DEGs))[0]  # Get indices of DEGs. 
    
        null_cellContact_DEGs_IDX = np.array([]) # Initialize DEG index array
        for gene in cellContact_DEGs:
            tmp_val = np.where(np.array(null_DEGs) == gene)[0]
            null_cellContact_DEGs_IDX = np.append(null_cellContact_DEGs_IDX, tmp_val) # Get indices of DEGs.
    
        null_cellContact_DEGs_IDX = null_cellContact_DEGs_IDX.astype(np.int64) # Convert float64 to int64 
        
        logRatio_null_cellContact = np.array(logRatio_null[null_cellContact_DEGs_IDX])
        fdr_null_cellContact = np.array(fdr_null[null_cellContact_DEGs_IDX])
    
        logRatio1_cellContact = logRatio1_total[0][cellContact_DEGs_IDX]
        logRatio1_total = logRatio1_total[0]
        pvalue1_cellContact = pvalue1_total[0][cellContact_DEGs_IDX]
        pvalue1_total = pvalue1_total[0]
        fdr1_cellContact = fdr1_total[0][cellContact_DEGs_IDX]
        fdr1_total = fdr1_total[0]
        
        logRatio2_cellContact = logRatio2_total[0][cellContact_DEGs_IDX]
        logRatio2_total= logRatio2_total[0]
        pvalue2_cellContact = pvalue2_total[0][cellContact_DEGs_IDX]
        pvalue2_total = pvalue2_total[0]
        fdr2_cellContact = fdr2_total[0][cellContact_DEGs_IDX]
        fdr2_total = fdr2_total[0]
        
        return cellContact_DEGs_IDX, cellContact_DEGs, logRatio1_cellContact, pvalue1_cellContact, fdr1_cellContact, logRatio2_cellContact, pvalue2_cellContact, fdr2_cellContact, logRatio_null_cellContact, fdr_null_cellContact, logRatio1_total, pvalue1_total, fdr1_total, logRatio2_total, pvalue2_total, fdr2_total
                
    else:
        
        print("Please choose a data type between Imga and NGS.")
     
    
def get_heatmap(data_type, center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, log_data_zvalue, gene_name, cell_id, cellContact_DEGs_IDX, logRatio1_cellContact, logRatio2_cellContact, folderName2):
    
    """
    Plot heatmaps for the given data.

    Parameters:
    - data_type (str): Type of data ('Image' or 'NGS').
    - center_celltype (str): Center cell type.
    - clusterSelect (list): List of cluster select values.
    - matchComb (numpy.ndarray): Array of match combinations.
    - neiCombUnique (list): List of unique neighbor combinations.
    - log_data (pandas.DataFrame): Log data.
    - log_data_zvalue (pandas.DataFrame): Z-value of log data.
    - gene_name (list): List of gene names.
    - cell_id (list): List of cell IDs.
    - cellContact_DEGs_IDX (numpy.ndarray): Indices of cell-contact DEGs.
    - logRatio1_cellContact (numpy.ndarray): Log ratio 1 of cell-contact.
    - logRatio2_cellContact (numpy.ndarray): Log ratio 2 of cell-contact.
    - folderName2 (str): Folder name to save the files.

    Returns:
    None
    """
    
    cellSelect1 = matchComb == clusterSelect[0]  # Selects cells based on the match combination with the first cluster
    cellSelect1 = np.array(cellSelect1)
    cellSelect2 = matchComb == clusterSelect[1]  # Selects cells based on the match combination with the second cluster
    cellSelect2 = np.array(cellSelect2)
    cellSelect3 = matchComb == clusterSelect[2]  # Selects cells based on the match combination with the third cluster
    cellSelect3 = np.array(cellSelect3)
        
    if (cellContact_DEGs_IDX.shape[0] > 0) and (data_type == 'Image') or (data_type == 'NGS'):  # Checks conditions for cell contact differentially expressed genes (DEGs)
        
        plt.close('all')  # Closes all existing plots
        
        fig = plt.figure(figsize=(20, 16))  # Creates a new figure with a size of 20x16 inches (TumorIII+Monocyte: 60x20)
        
        plt.subplot(3, 3, 1)  # Creates a subplot in a 3x3 grid at position 1
        
        cellSelect1_IDX = np.where(cellSelect1)  # Gets the indices of selected cells for cluster 1
        cellSelect1_array = np.array(cellSelect1_IDX[0], dtype=str)  # Converts the indices to an array of strings
        plt.imshow(log_data_zvalue.loc[cellContact_DEGs_IDX][cellSelect1_array], aspect='auto')  # Plots the heatmap of log-transformed data for selected cells
        
        a = plt.gca().get_yticklabels()  # Gets the y-axis tick labels
        plt.gca().set_yticklabels(a, fontsize=12)  # Sets the y-axis tick labels with fontsize 12
        plt.gca().tick_params(axis='y', which='both', labelsize=12)  # Sets tick parameters for the y-axis
        plt.yticks(np.arange(cellContact_DEGs_IDX.shape[0]))  # Sets the y-axis ticks
        plt.gca().set_yticklabels(gene_name[0][cellContact_DEGs_IDX])  # Sets the y-axis tick labels with gene names
        
        if cellContact_DEGs_IDX.shape[0] > 30:  # Checks if the number of DEGs is greater than 30
            ax = plt.gca()
            ax.yaxis.set_tick_params(labelsize=6)  # Sets smaller tick label font size for y-axis
        
        plt.xticks([])  # Disables x-axis ticks
        plt.clim(-3, 3)  # Sets the color limit for the colormap
        plt.set_cmap('bwr')  # Sets the colormap to 'blue-white-red'
        plt.title(f"{neiCombUnique[0]}\n(n={np.sum(cellSelect1)})", fontsize=12)  # Sets the title of the subplot
        
        # Iterate over the range of cellContact_DEGs_IDX shape
        for i in range(cellContact_DEGs_IDX.shape[0]):
            # Create the filename using neiCombUnique, gene_name, and cellContact_DEGs_IDX
            filename = f"{neiCombUnique[0]}_{gene_name[0][cellContact_DEGs_IDX[i]]}.txt"
            
            # Get log_values from log_data using cellContact_DEGs_IDX and cellSelect1_array
            log_values = log_data.loc[cellContact_DEGs_IDX[i]][cellSelect1_array]
            
            # Get zvalues from log_data_zvalue using cellContact_DEGs_IDX and cellSelect1_array
            zvalues = log_data_zvalue.loc[cellContact_DEGs_IDX[i]][cellSelect1_array]
            
            # Get barcodes from cell_id using cellSelect1_IDX and convert to numpy array
            barcodes = cell_id[0][np.array(cellSelect1_IDX[0])]
            
            # Create the integrated array by vertically stacking barcodes, log_values, and zvalues
            integrated = np.vstack((barcodes, log_values, zvalues)).T
            
            # Sort the integrated array based on the second column (index 1)
            integrated = integrated[integrated[:, 1].argsort()]
            
            # Convert the integrated array to a DataFrame
            integrated = pd.DataFrame(integrated)
            
            # Save the integrated DataFrame to a CSV file with the specified filename in folderName2
            integrated.to_csv(f"{folderName2}/{filename}", header=None, index=False)

        # Set up subplot configuration
        plt.subplot(3, 3, 2)

        # Get cellSelect2_IDX using np.where on cellSelect2
        cellSelect2_IDX = np.where(cellSelect2)

        # Convert cellSelect2_IDX to a numpy array of strings
        cellSelect2_array = np.array(cellSelect2_IDX[0], dtype=str)

        # Plot the image using log_data_zvalue, cellContact_DEGs_IDX, and cellSelect2_array
        plt.imshow(log_data_zvalue.loc[cellContact_DEGs_IDX][cellSelect2_array], aspect='auto')

        # Remove y-axis ticks
        plt.yticks([])

        # Remove x-axis ticks
        plt.xticks([])

        # Set the color limit for the image
        plt.clim(-3, 3)

        # Set the color map to 'bwr'
        plt.set_cmap('bwr')

        # Set the title for the subplot using neiCombUnique and the count of cellSelect2
        plt.title(f"{neiCombUnique[1]}\n(n={np.sum(cellSelect2)})", fontsize=12)
        
        # Iterate over the range of cellContact_DEGs_IDX shape
        for i in range(cellContact_DEGs_IDX.shape[0]):
            # Create the filename using neiCombUnique, gene_name, and cellContact_DEGs_IDX
            filename = f"{neiCombUnique[1]}_{gene_name[0][cellContact_DEGs_IDX[i]]}.txt"
            
            # Get log_values from log_data using cellContact_DEGs_IDX and cellSelect2_array
            log_values = log_data.loc[cellContact_DEGs_IDX[i]][cellSelect2_array]
            
            # Get zvalues from log_data_zvalue using cellContact_DEGs_IDX and cellSelect2_array
            zvalues = log_data_zvalue.loc[cellContact_DEGs_IDX[i]][cellSelect2_array]
            
            # Get barcodes from cell_id using cellSelect2_IDX and convert to numpy array
            barcodes = cell_id[0][np.array(cellSelect2_IDX[0])]
            
            # Create the integrated array by vertically stacking barcodes, log_values, and zvalues
            integrated = np.vstack((barcodes, log_values, zvalues)).T
            
            # Sort the integrated array based on the second column (index 1)
            integrated = integrated[integrated[:, 1].argsort()]
            
            # Convert the integrated array to a DataFrame
            integrated = pd.DataFrame(integrated)
            
            # Save the integrated DataFrame to a CSV file with the specified filename in folderName2
            integrated.to_csv(f"{folderName2}/{filename}", header=None, index=False)

        # Create an array of selected cells indices
        cellSelect3_IDX = np.where(cellSelect3)
        cellSelect3_array = np.array(cellSelect3_IDX[0], dtype=str)
        for i in range(cellContact_DEGs_IDX.shape[0]): # Loop over the indices of DEGs
                
                # Create a filename for saving the data
                filename = f"{neiCombUnique[2]}_{gene_name[0][cellContact_DEGs_IDX[i]]}.txt"
                
                # Extract relevant data for saving
                log_values = log_data.loc[cellContact_DEGs_IDX[i]][cellSelect3_array]
                zvalues = log_data_zvalue.loc[cellContact_DEGs_IDX[i]][cellSelect3_array]
                barcodes = cell_id[0][np.array(cellSelect3_IDX[0])]
                
                # Stack the data vertically
                integrated = np.vstack((barcodes, log_values, zvalues)).T
                
                # Sort the integrated data
                integrated = integrated[integrated[:, 1].argsort()]
                
                # Convert the integrated data to a DataFrame
                integrated = pd.DataFrame(integrated)
                
                # Save the integrated data to a CSV file
                integrated.to_csv(f"{folderName2}/{filename}", header=None, index=False)

        if data_type == 'NGS': # Checks if the data type is 'NGS'
            
            # Create a subplot at position (3,3,3)
            plt.subplot(3,3,3)
            
            # Create an array of selected cells indices
            cellSelect3_IDX = np.where(cellSelect3)
            cellSelect3_array = np.array(cellSelect3_IDX[0], dtype=str)
            
            # Plot the heatmap with selected data
            plt.imshow(log_data_zvalue.loc[cellContact_DEGs_IDX][cellSelect3_array], aspect= 'auto')
            
            # Remove y-axis and x-axis ticks
            plt.yticks([])
            plt.xticks([])
            
            # Set the colorbar limits
            plt.clim(-3, 3)
            
            # Set the colormap to 'bwr' (blue-white-red)
            plt.set_cmap('bwr')
            
            # Set the title of the plot
            plt.title(f"{neiCombUnique[2]}\n(n={np.sum(cellSelect3)})", fontsize=12)
            
            for i in range(cellContact_DEGs_IDX.shape[0]): # Loop over the indices of DEGs
                
                # Create a filename for saving the data
                filename = f"{neiCombUnique[2]}_{gene_name[0][cellContact_DEGs_IDX[i]]}.txt"
                
                # Extract relevant data for saving
                log_values = log_data.loc[cellContact_DEGs_IDX[i]][cellSelect3_array]
                zvalues = log_data_zvalue.loc[cellContact_DEGs_IDX[i]][cellSelect3_array]
                barcodes = cell_id[0][np.array(cellSelect3_IDX[0])]
                
                # Stack the data vertically
                integrated = np.vstack((barcodes, log_values, zvalues)).T
                
                # Sort the integrated data
                integrated = integrated[integrated[:, 1].argsort()]
                
                # Convert the integrated data to a DataFrame
                integrated = pd.DataFrame(integrated)
                
                # Save the integrated data to a CSV file
                integrated.to_csv(f"{folderName2}/{filename}", header=None, index=False)
                
        # Adjust the subplot layout    
        plt.subplots_adjust(left=0.12, right=0.94, bottom=0.08, top=0.94, wspace=0.2, hspace=0.2)
        
        # Add a colorbar to the plot
        plt.colorbar()
        
        # Adjust the plot layout to avoid overlapping elements
        plt.tight_layout()
        
        # Display the plot
        #plt.show()
        
        # Assign a filename for the heatmap plot
        filename_heatmap = center_celltype
        
        # Save the heatmap plot as a PDF file
        fig.savefig(f"{folderName2}/{filename_heatmap}_heatmap.pdf", dpi=300)


def get_volcano_plot(df, data_type, center_celltype, pCutoff, lrCutoff, top_genes, folderName2):
    """
    Generate a volcano plot based on the input DataFrame.

    Parameters:
        df (pandas.DataFrame): The input DataFrame containing the data.
        data_type (str): Choose between Image and NGS.
    """
    if data_type == "Image":
        x_name = 'logRatio1'
        y_name = 'fdr1'
        
    elif data_type == "NGS":
        final_logRatio = []
        final_pValue = []
        for i in range(len(df)):
            
            orgLogRatio = [df['logRatio1'][i], df['logRatio2'][i], df['logRatio_null'][i]]
            absLogRatio = [abs(df['logRatio1'][i]), abs(df['logRatio2'][i]), abs(df['logRatio_null'][i])]
            min_logRatio = min(absLogRatio)
            min_index = absLogRatio.index(min_logRatio)
            
            tempPval = [df['pvalue1'][i], df['pvalue2'][i], df['fdr_null'][i]]
            max_pVal = max(tempPval)
            max_index = tempPval.index(max_pVal)
         
            if min_logRatio <= lrCutoff:
                if min_index == 0:
                    final_logRatio.append(df['logRatio1'][i])
                    final_pValue.append(df['pvalue1'][i])
                elif min_index == 1:    
                    final_logRatio.append(df['logRatio2'][i])
                    final_pValue.append(df['pvalue2'][i])
                else:
                    final_logRatio.append(df['logRatio_null'][i])
                    final_pValue.append(df['fdr_null'][i])
                    
            elif max_pVal >= pCutoff:    
                  if max_index == 0:
                      final_logRatio.append(df['logRatio1'][i])
                      final_pValue.append(df['pvalue1'][i])
                  elif max_index == 1:    
                      final_logRatio.append(df['logRatio2'][i])
                      final_pValue.append(df['pvalue2'][i])
                  else:
                      final_logRatio.append(df['logRatio_null'][i])
                      final_pValue.append(df['fdr_null'][i])
            
            else:
                all_positive = all(num > 0 for num in orgLogRatio)
                all_negative = all(num < 0 for num in orgLogRatio)
                if (min_logRatio > 0.4) and (max_pVal < 0.01) and ((all_positive == True) or (all_negative == True)):
                    final_logRatio.append(df['logRatio_null'][i])
                    final_pValue.append(df['fdr_null'][i])
                else:    
                    final_logRatio.append(0)
                    final_pValue.append(0)
        
        df['smallest_logRatio'] = final_logRatio
        df['largest_pvalue'] = final_pValue      
            
        x_name = 'smallest_logRatio'
        y_name = 'largest_pvalue'
        
    else:
        raise ValueError("Please choose a data type between Imga and NGS.")
    
    plt.close('all')
    fig = plt.figure(figsize=(6, 4)) 

    plt.scatter(x=df[x_name],y=df[y_name].apply(lambda x:-np.log10(x)),s=4,label="Not significant")
 
    down = df[(df[x_name]<=-lrCutoff)&(df[y_name]<=pCutoff)]
    up = df[(df[x_name]>=lrCutoff)&(df[y_name]<=pCutoff)]

    plt.scatter(x=down[x_name],y=down[y_name].apply(lambda x:-np.log10(x)),s=4,label="Down-regulated",color="blue")
    plt.scatter(x=up[x_name],y=up[y_name].apply(lambda x:-np.log10(x)),s=4,label="Up-regulated",color="red")

    # Label top 10 up-regulated genes
    texts_up=[]
    up_top = up.nlargest(top_genes, x_name)
    for i, r in up_top.iterrows():
        texts_up.append(plt.text(x=r[x_name], y=-np.log10(r[y_name]), s=up_top['gene'][i]))
    adjust_text(texts_up,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    # Label top 10 down-regulated genes
    texts_down=[]
    down_top = down.nsmallest(top_genes, x_name)
    for i, r in down_top.iterrows():
        texts_down.append(plt.text(x=r[x_name], y=-np.log10(r[y_name]), s=down_top['gene'][i]))

    if texts_down != []:    
        adjust_text(texts_down,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    filename_volcano = center_celltype

    plt.title("Cell contact-dependent genes in "+f"{filename_volcano}")
    plt.xlabel("logFC")
    plt.ylabel("-logFDR") 
    plt.axvline(-lrCutoff,color="grey",linestyle="--")
    plt.axvline(lrCutoff,color="grey",linestyle="--")
    plt.axhline(-np.log10(pCutoff),color="grey",linestyle="--")
    #plt.legend(bbox_to_anchor=(1.0, 1.0))
    #plt.show()
    
    fig.savefig(f"{folderName2}/{filename_volcano}_volcano.pdf", dpi=300)


def delete_files_with_keyword(directory, keyword):
    
    file_list = os.listdir(directory)
    
    for filename in file_list:
        
        if keyword in filename:
            file_path = os.path.join(directory, filename)
            try:
                os.remove(file_path)
                print(f"Deleted file: {file_path}")
            except Exception as e:
                print(f"Error while deleting {file_path}: {e}")


def analyze_data(df_cell_id, df_gene_name, df_log_data, path_categorization, data_type, lrCutoff, pCutoff, pCutoff2, direction, normality_test, top_genes, save:bool, root ='DE_results/'):
    
    # Conducting neighbor-dependent gene expression analysis
    cellType_total_contact_DEGs = []
    total_contact_DEGs = []
    logRatio1_total_contact_DEGs = []
    pvalue1_total_contact_DEGS = []
    
    if data_type == "Image":
        
        files = os.listdir(path_categorization)  
        for file in files:
            if file == '.DS_Store':
                files.remove(file)
    
        fdr1_total_contact_DEGs = []
        
        org_cell_id = df_cell_id
        org_gene_name = df_gene_name
        org_log_data = df_log_data
        new_columns = range(len(org_log_data.columns))
        org_log_data = org_log_data.rename(columns=dict(zip(org_log_data.columns, new_columns)))
        
    elif data_type == "NGS":
        
        # Import data
        files = os.listdir(path_categorization)
        for file in files:
            if file == '.DS_Store':
                files.remove(file)
            
        logRatio2_total_contact_DEGs = []
        pvalue2_total_contact_DEGs = []
        logRatio_null_total_contact_DEGs = []
        fdr_null_total_contact_DEGs = []   
    
        org_cell_id = df_cell_id
        org_gene_name = df_gene_name
        org_log_data = df_log_data
        new_columns = range(len(org_log_data.columns))
        org_log_data = org_log_data.rename(columns=dict(zip(org_log_data.columns, new_columns)))     
                
    else:
        print("Please choose a data type between Imga and NGS.")   
    
             
    center_celltype_set = set()  # Set to store unique third elements
    for file_name in files:
        
        # Remove the '.csv' extension
        file_name = file_name.replace('.csv', '')
    
        # Split the name by '_' to get individual elements
        elements = file_name.split('_')
    
        # Add the third element to the set
        if len(elements) >= 3:
            center_celltype_set.add(elements[2])
    
            
    for k in tqdm(range(len(center_celltype_set)), desc = "neighbor-dependent gene expression analysis"):
        
        print(k)
        center_celltype = list(center_celltype_set)[k]
        cell_id = org_cell_id.copy()
        gene_name = org_gene_name.copy()
        log_data = org_log_data.copy()
        
        # Importing data
        index = pd.read_csv(path_categorization + 'index_categorized_' + center_celltype + '.csv', header=None)
        matchComb = pd.read_csv(path_categorization + 'matchComb_categorized_' + center_celltype + '.csv', header=None)
        matchComb = list(matchComb[0])
        neiCombUnique = pd.read_csv(path_categorization + 'neiCombUnique_categorized_' + center_celltype + '.csv', header=None)
        neiCombUnique = list(neiCombUnique[0])
        prop = pd.read_csv(path_categorization + 'prop_categorized_' + center_celltype + '.csv', delimiter=",", names=['prop1', 'prop2'])
        
        cell_id_total = cell_id[0]
        selected_row = index[index[0]==1]
        cell_id = cell_id_total[cell_id_total.index.isin(selected_row.index)]
        cell_id = cell_id.reset_index()
        
        log_data_total = log_data
        log_data = log_data_total.iloc[:, selected_row.index]
        column_names = [str(i) for i in range(len(matchComb))]
        log_data.columns = column_names
        
        mean_values = log_data_total.mean(axis=1)
        std_values = log_data_total.std(axis=1)
        log_data_zvalue = (log_data_total - mean_values[:, np.newaxis]) / std_values[:, np.newaxis]
        log_data_zvalue = log_data_zvalue.iloc[:, selected_row.index]
        log_data_zvalue.columns = column_names
        
        if data_type == "Image":
            clusterSelect = list(np.unique(matchComb))
            
        elif data_type == "NGS":
            
            seedNumber = 1
            clusterSelect = list(np.unique(matchComb))
    
            heteroSpotNum = np.sum(matchComb == clusterSelect[0])
            if heteroSpotNum > 100:
                randSize = 30
            else:
                randSize = 100
    
            log_data_artificialHeteroSpots, normalized_props = create_nullmodel(seedNumber, randSize, prop, matchComb, clusterSelect, log_data)
    
        else:
            print("Please choose a data type between Imga and NGS.")
            
        folderName2 = os.path.join(root, center_celltype)
        df = pd.DataFrame()
        df2 = pd.DataFrame()
        df_null2 = pd.DataFrame()
        if data_type == "Image":
            
            null_DEGs = None
            fdr_null = None
            logRatio_null = None
            
            cellContact_DEGs_IDX, cellContact_DEGs, logRatio1_cellContact, pvalue1_cellContact, fdr1_cellContact, logRatio1_total, pvalue1_total, fdr1_total = find_contactDEGs(data_type, center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, gene_name, pCutoff, pCutoff2, lrCutoff, null_DEGs, fdr_null, logRatio_null, direction, normality_test)
            if len(cellContact_DEGs) > 0:
                
                os.makedirs(folderName2, exist_ok=True)
                save_path = os.path.join(folderName2, center_celltype +'_stat_contact_DEGs.csv')
                save_path2 = os.path.join(folderName2, center_celltype +'_stat_total_genes.csv')
                df = pd.DataFrame({'cellContact_DEGs_IDX': cellContact_DEGs_IDX, 'cellContact_DEGs': cellContact_DEGs, 'logRatio1': logRatio1_cellContact, 'pvalue1': pvalue1_cellContact, 'fdr1': fdr1_cellContact})
                df.to_csv(save_path, index=False)
                df2 = pd.DataFrame({'gene':list(gene_name[0]),'logRatio1': logRatio1_total,'pvalue1': pvalue1_total, 'fdr1': fdr1_total})
                df2.to_csv(save_path2, index=False)
                
                tmp_cellType=[]
                for i in range(len(cellContact_DEGs)):
                    tmp_cellType.append(center_celltype)
                
                cellType_total_contact_DEGs.extend(tmp_cellType)    
                total_contact_DEGs.extend(cellContact_DEGs)
                logRatio1_total_contact_DEGs.extend(logRatio1_cellContact)
                pvalue1_total_contact_DEGS.extend(pvalue1_cellContact)
                fdr1_total_contact_DEGs.extend(fdr1_cellContact)
    
        elif data_type == "NGS":
            
            null_DEGs, logRatio_null, pvalue_null, fdr_null, logRatio_null_total, fdr_null_total  = find_nullDEGs(center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, log_data_artificialHeteroSpots, gene_name, pCutoff, pCutoff2, lrCutoff, direction, normality_test)
                                                                            
            if len(null_DEGs) > 0:
                os.makedirs(folderName2, exist_ok=True)
                save_path = os.path.join(folderName2, center_celltype +'_stat_null_cotact_genes.csv')
                df_null = pd.DataFrame({'null_DEGs': null_DEGs, 'pvalue_null': pvalue_null, 'fdr_null': fdr_null, 'logRatio_null': logRatio_null})
                #df_null.to_csv(save_path, index=False)
                save_path2 = os.path.join(folderName2, center_celltype +'_stat_null_total_genes.csv')
                df_null2 = pd.DataFrame({'gene': list(gene_name[0]),'logRatio_null': logRatio_null_total, 'fdr_null': fdr_null_total})
                #df_null2.to_csv(save_path2, index=False)
    
                #load_path = os.path.join(folderName2, center_celltype +'_stat_null_cotact_genes.csv')
                #loaded_data = pd.read_csv(load_path)
                
                cellContact_DEGs_IDX, cellContact_DEGs, logRatio1_cellContact, pvalue1_cellContact, fdr1_cellContact, logRatio2_cellContact, pvalue2_cellContact, fdr2_cellContact, logRatio_null_cellContact, fdr_null_cellContact, logRatio1_total, pvalue1_total, fdr1_total, logRatio2_total, pvalue2_total, fdr2_total = find_contactDEGs(data_type, center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, gene_name, pCutoff, pCutoff2, lrCutoff, null_DEGs, df_null['fdr_null'], df_null['logRatio_null'], direction, normality_test)
                
                if len(cellContact_DEGs) > 0:
                    save_path = os.path.join(folderName2, center_celltype +'_stat_contact_DEGs.csv')
                    df = pd.DataFrame({'cellContact_DEGs_IDX': cellContact_DEGs_IDX, 'cellContact_DEGs': cellContact_DEGs, 'pvalue1_cellContact': pvalue1_cellContact, 'fdr1_cellContact': fdr1_cellContact, 'logRatio1_cellContact': logRatio1_cellContact, 'pvalue2_cellContact': pvalue2_cellContact, 'fdr2_cellContact': fdr2_cellContact, 'logRatio2_cellContact': logRatio2_cellContact})
                    df.to_csv(save_path, index=False)
                
                    save_path2 = os.path.join(folderName2, center_celltype +'_stat_total_genes.csv')
                    df2 = pd.DataFrame({'gene': list(gene_name[0]), 'logRatio1': logRatio1_total, 'pvalue1': pvalue1_total,'logRatio2':logRatio2_total, 'pvalue2': pvalue2_total, 'logRatio_null': df_null2['logRatio_null'], 'fdr_null': df_null2['fdr_null']})
                    df2.to_csv(save_path2, index=False)
                    
                    tmp_cellType = np.full((len(cellContact_DEGs),), center_celltype)
                    cellType_total_contact_DEGs.extend(tmp_cellType)
                    total_contact_DEGs.extend(cellContact_DEGs)
                    logRatio1_total_contact_DEGs.extend(logRatio1_cellContact)
                    pvalue1_total_contact_DEGS.extend(pvalue1_cellContact)
                    logRatio2_total_contact_DEGs.extend(logRatio2_cellContact)
                    pvalue2_total_contact_DEGs.extend(pvalue2_cellContact)
                    logRatio_null_total_contact_DEGs.extend(logRatio_null_cellContact)
                    fdr_null_total_contact_DEGs.extend(fdr_null_cellContact)
            else:
                cellContact_DEGs = []
        else:
            print("Please choose a data type between Imga and NGS.")
            
        if len(cellContact_DEGs) > 0:
            if data_type == "Image":
                logRatio2_total_contact_DEGs = None
                get_volcano_plot(df2, data_type, center_celltype, pCutoff2, lrCutoff, top_genes, folderName2)
                get_heatmap(data_type, center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, log_data_zvalue, gene_name, cell_id, cellContact_DEGs_IDX, logRatio1_total_contact_DEGs, logRatio2_total_contact_DEGs, folderName2)
            
            elif data_type == "NGS": 
                get_volcano_plot(df2, data_type, center_celltype, pCutoff2, lrCutoff, top_genes, folderName2)
                get_heatmap(data_type, center_celltype, clusterSelect, matchComb, neiCombUnique, log_data, log_data_zvalue, gene_name, cell_id, cellContact_DEGs_IDX, logRatio1_total_contact_DEGs, logRatio2_total_contact_DEGs, folderName2)
    
    if data_type == "Image" and len(total_contact_DEGs) > 0:
        total_df = pd.DataFrame({'celltype': cellType_total_contact_DEGs, 'DEG': total_contact_DEGs, 'logRatio1': logRatio1_total_contact_DEGs, 'pvalue1': pvalue1_total_contact_DEGS, 'fdr1': fdr1_total_contact_DEGs})
        
    elif data_type == "NGS" and len(total_contact_DEGs) > 0:
        total_df = pd.DataFrame({'celltype': cellType_total_contact_DEGs, 'DEG': total_contact_DEGs, 'logRatio1': logRatio1_total_contact_DEGs, 'pvalue1': pvalue1_total_contact_DEGS, 'logRatio2': logRatio2_total_contact_DEGs, 'pvalue2': pvalue2_total_contact_DEGs, 'logRatio_null': logRatio_null_total_contact_DEGs, 'fdr_null': fdr_null_total_contact_DEGs})
    
    # Check if saving is requested
    if save == True:
        if not os.path.exists(root):
            os.makedirs(root)  
        total_df.to_csv(root + 'DEG_list.csv', index = None)
    
    return total_df