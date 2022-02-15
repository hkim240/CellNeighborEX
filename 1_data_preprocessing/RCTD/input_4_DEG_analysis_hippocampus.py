####### input for DEG analysis without neighbor analysis #######
import pandas as pd
from collections import Counter
import seaborn as sns
import copy

### Importing data
# All beads with unassinged cell type were removed
path_hippo = '/Users/hyobinkim/RCTD_Plots/mouse_hippocampus/'
df_RCTD_hippo = pd.read_csv(path_hippo + 'class_unassigned_removed.csv', sep = ';', header = 0)
df_name_hippo = pd.read_csv(path_hippo + 'hippocampus_cell_types.csv', sep = ';', header = 0)
df_prop_hippo = pd.read_csv(path_hippo + 'composition.csv', sep = ';', header = 0)
df_prop_hippo = df_prop_hippo.rename(columns={'Unnamed: 0': 'barcode'})
df_prop_hippo_rev = df_prop_hippo.stack().str.replace(',','.').unstack()  
df_location_hippo = pd.read_csv(path_hippo + 'data/Hippocampus_BeadLocationsForR.csv', sep = ',', header = 0)

### Changing the value of RCTD['scond_type'] when the bead is a singlet
for i in range(len(df_RCTD_hippo['barcode'])):
    if df_RCTD_hippo['spot_class'][i] == 'singlet':
        df_RCTD_hippo['second_type'][i] = df_RCTD_hippo['first_type'][i]

### Adding the info of cell type names and abbrevs to RCTD
df_RCTD_hippo['celltype1'] = 'na'
df_RCTD_hippo['celltype2'] = 'na'
for i in range(len(df_RCTD_hippo['barcode'])):
    print(i)
    type1 = df_RCTD_hippo['first_type'][i]
    type2 = df_RCTD_hippo['second_type'][i]
    for j in range(len(df_name_hippo['name'])):   
        
        if type1 == df_name_hippo['cluster_number'][j]:
            df_RCTD_hippo['celltype1'][i] = df_name_hippo['name'][j]
        
        if type2 == df_name_hippo['cluster_number'][j]:
            df_RCTD_hippo['celltype2'][i] = df_name_hippo['name'][j]              

###################
# comb_celltype_hippcampus.csv

comb_celltype = []
for i in range(len(df_RCTD_hippo['barcode'])):
    print(i)
    ct1 = df_RCTD_hippo['celltype1'][i]
    ct2 = df_RCTD_hippo['celltype2'][i]
    comb_celltype.append(ct1 + '+' + ct2)

comb_celltype_hippo = pd.DataFrame({'barcode':df_RCTD_hippo['barcode'].tolist(), 'comb_celltype': comb_celltype})
comb_celltype_hippo.to_csv('/Users/hyobinkim/Documents/MATLAB/FindDEGs_artificialDoublets/final/11_NicheNet/hippocampus/Slide-seq/input_data/comb_celltype_hippocampus.csv', index = False) 

###################

df_RCTD_hippo['abbrev1'] = 'na'
df_RCTD_hippo['abbrev2'] = 'na'
for i in range(len(df_RCTD_hippo['barcode'])):
    print(i)
    type1 = df_RCTD_hippo['first_type'][i]
    type2 = df_RCTD_hippo['second_type'][i]
    for j in range(len(df_name_hippo['abbrev'])):   
        
        if type1 == df_name_hippo['cluster_number'][j]:
            df_RCTD_hippo['abbrev1'][i] = df_name_hippo['abbrev'][j]
        
        if type2 == df_name_hippo['cluster_number'][j]:
            df_RCTD_hippo['abbrev2'][i] = df_name_hippo['abbrev'][j] 
            
### Adding the info of proportions to RCTD
# extracting the same barcodes from df_prop_hippo_rev
list_barcode = df_RCTD_hippo['barcode'].tolist()
df_prop_hippo_rev2 = copy.deepcopy(df_prop_hippo_rev)
for i in range(len(df_prop_hippo_rev['barcode'])):
    print(i)
    if df_prop_hippo_rev['barcode'][i] not in list_barcode:
        df_prop_hippo_rev2 = df_prop_hippo_rev2.drop([i])
df_prop_hippo_rev2 = df_prop_hippo_rev2.reset_index()

# checking if df_prop_hippo_rev2['barcode'] is equal to df_RCTD_hippo['barcode']
tmp=[]
for i in range(len(df_prop_hippo_rev2['barcode'])):
    if df_prop_hippo_rev2['barcode'][i] == df_RCTD_hippo['barcode'][i]:
        tmp.append(1)

# Adding props to df_RCTD_hippo['prop1'] & ['prop2']
df_RCTD_hippo['prop1'] = 0
df_RCTD_hippo['prop2'] = 0
for i in range(len(df_RCTD_hippo['barcode'])):
    print(i)
    type1 = str(df_RCTD_hippo['first_type'][i])
    type2 = str(df_RCTD_hippo['second_type'][i])
    if type1 != type2:
        df_RCTD_hippo['prop1'][i] =  df_prop_hippo_rev2[type1][i]
        df_RCTD_hippo['prop2'][i] =  df_prop_hippo_rev2[type2][i]
    else:
        df_RCTD_hippo['prop1'][i] =  df_prop_hippo_rev2[type1][i]

### Adding the info of bead location to RCTD
# extracting the same barcodes from df_location_hippo
df_location_hippo2 = copy.deepcopy(df_location_hippo)
for i in range(len(df_location_hippo['barcodes'])):
    print(i)
    if df_location_hippo['barcodes'][i] not in list_barcode:
        df_location_hippo2 = df_location_hippo2.drop([i])
df_location_hippo2 = df_location_hippo2.reset_index()

# checking if df_prop_hippo_rev2['barcode'] is equal to df_RCTD_hippo['barcode']
tmp2=[]
for i in range(len(df_location_hippo2['barcodes'])):
    if df_location_hippo2['barcodes'][i] == df_RCTD_hippo['barcode'][i]:
        tmp2.append(1)

# Adding location to df_RCTD_hippo['x'] & ['y']
df_RCTD_hippo['x'] = -1
df_RCTD_hippo['y'] = -1
for i in range(len(df_RCTD_hippo['barcode'])):
    print(i)
    df_RCTD_hippo['x'][i] =  df_location_hippo2['xcoord'][i]
    df_RCTD_hippo['y'][i] =  df_location_hippo2['ycoord'][i]

df_RCTD_hippo.to_csv(path_hippo + 'Slide-seq_hippocampus.csv', index = False)




#%%%

### Importing data
path_hippo = '/Users/hyobinkim/RCTD_Plots/mouse_hippocampus/'
df_RCTD_hippo = pd.read_csv(path_hippo + 'Slide-seq_hippocampus.csv', sep = ',', header = 0)

### Assigning a color into each bead
palette = sns.color_palette("hls", 17)
list_color_type1 = []
list_color_type2 = []
for i in range(len(df_RCTD_hippo['barcode'])):
    num1 = df_RCTD_hippo['first_type'][i] - 1
    list_color_type1.append(palette[num1])
    num2 = df_RCTD_hippo['second_type'][i]- 1
    list_color_type2.append(palette[num2])

df_RCTD_hippo['color1'] = list_color_type1
df_RCTD_hippo['color2'] = list_color_type2

### Adding the info of doublet types
df_RCTD_hippo['value'] = df_RCTD_hippo.apply(lambda x: hash(frozenset([x["celltype1"],x["celltype2"]])),axis=1)

### Sorting the doublet types along the sample size
counter_doublets = Counter(df_RCTD_hippo['value'])
sorted_counter = sorted(counter_doublets.items(), key=lambda x: x[1], reverse=True)

### Generating input files for DEG analysis
sample_size = 30
min_sample_size = 1
doublet_type = []
bg_length = []
for i in range(len(sorted_counter)):
    
    print(i)
    if sorted_counter[i][1] >= sample_size:
        
        val1 = sorted_counter[i][0]
        df_RCTD_sub = df_RCTD_hippo[df_RCTD_hippo['value'] == val1]
        # selecting doublet_certain
        #df_RCTD_sub = df_RCTD_sub[df_RCTD_sub['spot_class'] == 'doublet_certain']
        list_index1 = df_RCTD_sub.index.tolist()
        
        if len(list_index1) >= sample_size: # for the number of heterotypic doublets with doublet certain
            idx = list_index1[0];
            if df_RCTD_sub['celltype1'][idx] != df_RCTD_sub['celltype2'][idx]:
                
                abbrev = df_RCTD_sub['abbrev1'][idx]
                abbrev2 = df_RCTD_sub['abbrev2'][idx]
                cell_type_abbrev = abbrev + '+' + abbrev2
                               
                cellTypeNum1 = df_RCTD_sub['first_type'][idx]
                cellTypeNum2 = df_RCTD_sub['second_type'][idx]
                
                # input: prop
                list_prop1 = []
                list_prop2 = []
                for i in range(len(df_RCTD_sub)):
                    index = list_index1[i]
                    if df_RCTD_sub['first_type'][index] == cellTypeNum1:
                        list_prop1.append(df_RCTD_sub['prop1'][index])
                        list_prop2.append(df_RCTD_sub['prop2'][index])
                    elif df_RCTD_sub['first_type'][index] == cellTypeNum2:
                        list_prop1.append(df_RCTD_sub['prop2'][index])
                        list_prop2.append(df_RCTD_sub['prop1'][index])
                    
                # input: index, neiCombUnique, matchComb    
                df_RCTD_hippo['index_boolean'] = 0
                for j in range(len(list_index1)):
                    df_RCTD_hippo['index_boolean'][list_index1[j]] = 1
    
                df_RCTD_tmp1 = df_RCTD_hippo[df_RCTD_hippo['celltype1'] == df_RCTD_hippo['celltype1'][idx]]
                df_RCTD_tmp1 = df_RCTD_tmp1[df_RCTD_tmp1['celltype2'] == df_RCTD_hippo['celltype1'][idx]]
                # selecting singlet
                #df_RCTD_tmp1 = df_RCTD_tmp1[df_RCTD_tmp1['spot_class'] == 'singlet']
                
                df_RCTD_tmp2 = df_RCTD_hippo[df_RCTD_hippo['celltype1'] == df_RCTD_hippo['celltype2'][idx]]
                df_RCTD_tmp2 = df_RCTD_tmp2[df_RCTD_tmp2['celltype2'] == df_RCTD_hippo['celltype2'][idx]]
                # selecting singlet 
                #df_RCTD_tmp2 = df_RCTD_tmp2[df_RCTD_tmp2['spot_class'] == 'singlet']
                
                list_neiCombUnique = []
                if len(df_RCTD_tmp1) >= min_sample_size and len(df_RCTD_tmp2) >= min_sample_size:
                    
                    list_index2 = df_RCTD_tmp1.index.tolist()
                    df_RCTD_tmp1 = df_RCTD_tmp1.reset_index()
                    val2 = df_RCTD_tmp1['value'][0]
                    
                    list_index3 = df_RCTD_tmp2.index.tolist()
                    df_RCTD_tmp2 = df_RCTD_tmp2.reset_index()
                    val3 = df_RCTD_tmp2['value'][0]
                
                    for j in range(len(list_index2)):
                        df_RCTD_hippo['index_boolean'][list_index2[j]] = 1
                    
                    list_neiCombUnique.append(cell_type_abbrev)
                    list_neiCombUnique.append(abbrev+'+'+abbrev)
                    list_neiCombUnique.append(abbrev2+'+'+abbrev2)
                    
                    for j in range(len(list_index3)):
                        df_RCTD_hippo['index_boolean'][list_index3[j]] = 1    
                    
                    list_matchComb = []
                    for rctdIDX in range(len(df_RCTD_hippo)):
                        if df_RCTD_hippo['index_boolean'][rctdIDX] > 0:
                            val = df_RCTD_hippo['value'][rctdIDX]
                            if val == val1:
                                list_matchComb.append(1)
                            elif val == val2:    
                                list_matchComb.append(2)
                            elif val == val3:    
                                list_matchComb.append(3)
                            
                # elif len(df_RCTD_tmp1) >= min_sample_size and len(df_RCTD_tmp2) < min_sample_size:
                    
                #     list_index2 = df_RCTD_tmp1.index.tolist()
                #     df_RCTD_tmp1 = df_RCTD_tmp1.reset_index()
                #     val2 = df_RCTD_tmp1['value'][0]
                    
                #     for j in range(len(list_index2)):
                #         df_RCTD['index_boolean'][list_index2[j]] = 1
                    
                #     list_neiCombUnique.append(cell_type_abbrev)
                #     list_neiCombUnique.append(abbrev+'+'+abbrev)
                    
                #     list_matchComb = []
                #     for rctdIDX in range(len(df_RCTD)):
                #         if df_RCTD['index_boolean'][rctdIDX] > 0:
                #             val = df_RCTD['value'][rctdIDX]
                #             if val == val1:
                #                 list_matchComb.append(1)
                #             elif val == val2:    
                #                 list_matchComb.append(2)
                                
                # elif len(df_RCTD_tmp1) < min_sample_size and len(df_RCTD_tmp2) >= min_sample_size:    
                    
                #     list_index3 = df_RCTD_tmp2.index.tolist()
                #     df_RCTD_tmp2 = df_RCTD_tmp2.reset_index()
                #     val3 = df_RCTD_tmp2['value'][0]
                    
                #     for j in range(len(list_index3)):
                #         df_RCTD['index_boolean'][list_index3[j]] = 1
                    
                #     list_neiCombUnique.append(cell_type_abbrev)
                #     list_neiCombUnique.append(abbrev2+'+'+abbrev2)    
                    
                #     list_matchComb = []
                #     for rctdIDX in range(len(df_RCTD)):
                #         if df_RCTD['index_boolean'][rctdIDX] > 0:
                #             val = df_RCTD['value'][rctdIDX]
                #             if val == val1:
                #                 list_matchComb.append(1)
                #             elif val == val3:    
                #                 list_matchComb.append(2)   
                
                if len(list_neiCombUnique) > 0:
                    
                    doublet_type.append(cell_type_abbrev)
                    bg_length.append(len(list_neiCombUnique)-1)
                    
                    index_boolean = pd.DataFrame({'index': df_RCTD_hippo['index_boolean'].tolist()})
                    index_boolean.to_csv(path_hippo + 'input_wo_NA/index_hippo_%s.csv' %cell_type_abbrev, index = False, header = None) 
                    
                    neiCombUnique = pd.DataFrame({'neiCombUnique':list_neiCombUnique})
                    neiCombUnique.to_csv(path_hippo + 'input_wo_NA/neiCombUnique_hippo_%s.csv' %cell_type_abbrev, index = False, header = None)
                    
                    matchComb = pd.DataFrame({'matchComb': list_matchComb})
                    matchComb.to_csv(path_hippo + 'input_wo_NA/matchComb_hippo_%s.csv' %cell_type_abbrev, index = False, header = None) 
                    
                    prop = pd.DataFrame({'prop1': list_prop1, 'prop2': list_prop2})
                    prop.to_csv(path_hippo + 'input_wo_NA/prop_hippo_%s.csv' %cell_type_abbrev, index = False, header = None)
                    
# bg_info = pd.DataFrame({'type':doublet_type, 'bg_len':bg_length})
# bg_info.to_csv(path_hippo + 'doublet_types_largerthan30_bg2_certain.csv', index = False, header = None) 
