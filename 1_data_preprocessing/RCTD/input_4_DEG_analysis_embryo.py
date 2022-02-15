####### input for DEG analysis without neighbor analysis #######
import pandas as pd
from collections import Counter
import seaborn as sns
import copy


### Importing data
path_RCTD = '/Users/hyobinkim/tmp/cecilia_data_files/3_celltype_composition/comparison/'
df_RCTD = pd.read_csv(path_RCTD + 'Puck190926_03_RCTD_class.csv', sep = ',', header = 0)
df_name = pd.read_csv(path_RCTD + 'celltype_number_name.csv', sep = ';', header = 0)
df_abbrev =  pd.read_csv(path_RCTD + 'abbreviation.csv', sep = ';', header = 0)

path_prop = '/Users/hyobinkim/RCTD_Plots/mouse_embryo/our_data_whole_genes/'
df_prop = pd.read_csv(path_prop + 'composition_whole.csv', sep = ';', header = 0)
df_prop = df_prop.rename(columns={'Unnamed: 0': 'barcode'})

### Extracting 42362 spots
list_barcode = df_RCTD['barcode'].tolist()
df_prop_rev = copy.deepcopy(df_prop)
for i in range(len(df_prop['barcode'])):
    print(i)
    if df_prop['barcode'][i] not in list_barcode:
        df_prop_rev = df_prop_rev.drop([i])
df_prop_rev = df_prop_rev.reset_index()

# just check if df_prop_rev['barcode'] is equal to df_RCTD['barcode']
tmp=[]
for i in range(len(df_prop_rev['barcode'])):
    if df_prop_rev['barcode'][i] == df_RCTD['barcode'][i]:
        tmp.append(1)
        
# change comma to dot in pandas
df_prop_rev = df_prop_rev.stack().str.replace(',','.').unstack()        

### Changing the value of RCTD['scond_type'] when the bead is a singlet
for i in range(len(df_RCTD['barcode'])):
    if df_RCTD['spot_class'][i] == 'singlet':
        df_RCTD['second_type'][i] = df_RCTD['first_type'][i]

### Adding the proportions of first_type and second_type to RCTD
df_RCTD['prop1'] = 0
df_RCTD['prop2'] = 0
for i in range(len(df_RCTD['barcode'])):
    print(i)
    type1 = str(df_RCTD['first_type'][i])
    type2 = str(df_RCTD['second_type'][i])
    if type1 != type2:
        df_RCTD['prop1'][i] =  df_prop_rev[type1][i]
        df_RCTD['prop2'][i] =  df_prop_rev[type2][i]
    else:
        df_RCTD['prop1'][i] =  df_prop_rev[type1][i]
    
### Assigning a color into each bead
palette = sns.color_palette("hls", 37)
list_color_type1 = []
list_color_type2 = []
for i in range(len(df_RCTD['barcode'])):
    num1 = df_RCTD['first_type'][i] - 1
    list_color_type1.append(palette[num1])
    num2 = df_RCTD['second_type'][i]- 1
    list_color_type2.append(palette[num2])

df_RCTD['color1'] = list_color_type1
df_RCTD['color2'] = list_color_type2

### Adding the info of doublet types
df_RCTD['value'] = df_RCTD.apply(lambda x: hash(frozenset([x["celltype"],x["celltype2"]])),axis=1)

### Sorting the doublet types along the sample size
counter_doublets = Counter(df_RCTD['value'])
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
        df_RCTD_sub = df_RCTD[df_RCTD['value'] == val1]
        # selecting doublet_certain
        #df_RCTD_sub = df_RCTD_sub[df_RCTD_sub['spot_class'] == 'doublet_certain']
        list_index1 = df_RCTD_sub.index.tolist()
        
        if len(list_index1) >= sample_size: # for the number of heterotypic doublets with doublet certain
            idx = list_index1[0];
            if df_RCTD_sub['celltype'][idx] != df_RCTD_sub['celltype2'][idx]:
                
                abbrev = 'na'
                abbrev2 = 'na'
                cellTypeNum1 = 0
                cellTypeNum2 = 0
                for celltypeIDX in range(len(df_abbrev)):
                    if df_RCTD_sub['celltype'][idx] == df_abbrev['cell_type'][celltypeIDX]:
                        abbrev = df_abbrev['abbreviation'][celltypeIDX]
                        cellTypeNum1 = df_RCTD_sub['first_type'][idx]
                    if df_RCTD_sub['celltype2'][idx] == df_abbrev['cell_type'][celltypeIDX]:
                        abbrev2 = df_abbrev['abbreviation'][celltypeIDX]
                        cellTypeNum2 = df_RCTD_sub['second_type'][idx]
                cell_type_abbrev = abbrev + '+' + abbrev2
                
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
                df_RCTD['index_boolean'] = 0
                for j in range(len(list_index1)):
                    df_RCTD['index_boolean'][list_index1[j]] = 1
    
                df_RCTD_tmp1 = df_RCTD[df_RCTD['celltype'] == df_RCTD['celltype'][idx]]
                df_RCTD_tmp1 = df_RCTD_tmp1[df_RCTD_tmp1['celltype2'] == df_RCTD['celltype'][idx]]
                # selecting singlet
                #df_RCTD_tmp1 = df_RCTD_tmp1[df_RCTD_tmp1['spot_class'] == 'singlet']
                
                df_RCTD_tmp2 = df_RCTD[df_RCTD['celltype'] == df_RCTD['celltype2'][idx]]
                df_RCTD_tmp2 = df_RCTD_tmp2[df_RCTD_tmp2['celltype2'] == df_RCTD['celltype2'][idx]]
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
                        df_RCTD['index_boolean'][list_index2[j]] = 1
                    
                    list_neiCombUnique.append(cell_type_abbrev)
                    list_neiCombUnique.append(abbrev+'+'+abbrev)
                    list_neiCombUnique.append(abbrev2+'+'+abbrev2)
                    
                    for j in range(len(list_index3)):
                        df_RCTD['index_boolean'][list_index3[j]] = 1    
                    
                    list_matchComb = []
                    for rctdIDX in range(len(df_RCTD)):
                        if df_RCTD['index_boolean'][rctdIDX] > 0:
                            val = df_RCTD['value'][rctdIDX]
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
                    
                    index_boolean = pd.DataFrame({'index': df_RCTD['index_boolean'].tolist()})
                    index_boolean.to_csv(path_RCTD + 'input_wo_NA/index_embryo_%s.csv' %cell_type_abbrev, index = False, header = None) 
                    
                    neiCombUnique = pd.DataFrame({'neiCombUnique':list_neiCombUnique})
                    neiCombUnique.to_csv(path_RCTD + 'input_wo_NA/neiCombUnique_embryo_%s.csv' %cell_type_abbrev, index = False, header = None)
                    
                    matchComb = pd.DataFrame({'matchComb': list_matchComb})
                    matchComb.to_csv(path_RCTD + 'input_wo_NA/matchComb_embryo_%s.csv' %cell_type_abbrev, index = False, header = None) 
                    
                    prop = pd.DataFrame({'prop1': list_prop1, 'prop2': list_prop2})
                    prop.to_csv(path_RCTD + 'input_wo_NA/prop_embryo_%s.csv' %cell_type_abbrev, index = False, header = None)
                    
bg_info = pd.DataFrame({'type':doublet_type, 'bg_len':bg_length})
bg_info.to_csv(path_RCTD + 'doublet_types_largerthan30_bg2_certain.csv', index = False, header = None) 

#%%
####### Plots for spatial maps #######   
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
  
path_input = '/Users/hyobinkim/Documents/MATLAB/FindDEGs/input_data/singlets_no_isolated/csv_files/'
df_index = pd.read_csv(path_input + 'index_NA_singlets_13.csv', sep = ',', header = None, names=['index'])
df_matchComb = pd.read_csv(path_input + 'matchComb_NA_singlets_13.csv', sep = ',', header = None, names=['matchComb'])
df_neiCombUnique = pd.read_csv(path_input + 'neiCombUnique_NA_singlets_13.csv', sep = ',', header = None, names=['neiCombUnique'])
           
### Visualization  
list_matchComb = []
flag = 0
for i in range(len(df_index)):
    
    idx_boolean = df_index['index'][i]    
    if idx_boolean == 1:
        list_matchComb.append(df_matchComb['matchComb'][flag])
        flag += 1
    else:
        list_matchComb.append(0)

df_RCTD['matchComb']= list_matchComb  

path_ex = '/Users/hyobinkim/Documents/MATLAB/FindDEGs/results/heatmap/En+Sc/'
df_LSc = pd.read_csv(path_ex + 'L+Sc_Hbb-y.txt', sep = ',', header = None, names=['barcode','log-val','z-val'])
df_EnSc = pd.read_csv(path_ex + 'En+Sc_Hbb-y.txt', sep = ',', header = None, names=['barcode','log-val','z-val'])

df_LSc = df_LSc.sort_values(by='log-val', ascending=False)
df_LSc = df_LSc.reset_index() 
list_LSc_barcode = df_LSc['barcode'].tolist()
list_LSc_barcode = list_LSc_barcode[:74]

df_RCTD['LSc'] = 0
for i in range(len(df_RCTD)):   
    if df_RCTD['barcode'][i] in list_LSc_barcode:
        df_RCTD['LSc'][i] = 1
            
comb = 10
df_RCTD_sub = df_RCTD[df_RCTD['LSc']==1]     
df_RCTD_sub = df_RCTD_sub.reset_index() 

fig, ax = plt.subplots()
fig.set_size_inches(30.,18.)

# c=df_RCTD_sub['color1']
# c=df_RCTD_sub['color2']
ax.scatter(df_RCTD_sub['x'], df_RCTD_sub['y'], s=200, c=df_RCTD_sub['matchComb'], edgecolor="black", marker=MarkerStyle("o", fillstyle="right"))
ax.scatter(df_RCTD_sub['x'], df_RCTD_sub['y'], s=200, c=df_RCTD_sub['matchComb'], edgecolor="black", marker=MarkerStyle("o", fillstyle="left"))
celltype = df_RCTD_sub['celltype'][0]+'+'+df_RCTD_sub['celltype2'][0]
matchCombName = df_neiCombUnique['neiCombUnique'][comb-1]
ax.set_title(celltype, fontsize=40)
fig.savefig('/Users/hyobinkim/tmp/cecilia_data_files/3_celltype_composition/comparison/RCTD_spatial_maps/%d_%s_L+Sc74.png' %(comb,matchCombName), dpi=300)   
plt.show()

# df_RCTD_sub = df_RCTD[df_RCTD['matchComb'] != 0]        
# df_RCTD_sub = df_RCTD_sub.reset_index()

# fig, ax = plt.subplots()
# fig.set_size_inches(30.,18.)

# ax.scatter(df_RCTD_sub['x'], df_RCTD_sub['y'], s=200, c=df_RCTD_sub['matchComb'], edgecolor="black", marker=MarkerStyle("o", fillstyle="right"))
# ax.scatter(df_RCTD_sub['x'], df_RCTD_sub['y'], s=200, c=df_RCTD_sub['matchComb'], edgecolor="black", marker=MarkerStyle("o", fillstyle="left"))
# celltype = df_RCTD_sub['celltype'][0]+'+'+df_RCTD_sub['celltype2'][0]
# ax.set_title(celltype, fontsize=40)
# fig.savefig('/Users/hyobinkim/tmp/cecilia_data_files/3_celltype_composition/comparison/RCTD_spatial_maps/all.png', dpi=300) 
# plt.show()



        