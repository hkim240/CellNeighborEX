import pandas as pd
from collections import Counter


### Importing data
path = '/Users/kimh15/Downloads/github_scripts/Python/categorization/'

# Choose data: seqFISH_embryo.csv, Slide-seq_embryo.csv, Slide-seq_hippocampus.csv, Slide-seq_liver_cancer.csv
df = pd.read_csv(path + 'input/Slide-seq_hippocampus.csv', header = 0) 


### Coverting doublet types to unique numbers and adding the nubmers to df

# For seqFISH dataset
# df['integrated'] = 'NA'
# for idx in range(len(df['integrated'])):
#     df['integrated'][idx] = df['celltype1'][idx] + '+' + df['celltype2'][idx]

# df['value'] = df.apply(lambda x: hash(frozenset([x["integrated"]])),axis=1)
  
# For Slide-seq datasets
df['value'] = df.apply(lambda x: hash(frozenset([x["celltype1"],x["celltype2"]])),axis=1)


### Sorting the doublet types along the sample size
counter_doublets = Counter(df['value'])
sorted_counter = sorted(counter_doublets.items(), key=lambda x: x[1], reverse=True)


### Generating input files for DEG analysis
sample_size = 30 # cut-off for the number of heteorypic spots
min_sample_size = 1
doublet_type = []
for i in range(len(sorted_counter)):
    
    print(i)
    if sorted_counter[i][1] >= sample_size:
        
        val1 = sorted_counter[i][0]
        df_sub = df[df['value'] == val1]
        list_index1 = df_sub.index.tolist()
        
        if len(list_index1) >= sample_size: 
            
            idx = list_index1[0];
            if df_sub['celltype1'][idx] != df_sub['celltype2'][idx]:
                
                abbrev = df_sub['celltype1'][idx]
                abbrev2 = df_sub['celltype2'][idx]
                cellTypeNum1 = df_sub['first_type'][idx]
                cellTypeNum2 = df_sub['second_type'][idx]
                cell_type_abbrev = abbrev + '+' + abbrev2
                
                # input: prop
                list_prop1 = []
                list_prop2 = []
                for xx in range(len(df_sub)):
                    index = list_index1[xx]
                    if df_sub['first_type'][index] == cellTypeNum1:
                        list_prop1.append(df_sub['prop1'][index])
                        list_prop2.append(df_sub['prop2'][index])
                    elif df_sub['first_type'][index] == cellTypeNum2:
                        list_prop1.append(df_sub['prop2'][index])
                        list_prop2.append(df_sub['prop1'][index])
                    
                # input: index, neiCombUnique, matchComb    
                df['index_boolean'] = 0
                for yy in range(len(list_index1)):
                    df['index_boolean'][list_index1[yy]] = 1
    
                df_tmp1 = df[df['celltype1'] == df['celltype1'][idx]]
                df_tmp1 = df_tmp1[df_tmp1['celltype2'] == df['celltype1'][idx]]
                
                df_tmp2 = df[df['celltype1'] == df['celltype2'][idx]]
                df_tmp2 = df_tmp2[df_tmp2['celltype2'] == df['celltype2'][idx]]
                
                list_neiCombUnique = []
                if len(df_tmp1) >= min_sample_size and len(df_tmp2) >= min_sample_size:
                    
                    list_index2 = df_tmp1.index.tolist()
                    df_tmp1 = df_tmp1.reset_index()
                    val2 = df_tmp1['value'][0]
                    
                    list_index3 = df_tmp2.index.tolist()
                    df_tmp2 = df_tmp2.reset_index()
                    val3 = df_tmp2['value'][0]
                
                    for zz in range(len(list_index2)):
                        df['index_boolean'][list_index2[zz]] = 1
                    
                    list_neiCombUnique.append(cell_type_abbrev)
                    list_neiCombUnique.append(abbrev+'+'+abbrev)
                    list_neiCombUnique.append(abbrev2+'+'+abbrev2)
                    
                    for ww in range(len(list_index3)):
                        df['index_boolean'][list_index3[ww]] = 1    
                    
                    list_matchComb = []
                    for rctdIDX in range(len(df)):
                        if df['index_boolean'][rctdIDX] > 0:
                            val = df['value'][rctdIDX]
                            if val == val1:
                                list_matchComb.append(1)
                            elif val == val2:    
                                list_matchComb.append(2)
                            elif val == val3:    
                                list_matchComb.append(3)
                            
                if len(list_neiCombUnique) > 0:
                    
                    doublet_type.append(cell_type_abbrev)
                    
                    index_boolean = pd.DataFrame({'index': df['index_boolean'].tolist()})
                    index_boolean.to_csv(path + 'output/hippocampus/index_hippocampus_%s.csv' %cell_type_abbrev, index = False, header = None) 
                    
                    neiCombUnique = pd.DataFrame({'neiCombUnique':list_neiCombUnique})
                    neiCombUnique.to_csv(path + 'output/hippocampus/neiCombUnique_hippocampus_%s.csv' %cell_type_abbrev, index = False, header = None)
                    
                    matchComb = pd.DataFrame({'matchComb': list_matchComb})
                    matchComb.to_csv(path + 'output/hippocampus/matchComb_hippocampus_%s.csv' %cell_type_abbrev, index = False, header = None) 
                    
                    prop = pd.DataFrame({'prop1': list_prop1, 'prop2': list_prop2})
                    prop.to_csv(path + 'output/hippocampus/prop_hippocampus_%s.csv' %cell_type_abbrev, index = False, header = None)
                    
