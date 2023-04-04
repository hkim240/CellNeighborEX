import pandas as pd

path = '/Users/kimh15/Downloads/github_scripts/Matlab/'
df_top2000 = pd.read_csv(path + 'input_logdata/gene_name_top_plus.txt', sep = ' ', names=['gene'], header=None)
df_markers = pd.read_csv(path + 'input_markers/liver_cell_type_markers.csv', sep = ',', header = 0)
    

total_markers = []
for j in range(len(df_markers)):    
    Markers = []
    for w in range(2): # number of markers for each cell type
        col = 'Marker'+str(w+1)
        #marker = ''.join(df_markers[col][j].split()) # embryo (because of space '')
        marker = df_markers[col][j] # hippocampus
        if marker != '-':
            Markers.append(marker)
    
    total_markers.append(Markers)
        
        
list_marker_check = []
for i in range(len(df_top2000)):
    geneCellType = []
    # gene = ''.join(df_top2000['gene'][i].split()) # embryo(in the case that gene name has space)
    gene = df_top2000['gene'][i] 
    #print(i)
    flag = 0
    
    for j in range(len(total_markers)):
        if gene in total_markers[j]:
            cellType = df_markers['Celltype'][j]+'+'+df_markers['Celltype'][j]
            #print(cellType)
            geneCellType.append(cellType)
            flag = 1
                
    if flag == 0:
        geneCellType.append('None')
    
    list_marker_check.append(geneCellType)

max_len = 1
for x in range(len(list_marker_check)):
    if len(list_marker_check[x]) > max_len:
        max_len = len(list_marker_check[x])

total_matchMarker = []
for i in range(max_len):
    tmp = []
    for j in range(len(list_marker_check)):
        if i < len(list_marker_check[j]):
            tmp.append(list_marker_check[j][i])
        else:
            tmp.append('None')
        
    total_matchMarker.append(tmp)           

                
matchMarker = pd.DataFrame({'list1': total_matchMarker[0], 'list2': total_matchMarker[1]})
matchMarker.to_csv(path + 'input_markers/output/liver_matchMarker_plus_selected.csv', index = False, header = None)
                                
            
