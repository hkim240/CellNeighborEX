import squidpy as sq
import pandas as pd
from collections import Counter


#### Load the pre-processed dataset: seqFISH in Slide mouse embryo 1
path = '/Users/kimh15/Downloads/github_scripts/Python/seqFISH_neighbor/'
adata = sq.datasets.seqfish() # 19416 cells x 351 genes


#### Create dataframe
# barcode
df_barcode = pd.DataFrame(adata.obs.index, columns=["barcode"])

# cell type
df_celltype = adata.obs.set_index(df_barcode.index)
df_celltype = df_celltype.drop(columns = "Area")
df_celltype = df_celltype.rename(columns={"celltype_mapped_refined":"celltype1"})

# x, y coordinates
xVal_list = []
yVal_list = []
for idx in range(len(adata.obsm["spatial"])):
    xVal = adata.obsm["spatial"][idx][0]
    yVal = adata.obsm["spatial"][idx][1]
    xVal_list.append(xVal)
    yVal_list.append(yVal)

df_coord = pd.DataFrame(list(zip(xVal_list, yVal_list)), columns=["x", "y"])

# merge sub-dataframes
frames = [df_barcode, df_celltype, df_coord]
df = pd.concat(frames, axis=1)

# add dataframe of "first_type", "second_type","celltype2","prop1","prop2"
df['first_type'] = -1
df['second_type'] = -1
df['celltype2'] = 'NA'
df['prop1'] = 1
df['prop2'] = 0


#### Find immediate neighbors
# detect neighbors using delaunay triangulation
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True) 

# calculate number of neighbors
matrix = adata.obsp["spatial_connectivities"]
neiNum = Counter(matrix.tocoo().row)
df['neiNum'] = 0
for idx in range(len(df)):
    df['neiNum'][idx] = neiNum[idx]
    
# identify neighboring cell type
df['neiType'] = 'NA'
df['catNum'] = 0
cumulativeIdx = 0
for idx in range(len(df)):
    print(idx)
    
    temp_neiType = []
    for idx2 in range(df['neiNum'][idx]):
        neiIdx = matrix.tocoo().col[cumulativeIdx+idx2]    
        temp_neiType.append(df['celltype1'][neiIdx]) 
    df['neiType'][idx] = temp_neiType 
    
    cat_num = 0
    categoryNumber = len(Counter(df['neiType'][idx]).keys()) 
    if df['celltype1'][idx] in Counter(df['neiType'][idx]).keys():
        cat_num = categoryNumber - 1
    else: 
        cat_num = categoryNumber
    df['catNum'][idx] = cat_num  

    cumulativeIdx = cumulativeIdx + df['neiNum'][idx]

# export dataframe
df.to_csv(path + 'seqFISH_neiInfo.csv', index=False)


#### Import seqFISH_neiInfo data 
df = pd.read_csv(path +'output/seqFISH_neiInfo.csv') 

      
#### Define a function of annotating cluster numbers
#df["celltype1"] = pd.Categorical(df["celltype1"])
#values = df["celltype1"].cat.codes
#values.unique()
#array([13, 12,  7,  0,  9,  6, 10, 11, 14, 19, 16, 20,  8,  3, 21,  4, 15, 1, 17,  5,  2, 18], dtype=int8)
#celltypes = df["celltype1"].cat.categories
#celltypes
# Index(['Allantois', 'Anterior somitic tissues', 'Cardiomyocytes', 'Cranial mesoderm', 'Definitive endoderm', 
# 'Dermomyotome', 'Endothelium', 'Erythroid', 'Forebrain/Midbrain/Hindbrain', 'Gut tube', 'Haematoendothelial progenitors', 
# 'Intermediate mesoderm', 'Lateral plate mesoderm', 'Low quality', 'Mixed mesenchymal mesoderm', 'NMP', 'Neural crest', 
# 'Presomitic mesoderm', 'Sclerotome', 'Spinal cord', 'Splanchnic mesoderm', 'Surface ectoderm'], dtype='object')

def annotate_num(ct):
    
    val = -1
    if ct == 'Allantois':
        val = 13
    elif ct == 'Anterior somitic tissues':
        val = 12
    elif ct == 'Cardiomyocytes':
        val = 7
    elif ct == 'Cranial mesoderm':
        val = 0
    elif ct == 'Definitive endoderm':
        val = 9
    elif ct == 'Dermomyotome':
        val = 6
    elif ct == 'Endothelium':
        val = 10
    elif ct == 'Erythroid':
        val = 11
    elif ct == 'Forebrain/Midbrain/Hindbrain':
        val = 14
    elif ct == 'Gut tube':
        val = 19
    elif ct == 'Haematoendothelial progenitors': 
        val = 16
    elif ct == 'Intermediate mesoderm':
        val = 20
    elif ct == 'Lateral plate mesoderm':
        val = 8
    elif ct == 'Low quality':
        val = 3
    elif ct == 'Mixed mesenchymal mesoderm':
        val = 21
    elif ct == 'NMP':
        val = 4
    elif ct == 'Neural crest':
        val = 15
    elif ct == 'Presomitic mesoderm':
        val = 1
    elif ct == 'Sclerotome':
        val = 17
    elif ct == 'Spinal cord':
        val = 5
    elif ct == 'Splanchnic mesoderm':
        val = 2
    elif ct == 'Surface ectoderm':
        val = 18                                
    
    return val


#### Extract cells with the only one neighboring cell type: catNum==0: itself, catNum==1: another cell type
df_filtered = df[df['neiNum'] > 0]
df_filtered = df_filtered[(df_filtered['catNum'] == 0) | (df_filtered['catNum'] == 1)]
df_filtered = df_filtered.reset_index(drop=True)


#### Fill out 'celltype2', 'first_type', 'second_type'
for idx in range(len(df_filtered)):
    
    print(idx)

    df_filtered['first_type'][idx] = annotate_num(df_filtered['celltype1'][idx])
    if df_filtered['catNum'][idx] == 0: # cell surrounded by itself cell type 
        df_filtered['celltype2'][idx] = df_filtered['celltype1'][idx]
    
    else:
        # neighboring cell type
       temp = df_filtered['neiType'][idx].split("'")
       idxFlag = 0
       while(True):
            if temp[idxFlag] != df_filtered['celltype1'][idx] and temp[idxFlag] != ', ' and temp[idxFlag] != '[' and temp[idxFlag] != ']':
                df_filtered['celltype2'][idx] = temp[idxFlag]
                break
            idxFlag = idxFlag + 1
        
    df_filtered['second_type'][idx] = annotate_num(df_filtered['celltype2'][idx])


#### remove space and modify the names of cell types
df_filtered['celltype1'] = df_filtered['celltype1'].str.replace(' ', '-')
df_filtered['celltype1'] = df_filtered['celltype1'].str.replace('/', '-')
df_filtered['celltype2'] = df_filtered['celltype2'].str.replace(' ', '-')
df_filtered['celltype2'] = df_filtered['celltype2'].str.replace('/', '-')


#### Export output
df_filtered.to_csv(path +'output/seqFISH_embryo.csv', index=False) 





