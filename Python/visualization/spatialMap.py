import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


#### Visualize spatial gene expression 
# example in seqFISH data:
heteroR = 'Gut-tube+Neural-crest'
heteroB = 'Gut-tube+Gut-tube'
homoY = 'Neural-crest' 
geneR = 'Pitx1' 

# example in liver cancer Slide-seq data:
# heteroR = 'TumorIII+Monocyte'
# heteroB = 'TumorIII+TumorIII'
# homoY = 'Monocyte+Monocyte' 
# geneR = 'F13a1' 

# example in hippocampus Slide-seq data:
# heteroR = 'EnT+A'
# heteroB = 'EnT+EnT'
# homoY = 'A+A'
# geneR = 'Fabp7'


#### Import data
# expression data
path_exp = '../Python/visualization/'
dfolder = 'expression_data/Pitx1/'

df_spatial_total = pd.DataFrame()
for filename in os.listdir(path_exp + dfolder):
    
    df_spatial = pd.DataFrame()
    if filename != '.DS_Store':
        
        print(filename)
        df_temp = pd.read_csv(path_exp + dfolder + filename, sep = ',', header = 0)
        df_spatial['barcode'] = df_temp['barcode']
        df_spatial['type'] = filename[:-4] # remove ".csv"
        df_spatial['zscore'] = df_temp['zscore']
        
    df_spatial_total = df_spatial_total.append(df_spatial, ignore_index=True)

# cell or bead location data
path_df = '../Python/categorization/input/'

# seqFISH dataset 
df = pd.read_csv(path_df + 'seqFISH_embryo.csv', sep = ',', header = 0)

# Slide-seq in liver cancer
#df = pd.read_csv(path_df + 'Slide-seq_liver_cancer.csv', sep = ',', header = 0)

# Slide-seq in hippocampus
#df = pd.read_csv(path_df + 'Slide-seq_hippocampus.csv', sep = ',', header = 0)


#### Set up the default colors of beads 
df['beadsize'] = 5 #seqFISH:5, Slide-seq:10
df['edgecolor'] = 'NA'
df['beadcolor'] = 'NA'
df['type'] = 'NA'
df['expression'] = 1.2 # The colors of cells or beads (background) with gray edges are set to white
                       # seqFISH:1.2, Slide-seq liver cancer:1.5, Slide-seq hippocampus:2.7
for i in range(len(df)):
    df['edgecolor'][i] = (0.85,0.85,0.85) # gray
    df['beadcolor'][i] = (0.85,0.85,0.85)
     
                 
#### Set up the colors of neighbor specific cells or beads
for idx in range(len(df)): 
    print(idx)
    for idx2 in range(len(df_spatial_total)):
        
        if df['barcode'][idx] == df_spatial_total['barcode'][idx2]:
            
            df['expression'][idx] = df_spatial_total['zscore'][idx2]
            
            if df_spatial_total['type'][idx2] == heteroR:
                df['edgecolor'][idx] = (1,0,0) # red
                df['type'][idx] = heteroR
                df['beadsize'][idx] = 60 #seqFISH:60, Slide-seq:600
                
            elif df_spatial_total['type'][idx2] == heteroB:
                df['edgecolor'][idx] = (0,0,1) # blue
                df['type'][idx] = heteroB
                df['beadsize'][idx] = 30 #seqFISH:30, Slide-seq:200
                
            elif df_spatial_total['type'][idx2] == homoY:
                df['edgecolor'][idx] = (0,0,0) # black
                df['type'][idx] = homoY
                df['beadsize'][idx] = 30 #seqFISH:30, Slide-esq:200


#### Visualization 
df_heteroR = pd.DataFrame()
df_heteroB = pd.DataFrame()   
df_homoY = pd.DataFrame()

df_heteroR = df[df['type'] == heteroR]
df_heteroR = df_heteroR.reset_index()
df_heteroB = df[df['type'] == heteroB]
df_heteroB = df_heteroB.reset_index()
df_homoY = df[df['type'] == homoY]
df_homoY = df_homoY.reset_index()


#### Plot spatial maps
fig, ax = plt.subplots()
fig.set_size_inches(20,28) # seqFISH:(20,28), Slide-seq:(28,28)

# seqFISH 
labels = ['Gut tube/Neural crest', 'Gut tube/Gut tube', 'Neural crest']

# Slide-seq liver cancer
#labels = ['TumorIII+Monocyte', 'TumorIII', 'Monocyte']

# Slide-seq hippocampus
#labels = ['Endothelial_Tip+Astrocyte', 'Endothelial_Tip', 'Astrocyte']

# Gray - background
# For seqFISH dataset, add '-' to the y coordinates
# Adjust the values of linewidths and zorder depending on datasets
cbar_flag=ax.scatter(df['x'], -df['y'], s=df['beadsize'], c=df['expression'].to_numpy(), cmap='bwr', edgecolor=df['edgecolor'], linewidths=2, zorder=2)

# Red
ax.scatter(df_heteroR['x'], -df_heteroR['y'], s=df_heteroR['beadsize'], c=df_heteroR['expression'].to_numpy(), cmap='bwr', edgecolor=df_heteroR['edgecolor'], linewidths=2, zorder=4.5)
ax.scatter(df_heteroR['x'], -df_heteroR['y'], s=df_heteroR['beadsize'], c='white', edgecolor=df_heteroR['edgecolor'], linewidths=2, label=labels[0])

# Blue
ax.scatter(df_heteroB['x'], -df_heteroB['y'], s=df_heteroB['beadsize'], c=df_heteroB['expression'].to_numpy(), cmap='bwr', edgecolor=df_heteroB['edgecolor'], linewidths=2, zorder=3.5)
ax.scatter(df_heteroB['x'], -df_heteroB['y'], s=df_heteroB['beadsize'], c='white', edgecolor=df_heteroB['edgecolor'], linewidths=2, label=labels[1]) 

# Black
ax.scatter(df_homoY['x'], -df_homoY['y'], s=df_homoY['beadsize'], c=df_homoY['expression'].to_numpy(), cmap='bwr', edgecolor=df_homoY['edgecolor'], linewidths=2, zorder=4.0)
ax.scatter(df_homoY['x'], -df_homoY['y'], s=df_homoY['beadsize'], c='white', edgecolor=df_homoY['edgecolor'], linewidths=2, label=labels[2])

# Legend for Slide-seq data
ax.legend(bbox_to_anchor=(0.98,0.42), prop={'size': 30}, frameon=False, fontsize=30) # seqFISH:30, Slide-seq:60
legend = ax.get_legend()
legend.set_title('Cell type', prop={'size':30}) # seqFISH:30, Slide-seq:60
legend._legend_box.align = "left"
legend.legendHandles[0]._sizes = [80] # seqFISH:80, Slide-seq:1600
legend.legendHandles[1]._sizes = [40] # seqFISH:40, Slide-seq:800
legend.legendHandles[2]._sizes = [40] # seqFISH:40, Slide-seq:800
   
# Spatial plot
ax.axis('off')
ax.xaxis.set_major_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.NullLocator())

cbar = plt.colorbar(cbar_flag, cax = fig.add_axes([0.91, 0.45, 0.03, 0.2]),ticks=[-0.5,0,1,2,3,4,5]) 
cbar.ax.tick_params(labelsize=30) # seqFISH:30, Slide-seq:60
cbar.ax.set_title('%s \n expression' %geneR, fontsize=30) # seqFISH:30, Slide-seq:60

plt.show()
#fig.savefig(path_exp + '/output/%s.png' %geneR, dpi=300)





