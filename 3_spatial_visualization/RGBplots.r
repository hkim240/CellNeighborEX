
library(ggplot2)
library(Seurat)
library(dplyr)
library(grid)

hippocampus=TRUE

source("plot_functions.r")

if (hippocampus){
  this_folder=""
  coords_file=""
  save_folder=""
  DGE_file=""
  count_file <- readRDS("")
  
}else{
  
  DGE_file=""
  this_folder=""
  save_folder=""
  coords_file=""
  count_file<- readRDS("")
  
}


this_value_type="z"

DGE=read.csv(DGE_file,sep=";")

DGE$comb=paste(DGE$Cell1,DGE$Cell2,sep="+")

DGEu=distinct(DGE,comb,.keep_all = TRUE)

DGE_table=as.data.frame(table(c(DGEu$Cell1,DGEu$Cell2)))

print(as.character(DGE_table$Var1[DGE_table$Freq>2]))

cells_use=as.character(DGE_table$Var1[DGE_table$Freq>2])



coords=readRDS(coords_file)

if (hippocampus){
  coords$cells=rownames(coords)
}


setwd(this_folder)

this_folders=list.files(this_folder)

p_theme=pres_theme("Arial",16)

cell_number=1

for (this_ct in cells_use){
  

  
  cell_number=cell_number+1
  
  curr_folders=c(this_folders[grep(paste("+",this_ct,sep=""),this_folders,fixed=TRUE)],this_folders[grep(paste(this_ct,"+",sep=""),this_folders,fixed=TRUE)])
  
  curr_DGE=subset(DGE,comb %in% curr_folders)
  curr_DGE=subset(curr_DGE,(Cell1 %in% this_ct) | (Cell2 %in% this_ct) )
  
  unique_pairs=unique(curr_DGE$comb)
  
  if (length(unique_pairs)>3){
    unique_pairs=unique_pairs[1:3] #1:3
  }
  
  
  if(length(unique_pairs)>2){
    comb=t(combn(unique_pairs,3,simplify=TRUE))
    
    
    
    unique_cells=comb
  }else{
    unique_cells=unique_pairs
  }
  
  first=TRUE
  
  color_genes=c()
  
  for (folder in  unique_cells){
    
    other_cells=subset(curr_DGE,!(comb %in% folder))
    
    
    print(folder)
    
    curr_folder=paste(this_folder,folder,sep="")  
    
    setwd(curr_folder)
    
    
    files2=list.files(curr_folder)
    
    
    
    this_files=divide_files(files2)
    
    available_genes=names(this_files)
    available_genes=available_genes[which(!(available_genes=="NA"))]
    
    
    
    
    if (length(available_genes[!(available_genes %in% other_cells$Gene)])>0){
      available_genes=available_genes[!(available_genes %in% other_cells$Gene)]
    }
    
    
    gene_name=available_genes[1]
    
    
    gene_file=this_files[[gene_name]]
    
    
    if (length(files2)>0){  
      data_set=importDGE(gene_file,coords)
      
    }
    
    if (first){
      res_data=data_set
      comb_data=data_set
   
      
      curr_rows=rownames(comb_data)[rownames(comb_data) %in% colnames(count_file)]
      
      comb_data[curr_rows,gene_name]=count_file@assays$RNA@data[gene_name,curr_rows]
      
      comb_data[rownames(data_set),gene_name]=data_set$log
      color_genes=c(color_genes,gene_name)
      
      first=FALSE
    }else{
      
      
      
      res_data=rbind(res_data,data_set)
      
      
      extra_cols=colnames(comb_data)[!(colnames(comb_data) %in% colnames(data_set))]
      
      for (col in extra_cols){
        data_set[,col]=rep(0,nrow(data_set))
        
      }
      
      
      
      comb_data=rbind(comb_data,data_set[!(rownames(data_set) %in% rownames(comb_data)),])
      curr_rows=rownames(comb_data)[rownames(comb_data) %in% colnames(count_file)]
      
      comb_data[curr_rows,gene_name]=count_file@assays$RNA@data[gene_name,curr_rows]
      
   
      curr_rows=rownames(data_set)[rownames(data_set) %in% colnames(count_file)]
      
      comb_data[curr_rows,gene_name]=count_file@assays$RNA@data[gene_name,curr_rows]
      color_genes=c(color_genes,gene_name)
      
      
      
      
    }
    
  }
  

  comb_data=comb_data[!rownames(comb_data) %in% c("NA","NA1","NA.1"),]

  color_genes=unique(color_genes)
  
  if (length(color_genes)<3){
    comb_data$Rcolor=rgb(comb_data[,color_genes[1]]/max(comb_data[,color_genes[1]]),comb_data[,color_genes[2]]/max(comb_data[,color_genes[2]])
                         ,0)
    
  }else{
    
    comb_data$Rcolor=rgb(comb_data[,color_genes[1]]/max(comb_data[,color_genes[1]]),comb_data[,color_genes[2]]/max(comb_data[,color_genes[2]])
                         ,comb_data[,color_genes[3]]/max(comb_data[,color_genes[3]]))
    
    
  }
  
  
  
  
  res_data=subset(comb_data,db %in% c("heterotypic"))
  res_data=rbind(res_data,subset(comb_data,celltype %in% this_ct))

  
  this_labels=unique(res_data$celltype)
  
  
  res_data$celltype= factor(res_data$celltype,levels=this_labels)
  
  this_doub=which(grepl("[+]",this_labels)==TRUE)
  this_sing=which(grepl("[+]",this_labels)==FALSE)
  
  this_colors=rep("",length(this_labels))
  this_colors[this_doub]=c("red","green","blue")
  available_cols=c("darkorange","pink","gold2","brown","cyan","purple")
  this_colors[this_sing]=available_cols[seq(1,length(this_sing))]
  
  
  
  
 
  extra_data=coords[!(rownames(coords) %in% rownames(res_data)),]
  
  
  
  grob <- grobTree(textGrob(color_genes[1], x=0.0,  y=0.95, hjust=0,
                            gp=gpar(col="red",fontsize=14)),
                   textGrob(color_genes[2], x=0.0,  y=0.9, hjust=0,
                            gp=gpar(col="green",fontsize=14))
                    ,textGrob( color_genes[3], x=0.0,  y=0.85, hjust=0,
                      gp=gpar(col="blue",fontsize=14))
                   
  )
  
  this_labels=fix_label(this_labels,this_ct)
  
  
  gp1=ggplot(extra_data,aes(x,y))+geom_point(stroke=0.5,shape=21,colour="grey",fill="grey",size=1) +
    geom_point(data=res_data,aes(x,y,fill=Rcolor,colour=celltype,size=db),stroke=1,shape=21)+
    scale_size_manual(values=c(3,1.5))+p_theme+scale_fill_identity(guide = "legend",breaks=c("red","blue","green"))+
    scale_color_manual(values=this_colors,labels=this_labels)+
    labs(size="Composition", colour="Cell type",fill="Gene expression")+
    guides(colour=guide_legend(override.aes = list(size =3,stroke=2)),size="none")+
    annotation_custom(grob)
  
  
  
  print(gp1)
  
  homo_data=subset(comb_data,db %in% c("homotypic"))
  
  pre_name=paste(unique(homo_data$celltype),collapse = "_")
  pre_gene=paste(unique(homo_data$gene),collapse="_")
  file_name=paste(save_folder,pre_name,"_",pre_gene,".png",sep="")
  zoom_name=paste(save_folder,pre_name,"_",pre_gene,"_zoom.png",sep="")
  if(hippocampus){
       z_sp=gp1+p_theme +scale_x_continuous(limits=c(0,3000))+scale_y_continuous(limits=c(4000,6000))
  }else{ 
    z_sp=gp1+p_theme +scale_x_continuous(limits=c(4000,5500))+scale_y_continuous(limits=c(2000,4000))
  }
  
  #ggsave(file_name,plot=gp1)
  
  #ggsave(zoom_name,plot=z_sp)
}


