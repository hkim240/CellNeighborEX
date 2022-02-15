pres_theme<-function(this_family,this_size){
  res_theme=theme(
    legend.title = element_text(family=this_family,size=this_size),
    text=element_text(family=this_family,size=this_size),
    legend.text=element_text(family=this_family,size=this_size),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    legend.background = element_rect(fill = "white", colour = "white"),
    legend.key = element_rect(fill = "white", colour = "white")
    
  )
}


divide_files=function(files){
  
  this_list=list()
  
  for (file1 in files){
    
    curr_items=unlist(strsplit(file1,"_"))
    gene=gsub(".txt","",curr_items[2])
    
    if (!(grepl("pdf",gene))){
      this_list[[gene]]=c(this_list[[gene]],file1)
    }
  }
  return(this_list)
}


importDGE=function(this_files,this_coords){
  
  res_data=data.frame()
  
  for (file1 in this_files){
    
    data1=read.table(file1,sep=",")
    
    colnames(data1)=c("barcode","log","z")
    data1$log=as.numeric(data1$log)
    curr_extra=this_coords[data1$barcode,c("cells","x","y")]
    data1=cbind(data1,curr_extra)
    
    
    
    curr_items=unlist(strsplit(file1,"_"))
    celltypes=unlist(strsplit(curr_items[1],"[+]"))
    celltype1=celltypes[1]
    celltype2=celltypes[2]
    gene=gsub(".txt","",curr_items[2])
    
    
    data1$gene=rep(gene,nrow(data1))
    data1$celltype1=rep(celltype1,nrow(data1))
    data1$celltype2=rep(celltype2,nrow(data1))
    
    
    
    if (nrow(res_data)>0){
      
      res_data=rbind(res_data,data1)
    }else{
      res_data=data1
    }
    
  }
  
  res_data$db=rep("homotypic",nrow(res_data))
  
  ht=which((res_data$celltype1==res_data$celltype2)==FALSE)
  
  res_data$db[ht]=rep("heterotypic",length(ht))
  res_data$celltype=res_data$celltype1
  res_data$celltype[ht]=paste(res_data$celltype1[ht],"+",res_data$celltype2[ht],sep="")
  
  
  return(res_data)
}


fix_label=function(this_lab,this_ct){
  res_lab=c()
  for (lb in this_lab){
    
    curr_lb=as.character(unlist(strsplit(lb,"[+]")))
    curr_lb=unique(c(this_ct,curr_lb))
    curr_lb=paste(curr_lb,collapse = "+")
    res_lab=c(res_lab,curr_lb)
  }
  
  
  return(res_lab)
}


spatial_plot=function(plot_data,this_gene,this_value,this_coords,is_small,this_cols,this_dict){
  
  colnames(plot_data)[which(colnames(plot_data)==this_value)]="value"
  
  extra_data=this_coords[!(rownames(this_coords) %in% rownames(plot_data)),]
  
  this_mid=max(plot_data$value)-((abs(max(plot_data$value))+abs(min(plot_data$value)))/2)
  
  this_labels=unique(plot_data$celltype)
  
  this_doub=which(grepl("[+]",this_labels)==TRUE)
  this_sing=which(grepl("[+]",this_labels)==FALSE)
  
  sing_labels=mapvalues(this_labels[this_sing],from=change_abbr$to,to=change_abbr$from,warn_missing = FALSE)
  doub_labels=unlist(strsplit(gsub("[+]","_",this_labels[this_doub]),"_"))
  
  doub_labels=mapvalues(doub_labels,from=change_abbr$to,to=change_abbr$from,warn_missing = FALSE)
  doub_labels=mgsub(doub_labels,"\\s+","\n")
  doub_labels=paste(doub_labels,collapse = "+\n")
 
  old_labels=c(this_labels[this_doub],this_labels[this_sing])
  this_labels=c(doub_labels,sing_labels)
  
  
  this_celltype=unique(as.character(plot_data$celltype))
  
  if (is_small){
    
    gene.plot=ggplot(extra_data,aes(x,y))+geom_point(stroke=0.5,shape=21,colour="grey",fill="grey",size=1) 
    
    gene.plot=gene.plot+geom_point(data=plot_data,aes(x,y,colour=celltype,size=db,fill=celltype),stroke=0.5,shape=21)+
      labs(title=paste(this_gene,", highly expressed in ",paste(this_celltype,collapse=", "),sep=""))+
      scale_fill_manual(values=c("orange","red"))+
      scale_colour_manual(values=this_cols$color, breaks=this_cols$label)+
      scale_size_manual(values=c(1.5,3))+guides(color="none",fill="none")+labs(size="Composition")
    
    color_labels=data.frame()
    
    
  }else{
    
    color_labels=data.frame(label=this_labels,color=c("red","blue","black"))
    
    
    
    plot_data$celltype=mapvalues(plot_data$celltype,from=old_labels,to=this_labels)
    
    this_sizes= c(3,1.5,1.5)
    
    
    if (doub_labels=="Endothelial+\nLens"){
      
      curr_data=subset(plot_data,celltype %in% c(doub_labels,"Lens"))
      
      curr_data1=subset(curr_data,y<2700)
      curr_data2=subset(curr_data,x>5200)
      curr_data3=subset(curr_data,x<4500)
      
      curr_miss=curr_data[unique(c(rownames(curr_data1),rownames(curr_data2),rownames(curr_data3))),]
      
      print(nrow(curr_miss))
      
      plot_data[rownames(curr_miss),"celltype"]=rep("misannotated",nrow(curr_miss))
      
      color_labels=rbind(color_labels,data.frame(label=c("misannotated"),color=c("yellow")))
      
      this_labels=c(this_labels,"misannotated")
      this_sizes=c(this_sizes,1.5)
      
    }
    
    
    plot_data$celltype=factor(plot_data$celltype,levels=this_labels)
    
    gene.plot=ggplot(extra_data,aes(x,y))+geom_point(stroke=0.5,shape=21,colour="lightgrey",fill="NA",size=0.5) 
    
    gene.plot=gene.plot+geom_point(data=plot_data,aes(x,y,fill=value,colour=celltype,size=db),stroke=1,shape=21)+
      
      scale_fill_gradient2(low = "blue",
                           mid = "white",
                           high = "red",midpoint=this_mid)+
      scale_colour_manual(values=color_labels$color,breaks=color_labels$label)+
      scale_size_manual(values=this_sizes)+
      labs(title="",size="Composition", colour="Cell type",fill=paste(this_gene,"\nexpression",sep=""))+
      guides(colour=guide_legend(override.aes = list(size=this_sizes,stroke=2)),size="none")
    
    
    
  }
  
  
  colnames(plot_data)[which(colnames(plot_data)=="value")]=this_value 
  
  return(list(pl=gene.plot,d=plot_data,lb=color_labels))
}


violin_plot_fig2=function(this_data,this_value,this_lb,this_gene){
  
  colnames(this_data)[which(colnames(this_data)==this_value)]="value"
  
  p_theme=bar_theme("Arial",20)
  
  this_plot=ggplot(this_data,aes(x=celltype,y=value,fill=celltype))+
    geom_violin()+
    geom_boxplot(fill="white",width=0.05)+ geom_point(shape=21,size=1,fill="white")+
    stat_summary(fun=mean, geom="point", shape=21, size=4, color="white")+
    stat_n_text(size = 5,family = "Arial") +
    theme(legend.position = "none")+labs(x="",y=this_lb)+p_theme+
    scale_fill_manual(values=c("red","blue","black"))+labs(title=this_gene)
  
  
  colnames(this_data)[which(colnames(this_data)=="value")]=this_value
  
  
  return(this_plot)
  
  
}


bar_theme<-function(this_family,this_size){
  res_theme=theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    text=element_text(family=this_family,size=this_size),
    axis.text = element_text(family=this_family,size=this_size),
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text=element_text(family=this_family,size=this_size)
    
    
    
  )
  
  return(res_theme)
}




get_file_Name=function(files,folder_name,this_ext,this_type){
  genes=c()
  cells=c()
  
  for (file1 in files){
    
    curr_items=unlist(strsplit(file1,"_"))
    celltypes=unlist(strsplit(curr_items[1],"[+]"))
    celltype1=celltypes[1]
    celltype2=celltypes[2]
    gene=gsub(".txt","",curr_items[2])
    genes=c(genes,gene)
    cells=c(celltype1,celltype2,cells)
    
  }
  
  
  cells_out=paste(sort(unique(cells)),collapse = "+")
  this_gene=unique(genes)
  
  this_name=paste(folder_name,cells_out,"_",this_gene,this_type,this_ext,sep="")
  
  
  return(this_name)
}
