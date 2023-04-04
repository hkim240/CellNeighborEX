function [cellContact_DEGs_IDX,cellContact_DEGs,pvalue1_cellContact,fdr1_cellContact,logRatio1_cellContact,pvalue2_cellContact,fdr2_cellContact,logRatio2_cellContact,logRatio_null_cellContact,fdr_null_cellContact]=findCellContactDEGs(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,gene_name,pCutoff,lrCutoff,null_DEGs,fdr_null,logRatio_null,direction)


%%%% Getting the information of heterotypic beads
criterionClusterIndex=-1;
cellSelect1=[];

for idx=1:size(clusterSelect,2)
    
    clusterIndex=clusterSelect(idx);
    combiTemp=neiCombUnique(clusterIndex);
    
    if combiTemp == center_celltype

        criterionClusterIndex=clusterIndex; % defining the index of heterotypic pair
        cellSelect1=matchComb==criterionClusterIndex; % collecting true heterotypic beads 

    end
    
end


%%%% Performing statistcal tests between heterotypic beads and homotypic beads
cellSelect2=matchComb==clusterSelect(2); % first cell type of homotypic beads
cellSelect3=matchComb==clusterSelect(3); % second cell type of homotypic beads

max_leng=1;
pvalue1_total=cell(max_leng,1);
fdr1_total=cell(max_leng,1);
logRatio1_total=cell(max_leng,1);
pvalue2_total=cell(max_leng,1);
fdr2_total=cell(max_leng,1);
logRatio2_total=cell(max_leng,1);

pvalue1=ones(size(log_data,1),1);
pvalue2=ones(size(log_data,1),1);

logRatio1=zeros(size(log_data,1),1);
logRatio2=zeros(size(log_data,1),1);

sample_size=30;

for i=1:size(log_data,1)
    
    % Comparing heterotypic beads with the first cell type of homotypic beads 
    if sum(cellSelect2) < sample_size

        % Wilcoxon-rank sum test   
        [p,h,stats] = ranksum(log_data(i,find(cellSelect2)),log_data(i,find(cellSelect1)));
        
    else

        % F-test: to test if variances are equal or not
        h_var = vartest2(log_data(i,find(cellSelect2)),log_data(i,find(cellSelect1)));

        if isnan(h_var)
            p = 1; 
      
        elseif h_var == 1 

            % Welch's t-test: uneqaul variance
            [h,p,ci,stats]=ttest2(log_data(i,find(cellSelect2)),log_data(i,find(cellSelect1)),'Vartype','unequal');
            
        elseif h_var == 0

            % Student's t-test: equal variance
            [h,p,ci,stats]=ttest2(log_data(i,find(cellSelect2)),log_data(i,find(cellSelect1)));
           
        end    

    end  

    pvalue1(i)=p;
    logRatio1(i)=mean(log_data(i,cellSelect1)+1)-mean(log_data(i,cellSelect2)+1);

    % Comparing heterotypic beads with the second cell type of homotypic beads
    if sum(cellSelect3) < sample_size

        % Wilcoxon-rank sum test   
        [p2,h2,stats2] = ranksum(log_data(i,find(cellSelect3)),log_data(i,find(cellSelect1)));
      
    else

        % F-test: to test if variances are equal or not
        h_var2 = vartest2(log_data(i,find(cellSelect3)),log_data(i,find(cellSelect1)));
        
        if isnan(h_var2)
            p2 = 1; 
       
        elseif h_var2 == 1 

            % Welch's t-test: uneqaul variance
            [h2,p2,ci,stats2]=ttest2(log_data(i,find(cellSelect3)),log_data(i,find(cellSelect1)),'Vartype','unequal');
            
        elseif h_var2 == 0

            % Student's t-test: equal variance
            [h2,p2,ci,stats2]=ttest2(log_data(i,find(cellSelect3)),log_data(i,find(cellSelect1)));
   
        end

    end 
   
    pvalue2(i)=p2;
    logRatio2(i)=mean(log_data(i,cellSelect1)+1)-mean(log_data(i,cellSelect3)+1); 
     
end

[fdr1,qvalue1] = mafdr(pvalue1,'lambda',0.15);
[fdr2,qvalue2] = mafdr(pvalue2,'lambda',0.15);

pvalue1_total{criterionClusterIndex}=pvalue1;
fdr1_total{criterionClusterIndex}=fdr1;
logRatio1_total{criterionClusterIndex}=logRatio1;
pvalue2_total{criterionClusterIndex}=pvalue2;
fdr2_total{criterionClusterIndex}=fdr2;
logRatio2_total{criterionClusterIndex}=logRatio2;

%%%% Finding DEGs: comparison between heterotypic beads and homotypic beads
max_leng=size(pvalue1_total,1);
DEGindex=zeros(size(gene_name,1),max_leng);

if direction==1 % up-regulated genes

    DEGindex(:,criterionClusterIndex)=pvalue1_total{criterionClusterIndex}<pCutoff...
        & logRatio1_total{criterionClusterIndex}>lrCutoff...
        & pvalue2_total{criterionClusterIndex}<pCutoff...
        & logRatio2_total{criterionClusterIndex}>lrCutoff;
 
else % down-regulated genes

    DEGindex(:,criterionClusterIndex)=pvalue1_total{criterionClusterIndex}<pCutoff...
        & logRatio1_total{criterionClusterIndex}<(-lrCutoff)...
        & pvalue2_total{criterionClusterIndex}<pCutoff...
        & logRatio2_total{criterionClusterIndex}<(-lrCutoff);

end

geneIndexTemp2=find(DEGindex(:,criterionClusterIndex) & sum(DEGindex(:,criterionClusterIndex),2)==1);
temp_DEGs=gene_name(geneIndexTemp2);


%%%% Confirming the significance of DEGs: intersection between DEGx and DEGy
% DEGx (null_DEGs): DEGs identified from the null model
% DEGy (temp_DEGs): DEGs identified from the comparison between heterotypic and homotypic beads

% final cell-contact DEGs
cellContact_DEGs=intersect(null_DEGs, temp_DEGs);   
cellContact_DEGs_IDX=find(ismember(gene_name,cellContact_DEGs)); 

null_cellContact_DEGs_IDX=[];
for idx=1:size(cellContact_DEGs,1)
    
    tmp_val=find(ismember(null_DEGs,gene_name(cellContact_DEGs_IDX(idx))));
    null_cellContact_DEGs_IDX=[null_cellContact_DEGs_IDX; tmp_val];

end

logRatio_null_cellContact=logRatio_null(null_cellContact_DEGs_IDX);
fdr_null_cellContact=fdr_null(null_cellContact_DEGs_IDX);

logRatio1_cellContact=logRatio1_total{criterionClusterIndex};
logRatio1_cellContact=logRatio1_cellContact(cellContact_DEGs_IDX);
pvalue1_cellContact=pvalue1_total{criterionClusterIndex};
pvalue1_cellContact=pvalue1_cellContact(cellContact_DEGs_IDX);
fdr1_cellContact=fdr1_total{criterionClusterIndex};
fdr1_cellContact=fdr1_cellContact(cellContact_DEGs_IDX);

logRatio2_cellContact=logRatio2_total{criterionClusterIndex};
logRatio2_cellContact=logRatio2_cellContact(cellContact_DEGs_IDX);
pvalue2_cellContact=pvalue2_total{criterionClusterIndex};
pvalue2_cellContact=pvalue2_cellContact(cellContact_DEGs_IDX);
fdr2_cellContact=fdr2_total{criterionClusterIndex};
fdr2_cellContact=fdr2_cellContact(cellContact_DEGs_IDX);


end
