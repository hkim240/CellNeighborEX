function [null_DEGs, pvalue_null,fdr_null,logRatio_null]=findNullDEGs(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,log_data_artificialHeteroSpots,gene_name,pCutoff,pCutoff2,lrCutoff)


%%%% Getting the information of heterotypic beads
cellSelect1=[];
criterionClusterIndex=-1;

for idx=1:size(clusterSelect,2)
    
    clusterIndex = clusterSelect(idx);
    combiTemp=neiCombUnique(clusterIndex);      
   
    if combiTemp == center_celltype

        criterionClusterIndex = clusterIndex; % defining the index of heterotypic pair
        cellSelect1=matchComb==clusterIndex; % collecting true heterotypic beads

    end

end


%%%% Performing statistcal tests between artificial heterotypic beads and true ones
max_leng=1;
pvalue_total=cell(max_leng,1);
fdr_total=cell(max_leng,1);
logRatio_total=cell(max_leng,1);

pvalue=ones(size(log_data,1),1);
fdr=ones(size(log_data,1),1);
logRatio=zeros(size(log_data,1),1);

for i=1:size(log_data,1)
    
    % F-test: to test if variances are equal or not
    h_var = vartest2(log_data_artificialHeteroSpots(i,:),log_data(i,find(cellSelect1)));   
    
    if h_var == 0 || h_var == 1
        
        if h_var == 1 

            % Welch's t-test: uneqaul variance
            [h,p,ci,stats]=ttest2(log_data_artificialHeteroSpots(i,:),log_data(i,find(cellSelect1)),'Vartype','unequal');
        
        elseif h_var == 0

            % Student's t-test: equal variance
            [h,p,ci,stats]=ttest2(log_data_artificialHeteroSpots(i,:),log_data(i,find(cellSelect1)));
        
        end

        pvalue(i)=p;
        logRatio(i)=mean(log_data(i,cellSelect1)+1)-mean(log_data_artificialHeteroSpots(i,:)+1);
        
    end

end

[fdr,qvalue] = mafdr(pvalue,'lambda',0.15);

pvalue_total{criterionClusterIndex}=pvalue;
fdr_total{criterionClusterIndex}=fdr;
logRatio_total{criterionClusterIndex}=logRatio;

%%%% Finding DEGs: comparison between heterotypic beads and artificial heterotypic beads
max_leng=size(pvalue_total,1);
DEGindex=zeros(size(gene_name,1),max_leng);
DEGindex(:,criterionClusterIndex)=pvalue_total{criterionClusterIndex}<pCutoff2 & fdr_total{criterionClusterIndex}<pCutoff & logRatio_total{criterionClusterIndex}>lrCutoff;
geneIndexTemp=find(DEGindex(:,criterionClusterIndex) & sum(DEGindex(:,criterionClusterIndex),2)==1);

[~,sortIndex]=sort(logRatio_total{criterionClusterIndex}(geneIndexTemp),'descend'); % sorted by log-ratio
geneIndex_sorted=geneIndexTemp(sortIndex); % DEGs identified from the null model

logRatio_null=logRatio_total{criterionClusterIndex};
logRatio_null=logRatio_null(geneIndex_sorted);
pvalue_null=pvalue_total{criterionClusterIndex};
pvalue_null=pvalue_null(geneIndex_sorted);
fdr_null=fdr_total{criterionClusterIndex};
fdr_null=fdr_null(geneIndex_sorted);

null_DEGs=gene_name(geneIndex_sorted);


end
