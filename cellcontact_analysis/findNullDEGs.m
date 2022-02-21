function [pvalue_total,fdr_total,logRatio_total, ci_total, tvalue_total ]=findNullDEGs(center_celltype,clusterSize,neiCombUnique,clusterSelect,matchComb,log_data,log_data_artificialHeteroSpots)

max_leng=1;
pvalue_total=cell(max_leng,1);
fdr_total=cell(max_leng,1);
logRatio_total=cell(max_leng,1);
ci_total=cell(max_leng,1);
tvalue_total=cell(max_leng,1);

cellSelect1=[];
criterionClusterIndex=-1;

for idx=1:clusterSize
    
    clusterIndex = clusterSelect(idx);
    combiTemp=neiCombUnique(clusterIndex);      
   
    if combiTemp == center_celltype
        cellSelect1=matchComb==clusterIndex;
        criterionClusterIndex = clusterIndex;
    end
end

pvalue=ones(size(log_data,1),1);
fdr=ones(size(log_data,1),1);
logRatio=zeros(size(log_data,1),1);
ci_lower=zeros(size(log_data,1),1);
ci_upper=zeros(size(log_data,1),1);
tvalue=zeros(size(log_data,1),1);

for i=1:size(log_data,1)
    h_var = vartest2(log_data_artificialHeteroSpots(i,:),log_data(i,find(cellSelect1)));   
    if h_var == 0 || h_var == 1
        if h_var == 1 
            %%%% Welch-t-test: uneqaul variance
            [h,p,ci,stats]=ttest2(log_data_artificialHeteroSpots(i,:),log_data(i,find(cellSelect1)),'Vartype','unequal');
        elseif h_var == 0
            %%%% Student's t-test: equal variance
            [h,p,ci,stats]=ttest2(log_data_artificialHeteroSpots(i,:),log_data(i,find(cellSelect1)));
        end

        pvalue(i)=p;
        logRatio(i)=mean(log_data(i,cellSelect1)+1)-mean(log_data_artificialHeteroSpots(i,:)+1);
        ci_lower=ci(1);
        ci_upper=ci(2);
        tvalue(i)=stats.tstat;
    end
end
% original script: fdr=mafdr(pvalue);
% Warning: The estimated PI0 is greater than 1. Please check the p-values are valid or try a different lambda method. 
% PI0 is set to 1. So the following modified script was added:
[fdr,qvalue] = mafdr(pvalue,'lambda',0.15);

pvalue_total{criterionClusterIndex}=pvalue;
fdr_total{criterionClusterIndex}=fdr;
logRatio_total{criterionClusterIndex}=logRatio;
ci_total{criterionClusterIndex}=[ci_lower, ci_upper]; 
tvalue_total{criterionClusterIndex}=tvalue;

end
