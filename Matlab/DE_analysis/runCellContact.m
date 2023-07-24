%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CellNeighborEX].
%
% CellNeighborEX is a framework that detects neighbor-dependent genes from
% spatial transcriptomics (ST) data. It identifies up-regulated and
% down-regulated genes when two different cell types contact each other.
% CellNeighborEX works for both image- and NGS-based ST data. 
% 
%
% [ST DATA]
%
% - Image-based ST data (single cell resolution): 
%   e.g., seqFISH data in a mouse embryo
% - NGS-based ST data (10 µm, near cellular resolution): 
%   e.g., Slide-seq data in a mouse embryo, hippocampus, and liver cancer 
% 
%
% [STEPS]
%
% The following four steps are performed in MATLAB:
% - STEP1: Converting log-normalized data to z-values
% - STEP2: Creating the null model of artificial heterotypic spots 
% - STEP3: Finding neighbor-dependent genes 
% - SETP4: Plotting heatmaps 
%
% * STEP2 is applied to NGS-based ST data only.
%
%
% [INPUT] 
%
% The details of the input files are described below corresponding command
% lines:
% 
% - seqFISH in a mouse embryo:
% 'input/seqFISH/seqFISH.mat'
% 'input/seqFISH/celltype1+celltype2.mat'
%
% - Slide-seq in a mouse embryo:
% 'input/embryo/embryo_matchMarker_plus_selected.csv'
% 'input/embryo/mouseEmbryo_Slideseq_RCTD_top2000_plus.mat'
% 'input/embryo/celltype1+celltype2.mat'
%
% - Slide-seq in mouse hippocampus:
% 'input/hippocampus/hippocampus_matchMarker_plus_selected.csv'
% 'input/hippocampus/mouseHippocampus_Slideseq_RCTD_top2000_plus.mat'
% 'input/hippocampus/celltype1+celltype2.mat'
%
% - Slide-seq in mouse liver cancer:
% 'input/liver/liver_matchMarker_plus_selected.csv'
% 'input/liver/mouseLiver_Slideseq_RCTD_top2000_plus.mat'
% 'input/liver/celltype1+celltype2.mat'
%
% * For Slide-seq datasets, the number of the heterotypic beads of 
% celltype1+celltype2 is equal to or larger than 30.
% 
%
% [OUTPUT]
%
% - Differentially expressed genes: DEG, log-ratio, p-value, fdr
% (i) DEGs detected from the null model (for only NGS-based ST data)
% 'celltype1+celltype2_stat_null.mat'
% (ii) > DEGs detected from the comparison of heterotypic and homotypic
%        neighbors (for image-based ST data)
%      > DEGs detected from the comparison of heterotypic and homotypic beads
%        (for NGS-based ST data)
% 'celltype1+celltype2_stat_cellContact.mat'
% 
% - Neighbor-dependent genes: bead barcode, log-value, z-value
% 'celltype1+celltype2_gene.txt'
% 'celltype1+celltype1_gene.txt'
% 'celltype2+celltype2_gene.txt'
%
% - Heatmaps: expression of neighbor-dependent genes 
%             (plus cell type markers for Slide-seq data)
% 'celltype1+celltype2.pdf'
% 
% * The output files above are saved in "output/celltype1+celltype2/"
%
% - DEG list for total cell type pairs
% > DEG list for image-based ST data: cell-type pair, DEG, 
%   homotypic1_log-ratio, homotypic_p-value, homotypic1_fdr
% > DEG list for NGS-based ST data: cell-type pair, DEG, 
%   homotypic1_log-ratio, homotypic1_p-value, homotypic2_log-ratio, 
%   homotypic2_p-value, artificial_log-ratio, artificial_fdr
% 'DEG_list.txt'
%
% *'DEG.list.txt' is saved in "output/"
%
%
% [PARAMETER SETTING]
%
% The details of the parameters are described below corresponding command
% lines:
%
% - Data type: dataType
% - Thresholds: pCutoff, pCutoff2, lrCutoff
% - Up/Down-regulated genes: direction
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % measuring execution time: start point  

dataType = 2; % 1 for image-based ST data, 2 for NGS-based ST data

if dataType == 1 

    files = dir('input/seqFISH/*.mat'); % embryo seqFISH

elseif dataType == 2

    %%%% Adding the informaton of cell type markers 
    % Importing data
    matchMarker = importdata('input/embryo/embryo_matchMarker_plus_selected.csv'); % choose embryo, hippocampus, or liver
    
    % Input: 'embryo_matchMarker_plus_selected.csv' where each row means a gene
    % It says which gene is which cell type marker among genes from 1st to 2089th
    
    % Input: 'hippocampus_matchMarker_plus_selected.csv' where each row means a gene
    % It says which gene is which cell type marker among genes from 1st to 2034th
    % If a gene works as a marker for two cell types, one cell type is written in the 1st column, 
    % and the other cell type is written in the 2nd column
    
    % Input: 'liver_matchMarker_plus_selected.csv' where each row means a gene
    % It says which gene is which cell type marker among genes from 1st to 2007th
    
    matchMarker = regexp(matchMarker, ',', 'split');
    
    for i = 1:size(matchMarker{1},2)
        
        arrayName = strcat('matchMarker',num2str(i));
        matchMarker_var.(arrayName) = cell(size(matchMarker,1),1);
        
        for j = 1:size(matchMarker,1)
    
            stringWords = strread(matchMarker{j,1}{i}, '%s');
            matchMarker_var.(arrayName)(j) = stringWords;
    
        end
        
        matchMarker_var.(arrayName) = string(matchMarker_var.(arrayName));
        
    end

    files = dir('input/embryo/*.mat'); % choose embryo, hippocampus, or liver

end


%%%% Conducting neighbor-dependent gene expression analysis
total_cellType=[];
total_cellContact_DEGs=[];
total_logRatio1_cellContact=[];
total_pvalue1_cellContact=[];

if dataType == 1

    total_fdr1_cellContact=[];

elseif dataType == 2

    total_logRatio2_cellContact=[];
    total_pvalue2_cellContact=[];
    total_logRatio_null_cellContact=[];
    total_fdr_null_cellContact=[];

end    

for k = 1:length(files)
    
    baseFileName = files(k).name;
    
    if ~strcmp(baseFileName,'seqFISH_data.mat') && ~strcmp(baseFileName,'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat') && ~strcmp(baseFileName,'mouseEmbryo_Slideseq_RCTD_top2000_plus.mat') && ~strcmp(baseFileName,'mouseLiver_Slideseq_RCTD_top2000_plus.mat')

        name = regexprep(baseFileName,'.mat','');
        
        center_celltype = name; % heterotypic cell type pair
        
        % Importing data
        if dataType == 1

            load 'input/seqFISH/seqFISH_data.mat'

        elseif dataType == 2    

            load 'input/embryo/mouseEmbryo_Slideseq_RCTD_top2000_plus.mat' % choose Embryo, Hippocampus, or Liver

        end

        % Input: 'seqFISH_data.mat'
        % 'cell_id' - 11797 cell IDs
        % 'gene_name' - 351 genes
        % 'log_data' - log-normalized data (351 genes x 11797 cells)

        % Input: 'mouseEmbryo_Slideseq_RCTD_top2000_plus.mat'
        % 'cell_id' - 42362 spot IDs
        % 'gene_name' - top 2000 plus cell type markers (2089 genes)
        % 'log_data' - log-normalized data (2089 genes x 42362 spots)
        
        % Input:'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat'
        % 'cell_id' - 41344 spot IDs
        % 'gene_name' - top 2000 plus cell type markers (2034 genes)
        % 'log_data' - log-normalized data (2034 genes x 41344 spots)

        % Input:'mouseLiver_Slideseq_RCTD_top2000_plus.mat'
        % 'cell_id' - 22841 spot IDs
        % 'gene_name' - top 2000 plus cell type markers (2007 genes)
        % 'log_data' - log-normalized data (2007 genes x 22841 spots)
    
        fullFileName = fullfile(files(k).folder, baseFileName);
        load(fullFileName);
    
        % Input:'celltype1(ct1)+celltype2(ct2).mat'
        % 'index' - Boolean membership for ct1+ct2, ct1+ct1, ct2+ct2 
        %           sum(index) means the total number of ct1+ct2, ct1+ct1, and ct2+ct2 cells/spots
        % 'neiCombUnique' - three types of pairs: ct1+ct2(1), ct1+ct1(2), ct2+ct2(3)
        % 'matchComb' - sum(index) spots classified into the three types
        % 'prop' - cell type proportions of true heterotypic spots(ct1+ct2) 
    
    
        %----------------------------------------------------------------------
        % STEP1: Converting log-normalized data to z-values 
        %----------------------------------------------------------------------
        cell_id_total = cell_id;
        cell_id = cell_id_total(index);
    
        log_data_total = log_data;
        log_data = log_data_total(:,index);

        log_data_zvalue = (log_data_total-repmat(mean(log_data_total,2),1,size(log_data_total,2)))./repmat(std(log_data_total')',1,size(log_data_total,2));
        log_data_zvalue(isnan(log_data_zvalue)) = 0;
        log_data_zvalue = log_data_zvalue(:,index);
        
        %----------------------------------------------------------------------
        % STEP2: Creating the null model of artificial heterotypic spots
        %----------------------------------------------------------------------
        if dataType == 2

            seedNumber = 1; % random generator
            clusterSelect=unique(matchComb);

            heteroSpotNum = sum(matchComb==clusterSelect(1)); % number of true heterotypic spots
            if heteroSpotNum > 100
                randSize = 30; % number of artificial heterotypic spots generated for each spot
            else    
                randSize = 100; % number of artificial heterotypic spots generated for each spot
            end
            
            [log_data_artificialHeteroSpots,normalized_props]=createNullModel(seedNumber,randSize,prop,matchComb,clusterSelect,log_data);
        
        end    
    
        %----------------------------------------------------------------------
        % STEP3: Finding neighbor-dependent genes
        %----------------------------------------------------------------------
        if dataType == 1

            clusterSelect=unique(matchComb);

        end    

        % Setting cut-off values
        pCutoff = 0.01; % p-value 
        pCutoff2 = 0.01; % fdr
        lrCutoff = 0.4;  % log-ratio
        
        % Finding up- or down-regulated genes
        direction = 1; % 1 for up-regulated genes, 0 for down-regulated genes % choose up or down

        % Comparing heterotypic beads with artificial heterotypic beads
        folderName='output/';  
        folderName2=[folderName, center_celltype];

        if dataType == 1
            
            % Comparing heterotypic beads with homotypic beads
            [cellContact_DEGs_IDX,cellContact_DEGs,pvalue1_cellContact,fdr1_cellContact,logRatio1_cellContact]=findCellContactDEGs_img(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,gene_name,pCutoff,pCutoff2,lrCutoff,direction);
        
            if size(cellContact_DEGs,1) > 0
                
                mkdir(folderName2);
                save([folderName2,'/',name,'_stat_cellContact.mat'],'cellContact_DEGs_IDX','cellContact_DEGs','pvalue1_cellContact','fdr1_cellContact','logRatio1_cellContact');
                
                tmp_cellType=[];
                for idx=1:size(cellContact_DEGs,1)
    
                    tmp_cellType = [tmp_cellType; convertCharsToStrings(center_celltype)];
    
                end   
    
                total_cellType=[total_cellType; tmp_cellType];
                total_cellContact_DEGs=[total_cellContact_DEGs; cellContact_DEGs];
                total_logRatio1_cellContact=[total_logRatio1_cellContact; logRatio1_cellContact];
                total_pvalue1_cellContact=[total_pvalue1_cellContact; pvalue1_cellContact];
                total_fdr1_cellContact=[total_fdr1_cellContact; fdr1_cellContact];
                
            end

        elseif dataType == 2

            [null_DEGs,pvalue_null,fdr_null,logRatio_null]=findNullDEGs(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,log_data_artificialHeteroSpots,gene_name,pCutoff,pCutoff2,lrCutoff,direction);
            
            if size(null_DEGs,1) > 0
        
                mkdir(folderName2);
                save([folderName2,'/',name,'_stat_null.mat'],'null_DEGs','pvalue_null','fdr_null','logRatio_null');
                
            end    
            
        
            % Comparing heterotypic beads with homotypic beads
            if size(null_DEGs,1) > 0
        
                load([folderName2,'/',name,'_stat_null.mat']);
                [cellContact_DEGs_IDX,cellContact_DEGs,pvalue1_cellContact,fdr1_cellContact,logRatio1_cellContact,pvalue2_cellContact,fdr2_cellContact,logRatio2_cellContact,logRatio_null_cellContact,fdr_null_cellContact]=findCellContactDEGs(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,gene_name,pCutoff,lrCutoff,null_DEGs,fdr_null,logRatio_null,direction);
                
                if size(cellContact_DEGs,1) > 0
                    
                    save([folderName2,'/',name,'_stat_cellContact.mat'],'cellContact_DEGs_IDX','cellContact_DEGs','pvalue1_cellContact','fdr1_cellContact','logRatio1_cellContact','pvalue2_cellContact','fdr2_cellContact','logRatio2_cellContact');
                    
                    tmp_cellType=[];
                    for idx=1:size(cellContact_DEGs,1)
    
                        tmp_cellType = [tmp_cellType; convertCharsToStrings(center_celltype)];
    
                    end   
    
                    total_cellType=[total_cellType; tmp_cellType];
                    total_cellContact_DEGs=[total_cellContact_DEGs; cellContact_DEGs];
                    total_logRatio1_cellContact=[total_logRatio1_cellContact; logRatio1_cellContact];
                    total_pvalue1_cellContact=[total_pvalue1_cellContact; pvalue1_cellContact];
                    total_logRatio2_cellContact=[total_logRatio2_cellContact; logRatio2_cellContact];
                    total_pvalue2_cellContact=[total_pvalue2_cellContact; pvalue2_cellContact];
                    total_logRatio_null_cellContact=[total_logRatio_null_cellContact; logRatio_null_cellContact];
                    total_fdr_null_cellContact=[total_fdr_null_cellContact; fdr_null_cellContact];
    
                end
                
            end
        
        end

        %----------------------------------------------------------------------
        % SETP4: Plotting heatmaps 
        %----------------------------------------------------------------------
        if dataType == 1  && size(cellContact_DEGs,1) > 0

            plotHeatmaps_img(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,log_data_zvalue,gene_name,cell_id,cellContact_DEGs_IDX,logRatio1_cellContact,folderName2);
       
        elseif dataType == 2 && size(null_DEGs,1) > 0 && size(cellContact_DEGs,1) > 0
            
            plotHeatmaps(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,log_data_zvalue,gene_name,cell_id,cellContact_DEGs_IDX,logRatio1_cellContact,logRatio2_cellContact,matchMarker_var,folderName2);
            
        end
    
    end

end   

if dataType == 1

    total=[total_cellType, total_cellContact_DEGs, total_logRatio1_cellContact, total_pvalue1_cellContact,total_fdr1_cellContact]; 

elseif dataType == 2  

    total=[total_cellType,total_cellContact_DEGs,total_logRatio1_cellContact,total_pvalue1_cellContact,total_logRatio2_cellContact,total_pvalue2_cellContact,total_logRatio_null_cellContact,total_fdr_null_cellContact];    

end

writematrix(total,[folderName,'/','DEG_list.txt']);
timeElapsed = toc; % measuring execution time: end point 