%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CellContact is a statistical framework to detect cell contact-dependent 
% gene expression. This is a new framework to analyze spatial transcriptomics 
% (ST) data using spatial beads composed of multiple cell types, especially 
% for Slide-seq data with 10 µm resolution. 
%
%
% [ST DATA]
%
% Slide-seq V2 data in mouse embryo and hippcampus
%
%
% [STEPS]
%
% The following four steps are performed:
% STEP1: Converting log-normalized data to z-values
% STEP2: Creating the null model of artificial heterotypic spots
% STEP3: Finding differentially expressed genes (DEGs)
% SETP4: Plotting heatmaps 
%
%
% [INPUT]
%
% Embryo:
% 'input/embryo/embryo_matchMarker_plus_selected.csv'
% 'input/embryo/mouseEmbryo_Slideseq_RCTD_top2000_plus.mat'
% 'input/embryo/celltype1+celltype2.mat'
%
% Hippocampus:
% 'input/hippocampus/hippocampus_matchMarker_plus_selected.csv'
% 'input/hippocampus/mouseHippocampus_Slideseq_RCTD_top2000_plus.mat'
% 'input/hippocampus/celltype1+celltype2.mat'
%
% The details of the input files are described below the command lines.
% "celltype1" and "celltype2" refer to respective cell types in mouse 
% embryo and hippocampus. For all heterotypic pairs, each sample size is 
% greater than or equal to 30 (i.e., # of the heterotypic beads of
% celltype1+celltype2 >= 30).
%
%
% [OUTPUT]
%
% Differentially expressed genes: DEGs, fdr, log-ratio, p-value
% (i) DEGs identified from the null model
% 'celltype1+celltype2_stat_null.mat'
% (ii) DEGs identified from the comparison between heterotypic and homotypic beads
% 'celltype1+celltype2_stat_cellContact.mat'
% 
% Cell contact-dependent genes: bead barcode, log-value, z-value
% 'celltype1+celltype2_gene.txt'
% 'celltype1+celltype1_gene.txt'
% 'celltype2+celltype2_gene.txt'
%
% Heatmaps: expression of cell contact-dependent genes & cell type markers in homotypic and heterotypic beads
% 'celltype1+celltype2.pdf'
% 
% The output files above are saved in "output/celltype1+celltype2/"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % measuring execution time: start point  


%%%% Adding the informaton of cell type markers 
% Importing data
matchMarker = importdata('input/embryo/embryo_matchMarker_plus_selected.csv'); % choose embryo or hippocampus

% Data: 'embryo_matchMarker_plus_selected.csv' where each row means a gene
% It says which gene is which cell type marker among genes from 1st to 2089th.

% Data: 'hippocampus_matchMarker_plus_selected.csv' where each row means a gene
% It says which gene is which cell type marker among genes from 1st to 2034th.
% If a gene is simultaneously a marker of two cell types, 
% one cell type is written in the 1st column, and the other cell type is written in the 2nd column. 

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


%%%% Conducting cell contact-dependent gene expression analysis
files = dir('input/embryo/*.mat'); % choose embryo or hippocampus

for k = 1:length(files)
    
    baseFileName = files(k).name;
    
    if ~strcmp(baseFileName,'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat') & ~strcmp(baseFileName,'mouseEmbryo_Slideseq_RCTD_top2000_plus.mat')

        name = regexprep(baseFileName,'.mat','');
        
        center_celltype = name; % heterotypic pair
           
        % Importing data
        load 'input/embryo/mouseEmbryo_Slideseq_RCTD_top2000_plus.mat' % choose Embryo or Hippocampus
        
        % Data: 'mouseEmbryo_Slideseq_RCTD_top2000_plus.mat'
        % 'cell_id' - 42362 spot IDs
        % 'gene_name' - top 2000 plus cell type markers (2089 genes)
        % 'log_data' - log-normalized data (2089 genes x 42362 spots)
        
        % Data:'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat'
        % 'cell_id' - 41344 spot IDs
        % 'gene_name' - top 2000 plus cell type markers (2034 genes)
        % 'log_data' - log-normalized data (2034 genes x 41344 spots)
    
        fullFileName = fullfile(files(k).folder, baseFileName);
        load(fullFileName);
    
        % Data:'celltype1(ct1)+celltype2(ct2).mat'
        % 'index' - Boolean membership for ct1+ct2, ct1+ct1, ct2+ct2 
        %           among 42362 spots for embryo/41344 spots for hippocampus
        %           sum(index) means the total number of ct1+ct2, ct1+ct1, and ct2+ct2 spots
        % 'neiCombUnique' - three types of pairs: ct1+ct2(1), ct1+ct1(2), ct2+ct2(3)
        % 'matchComb' - sum(index) spots classified into the three types
        % 'prop' - cell type proportions of true heterotypic spots (ct1+ct2)
    
    
        %----------------------------------------------------------------------
        % STEP1: Converting log-normalized data to z-values 
        %----------------------------------------------------------------------
        cell_id_total = cell_id;
        cell_id = cell_id_total(index);
    
        log_data_total=log_data;
        log_data=log_data_total(:,index);
        log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
        log_data_zvalue(isnan(log_data_zvalue))=0;
    
    
        %----------------------------------------------------------------------
        % STEP2: Creating the null model of artificial heterotypic spots
        %----------------------------------------------------------------------
        seedNumber=1; randSize=100;
        clusterSelect=unique(matchComb);
        log_data_artificialHeteroSpots=createNullModel(seedNumber,randSize,prop,matchComb,clusterSelect,log_data);
        
    
        %----------------------------------------------------------------------
        % STEP3: Finding differentially expressed genes (DEGs)
        %----------------------------------------------------------------------
        % Setting cut-off values
        pCutoff=0.01; pCutoff2=0.001; 
        lrCutoff=0.58; % set log-ratio: for embryo=0.58, for hippocampus=0.7
    
        % Comparing heterotypic beads with artificial heterotypic beads
        folderName='output/';  
        folderName2=[folderName, center_celltype];
        [null_DEGs, pvalue_null,fdr_null,logRatio_null]=findNullDEGs(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,log_data_artificialHeteroSpots,gene_name,pCutoff,pCutoff2,lrCutoff);
        
        if size(null_DEGs,1) > 0
    
            mkdir(folderName2);
            save([folderName2,'/',name,'_stat_null.mat'],'null_DEGs','pvalue_null','fdr_null','logRatio_null');
        
        end    
        
    
        % Comparing heterotypic beads with homotypic beads
        if size(null_DEGs,1) > 0
    
            load([folderName2,'/',name,'_stat_null.mat']);
            [cellContact_DEGs_IDX,cellContact_DEGs,pvalue1_cellContact,fdr1_cellContact,logRatio1_cellContact,pvalue2_cellContact,fdr2_cellContact,logRatio2_cellContact]=findCellContactDEGs(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,gene_name,pCutoff,lrCutoff,null_DEGs);
            
            if size(cellContact_DEGs,1) > 0
    
                save([folderName2,'/',name,'_stat_cellContact.mat'],'cellContact_DEGs_IDX','cellContact_DEGs','pvalue1_cellContact','fdr1_cellContact','logRatio1_cellContact','pvalue2_cellContact','fdr2_cellContact','logRatio2_cellContact');
            
            end
    
        end
        
    
        
        %----------------------------------------------------------------------
        % SETP4: Plotting heatmaps 
        %----------------------------------------------------------------------
        if size(null_DEGs,1) > 0 & size(cellContact_DEGs,1) > 0
    
            plotHeatmaps(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,log_data_zvalue,gene_name,cell_id,cellContact_DEGs_IDX,logRatio1_cellContact,logRatio2_cellContact,matchMarker_var,folderName2);
        
        end
    
    end

end    


timeElapsed = toc % measuring execution time: end point 