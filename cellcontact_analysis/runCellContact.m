%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CellContact: Detecting cell contact-dependent gene expression 
% 
% Input ...
% Output ...
%
% Created by Hyobin Kim 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % measure execution time: start point  




%%%% Retrieving informaton of cell type markers 
% Importing data
matchMarker = importdata('hippocampus_matchMarker_plus_selected.csv'); % embryo/hippocampus

% Data: 'embryo_matchMarker_plus_selected.csv'
% It says which gene is which cell type marker among genes from 1st to 2089th.

% Data: 'hippocampus_matchMarker_plus_selected.csv'
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




%%%% Cell contact-dependent gene expression analysis
files = dir('input_data/hippocampus/all/temp/*.mat'); % embryo/hippocampus
for k = 1:length(files)
    
    baseFileName = files(k).name;
    name = regexprep(baseFileName,'.mat','');
    
    center_celltype = name;
       
    %%%% Importing data
    load 'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat' % Embryo/Hippocampus
    
    %%%% Data: 'mouseEmbryo_Slideseq_RCTD_top2000_plus.mat'
    % 'cell_id' - 42362 spot IDs
    % 'gene_name' - top 2000 plus cell type markers (2089 genes)
    % 'log_data' - log-normalized data (2089 genes x 42362 spots)
    
    %%%% Data:'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat'
    % 'cell_id' - 41344 spot IDs
    % 'gene_name' - top 2000 plus cell type markers (2034 genes)
    % 'log_data' - log-normalized data (2034 genes x 41344 spots)

    fullFileName = fullfile(files(k).folder, baseFileName);
    load(fullFileName);

    %%%% Data:'celltype1(ct1)+celltype2(ct2).mat'
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
    clusterSize=length(clusterSelect);
    log_data_artificialHeteroSpots=createNull(seedNumber,randSize,prop,matchComb,clusterSelect,log_data);
    

    %----------------------------------------------------------------------
    % STEP3: Finding DEGs
    %----------------------------------------------------------------------
    folderName='results/';   
    [pvalue_total,fdr_total,logRatio_total,ci_total, tvalue_total]=findNullDEGs(center_celltype,clusterSize,neiCombUnique,clusterSelect,matchComb,log_data,log_data_artificialHeteroSpots);
    save([folderName,name,'_stat.mat'],'pvalue_total','fdr_total','logRatio_total','ci_total', 'tvalue_total');

    
    %----------------------------------------------------------------------
    % SETP4: Plotting heatmaps 
    %----------------------------------------------------------------------
    load([folderName,name,'_stat.mat']);
    folderName2='heatmap/';
    
    heatmapFolderName=[folderName, folderName2, center_celltype];
    mkdir(heatmapFolderName);
    
    %%%% cut_off parameter setting for finding DEGs
    pCutoff=0.01; 
    lrCutoff=0.7; % embryo=0.58, hippocampus=0.7
    % delete: DEGnumber=10;
    findCellContactDEGs(heatmapFolderName,center_celltype,cell_id,matchComb,neiCombUnique,matchMarker_var,clusterSelect,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data,log_data_zvalue,gene_name);

end    




timeElapsed = toc % measure execution time: end point 