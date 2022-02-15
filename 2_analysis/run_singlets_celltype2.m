tic % measure execution time: start point  

matchMarker = importdata('hippocampus_matchMarker_plus_selected.csv'); % embryo/hippocampus
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

files = dir('input_data/hippocampus/all/*.mat'); % embryo/hippocampus
for k = 1:length(files)
    
    baseFileName = files(k).name;
    name = regexprep(baseFileName,'.mat','');
    
    center_celltype = name;
       
    %%%%%%%%%%%%% Importing data %%%%%%%%%%%%%
    %load 'mouseEmbryo_Slideseq_RCTD_whole.mat' 
    load 'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat' % Embryo/Hippocampus
    
    % 'mouseEmbryo_Slideseq_RCTD_top2000_plus.mat'
    % 1. 'cell_id' - 42362 spots
    % 2. 'gene_name' - top 2000 plus cell type markers (# of the whole genes=23124)
    % 3. 'log_data' - 2089 by 42362
    
    % 'mouseHippocampus_Slideseq_RCTD_top2000_plus.mat'
    % 1. 'cell_id' - 41344 spots
    % 2. 'gene_name' - top 2000 plus cell type markers (# of the whole genes=23265)
    % 3. 'log_data' - 2034 by 41344

    fullFileName = fullfile(files(k).folder, baseFileName);
    load(fullFileName);
    % datafile 'ct1+ct2.mat' where 'ct' means cell type
    % 4. 'index' - Boolean membership for celltype1+celltype2 among 42362 beads
    % 5. 'neiCombUnique' - three types of compositions: ct1+ct2(heterotypic), ct1+ct1(homotypic), ct2+ct2(homotypic)
    % 6. 'matchComb' - three types matched to homotypic and heterotypic beads

    %%%%%%%%%%%%% Converting log data to z-values %%%%%%%%%%%%%
    cell_id_total = cell_id;
    cell_id = cell_id_total(index);

    log_data_total=log_data;
    log_data=log_data_total(:,index);
    log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
    log_data_zvalue(isnan(log_data_zvalue))=0;

    %%%%%%%%%%%%% Generating artificial heterotypic beads %%%%%%%%%%%%%
    seedNumber=1; randSize=100;
    clusterSelect=unique(matchComb);
    clusterSize=length(clusterSelect);
    [log_data_artificialDoublets, artificialDoubletsComb]=generateAD(seedNumber,randSize,prop,matchComb,clusterSelect,log_data);
    
    %%%%%%%%%%%%% Finding DEGs %%%%%%%%%%%%%
    folderName='results/';   
    [pvalue_total,fdr_total,logRatio_total,ci_total, tvalue_total]=DEG_ranksum4cluster2_comp(center_celltype,clusterSize,neiCombUnique,clusterSelect,matchComb,log_data,log_data_artificialDoublets);
    save([folderName,name,'_stat.mat'],'pvalue_total','fdr_total','logRatio_total','ci_total', 'tvalue_total');

    %%%%%%%%%%%%% Plotting heatmaps %%%%%%%%%%%%%
    load([folderName,name,'_stat.mat']);
    folderName2='heatmap/';
    
    heatmapFolderName = [folderName, folderName2, center_celltype];
    mkdir(heatmapFolderName);
    
    % cut_off parameter setting for finding DEGs
    pCutoff=0.01; lrCutoff=0.7; % embryo=0.58, hippocampus=0.7
    DEGnumber=10;
    drawHeatmap_cluster2_comp(heatmapFolderName,center_celltype,cell_id,matchComb,neiCombUnique,matchMarker_var,clusterSelect,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data,log_data_zvalue,gene_name,DEGnumber);

end    

timeElapsed = toc % measure execution time: end point 