files = dir('*.csv');

for k = 1:length(files)
    
    fileName = files(k).name;
    name = regexprep(fileName,'.csv','');
    
    center_celltype = strsplit(name,'_');
    center_celltype = center_celltype(3);
    center_celltype = center_celltype{1};
    
    indexFileName = ['index','_','embryo','_',center_celltype,'.csv'];
    index = importdata(indexFileName);
    index= logical(index);
    
    neiCombUniqueFileName = ['neiCombUnique','_','embryo','_',center_celltype,'.csv'];
    neiCombUnique = importdata(neiCombUniqueFileName);
    neiCombUnique = string(neiCombUnique');

    matchCombFileName = ['matchComb','_','embryo','_',center_celltype,'.csv'];
    matchComb = importdata(matchCombFileName);
    matchComb = matchComb';
    
    saveFileName = [center_celltype,'.mat'];
    save(saveFileName,'-v7.3','index','neiCombUnique','matchComb');
    
end    