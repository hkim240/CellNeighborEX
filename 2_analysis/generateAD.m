function [log_data_artificialDoublets, artificialDoubletsComb]=generateAD(seedNumber, randSize, prop, matchComb, clusterSelect, log_data)

rand('seed',seedNumber);
real_doublet_size=sum(matchComb==clusterSelect(1)); % real heterotypic beads
log_data_artificialDoublets=zeros(size(log_data,1),real_doublet_size*randSize);
artificialDoubletsComb=zeros(1,real_doublet_size*randSize);
for idx=1:real_doublet_size
    
    clusterIndex1=clusterSelect(2); % cluster number for singlet1
    clusterIndex2=clusterSelect(3); % cluster number for singlet2
    cellIndex1=find(matchComb==clusterIndex1);
    cellIndex2=find(matchComb==clusterIndex2);
    prop1=prop(idx,1); % first cell type proportion from a real heterotypic bead
    prop2=prop(idx,2); % second cell type proportion from a real heterotypic bead
    randIndex1=randi([1 size(cellIndex1,2)],1,randSize);
    randIndex2=randi([1 size(cellIndex2,2)],1,randSize); 
    % generate artificial heterotypic beads
    tmp_log_data=prop1*(log_data(:,cellIndex1(randIndex1)))+prop2*(log_data(:,cellIndex2(randIndex2)));
    log_data_artificialDoublets(:,(idx-1)*randSize+1:idx*randSize)=tmp_log_data;
    artificialDoubletsComb((idx-1)*randSize+1:idx*randSize)=idx;
    
end

end