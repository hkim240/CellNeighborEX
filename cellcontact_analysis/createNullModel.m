function [log_data_artificialHeteroSpots]=createNullModel(seedNumber, randSize, prop, matchComb, clusterSelect, log_data)


%%%% Creating artificial heterotypic beads
rand('seed',seedNumber);
true_heteroSpot_size=sum(matchComb==clusterSelect(1)); % The number of true heterotypic beads
log_data_artificialHeteroSpots=zeros(size(log_data,1),true_heteroSpot_size*randSize);

for idx=1:true_heteroSpot_size
    
    clusterIndex1=clusterSelect(2); % cluster number for the first cell type of homotypic beads 
    clusterIndex2=clusterSelect(3); % cluster number for the second cel type of homotypic beads
    cellIndex1=find(matchComb==clusterIndex1);
    cellIndex2=find(matchComb==clusterIndex2);
    prop1=prop(idx,1); % first cell type proportion from a true heterotypic bead
    prop2=prop(idx,2); % second cell type proportion from a true heterotypic bead
    randIndex1=randi([1 size(cellIndex1,2)],1,randSize);
    randIndex2=randi([1 size(cellIndex2,2)],1,randSize); 
    
    % How to create an artificial heterotypic bead:
    % (i) Respective homotypic spots are randomly selected.  
    % (ii) Their transcriptomes are mixed according to cell type proportions in a true heterotypic spot. 
    tmp_log_data=prop1*(log_data(:,cellIndex1(randIndex1)))+prop2*(log_data(:,cellIndex2(randIndex2)));
    log_data_artificialHeteroSpots(:,(idx-1)*randSize+1:idx*randSize)=tmp_log_data;
        
end


end