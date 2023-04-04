function plotHeatmaps(center_celltype,clusterSelect,matchComb,neiCombUnique,log_data,log_data_zvalue,gene_name,cell_id,cellContact_DEGs_IDX,logRatio1_cellContact,logRatio2_cellContact,matchMarker_var,folderName2)


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

%%%% Getting the heatmaps
cellSelect2=matchComb==clusterSelect(2);
cellSelect3=matchComb==clusterSelect(3);

if size(cellContact_DEGs_IDX,1) > 0
    
    close all
    
    %%%% Heatmaps for the expression of DEGs
    % Plotting subplot(3,3,1)
    subplot(3,3,1);
    imagesc(log_data_zvalue(cellContact_DEGs_IDX,find(cellSelect1)))
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',12)
    set(gca,'YTickLabelMode','auto')
    yticks(1:size(cellContact_DEGs_IDX,1))
    yticklabels(gene_name(cellContact_DEGs_IDX))
    
    if size(cellContact_DEGs_IDX,1) > 30
        
        ax = gca;
        ax.YAxis.FontSize = 6;

    end
    
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial','DefaultAxesTitleFontWeight','normal') 
    title({neiCombUnique(criterionClusterIndex);['(n=',num2str(length(find(cellSelect1))),')']},'FontSize',10); % font size: embryo & hippocampus=14, liver=10


    % Printing out: bead barcodes, log data, and z-values of heterotypic beads per DEG
    for i=1:size(cellContact_DEGs_IDX,1)
        
        filename=sprintf('%s_%s.txt',[neiCombUnique(criterionClusterIndex),gene_name(cellContact_DEGs_IDX(i))]);
        log_values=log_data(cellContact_DEGs_IDX(i),find(cellSelect1));
        zvalues=log_data_zvalue(cellContact_DEGs_IDX(i),find(cellSelect1));
        barcodes=cell_id(find(cellSelect1));
        integrated=[barcodes; log_values; zvalues];
        integrated=integrated';
        integrated=sortrows(integrated,2);
        writematrix(integrated,[folderName2,'/',filename]);

    end
    
  
    % Plotting subplot(3,3,2)
    subplot(3,3,2);
    imagesc(log_data_zvalue(cellContact_DEGs_IDX,find(cellSelect2)))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial','DefaultAxesTitleFontWeight','normal')
    title({neiCombUnique(clusterSelect(2));['(n=',num2str(length(find(cellSelect2))),')']},'FontSize',10); % font size: embryo & hippocampus=14, liver=10

    % Printing out: bead barcodes, log data, and z-values of homotypic beads per DEG
    for i=1:size(cellContact_DEGs_IDX,1)

        filename = sprintf('%s_%s.txt',[neiCombUnique(clusterSelect(2)),gene_name(cellContact_DEGs_IDX(i))]);
        log_values=log_data(cellContact_DEGs_IDX(i),find(cellSelect2));
        zvalues=log_data_zvalue(cellContact_DEGs_IDX(i),find(cellSelect2));
        barcodes=cell_id(find(cellSelect2));
        integrated=[barcodes; log_values; zvalues];
        integrated=integrated';
        integrated=sortrows(integrated,2);
        writematrix(integrated,[folderName2,'/',filename]);

    end
    
    % Plotting subplot(3,3,3)
    subplot(3,3,3);
    imagesc(log_data_zvalue(cellContact_DEGs_IDX,find(cellSelect3)))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial','DefaultAxesTitleFontWeight','normal')
    title({neiCombUnique(clusterSelect(3));['(n=',num2str(length(find(cellSelect3))),')']},'FontSize',10); % font size: embryo & hippocampus=14, liver=10

    % Printing out: bead barcodes, log data, and z-values of homotypic beads per DEG
    for i=1:size(cellContact_DEGs_IDX,1)

        filename = sprintf('%s_%s.txt',[neiCombUnique(clusterSelect(3)),gene_name(cellContact_DEGs_IDX(i))]);
        log_values=log_data(cellContact_DEGs_IDX(i),find(cellSelect3));
        zvalues=log_data_zvalue(cellContact_DEGs_IDX(i),find(cellSelect3));
        barcodes=cell_id(find(cellSelect3));
        integrated=[barcodes; log_values; zvalues];
        integrated=integrated';
        integrated=sortrows(integrated,2);
        writematrix(integrated,[folderName2,'/',filename]);

    end
    
  
    %%%% Heatmaps for the expression of cell type markers    
    % Getting the information of cell type marker index
    fns = fieldnames(matchMarker_var);
    cellSelect2FinalIDX=[];

    for var_i=1:size(fns,1)

        markerIDX = matchMarker_var.(fns{var_i})==neiCombUnique(2);
        tmp = find(markerIDX==1);
        cellSelect2FinalIDX=[cellSelect2FinalIDX; tmp];

    end
    
    % Plotting subplot(3,3,4)
    subplot(3,3,4); 
    imagesc(log_data_zvalue(cellSelect2FinalIDX,cellSelect1))
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',12)
    set(gca,'YTickLabelMode','auto')
    yticks(1:size(cellSelect2FinalIDX,1))
    yticklabels(gene_name(cellSelect2FinalIDX))
    
    if size(cellSelect2FinalIDX,1) > 30

        ax = gca;
        ax.YAxis.FontSize = 6;

    end

    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    
    % Plotting subplot(3,3,5)
    subplot(3,3,5);
    imagesc(log_data_zvalue(cellSelect2FinalIDX,cellSelect2))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
   
    % Plotting subplot(3,3,6)
    subplot(3,3,6);
    imagesc(log_data_zvalue(cellSelect2FinalIDX,cellSelect3))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    
    % Getting the information of cell type marker index
    cellSelect3FinalIDX=[];
    for var_i=1:size(fns,1)

        markerIDX = matchMarker_var.(fns{var_i})==neiCombUnique(3);
        tmp = find(markerIDX==1);
        cellSelect3FinalIDX = [cellSelect3FinalIDX; tmp];

    end
    
    % Plotting subplot(3,3,7)
    subplot(3,3,7);
    imagesc(log_data_zvalue(cellSelect3FinalIDX,cellSelect1))
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',12)
    set(gca,'YTickLabelMode','auto')
    yticks(1:size(cellSelect3FinalIDX,1))
    yticklabels(gene_name(cellSelect3FinalIDX))
    
    if size(cellSelect3FinalIDX,1) > 30

        ax = gca;
        ax.YAxis.FontSize = 6;

    end

    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    
    % Plotting subplot(3,3,8)
    subplot(3,3,8); 
    imagesc(log_data_zvalue(cellSelect3FinalIDX,cellSelect2))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    
    % Plotting subplot(3,3,9)
    subplot(3,3,9); 
    imagesc(log_data_zvalue(cellSelect3FinalIDX,cellSelect3))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    
    % Arranging the positions of heatmaps
    set(gcf, 'Position', [100, 100, 550, 650])

    for i=7:9

        subplot(3,3,i)
        p = get(gca, 'Position');

        if i == 7

            p(1) = p(1) - (p(4)*1.2) + 0.25; 

        elseif i == 8

            p(1) = p(1) - (p(4)*1.2) + 0.20;

        else   

            p(1) = p(1) - (p(4)*1.2) + 0.15; 

        end 

        p(4) = p(4) / 1.2;
        
        set(gca, 'Position', p);

    end

    for i=4:6

        subplot(3,3,i)
        p = get(gca, 'Position');

        if i == 4 

            p(1) = p(1) - (p(4)*1.2) + 0.25;

        elseif i == 5

            p(1) = p(1) - (p(4)*1.2) + 0.20;

        else   

            p(1) = p(1) - (p(4)*1.2) + 0.15; 

        end    

        p(2) = p(2) - p(4)*0.5;
        p(4) = p(4) / 1.2;
        set(gca, 'Position', p);

    end 

    for i=1:3

        subplot(3,3,i)
        p = get(gca, 'Position');

        if i == 1

            p(1) = p(1) - (p(4)*1.2) + 0.25; 

        elseif i == 2

            p(1) = p(1) - (p(4)*1.2) + 0.20;

        else   

            p(1) = p(1) - (p(4)*1.2) + 0.15; 

        end    
       
        p_diff = p(4) * 1;
        p(4) = p(4) + p_diff;
        p(2) = p(2) - p_diff;
        set(gca, 'Position', p);

    end

    colorbar('Position', [p(1)+p(3)+0.02  p(2)-0.38  0.05  p(3)*3.8]);
    
    % Saving the heatmaps
    filename_heatmap = char(center_celltype);
    saveas(gcf,[folderName2,'/',filename_heatmap,'.pdf'])

end


end