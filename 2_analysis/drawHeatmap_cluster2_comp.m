function drawHeatmap_cluster(outputFolder,center_celltype,cell_id,matchComb,neiCombUnique,matchMarker_var,clusterSelect,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data,log_data_zvalue,gene_name,DEGnumber)

%load 'colormap_2to19grey.mat'
     
criterionClusterIndex=-1;
cellSelect1=[];
for idx=1:size(clusterSelect,2)
    
    clusterIndex=clusterSelect(idx);
    combiTemp=neiCombUnique(clusterIndex);
    
    if combiTemp == center_celltype
        criterionClusterIndex=clusterIndex;
        cellSelect1=matchComb==criterionClusterIndex;
    end
    
end
    
max_leng=size(pvalue_total,1);
DEGindex=zeros(size(gene_name,1),max_leng);
DEGindex(:,criterionClusterIndex)=pvalue_total{criterionClusterIndex}<pCutoff & logRatio_total{criterionClusterIndex}>lrCutoff;

geneIndexTemp=find(DEGindex(:,criterionClusterIndex) & sum(DEGindex(:,criterionClusterIndex),2)==1);

%%%%% Artificial spots: printing out the values of log-ratio and p-val %%%%
% logRatio_artificial = logRatio_total{criterionClusterIndex};
% logRatio_artificial = logRatio_artificial(geneIndexTemp);
% pVal_artificial = pvalue_total{criterionClusterIndex};
% pVal_artificial = pVal_artificial(geneIndexTemp);
% 
% path_DEGs='results/heatmap/';
% path2_DEGs = [path_DEGs, center_celltype, '/'];
% fname = sprintf('%s_artificial.txt',neiCombUnique(criterionClusterIndex));
% writematrix([gene_name(geneIndexTemp), logRatio_artificial, pVal_artificial],[path2_DEGs,fname]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% printing out prospective DEGs
% prospectiveDEGs = gene_name(geneIndexTemp);
% celltypes = [];
% for i=1:size(prospectiveDEGs,1)
%     celltypes = [celltypes; neiCombUnique(criterionClusterIndex)];
% end    
% prospective = [celltypes, prospectiveDEGs];
% 
% path_prospective='results/heatmap/';
% path2_prospective = [path_prospective, center_celltype, '/'];
% fname = sprintf('%s_prospective_DEGs.txt',neiCombUnique(criterionClusterIndex));
% writematrix(prospective,[path2_prospective,fname]);
%%%%%

[~,sortIndex]=sort(logRatio_total{criterionClusterIndex}(geneIndexTemp),'descend');

% if size(geneIndexTemp,1)>DEGnumber
%     geneIndexTemp2=geneIndexTemp(sortIndex(1:DEGnumber));
% else
%     geneIndexTemp2=geneIndexTemp(sortIndex);
% end

%disp(size(geneIndexTemp2,1))
geneIndexTemp3=geneIndexTemp(sortIndex);

% investigating the whole genes without generating artificial doublets %%%%
% geneIndexTemp3=[];
% for i=1:size(log_data,1)
%     geneIndexTemp3=[geneIndexTemp3; i];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
cellSelect2=matchComb==clusterSelect(2);
cellSelect3=matchComb==clusterSelect(3);
geneIndexTemp2=[];
logRatio_DEGs=[];
logRatio2_DEGs=[];
pVal_DEGs=[];
pVal2_DEGs=[];
for i=1:size(geneIndexTemp3,1)
    
    gene_idx = geneIndexTemp3(i);
    if sum(cellSelect2) < 30
        % Wilcoxon-rank sum test   
        [p,h,stats] = ranksum(log_data(gene_idx,find(cellSelect2)),log_data(gene_idx,find(cellSelect1)));
        logRatio=mean(log_data(gene_idx,cellSelect1)+1)-mean(log_data(gene_idx,cellSelect2)+1);
    else
        if sum(cellSelect2) <= sum(cellSelect1)
            h_var = vartest2(log_data(gene_idx,find(cellSelect2)),log_data(gene_idx,find(cellSelect1)));
        else
            % top_sampled: control panel %-----------------------------
%             top_sampled = log_data(gene_idx,find(cellSelect2));
%             top_sampled = sort(top_sampled,'descend');
%             top_sampled = top_sampled(:, 1:sum(cellSelect1));
            
            % if you do not want to select top samples
            top_sampled = log_data(gene_idx,find(cellSelect2));
            
            h_var = vartest2(top_sampled,log_data(gene_idx,find(cellSelect1)));
            %---------------------------------------------
        end
        
        if h_var == 0 || h_var == 1
            if h_var == 1 
                % Welch-t-test: uneqaul variance
                if sum(cellSelect2) <= sum(cellSelect1)
                    [h,p,ci,stats]=ttest2(log_data(gene_idx,find(cellSelect2)),log_data(gene_idx,find(cellSelect1)),'Vartype','unequal');
                else
                    % top_sampled %-----------------------------
                    [h,p,ci,stats2]=ttest2(top_sampled,log_data(gene_idx,find(cellSelect1)),'Vartype','unequal');
                    % ------------------------------------------
                end
                
            elseif h_var == 0
                % Student's t-test: equal variance
                if sum(cellSelect2) <= sum(cellSelect1)
                    [h,p,ci,stats]=ttest2(log_data(gene_idx,find(cellSelect2)),log_data(gene_idx,find(cellSelect1)));
                else
                    % top_sampled %-----------------------------
                    [h,p,ci,stats2]=ttest2(top_sampled,log_data(gene_idx,find(cellSelect1)));
                    % ------------------------------------------
                end
                
            end    
        end
        
        if sum(cellSelect2) <= sum(cellSelect1)
            logRatio=mean(log_data(gene_idx,cellSelect1)+1)-mean(log_data(gene_idx,cellSelect2)+1);
        else
            logRatio=mean(log_data(gene_idx,cellSelect1)+1)-mean(top_sampled+1); 
        end    
    end  
    
    
    if sum(cellSelect3) < 30
        % Wilcoxon-rank sum test   
        [p2,h2,stats2] = ranksum(log_data(gene_idx,find(cellSelect3)),log_data(gene_idx,find(cellSelect1)));
        logRatio2=mean(log_data(gene_idx,cellSelect1)+1)-mean(log_data(gene_idx,cellSelect3)+1);  
    else
        
        if sum(cellSelect3) <= sum(cellSelect1)
            h_var2 = vartest2(log_data(gene_idx,find(cellSelect3)),log_data(gene_idx,find(cellSelect1)));
        else
            % top_sampled2: control panel %-----------------------------
%             top_sampled2 = log_data(gene_idx,find(cellSelect3));
%             top_sampled2 = sort(top_sampled2,'descend');
%             top_sampled2 = top_sampled2(:, 1:sum(cellSelect1));
            
            % if you do not want to select top samples
            top_sampled2 = log_data(gene_idx,find(cellSelect3));
            
            h_var2 = vartest2(top_sampled2,log_data(gene_idx,find(cellSelect1)));
            %---------------------------------------------
        end
        
        if h_var2 == 0 || h_var2 == 1
            if h_var2 == 1 
                % Welch-t-test: uneqaul variance
                if sum(cellSelect3) <= sum(cellSelect1)
                    [h2,p2,ci,stats2]=ttest2(log_data(gene_idx,find(cellSelect3)),log_data(gene_idx,find(cellSelect1)),'Vartype','unequal');
                else
                    % top_sampled %-----------------------------
                    [h2,p2,ci,stats2]=ttest2(top_sampled2,log_data(gene_idx,find(cellSelect1)),'Vartype','unequal');
                    % ------------------------------------------
                end
                
            elseif h_var2 == 0
                % Student's t-test: equal variance
                if sum(cellSelect3) <= sum(cellSelect1)
                    [h2,p2,ci,stats2]=ttest2(log_data(gene_idx,find(cellSelect3)),log_data(gene_idx,find(cellSelect1)));
                else
                    % top_sampled %-----------------------------
                    [h2,p2,ci,stats2]=ttest2(top_sampled2,log_data(gene_idx,find(cellSelect1)));
                    % ------------------------------------------
                end    
            end    
        end
        
        if sum(cellSelect3) <= sum(cellSelect1)
            logRatio2=mean(log_data(gene_idx,cellSelect1)+1)-mean(log_data(gene_idx,cellSelect3)+1);  
        else 
            logRatio2=mean(log_data(gene_idx,cellSelect1)+1)-mean(top_sampled2+1);
        end
        
    end 
     
    if p < pCutoff & logRatio > lrCutoff & p2 < pCutoff & logRatio2 > lrCutoff
        print_out = [neiCombUnique(criterionClusterIndex),gene_name(gene_idx)];
        disp(print_out)
        geneIndexTemp2 = [geneIndexTemp2, gene_idx];
        logRatio_DEGs = [logRatio_DEGs, logRatio];
        logRatio2_DEGs = [logRatio2_DEGs, logRatio2];
        pVal_DEGs = [pVal_DEGs, p];
        pVal2_DEGs = [pVal2_DEGs, p2];
    end
    
end

%%% homotypic1, homotypic2, heterotypic data files: barcodes, log-data, z-scores
path='results/heatmap/';
path2 = [path, center_celltype, '/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Homotypic spots: printing out the values of log-ratio and p-val %%%%%
% fname2 = sprintf('%s_homotypic.txt',neiCombUnique(criterionClusterIndex));
% writematrix([gene_name(geneIndexTemp2), logRatio_DEGs', pVal_DEGs', logRatio2_DEGs', pVal2_DEGs'],[path2_DEGs,fname2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geneIndexTemp2 = geneIndexTemp2';
if size(geneIndexTemp2,1) > 0
    %cellSelect2=matchComb==clusterSelect(2);
    
    logRatioMultip = [];
    for x=1:size(logRatio_DEGs,2)
        val = logRatio_DEGs(x)*logRatio2_DEGs(x);
        logRatioMultip = [logRatioMultip, val];
    end
    [~,sortIndex]=sort(logRatioMultip,'descend');
    geneIndexTemp2=geneIndexTemp2(sortIndex);
    
%     if size(geneIndexTemp2,1) > 4
%          geneIndexTemp2=geneIndexTemp2(1:4);
%     end
    
    % fig1
    %geneIndexTemp2 = [2077]; % embryo
    %geneIndexTemp2 = [2004]; %hippocampus
     
    close all
    figure(1)
    
    
    
    subplot(3,3,1);
    imagesc(log_data_zvalue(geneIndexTemp2,find(cellSelect1)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:size(geneIndexTemp2)
        filename = sprintf('%s_%s.txt',[neiCombUnique(criterionClusterIndex),gene_name(geneIndexTemp2(i))]);
        log_values=log_data(geneIndexTemp2(i),find(cellSelect1));
        zvalues=log_data_zvalue(geneIndexTemp2(i),find(cellSelect1));
        barcodes=cell_id(find(cellSelect1));
        integrated=[barcodes; log_values; zvalues];
        integrated=integrated';
        integrated=sortrows(integrated,2);
        writematrix(integrated,[path2,filename]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',12)
    set(gca,'YTickLabelMode','auto')
    yticks(1:size(geneIndexTemp2,1))
    yticklabels(gene_name(geneIndexTemp2))
    
    if size(geneIndexTemp2,1) > 30
        ax = gca;
        ax.YAxis.FontSize = 6;
    end
    
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    
    %colormap 'jet'
    %colorbar
    set(gca, 'FontName', 'Arial','DefaultAxesTitleFontWeight','normal') %Sans Serif
    t1=sprintf('Endothelial &\nLens cells');
    title({neiCombUnique(criterionClusterIndex);['(n=',num2str(length(find(cellSelect1))),')']},'FontSize',18); %t1
    %get(gca,'FontName')
    %title([neiCombUnique(criterionClusterIndex),length(find(cellSelect1))],'FontSize',12)
    
    
    subplot(3,3,2);
    imagesc(log_data_zvalue(geneIndexTemp2,find(cellSelect2)))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:size(geneIndexTemp2)
        filename = sprintf('%s_%s.txt',[neiCombUnique(clusterSelect(2)),gene_name(geneIndexTemp2(i))]);
        log_values=log_data(geneIndexTemp2(i),find(cellSelect2));
        zvalues=log_data_zvalue(geneIndexTemp2(i),find(cellSelect2));
        barcodes=cell_id(find(cellSelect2));
        integrated=[barcodes; log_values; zvalues];
        integrated=integrated';
        integrated=sortrows(integrated,2);
        writematrix(integrated,[path2,filename]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    %colormap 'jet'
    %colorbar
    set(gca, 'FontName', 'Arial','DefaultAxesTitleFontWeight','normal')
    t2=sprintf('Endothelial\ncells');
    title({neiCombUnique(clusterSelect(2));['(n=',num2str(length(find(cellSelect2))),')']},'FontSize',18);%t2
    %get(gca,'FontName')
    %title([neiCombUnique(clusterSelect(2)),length(find(cellSelect2))],'FontSize',12)

    %cellSelect3=matchComb==clusterSelect(3);
    subplot(3,3,3);
    imagesc(log_data_zvalue(geneIndexTemp2,find(cellSelect3)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:size(geneIndexTemp2)
        filename = sprintf('%s_%s.txt',[neiCombUnique(clusterSelect(3)),gene_name(geneIndexTemp2(i))]);
        log_values=log_data(geneIndexTemp2(i),find(cellSelect3));
        zvalues=log_data_zvalue(geneIndexTemp2(i),find(cellSelect3));
        barcodes=cell_id(find(cellSelect3));
        integrated=[barcodes; log_values; zvalues];
        integrated=integrated';
        integrated=sortrows(integrated,2);
        writematrix(integrated,[path2,filename]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    %colormap 'jet'
    %colorbar
    set(gca, 'FontName', 'Arial','DefaultAxesTitleFontWeight','normal')
    t3=sprintf('Lens\ncells');
    title({neiCombUnique(clusterSelect(3));['(n=',num2str(length(find(cellSelect3))),')']},'FontSize',18);%t3
    %get(gca,'FontName')
    %title([neiCombUnique(clusterSelect(3)),length(find(cellSelect3))],'FontSize',12)
 


    %=====================
    % cell type markers    
    fns = fieldnames(matchMarker_var);
    cellSelect2FinalIDX=[];
    for var_i=1:size(fns,1)
        markerIDX = matchMarker_var.(fns{var_i})==neiCombUnique(2);
        tmp = find(markerIDX==1);
        % cellSelect2FinalIDX = union(cellSelect2FinalIDX, tmp);
        cellSelect2FinalIDX=[cellSelect2FinalIDX; tmp];
    end
    
    % when matchMarker has only one column
    %markerIDX = matchMarker==neiCombUnique(2);
    %cellSelect2FinalIDX = find(markerIDX==1);
    
    subplot(3,3,4); % heatmap with only cell type markers: subplot(3,3,1);
    imagesc(log_data_zvalue(cellSelect2FinalIDX,cellSelect1))
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',18)
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
    %title({[neiCombUnique(criterionClusterIndex)];['(n=',num2str(length(find(cellSelect1))),')']},'FontSize',12); % heatmap with only cell type markers
    %get(gca,'FontName')
    %colormap 'jet'
    
    subplot(3,3,5);  % heatmap with only cell type markers: subplot(3,3,2);
    imagesc(log_data_zvalue(cellSelect2FinalIDX,cellSelect2))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    %title({[neiCombUnique(clusterSelect(2))];['(n=',num2str(length(find(cellSelect2))),')']},'FontSize',12); % heatmap with only cell type markers
    %get(gca,'FontName')
    %colormap 'jet'
    
    subplot(3,3,6);  % heatmap with only cell type markers: subplot(3,3,3);
    imagesc(log_data_zvalue(cellSelect2FinalIDX,cellSelect3))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    %title({[neiCombUnique(clusterSelect(3))];['(n=',num2str(length(find(cellSelect3))),')']},'FontSize',12); % heatmap with only cell type markers
    %get(gca,'FontName')
    %colormap 'jet'
    
  
    cellSelect3FinalIDX=[];
    for var_i=1:size(fns,1)
        markerIDX = matchMarker_var.(fns{var_i})==neiCombUnique(3);
        tmp = find(markerIDX==1);
        %cellSelect3FinalIDX = union(cellSelect3FinalIDX, tmp);
        cellSelect3FinalIDX = [cellSelect3FinalIDX; tmp];
    end
    % when matchMarker has only one column
    %markerIDX = matchMarker==neiCombUnique(3);
    %cellSelect3FinalIDX = find(markerIDX==1);
    
%     cellSelect3LogVal = sum(log_data(:, find(cellSelect3)),2);
%     [~,cellSelect3SortedIDX] = sort(cellSelect3LogVal,'descend');
%     cellSelect3FinalIDX = cellSelect3SortedIDX(1:DEGnumber-9);
    
    subplot(3,3,7); % heatmap with only cell type markers: subplot(3,3,4);
    imagesc(log_data_zvalue(cellSelect3FinalIDX,cellSelect1))
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',18)
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
    %get(gca,'FontName')
    %colormap 'jet'
    
    subplot(3,3,8); % heatmap with only cell type markers: subplot(3,3,5);
    imagesc(log_data_zvalue(cellSelect3FinalIDX,cellSelect2))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    %get(gca,'FontName')
    %colormap 'jet'
   
    subplot(3,3,9); % heatmap with only cell type markers: subplot(3,3,6);
    imagesc(log_data_zvalue(cellSelect3FinalIDX,cellSelect3))
    yticks([])
    xticks([])
    caxis([-3 3])
    colormap(redblue(256));
    set(gca, 'FontName', 'Arial')
    %get(gca,'FontName')
    %colormap 'jet'
    
    set(gcf, 'Position', [100, 100, 550, 650])
    for i=7:9 % heatmap with only cell type markers: i=4:6
        subplot(3,3,i)
        p = get(gca, 'Position');
        if i == 7 % heatmap with only cell type markers: i=4
            p(1) = p(1) - (p(4)*1.2) + 0.25; 
        elseif i == 8 % heatmap with only cell type markers: i=5
            p(1) = p(1) - (p(4)*1.2) + 0.20;
        else   
            p(1) = p(1) - (p(4)*1.2) + 0.15; 
        end    
        p(4) = p(4) / 1.2;
        
        set(gca, 'Position', p);
    end
    for i=4:6 % heatmap with only cell type markers: i=1:3
        subplot(3,3,i)
        p = get(gca, 'Position');
        if i == 4 % heatmap with only cell type markers: i=1
            p(1) = p(1) - (p(4)*1.2) + 0.25; 
        elseif i == 5 % heatmap with only cell type markers: i=2
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
        
%         p(2) =  p(2) - p(4); % fig1_heatmap
%         p(4) = p(4) / 1.2; % fig1_heatmap
        p_diff = p(4) * 1;
        p(4) = p(4) + p_diff;
        p(2) = p(2) - p_diff;
        set(gca, 'Position', p);
    end

    % [x y width height] [p(1)+p(3)+0.02  p(2)+0.005  0.05  p(3)*2]
    %colorbar('Position', [p(1)+p(3)+0.02  p(2)-0.191  0.02  p(3)*1.74]); % heatmap with only cell type markers
    %colorbar('Position', [p(1)+p(3)+0.02  p(2)-0.383  0.05  p(3)*2.641]); % fig1_heatmap
    colorbar('Position', [p(1)+p(3)+0.02  p(2)-0.38  0.05  p(3)*3.8]);
    
    %=====================
    
    
%     set(gcf, 'Position', [100, 100, 550, 650])
%     for i=1:3
%         subplot(3,3,i)
%         p = get(gca, 'Position');
%         p_diff = p(4) * 1;
%         p(4) = p(4) + p_diff;
%         p(2) = p(2) - p_diff;
%         set(gca, 'Position', p);
%     end

    % filename_heatmap = char(center_celltype);
    % saveas(gcf,[outputFolder,'/',filename_heatmap,'.pdf'])
end
    
end
