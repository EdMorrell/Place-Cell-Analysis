close all
clear

%   --- Calculate rate remapping index

fn_pos1 = 'D:\Ed\Data\Matlab Outputs\Single Units\Place Cell Stats\Novel Track\Pos1';
fn_pos2 = 'D:\Ed\Data\Matlab Outputs\Single Units\Place Cell Stats\Novel Track\Pos2';

addpath('D:\Ed\Scripts\Tools')

load([fn_pos1 filesep 'L_Track_Stats.mat']);
L_Track_Stats1 = L_Track_Stats;
clear L_Track_Stats

load([fn_pos2 filesep 'L_Track_Stats.mat']);
L_Track_Stats2 = L_Track_Stats;
clear L_Track_Stats

genotype = {'WT';'Het'};
Dir = {'LR';'RL'};

%Only incl. verified PCs
PC_Only = true;

%Constrain Firing Rate - excludes cells with FR's lower than threshold (in
%both environments)
constr_FR = true;
FR_Th = 0.27;

%% Compute Rate Remapping Index for each cell
WT_R_RM_Index = [];
Het_R_RM_Index = [];

WT_Mean_FR_Pos1 = [];
WT_Mean_FR_Pos2 = [];
Het_Mean_FR_Pos1 = [];
Het_Mean_FR_Pos2 = [];
for iGene = 1:size(genotype,1)
    
    eval(sprintf('%s_R_RM_Index = [];',genotype{iGene}))
    
    eval(sprintf('An_Num = numel(fieldnames(L_Track_Stats1.%s));',...
        genotype{iGene}))
    for iAnimal = 1:An_Num    
        eval(sprintf('f_names = fieldnames(L_Track_Stats1.%s);',...
            genotype{iGene}))
        An_Name = f_names{iAnimal,1};
        
        Pos1_mean_frs = [];
        Pos2_mean_frs = [];
        for iDir = 1:size(Dir,1)
            eval(sprintf('Pos1_mean_frs = L_Track_Stats1.%s.%s.%s.mean_frs;',...
                genotype{iGene},An_Name,Dir{iDir}))
            eval(sprintf('Pos2_mean_frs = L_Track_Stats2.%s.%s.%s.mean_frs;',...
                genotype{iGene},An_Name,Dir{iDir}))
            
            %If runs in one direction missing removes from analysis
            if size(Pos1_mean_frs,1) ~= size(Pos2_mean_frs,1)
                continue
            end

            %Removes non-PCs if true
            if PC_Only
                eval(sprintf('pvPos1 = L_Track_Stats1.%s.%s.%s.PC_Ver;',...
                    genotype{iGene},An_Name,Dir{iDir}))
                eval(sprintf('pvPos2 = L_Track_Stats2.%s.%s.%s.PC_Ver;',...
                    genotype{iGene},An_Name,Dir{iDir}))
                Pos1_PC = find(pvPos1);
                Pos2_PC = find(pvPos2);
                all_PC = [Pos1_PC;Pos2_PC];
                all_PC = sortrows(all_PC);
                pv = unique(all_PC);
                Pos1_mean_frs = Pos1_mean_frs(find(pv<=size(Pos1_mean_frs,1)),:);
                Pos2_mean_frs = Pos2_mean_frs(find(pv<=size(Pos2_mean_frs,1)),:);
            end
            
            if constr_FR               
                Pos1_A_TH = find(nanmean(Pos1_mean_frs,2)>FR_Th);
                Pos2_A_TH = find(nanmean(Pos2_mean_frs,2)>FR_Th);
                all_AT = [Pos1_A_TH;Pos2_A_TH];
                all_AT = sortrows(all_AT);
                A_Th_Ind = unique(all_AT);
                Pos1_mean_frs = Pos1_mean_frs(A_Th_Ind,:);
                Pos2_mean_frs = Pos2_mean_frs(A_Th_Ind,:);
            end
                       
            F1 = nanmean(Pos1_mean_frs,2);
          
            F2 = nanmean(Pos2_mean_frs,2);

            RR_Scores = abs((F1-F2)./(F1+F2));
            
            eval(sprintf('%s_Mean_FR_Pos1 = [%s_Mean_FR_Pos1 F1''];',...
                genotype{iGene},genotype{iGene}))
            eval(sprintf('%s_Mean_FR_Pos2 = [%s_Mean_FR_Pos2 F2''];',...
                genotype{iGene},genotype{iGene}))

            eval(sprintf('%s_R_RM_Index = [%s_R_RM_Index RR_Scores''];',...
                genotype{iGene},genotype{iGene}));
            
            if PC_Only
                clear pv pv_Pos1 pv_Pos2 Pos1_PC Pos2_PC all_PC F1 F2
            end
                    
        end
    end
end
%% Plotter 
% 
[mean_array,SEM_array] = Plot_SampleMeans(WT_R_RM_Index,...
    Het_R_RM_Index,1,'WT','Het','Rate Remapping Index',1,0);

%Stats
[hR_RM, pR_RM, nR_RM] = MeanStats(WT_R_RM_Index, ...
    Het_R_RM_Index);

%% Rate Change Histogram
num_bins = 20;

WT_R_RM_Index(WT_R_RM_Index == 0) = [];
Het_R_RM_Index(Het_R_RM_Index == 0) = [];

sb_RR = min(max(WT_R_RM_Index,max(Het_R_RM_Index)));
edges = linspace(0,sb_RR,num_bins);

[W_NCells] = histcounts(WT_R_RM_Index,edges);
[H_NCells] = histcounts(Het_R_RM_Index,edges);

% WT_PCells = W_NCells / size(WT_R_RM_Index,2);
% Het_PCells = H_NCells / size(Het_R_RM_Index,2);

WT_PCells = W_NCells / sum(W_NCells);
Het_PCells = H_NCells / sum(H_NCells);

%Cumulative Values
WT_CSum = cumsum(WT_PCells);
Het_CSum = cumsum(Het_PCells);

%KS - Test
[hKS,pKS] = kstest2(WT_PCells,Het_PCells);

%Smooth
[WT_PCells] = smooth_array(WT_PCells,1.5);
[Het_PCells] = smooth_array(Het_PCells,1.5);

edges = linspace(0,sb_RR,size(WT_PCells,2)+1);

x = edges(1:(size(edges,2)-1)) + (mean(diff(edges))/2);
figure; hold on
p = bar(x,WT_PCells, 'k');
g = bar(x,Het_PCells, 'c');
p.LineWidth = 2;
g.LineWidth = 2;
alpha(0.25)

[ax] = plot_prop();

%Cumulative Distribution Plotter
figure; hold on
plot(WT_CSum,'k', 'LineWidth', 1.5)
plot(Het_CSum,'c', 'LineWidth', 1.5)
ax = plot_prop();

%% FR Plotter
for iGene = 1:size(genotype,1)
    eval(sprintf('mean_FR_Array(1,iGene) = nanmean(%s_Mean_FR_Pos1);',...
        genotype{iGene}))
    eval(sprintf(['SEM_array(1,iGene) = nanstd(%s_Mean_FR_Pos1) ',...
        '/ sqrt(size(%s_Mean_FR_Pos1,2));'],genotype{iGene},...
        genotype{iGene}))
    eval(sprintf('mean_FR_Array(2,iGene) = nanmean(%s_Mean_FR_Pos2);',...
        genotype{iGene})) 
    eval(sprintf(['SEM_array(2,iGene) = nanstd(%s_Mean_FR_Pos2) ',...
        '/ sqrt(size(%s_Mean_FR_Pos2,2));'],genotype{iGene},...
        genotype{iGene}))
end

figure;
name = {'Position 1';'Position 2'};
edges = 1:2;
b = bar(edges,mean_FR_Array);
b(1,1).LineWidth = 1.5;
b(1,2).LineWidth = 1.5;
set(gca,'xticklabel',name)
hold on
title('Mean Firing Rate')

% Adds error bars
ngroups = size(mean_FR_Array,1);
nbars = size(mean_FR_Array,2);
% Calculates width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Sets the position of each error bar in the centre of the main bar
for i = 1:nbars
    % Calculates centre of each bar
    x_er = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x_er, mean_FR_Array(:,i),SEM_array(:,i),SEM_array(:,i),...
        'k', ...
        'linestyle', 'none',...
        'linewidth', 1.5);
end
%Plot Properties
ax = plot_prop();

%% Firing Rate Stats

Pos = {'Pos1';'Pos2'};
%Formats data for ANOVA
[FRs,FR_G,FR_P] = Anova_Dat_Struct(WT_Mean_FR_Pos1,...
    WT_Mean_FR_Pos2,Het_Mean_FR_Pos1,...
    Het_Mean_FR_Pos2,genotype,Pos);

[~,~,stats] = anovan(FRs,{FR_G FR_P},'model',2,'varnames',{'Genotype','Environment'});