close all
clear

%   --- Calculates Place Field Correlations between 2 conditions

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

%% Compute PF Spatial Correlation
WT_SpatialCorr = [];
Het_SpatialCorr = [];

for iGene = 1:size(genotype,1)
    
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
            
            %Computes spatial correlation between each firing rate curve
            xcorr_array = [];
            for iCell = 1:size(Pos1_mean_frs,1)
                P1_FRz = zscore(Pos1_mean_frs(iCell,:));
                P2_FRz = zscore(Pos2_mean_frs(iCell,:));
                xcorr_array(iCell) = corr2(P1_FRz,P2_FRz);
                clear P1_FRz P2_FRz
            end
            
            eval(sprintf('%s_SpatialCorr = [%s_SpatialCorr xcorr_array];',...
                genotype{iGene},genotype{iGene}))
            
        end
    end
end

%% Plotter 
% Mean Field Width
[mean_array,SEM_array] = Plot_SampleMeans(WT_SpatialCorr,...
    Het_SpatialCorr,1,'WT','Het','Place Field Spatial Correlation',1,0);

%Stats
[hXC, pXC, nXC] = MeanStats(WT_SpatialCorr, ...
    Het_SpatialCorr);

%% Rate Change Histogram
num_bins = 20;

WT_SpatialCorr(WT_SpatialCorr == 0) = [];
Het_SpatialCorr(Het_SpatialCorr == 0) = [];

sb_XC = min(max(WT_SpatialCorr,max(Het_SpatialCorr)));
edges = linspace(0,sb_XC,num_bins);

[W_NCells] = histcounts(WT_SpatialCorr,edges);
[H_NCells] = histcounts(Het_SpatialCorr,edges);

% WT_PCells = W_NCells / size(WT_SpatialCorr,2);
% Het_PCells = H_NCells / size(Het_SpatialCorr,2);

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

edges = linspace(0,sb_XC,size(WT_PCells,2)+1);

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