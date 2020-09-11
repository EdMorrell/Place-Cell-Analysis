close all
clear

%   --- Calculates change in centre of mass of PFs between conditions

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
WT_CoM = [];
Het_CoM = [];

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
            
            %Computes Centre of Mass for each cell
            COM_array_Pos1 = [];
            for iCell = 1:size(Pos1_mean_frs,1)
                pos_fr = 0;
                for iBin = 1:size(Pos1_mean_frs(iCell,:),2)
                    ipfr = iBin*Pos1_mean_frs(iCell,iBin);
                    pos_fr = pos_fr + ipfr;
                end
                COM_array_Pos1(iCell) = pos_fr/nansum(Pos1_mean_frs(iCell,:));
            end
            %Convert to cm
            COM_array_Pos1 = (COM_array_Pos1/size(Pos1_mean_frs,2))*175;
            
            COM_array_Pos2 = [];
            for iCell = 1:size(Pos2_mean_frs,1)
                pos_fr = 0;
                for iBin = 1:size(Pos2_mean_frs(iCell,:),2)
                    ipfr = iBin*Pos2_mean_frs(iCell,iBin);
                    pos_fr = pos_fr + ipfr;
                end
                COM_array_Pos2(iCell) = pos_fr/nansum(Pos2_mean_frs(iCell,:));
            end
            %Convert to cm
            COM_array_Pos2 = (COM_array_Pos2/size(Pos2_mean_frs,2))*175;
            
            COM_Change_Array = abs(COM_array_Pos2 - COM_array_Pos1);
            
            eval(sprintf('%s_CoM = [%s_CoM COM_Change_Array];',...
                genotype{iGene},genotype{iGene}))
            
        end
    end
end

%% Plotter 
% Mean Field Width
[mean_array,SEM_array] = Plot_SampleMeans(WT_CoM,...
    Het_CoM,1,'WT','Het','Change in Centre of Mass',1,0);

%Stats
[hCoM, pCoM, n] = MeanStats(WT_CoM, ...
    Het_CoM);

%% Rate Change Histogram
num_bins = 20;

WT_CoM(WT_CoM == 0) = [];
Het_CoM(Het_CoM == 0) = [];

sb_CoM = min(max(WT_CoM,max(Het_CoM)));
edges = linspace(0,sb_CoM,num_bins);

[W_NCells] = histcounts(WT_CoM,edges);
[H_NCells] = histcounts(Het_CoM,edges);

% WT_PCells = W_NCells / size(WT_CoM,2);
% Het_PCells = H_NCells / size(Het_CoM,2);

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

edges = linspace(0,sb_CoM,size(WT_PCells,2)+1);

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
