close all 
clear

% --- Mean Directionality of Each Place Cell
%           1 - Indicates complete bidrectionality (only active in one
%               direction)
%           2 - Indicates complete unidirectionality (Equally active in
%               both directions)

fn_TStats = 'D:\Ed\Data\Matlab Outputs\Single Units\Place Cell Stats\Familiar Track\L_Track_Stats.mat';

load(fn_TStats)

genotype = {'WT';'Het'};
dir = {'LR';'RL'};

%Only incl. verified PCs
PC_Only = true;

%% Compute BiDir Index for each cell
for iGene = 1:size(genotype,1)
    
    eval(sprintf('%s_BD_Index = [];',genotype{iGene}))
    
    eval(sprintf('An_Num = numel(fieldnames(L_Track_Stats.%s));',...
        genotype{iGene}))
    for iAnimal = 1:An_Num    
        eval(sprintf('f_names = fieldnames(L_Track_Stats.%s);',...
            genotype{iGene}))
        An_Name = f_names{iAnimal,1};
        
        eval(sprintf('LR_mean_frs = L_Track_Stats.%s.%s.LR.mean_frs;',...
            genotype{iGene},An_Name))
        eval(sprintf('RL_mean_frs = L_Track_Stats.%s.%s.RL.mean_frs;',...
            genotype{iGene},An_Name))
        
        %Removes non-PCs if true
        if PC_Only
            eval(sprintf('pvLR = L_Track_Stats.%s.%s.LR.PC_Ver;',...
                genotype{iGene},An_Name))
            eval(sprintf('pvRL = L_Track_Stats.%s.%s.RL.PC_Ver;',...
                genotype{iGene},An_Name))
            if ~isnan(pvLR) & ~isnan(pvRL)
                LR_PC = find(pvLR);
                RL_PC = find(pvRL);
                all_PC = [LR_PC;RL_PC];
                all_PC = sortrows(all_PC);
                pv = unique(all_PC);
                LR_mean_frs = LR_mean_frs(find(pv),:);
                RL_mean_frs = RL_mean_frs(find(pv),:);
            end
        end
        
        F1 = mean(LR_mean_frs,2);
          
        F2 = mean(RL_mean_frs,2);
        
        if ~isnan(F1) & ~isnan(F2)
        
            BD_Scores = abs((F1-F2)./(F1+F2));

            eval(sprintf('%s_BD_Index = [%s_BD_Index BD_Scores''];',...
                genotype{iGene},genotype{iGene}));
            
        end
        
        if PC_Only
            clear pv
        end
    end
end
%% Plotter 
% Mean Field Width
[mean_array,SEM_array] = Plot_SampleMeans(WT_BD_Index,...
    Het_BD_Index,1,'WT','Het','Bi-Directionality Index',1,0);

%Stats
[hBD, pBD, nBD] = MeanStats(WT_BD_Index, ...
    Het_BD_Index);