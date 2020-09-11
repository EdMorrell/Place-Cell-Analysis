close all hidden
clear

%   --- Looks at changes in key PF stats (previously generated)
%       between novel and familiar track

fn_pos1 = 'D:\Ed\Data\Matlab Outputs\Single Units\Place Cell Stats\Familiar Track';
fn_pos2 = 'D:\Ed\Data\Matlab Outputs\Single Units\Place Cell Stats\Novel Track\Pos2';

addpath('D:\Ed\Scripts\Tools')

load([fn_pos1 filesep 'L_Track_Stats.mat']);
L_Track_Stats_Pos1 = L_Track_Stats;
clear L_Track_Stats

load([fn_pos2 filesep 'L_Track_Stats.mat']);
L_Track_Stats_Pos2 = L_Track_Stats;
clear L_Track_Stats

genotype = {'WT';'Het'};
dir = {'LR';'RL'};
Pos = {'Pos1';'Pos2'};

%Only incl. verified PCs
PC_Only = true;

%Constrain to cells with a min Firing Rate 
constr_FR = true;
FR_Th = 0.2;


%%
%Loops through L_Track_Stats to extract field widths and number of place
%fields and places into a super array
for iGene = 1:size(genotype,1)
    for iPos = 1:size(Pos,1)
        eval(sprintf('%s_field_width_array_%s = [];',genotype{iGene},...
            Pos{iPos}))
        eval(sprintf('%s_mean_field_FR_%s = [];',genotype{iGene},...
            Pos{iPos}))
        eval(sprintf('%s_pfield_fr_%s = [];',genotype{iGene},Pos{iPos}))
        eval(sprintf('%s_spat_info_array_%s = [];',genotype{iGene},...
            Pos{iPos}))

        eval(sprintf('An_Num = numel(fieldnames(L_Track_Stats_%s.%s));',...
            Pos{iPos},genotype{iGene}))
        for iAnimal = 1:An_Num    
            eval(sprintf('f_names = fieldnames(L_Track_Stats_%s.%s);',...
                Pos{iPos},genotype{iGene}))
            An_Name = f_names{iAnimal,1};
            for idir = 1:size(dir)
                %Extracts fieldwidths
                eval(sprintf(['fws = L_Track_Stats_%s.%s.%s.',...
                    '%s.field_widths;'],Pos{iPos},genotype{iGene},...
                    An_Name,dir{idir}))

                %Extracts mean firing rates
                eval(sprintf(['m_frs = L_Track_Stats_%s.%s.%s.',...
                    '%s.mean_frs;'],Pos{iPos},genotype{iGene},...
                    An_Name,dir{idir}))

                %Extracts spatial info
                eval(sprintf('si = L_Track_Stats_%s.%s.%s.%s.spatial_info;',...
                    Pos{iPos},genotype{iGene},An_Name,dir{idir}))

                %Removes non-PCs if true (Only removes if PC in neither condition)                
%                 if PC_Only
%                     eval(sprintf('pvPos1 = L_Track_Stats_Pos1.%s.%s.%s.PC_Ver;',...
%                         genotype{iGene},An_Name,dir{idir}))
%                     eval(sprintf('pvPos2 = L_Track_Stats_Pos2.%s.%s.%s.PC_Ver;',...
%                         genotype{iGene},An_Name,dir{idir}))
%                     Pos1_PC = find(pvPos1);
%                     Pos2_PC = find(pvPos2);
%                     all_PC = [Pos1_PC;Pos2_PC];
%                     all_PC = sortrows(all_PC);
%                     pv = unique(all_PC);                    
%                     fws = fws(find(pv<=size(fws,1)));
%                     m_frs = m_frs(find(pv<=size(m_frs,1)),:);
%                     si = si(find(pv<=size(si,1)));
%                 end
                
                %Removes non-PCs if true
                if PC_Only
                    eval(sprintf('pv = L_Track_Stats_%s.%s.%s.%s.PC_Ver;',...
                        Pos{iPos},genotype{iGene},An_Name,dir{idir}))
                    fws = fws(find(pv));
                    m_frs = m_frs(find(pv),:);
                    si = si(find(pv));
                end
                
                if constr_FR
                    %Above threshold index
                    a_th_ind = find(nanmean(m_frs,2)>FR_Th);
                    m_frs = m_frs(a_th_ind,:);
                    fws = fws(a_th_ind);
                    si = si(a_th_ind,:);
                end
                
                %Adds spatial to new array
                eval(sprintf(['%s_spat_info_array_%s = ',...
                    '[%s_spat_info_array_%s si''];'],...
                    genotype{iGene},Pos{iPos},genotype{iGene},...
                    Pos{iPos}))
                    
                for iCell = 1:size(fws,1)
                    
                    %Adds field widths to new array
                    if isempty(fws{iCell,1}) | isnan(fws{iCell,1})
                        fws{iCell,1}(1,1:3) = NaN;
                    end
                    eval(sprintf(['%s_field_width_array_%s = ',...
                        '[%s_field_width_array_%s fws{iCell,1}(:,1)''];'],...
                        genotype{iGene},Pos{iPos},genotype{iGene},...
                        Pos{iPos}))

                    %Adds mean in-field FR to new array
                    eval(sprintf(['%s_mean_field_FR_%s = ',...
                        '[%s_mean_field_FR_%s fws{iCell,1}(:,3)''];'],...
                        genotype{iGene},Pos{iPos},genotype{iGene},...
                        Pos{iPos}))

                    %Find peaks firing rate of each PF and adds to new
                    %array
                    if isnan(fws{iCell,1})
                        cell_pfr = NaN;
                    else
                        cell_pfr = m_frs(iCell,fws{iCell,1}(:,2));
                        eval(sprintf(['%s_pfield_fr_%s = [%s_pfield_fr_%s ',...
                            'cell_pfr];'],genotype{iGene},Pos{iPos},...
                            genotype{iGene},Pos{iPos}))
                    end
                    clear cell_pfr 

                end
                clear fws m_frs si 
                if PC_Only
                    clear pv pvPos1 pvPos2 Pos1_PC Pos2_PC all_PC
                end
            end
        end
    end
end

%% Plotter
names = {'Pos1';'Pos2'};
%% Field Width
for iPos = 1:size(Pos,1)
    for iGene = 1:size(genotype,1)
        eval(sprintf('mean_FW_Array(iPos,iGene) = nanmean(%s_field_width_array_%s);',...
            genotype{iGene},Pos{iPos}))
        eval(sprintf(['SEM_array(iPos,iGene) = nanstd(%s_field_width_array_%s) ',...
            '/ sqrt(size(%s_field_width_array_%s,2));'],genotype{iGene},...
            Pos{iPos},genotype{iGene},Pos{iPos}))
    end
end
[b_hand] = multi_bar_plot(mean_FW_Array,1,SEM_array,'Field Width Change',...
    names);

%Formats data for ANOVA
[FWs,FW_G,FW_P] = Anova_Dat_Struct(WT_field_width_array_Pos1,...
    WT_field_width_array_Pos2,Het_field_width_array_Pos1,...
    Het_field_width_array_Pos2,genotype,Pos);

[~,~,stats] = anovan(FWs,{FW_G FW_P},'model',2,'varnames',{'Genotype','Environment'});
%% Peak In-Field Firing Rate
for iPos = 1:size(Pos,1)
    for iGene = 1:size(genotype,1)
        eval(sprintf('mean_PF_Array(iPos,iGene) = nanmean(%s_pfield_fr_%s);',...
            genotype{iGene},Pos{iPos}))
        eval(sprintf(['SEM_array(iPos,iGene) = nanstd(%s_pfield_fr_%s) ',...
            '/ sqrt(size(%s_pfield_fr_%s,2));'],genotype{iGene},...
            Pos{iPos},genotype{iGene},Pos{iPos}))
    end
end
[b_hand] = multi_bar_plot(mean_PF_Array,1,SEM_array,'Peak Field FR Change',...
    names);

%Formats data for ANOVA
[PFfrs,PF_G,PF_P] = Anova_Dat_Struct(WT_pfield_fr_Pos1,...
    WT_pfield_fr_Pos2,Het_pfield_fr_Pos1,...
    Het_pfield_fr_Pos2,genotype,Pos);

[~,~,stats] = anovan(PFfrs,{PF_G PF_P},'model',2,'varnames',{'Genotype','Environment'});
%% Mean In-Field Firing Rate
for iPos = 1:size(Pos,1)
    for iGene = 1:size(genotype,1)
        eval(sprintf('mean_MF_Array(iPos,iGene) = nanmean(%s_mean_field_FR_%s);',...
            genotype{iGene},Pos{iPos}))
        eval(sprintf(['SEM_array(iPos,iGene) = nanstd(%s_mean_field_FR_%s) ',...
            '/ sqrt(size(%s_mean_field_FR_%s,2));'],genotype{iGene},...
            Pos{iPos},genotype{iGene},Pos{iPos}))
    end
end
[b_hand] = multi_bar_plot(mean_MF_Array,1,SEM_array,'Mean Field FR Change',...
    names);

%Formats data for ANOVA
[FRs,FR_G,FR_P] = Anova_Dat_Struct(WT_mean_field_FR_Pos1,...
    WT_mean_field_FR_Pos2,Het_mean_field_FR_Pos1,...
    Het_mean_field_FR_Pos2,genotype,Pos);

[~,~,stats] = anovan(FRs,{FR_G FR_P},'model',2,'varnames',{'Genotype','Environment'});
%% Mean Spatial Info
for iPos = 1:size(Pos,1)
    for iGene = 1:size(genotype,1)
        eval(sprintf('mean_SI_Array(iPos,iGene) = nanmean(%s_spat_info_array_%s);',...
            genotype{iGene},Pos{iPos}))
        eval(sprintf(['SEM_array(iPos,iGene) = nanstd(%s_spat_info_array_%s) ',...
            '/ sqrt(size(%s_spat_info_array_%s,2));'],genotype{iGene},...
            Pos{iPos},genotype{iGene},Pos{iPos}))
    end
end
[b_hand] = multi_bar_plot(mean_SI_Array,1,SEM_array,'Mean Spatial Info Change',...
    names);

%Formats data for ANOVA
[SIs,SI_G,SI_P] = Anova_Dat_Struct(WT_spat_info_array_Pos1,...
    WT_spat_info_array_Pos2,Het_spat_info_array_Pos1,...
    Het_spat_info_array_Pos2,genotype,Pos);

[~,~,stats] = anovan(SIs,{SI_G SI_P},'model',2,'varnames',{'Genotype','Environment'});
%results = multcompare(stats,'Dimension',[1 2]);
%%

% for iPos = 1:size(Pos,1)
%     eval(sprintf('WT_mean_field_FR = WT_mean_field_FR_%s;',Pos{iPos}))
%     eval(sprintf('Het_mean_field_FR = Het_mean_field_FR_%s;',Pos{iPos}))
%     
%     WT_mean_field_FR(WT_mean_field_FR == 0) = [];
%     Het_mean_field_FR(Het_mean_field_FR == 0) = [];
% 
%     sb_SI = min(max(WT_mean_field_FR,max(Het_mean_field_FR)));
%     edges = linspace(0,sb_SI,20);
% 
%     [W_NCells] = histcounts(WT_mean_field_FR,edges);
%     [H_NCells] = histcounts(Het_mean_field_FR,edges);
% 
%     WT_PCells = W_NCells / size(WT_mean_field_FR,2);
%     Het_PCells = H_NCells / size(Het_mean_field_FR,2);
% 
%     %Smooth
%     [WT_PCells] = smooth_array(WT_PCells,1.5);
%     [Het_PCells] = smooth_array(Het_PCells,1.5);
% 
%     edges = linspace(0,sb_SI,size(WT_PCells,2)+1);
% 
%     x = edges(1:(size(edges,2)-1)) + (mean(diff(edges))/2);
%     figure; hold on
%     p = plot(x,WT_PCells, 'k');
%     g = plot(x,Het_PCells, 'c');
%     p.LineWidth = 2;
%     g.LineWidth = 2;
% 
%     [ax] = plot_prop();
%     
%     clear WT_mean_field_FR Het_mean_field_FR
% end