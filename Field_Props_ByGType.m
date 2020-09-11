close all 
clear

fn_TStats = 'D:\Ed\Data\Matlab Outputs\Single Units\Place Cell Stats\Familiar Track\L_Track_Stats.mat';

load(fn_TStats)

genotype = {'WT';'Het'};
dir = {'LR';'RL'};

%Only incl. verified PCs
PC_Only = true;

%Remove outliers?
rm_out = false;

%Constrain to cells with a min Firing Rate 
constr_FR = true;
FR_Th = 0.2;

%%
%Loops through L_Track_Stats to extract field widths and number of place
%fields and places into a super array
for iGene = 1:size(genotype,1)
    
    eval(sprintf('%s_field_width_array = [];',genotype{iGene}))
    eval(sprintf('%s_mean_field_FR = [];',genotype{iGene}))
    eval(sprintf('%s_field_num_array = zeros(1,10);',genotype{iGene}))
    eval(sprintf('%s_pfield_fr = [];',genotype{iGene}))
    eval(sprintf('%s_spat_info_array = [];',genotype{iGene}))
    eval(sprintf('An_Num = numel(fieldnames(L_Track_Stats.%s));',...
        genotype{iGene}))
    for iAnimal = 1:An_Num    
        eval(sprintf('f_names = fieldnames(L_Track_Stats.%s);',...
            genotype{iGene}))
        An_Name = f_names{iAnimal,1};
        for iDir = 1:size(dir)
            %Extracts fieldwidths
            eval(sprintf(['fws = L_Track_Stats.%s.%s.',...
                '%s.field_widths;'],genotype{iGene},An_Name,dir{iDir}))
            
            %Extracts mean firing rates
            eval(sprintf(['m_frs = L_Track_Stats.%s.%s.',...
                '%s.mean_frs;'],genotype{iGene},An_Name,dir{iDir}))
            
            %Extracts spatial info
            eval(sprintf('si = L_Track_Stats.%s.%s.%s.spatial_info;',...
                genotype{iGene},An_Name,dir{iDir}))
            
            %Removes non-PCs if true
            if PC_Only
                eval(sprintf('pv = L_Track_Stats.%s.%s.%s.PC_Ver;',...
                    genotype{iGene},An_Name,dir{iDir}))
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
            
            %Adds spatial info to new array
            eval(sprintf(['%s_spat_info_array = ',...
                '[%s_spat_info_array si''];'],...
                genotype{iGene},genotype{iGene}))
            
            for iCell = 1:size(fws,1)
                if ~isempty(fws{iCell,1})
                    %Adds field widths to new array
                    eval(sprintf(['%s_field_width_array = ',...
                        '[%s_field_width_array fws{iCell,1}(:,1)''];'],...
                        genotype{iGene},genotype{iGene}))
                    
                    %Finds and adds number of fields to new array
                    num_fields = size(fws{iCell,1},1);
                    eval(sprintf(['%s_field_num_array(num_fields) ',...
                        '= %s_field_num_array(num_fields) + 1;'],...
                        genotype{iGene},genotype{iGene}))
                    
                    %Adds mean in-field FR to new array
                    if ~isnan(fws{iCell,1})
                        eval(sprintf(['%s_mean_field_FR = ',...
                            '[%s_mean_field_FR fws{iCell,1}(:,3)''];'],...
                            genotype{iGene},genotype{iGene}))
                    end
                                        
                    %Find peaks firing rate of each PF and adds to new
                    %array
                    if ~isnan(m_frs(iCell,:))
                        cell_pfr = m_frs(iCell,fws{iCell,1}(:,2));
                        eval(sprintf(['%s_pfield_fr = [%s_pfield_fr ',...
                            'cell_pfr];'],genotype{iGene},genotype{iGene}))
                        clear cell_pfr 
                    end
                end
            end
            clear fws m_frs si 
            if PC_Only
                clear pv
            end
        end
    end
    if rm_out
        %Remove Outliers
        %Inter-quartile range of FRs
        eval(sprintf('IQR_FR = iqr(nonzeros(%s_pfield_fr));',genotype{iGene}))
        %75th percentile
        eval(sprintf('u_perc_FR = prctile(nonzeros(%s_pfield_fr),75);',...
            genotype{iGene}))
        Max_FR = u_perc_FR + (1.5*IQR_FR);

        eval(sprintf('in_ind = find(%s_pfield_fr < Max_FR);',...
            genotype{iGene}))
        eval(sprintf('%s_pfield_fr = %s_pfield_fr(1,in_ind);',...
            genotype{iGene},genotype{iGene}))
        eval(sprintf('%s_field_width_array = %s_field_width_array(1,in_ind);',...
            genotype{iGene},genotype{iGene}))
    end
end


%% Plotter
%% Field Width
% Mean Field Width
[mean_array,SEM_array] = Plot_SampleMeans(WT_field_width_array,...
    Het_field_width_array,1,'WT','Het','Place Field Width',1,0);

%Stats
[hPF_Width, pPF_Width, nPF_Width] = MeanStats(WT_field_width_array, ...
    Het_field_width_array);

%% Peak Within-Field Firing Rate
[median_array] = Plot_SampleMedians(WT_pfield_fr,...
    Het_pfield_fr,1,'WT','Het','Peak Within Field Firing Rate');

%Stats
[hPF_Fr, pPF_FR, nPF_FR] = CompareMeans_Stats(WT_pfield_fr,Het_pfield_fr,1);

%% Mean Within-Field Firing Rate
[mean_array,SEM_array] = Plot_SampleMeans(WT_mean_field_FR,...
    Het_mean_field_FR,1,'WT','Het','Mean In-Field Firing Rate',1,0);

%Stats
[hPF_mFR, pPF_mFR, nPF_mFR] = MeanStats(WT_mean_field_FR, ...
    Het_mean_field_FR);

%% Mean Spatial Information
[mean_array,SEM_array] = Plot_SampleMeans(WT_spat_info_array,...
    Het_spat_info_array,1,'WT','Het','Spatial Information',1,0);

%Stats
[hSI, pSI, nSI] = MeanStats(WT_spat_info_array, ...
    Het_spat_info_array);

%% Spatial Information Histogram
WT_spat_info_array(WT_spat_info_array == 0) = [];
Het_spat_info_array(Het_spat_info_array == 0) = [];

sb_SI = min(max(WT_spat_info_array,max(Het_spat_info_array)));
edges = linspace(0,sb_SI,20);

[W_NCells] = histcounts(WT_spat_info_array,edges);
[H_NCells] = histcounts(Het_spat_info_array,edges);

WT_PCells = W_NCells / size(WT_spat_info_array,2);
Het_PCells = H_NCells / size(Het_spat_info_array,2);

%Smooth
WT_PCells = smoothdata(WT_PCells,'gaussian',8);
Het_PCells = smoothdata(Het_PCells,'gaussian',8);

edges = linspace(0,sb_SI,size(WT_PCells,2)+1);

x = edges(1:(size(edges,2)-1)) + (mean(diff(edges))/2);
figure; hold on
p = plot(x,WT_PCells, 'k');
g = plot(x,Het_PCells, 'c');
p.LineWidth = 2;
g.LineWidth = 2;

[ax] = plot_prop();

%% Number of place fields
for iGene = 1:size(genotype,1)
    %Normalize by number of cells
    eval(sprintf(['%s_field_num_array = %s_field_num_array / ',...
        'sum(%s_field_num_array);'],genotype{iGene},genotype{iGene},...
        genotype{iGene}))
end

%Highest number of place fields
max_numpf = max(max(find(WT_field_num_array)),...
    max(find(Het_field_num_array)));

for iField = 1:max_numpf
    pf_per_cell(iField,1) = WT_field_num_array(iField);
    pf_per_cell(iField,2) = Het_field_num_array(iField);
end

figure;
%name = {'1';'2'};
edges = 1:max_numpf;
b = bar(edges,pf_per_cell);
b(1,1).LineWidth = 1.5;
b(1,2).LineWidth = 1.5;
%set(gca,'xticklabel',name)
hold on
title('Place Fields Per Place Cell')
%Plot Properties
ax = gca;
ax.FontName = 'Arial';
ax.FontWeight = 'bold';
ax.Box = 'off';
ax.LineWidth = 1.5;