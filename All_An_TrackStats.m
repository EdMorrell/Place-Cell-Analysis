close all
clear

% --- Creates a stucture detailing key stats for each cell in each animal

fn = 'D:\Ed\Data\Matlab Outputs\Filenames\Familiar Track';

save_path = 'D:\Ed\Data\Matlab Outputs\Single Units\Place Cell Stats\Familiar Track';

%Path
addpath('D:\Ed\Scripts\Tools')
addpath('D:\Ed\Scripts\Load Data')
addpath('D:\Ed\Toolboxes\MatlabImportExport_v6.0.0')
addpath(genpath('D:\Ed\Toolboxes\jolato'))
addpath(genpath('D:\Ed\Scripts\Single Units'))
rmpath(genpath('D:\Ed\Toolboxes\chronux_2_12'))

fn_files = [fn filesep 'fns.mat'];
load(fn_files)

fn_csc = [fn filesep 'csc_Run.mat'];
load(fn_csc)

%Spike position name (_Pos2 on end if multiple posns) 
s_fName = 'spike_pos';

%Position (if more than one - otherwise leave as '')
posn = '';

genotype = {'WT';'Het'};

WT = {'Durga';'Elijah';'Ganesh';'Ben';'Felix'};
Het = {'Abraham';'Buddha';'Fabulinus';'Avalon';'Craig';'Daniel';'Eddie'};
% WT = {};
% Het = {'Eddie'};

Dir = {'LR';'RL'};

%Save
save_file = false;

%Constrain Firing Rates?
constr = 0; %If 1, don't include cells with excessively high peak firing rates
fr_lim = 50;

%Plotting
plot_raster = true; %Plots a raster for every cell
%%
num_breaks = 0;
num_dos = 0;
num_runs = 0;
for iGene = 1:size(genotype,1)
    
    eval(sprintf('an_num = size(%s,1);',genotype{iGene}));
    for iAnimal = 1:an_num
        
        %Animal Name
        eval(sprintf('An_Name = %s{iAnimal};',genotype{iGene})) 
        
        %LFP for first timestamps
        eval(sprintf(['fLFP = [fns.%s filesep ''CSC'' ',...
                'csc.%s ''.ncs''];'],An_Name,An_Name))
        [LFP_Samples, LFP_Timestamps, Fs] = load_csc(fLFP, 'all');
        floor = LFP_Timestamps(1,1);
        
        %Load Tracking
        eval(sprintf('f_track = [fns.%s filesep ''Tracking'' posn ''.mat''];',...
            An_Name))
        load(f_track)
        
        %Load Spike Positions
        eval(sprintf('f_spos = [fns.%s filesep s_fName ''.mat''];',...
            An_Name))
        load(f_spos)
        spike_pos = spike_pos';
        
        for iDir = 1:size(Dir,1)
            
            %Load Run Timestamps
            eval(sprintf('f_run_ts = [fns.%s filesep ''%s_TS'' posn ''.mat''];',...
                An_Name,Dir{iDir}))
            load(f_run_ts)
            
            %Converts timestamps to seconds
            timestamps = timestamps / 1e+06;
            timestamps = timestamps - floor;
  
            % Constrain RunSpeed
            [timestamps] = Constrain_Run_Speed(timestamps, querycoords, 5);
                                   
            %Add n seconds before and after each timestamp to include 
            %start segment of track
            ts_pad = 4;
            timestamps(:,1) = timestamps(:,1) - ts_pad;   
            timestamps(:,2) = timestamps(:,2) + ts_pad; 
            
            %% Check for tracking dropouts and remove runs with too many
            dl_allowed = 0.5; %Allowed max data loss of 0.25s;
            pre_ts = size(timestamps,1);
            num_runs = num_runs + pre_ts;
            [timestamps] = Check_Run_Dropouts(Tracking,timestamps,...
                dl_allowed); 
            post_ts = size(timestamps,1);
            num_dos = num_dos + (pre_ts - post_ts);
            if isempty(timestamps)
                num_breaks = num_breaks + 1;
                eval(sprintf('L_Track_Stats.%s.%s.%s.mean_frs = NaN(1,30);',...
                    genotype{iGene},An_Name,Dir{iDir}))
                eval(sprintf('L_Track_Stats.%s.%s.%s.mean_occ = NaN(1,29);',...
                    genotype{iGene},An_Name,Dir{iDir}))
                eval(sprintf('L_Track_Stats.%s.%s.%s.spatial_info = NaN;',...
                    genotype{iGene},An_Name,Dir{iDir}))
                eval(sprintf('L_Track_Stats.%s.%s.%s.field_widths = {NaN};',...
                    genotype{iGene},An_Name,Dir{iDir}))
                eval(sprintf('L_Track_Stats.%s.%s.%s.PC_Ver = NaN;',...
                    genotype{iGene},An_Name,Dir{iDir}))
                continue
            end
            %% Get Spike positions for each run 
            [spike_pos_iRun] = get_SpikePos_IndRuns(spike_pos,timestamps); 
            
            eval(sprintf('spike_posns.%s.%s.%s = spike_pos_iRun;',...
                genotype{iGene},An_Name,Dir{iDir}))                  
            %% Generate a firing rate curve
            [mean_frs,s_mean_frs,mean_occ,spike_pos_iRun] = ...
                Mean_Firing_Rate(spike_pos_iRun,timestamps,Tracking,1);
            
            eval(sprintf('L_Track_Stats.%s.%s.%s.mean_frs = s_mean_frs;',...
                genotype{iGene},An_Name,Dir{iDir}))
            eval(sprintf('L_Track_Stats.%s.%s.%s.mean_occ = mean_occ;',...
                genotype{iGene},An_Name,Dir{iDir}))
            %% Verify whether cells are place cells
            [PC_Ver] = Place_Cell_Verify(mean_frs, 1);
            eval(sprintf('L_Track_Stats.%s.%s.%s.PC_Ver = PC_Ver;',...
                genotype{iGene},An_Name,Dir{iDir}))
            %% Plotting
            if plot_raster
                RasterPlotter(spike_pos_iRun,timestamps,s_mean_frs,PC_Ver)
            end     
            %% Get Spatial Information
            [spatial_info_array] = get_SpatialInfo(mean_frs,mean_occ);
            
            eval(sprintf(['L_Track_Stats.%s.%s.%s.spatial_info = ',...
                'spatial_info_array;'],genotype{iGene},An_Name,Dir{iDir}))
            %% Get Place Field (width and number) info
            [field_widths] = get_FieldWidths(s_mean_frs);
            
            eval(sprintf(['L_Track_Stats.%s.%s.%s.field_widths = ',...
                'field_widths;'],genotype{iGene},An_Name,Dir{iDir}))
            
            clear spatial_info_array field_widths mean_frs s_mean_frs ...
                mean_occ spike_pos_iRun PC_Ver
            
        end
    end
end

if save_file
    cd(save_path)
    save L_Track_Stats L_Track_Stats
    save spike_posns spike_posns
end