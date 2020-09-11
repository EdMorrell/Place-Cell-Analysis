close all;
clear;

% This script takes a cell array of spike times (1 cell per cell) and gets
% spike positions

fn_LFP = 'D:\Ed\Data\Matlab Outputs\Filenames\Familiar Track\fns.mat';
fn_LFP_tail = 'D:\Ed\Data\Matlab Outputs\Filenames\Familiar Track\csc_Run.mat';

addpath('D:\Ed\Scripts\Tools')
addpath('D:\Ed\Scripts\Single Units')

%Spike file name (Run_spike_times for familiar, spike_times_(posn) for
% novel)
s_fName = 'Run_spike_times';

%Position (if more than one - otherwise leave as '')
posn = '';

% Parameters

%Thresholding firing rate
fr_th = 0.2;
rem_cells = true; %If true removes cells above a threshold

%Save spike position array
save_pos_array = false; %If true then the position array is saved in the current path
spike_pos_fname = 'spike_pos_Pos2'; %Name to save spike position file as

%% Load filename structure
load(fn_LFP); load(fn_LFP_tail);

genotype = {'WT';'Het'};

WT = {'Durga';'Elijah';'Ganesh';'Ben';'Felix'};
Het = {'Abraham';'Buddha';'Fabulinus';'Avalon';'Craig';'Daniel';'Eddie'};
% WT = {'Felix'};
% Het = {};

%% Runs through each animal's spike times file to get spike positions
for iGene = 1:size(genotype,1)
    
    eval(sprintf('%s_cell_count = 1;',genotype{iGene}))
    cell_count = 1;
    eval(sprintf('an_num = size(%s,1);',genotype{iGene}));
    for iAnimal = 1:an_num
        
        %Animal Name
        eval(sprintf('An_Name = %s{iAnimal};',genotype{iGene})) 
        
        %% Load spikes and tracking
        eval(sprintf(['fn_spikes = [fns.%s filesep ',...
            's_fName ''.mat''];'],An_Name))
        eval(sprintf(['fn_tracking = [fns.%s filesep ',...
            '''Tracking'' posn ''.mat''];'],An_Name))
        eval(sprintf(['fn_events = [fns.%s filesep ',...
            '''Event_Timestamps.mat''];'],An_Name))
        
        load(fn_spikes); load(fn_tracking); load(fn_events);

        %% Remove Cells below a pre-determined firing rate?        
        if rem_cells
            [spike_times, Cell_ID] = FR_THold(spike_times, Cell_ID, fr_th);
        end
        %% Remove cells above Firing Rate
%         max_fr = 20; %10Hz tends to be peak FR for most pyr cells
%         [spike_times, Cell_ID] = FR_THold(spike_times,Cell_ID,max_fr,1);
        %% Change spike times to start from 0
        for iCell = 1:size(spike_times,1)
            
            spike_times{iCell,1} = spike_times{iCell,1} - ...
                Event_Timestamps(2,1);
        end
        %% Get spike positions     
        num_cells = size(spike_times,1);
        for iCell = 1:num_cells 
            %Makes a cell array of individual spike positions and times
            [pos_mat] = get_spike_positions(spike_times{iCell,1},[Tracking(:,2)...
                Tracking(:,3) Tracking(:,4)]);   
            spike_pos_all{cell_count,1} = pos_mat;
            spike_pos{iCell} = pos_mat;
            cell_count = cell_count + 1;
        end
        
        if save_pos_array
            eval(sprintf('cd(fns.%s)',An_Name))
            save(spike_pos_fname, 'spike_pos')
        end
        
        eval(sprintf('%s_cell_count = %s_cell_count + size(spike_times,1);',...
            genotype{iGene},genotype{iGene}))
        
        eval(sprintf('%s_spike_pos = spike_pos_all;',genotype{iGene}))
    
        clear spike_pos spike_pos_all Tracking spike_times Cell_ID
    end
end