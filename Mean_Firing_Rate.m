function [mean_frs,s_mean_frs,mean_occ,spike_pos_iRun] = Mean_Firing_Rate(spike_pos_iRun,...
    timestamps,Tracking,rm_NaN,constr_fr,fr_lim)

% --- Mean_Firing_Rate 
%      - Computes a firing rate curve at every point on the track for each
%        cell in the array
%
%   Inputs:
%          - spike_pos_iRun: Spike positions of individual runs, generated
%            using get_SpikePos_IndRuns.m
%          - timestamps: each run timestamps, form: [start end]
%          - Tracking: in form [raw_ts x_pos y_pos ts(in secs)]
%                      (only column 2 and 4 important to this function)
%          - rm_NaN: Remove spikes corresponding to NaN occupancy bins 
%                    1 - remove spikes 0 - no removal (default = 0)
%          - constr_fr: 1 - Don't incl. cells with a peak FR above that
%                           specified by fr_lim (or default = 40Hz)
%                       0 - Incl. all cells
%          - fr_lim: Maximum firing rate allowed for inclusion
%   Outputs
%          - mean_frs: Mean firing rate in form: row: cells, cols: firing
%                      rate at each point on track
%          - s_mean_frs: As mean_frs but smoothed
%          - mean_occ: Mean occupancy at every point along track
%          - spike_pos_iRuns: Inputted cell array with NaN spikes removed

%% Calculates mean firing rates at each point based on occupancy

if nargin < 4
    rm_NaN = 0;
end
if nargin < 5
    constr_fr = 0;
end
if nargin < 6
    fr_lim = 40;
end

edges = (0:25:725);

%Smoothing params
sigma = 1.5;
bw = 1;

%Calculates occupancy in every bin and creates an occ_map
occ_array = [];
for iRun = 1:size(timestamps,1)
    itrack_seg = find(Tracking(:,4) >= timestamps(iRun,1) & Tracking(:,4) ...
        <= timestamps(iRun,2));
    track_seg = Tracking(itrack_seg,2);
    occ_array(iRun,:) = histcounts(track_seg,edges);
end

%Convert time in occ_map to seconds
dt = mean(diff(Tracking(:,4)));
occ_array = occ_array * dt;

%Outlier removal
%Finds bin in which occupancy is greater than 1.5iqrs from 75th percentile
%(used due to non normal distribution of occupancy matrix)
Med_dwell = median(nonzeros(occ_array)); %Median dwell time per bin
IQR_dwell = iqr(nonzeros(occ_array)); %Inter-quartile range of dwell time
u_perc_dwell = prctile(nonzeros(occ_array),75); %Dwell time 75th percentile
Max_dwell = u_perc_dwell + (1.5*IQR_dwell); %Outlier threshold based on a distance of 1.5xIQR from 75th percentile

occ_array(occ_array > Max_dwell) = NaN;
occ_array(occ_array == 0) = NaN;

%Creates mean occupancy map
mean_occ = nanmean(occ_array,1);
% mean_occ = mean(occ_array,1);

%Removes spikes during NaN occupancy periods
if rm_NaN == 1
    [spike_pos_iRun] = Remove_NaN_Occ_Spikes(spike_pos_iRun,occ_array);
end

%Creates a cell array of firing rates for every run
fr_array = {};
for iCell = 1:size(spike_pos_iRun,1)
    for iRun = 1:size(timestamps,1)
        if isempty(spike_pos_iRun{iCell,1}{iRun,1})
            spike_array = zeros(1,size(edges,2)-1);
        else
            spike_array = ...
                histcounts(spike_pos_iRun{iCell,1}{iRun,1}(:,1),edges);
        end
        fr_array{iCell,1}(iRun,:) = spike_array ./ occ_array(iRun,:);
    end
end

%Calculates mean firing rate across all runs
x = (0:25:700);
mean_frs = zeros(size(spike_pos_iRun,1),size(x,2));
for iCell = 1:size(fr_array,1)
    mean_fr = nanmean(fr_array{iCell,1},1);
    i_bt = find(isnan(mean_fr) | isinf(mean_fr));
    mean_fr(1,i_bt) = 0;
    if constr_fr == 1
        if max(mean_fr) > fr_lim
            mean_fr = zeros(1,size(x,2));
        end
    end
    mean_frs(iCell,:) = mean_fr;
end

%Generate a smoothed mean_fr array
for iCell = 1:size(mean_frs,1)
    [s_mean_fr] = smooth_array(mean_frs(iCell,:),sigma,bw);
    s_mean_frs(iCell,:) = s_mean_fr;
end