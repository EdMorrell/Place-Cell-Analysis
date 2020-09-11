function [spike_times,Cell_ID] = check_refractory_period(spike_times,Cell_ID,p_viol)

%   ---  check_refractory_period
%           - Function to remove cells with an above threshold proportion
%             of cells with ISI's less than (2ms)
%           
% Inputs:
%        - spike_times: spike-times cell array
%        - cell_ID: Array of cell names
%        - p_viol: Allowed proportion of cells violating refractory period
%                  (default 0.01 ie less than 1%)
% Outputs:
%         - spike_times: Updated spike-time file

if nargin < 2
    p_viol = 0.01;
end

iCheck = 1;
for iCell = 1:size(spike_times,1)
    t_ds = diff(spike_times{iCell,1});
    pV = size(find(t_ds < 0.002),2) / size(t_ds,2);
    
    if pV <= p_viol
        new_spike_times{iCheck,1} = spike_times{iCell,1};
        new_Cell_ID{iCheck,1} = Cell_ID{iCell,1};
        iCheck = iCheck + 1;
    end
end

if exist('new_spike_times')
    clear spike_times
    spike_times = new_spike_times;
end

if exist('new_Cell_ID')
    clear Cell_ID
    Cell_ID = new_Cell_ID;
end
end