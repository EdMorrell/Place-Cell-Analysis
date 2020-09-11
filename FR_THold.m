function [spike_times, Cell_ID] = FR_THold(spike_times, Cell_ID, fr_th,...
    dir)

% --- FR_THold
%        - Function to remove units with a mean firing rate above/below a
%          pre-assigned firing rate threshold
%    Inputs:
%           - spike_times: cell array of spike times
%           - Cell_ID: cell array of cell IDs
%           - fr_th: Firing rate threshold (Hz)
%           - dir: 0 - Remove cells below this threshold (default = 0)
%                  1 - Remove cell above this threshold
%    Outputs:
%           - spike_times & Cell_ID with low FR cells removed

if nargin < 4
    dir = 0;
end

%Calculates number of cells
num_cells = size(spike_times,1);

%Calculates and creates an array of mean firing rate for every cell
for iCell = 1:num_cells
    if isempty(spike_times{iCell,1})
        fr = NaN;
    else
        fr = 1/(mean(diff(spike_times{iCell,1})));
    end
    fr_array(1,iCell) = fr;
end

%Removes cells below a pre-determined threshold value
new_spike_times = {};
cell_count = 1;
for iCell = 1:size(fr_array,2)
    if dir == 0 & fr_array(1,iCell) < fr_th || isnan(fr_array(1,iCell))
        continue
    elseif dir == 1 & fr_array(1,iCell) > fr_th || isnan(fr_array(1,iCell))
        continue
    else
        new_spike_times{cell_count,1} = spike_times{iCell,1};
        New_Cell_ID(cell_count) = Cell_ID(iCell);
        cell_count = cell_count + 1;
    end
end

clear spike_times

spike_times = new_spike_times;
Cell_ID = New_Cell_ID;
    
end