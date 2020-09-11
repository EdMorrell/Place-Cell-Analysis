function [spike_pos_iRun] = get_SpikePos_IndRuns(spike_pos,timestamps)

% --- get_SpikePos_IndRuns
%       - Function to get the spike positions for each individual run
%    
%    Inputs: 
%           - spike_pos: Cell array os spike positions and timestamps with
%             each cell in form: [x y time(secs)]
%           - timestamps: Array of timestamps corresponding to start and
%             end of each run in form: [start end]
%    Outputs:
%            - spike_pos_iRun: Cell array of spike positions for each run

%% Creates a cell array with the spike positions in each indivdual run
spike_pos_iRun = {};

for iCell = 1:size(spike_pos,1)
    for iRun = 1:size(timestamps,1)
        spike_index = find(spike_pos{iCell,1}(:,3) >= timestamps(iRun,1) & ...
            spike_pos{iCell,1}(:,3) <= timestamps(iRun,2));
        if isempty(spike_index)
            spike_pos_iRun{iCell,1}{iRun,1} = [];
        else
            spikes(:,1) = spike_pos{iCell,1}(spike_index(:,1),1);
            spikes(:,2) = spike_pos{iCell,1}(spike_index(:,1),2);
            spikes(:,3) = spike_pos{iCell,1}(spike_index(:,1),3);
            spike_pos_iRun{iCell,1}{iRun,1} = spikes;
            clear spikes      
        end
    end
end