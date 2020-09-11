function [spike_pos_iRun] = Remove_NaN_Occ_Spikes(spike_pos_iRun,occ_array)

% --- Remove_NaN_Occ_Spikes
%        - Function removes spikes corresponding to NaN occupancy on track
%          (ie periods where animal spent too long)
%
% Inputs: 
%        - spike_pos_iRun: Spike positions for each run
%        - occ_array: Occupancy at each bin
% Outputs:
%         - spike_pos_iRun: With NaN spikes removed

%%
for iCell = 1:size(spike_pos_iRun,1)
    for iRun = 1:size(spike_pos_iRun{1,1},1)
        if isempty(spike_pos_iRun{iCell,1}{iRun,1})
            continue
        else
            [~,nan_x_pos] = find(isnan(occ_array(iRun,:)));

            bounds(1,:) = (nan_x_pos - 1)*25;
            bounds(2,:) = nan_x_pos *25;

            in_bound_array = [];
            for iBound = 1:size(bounds,2)
                x_ind = find(spike_pos_iRun{iCell,1}{iRun,1}(:,1) >= ...
                    bounds(1,iBound) & ...
                    spike_pos_iRun{iCell,1}{iRun,1}(:,1) < ...
                    bounds(2,iBound));
                if isempty(x_ind)
                    continue
                else
                    in_bound_array = [in_bound_array x_ind'];
                end
            end
            spike_pos_iRun{iCell,1}{iRun,1}(x_ind,:) = [];
            clear x_ind in_bound_array bounds
    end
end

end