function [spike_pos_mat] = get_spike_positions(spike_times, tracking)
%%
% - This function creates a Nx3 matrix of x_positions, y-positions and
% timestamps for an individual cell

% - Add findClosestValue to the path

%Tracking is a Nx3 matrix in the form: [x-positions y-positions timestamps]

spike_pos_mat = zeros(length(spike_times),3);
for iSpike = 1:length(spike_times)
    [i,v] = findClosestValue(tracking(:,3),spike_times(iSpike));
    spike_pos_mat(iSpike,1) = tracking(i,1);
    spike_pos_mat(iSpike,2) = tracking(i,2);
    spike_pos_mat(iSpike,3) = spike_times(iSpike);
end