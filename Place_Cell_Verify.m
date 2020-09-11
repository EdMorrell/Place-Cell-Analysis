function [PC_Ver] = Place_Cell_Verify(mean_frs, measure)

%  ---  Place_Cell_Verify
%         - Function to determine whether a cell is a place cell based on
%           chosen measure
%
%  Inputs:
%         - mean_frs: Mean Firing Rate array
%         - measure: 1 - Firing Rate Measure: Mean firing rate of at least
%                                             0.2 Hz on the track and a
%                                             peak firing rate of at least
%                                             5HZ
%         - measure: 2 - Spatial Information Measure
%  Outputs:
%          - PC_Ver: Array with 1 for PC and 0 for not PC

if measure == 1
    %Place cell if mean on track firing rate above 0.1Hz and peak firing
    %rate above 5Hz
    PC_Ver = zeros(size(mean_frs,1),1);
    for iCell = 1:size(mean_frs,1)
        if mean(mean_frs(iCell,:)) > 0.2 && max(mean_frs(iCell,:)) > 5
            PC_Ver(iCell,1) = 1;
        else
            PC_Ver(iCell,1) = 0;
        end
    end
end