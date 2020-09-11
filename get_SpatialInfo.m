function [spatial_info_array] = get_SpatialInfo(mean_frs,mean_occ)

% --- get_SpatialInfo
%       - Function to compute spatial information for every cell based on
%         Skaggs (1993)
%     
%   Inputs:
%        - mean_frs: Mean firing rate at each point on the track for each cell
%                    computed using Mean_Firing_Rate.m
%        - mean_occ: Mean occupancy at each point on track computed using 
%                    Mean_Firing_Rate.m
%   Outputs:
%           - spatial_info_array: Array of spatial information for each
%                                 cell

%%
%Calculates spatial info
for iCell = 1:size(mean_frs,1)
    spatial_info = 0;
    ov_mean_fr = 0;
    %This loop computes mean overall firing rate (= ?pi.ri)
    for iBin = 1:size(mean_occ,2)    
        if isnan(mean_occ(iBin))
            pi = 0;
        else
            pi = mean_occ(iBin)/nansum(mean_occ); %Probability of being in that bin
        end    
        ri = mean_frs(iCell,iBin); %Firing rate within that bin   
        ov_mean_fr = ov_mean_fr + (pi * ri);    
    end
    %Loop computes ri and pi and bits per bin (using ov_mean_fr from
    %previous loop)
    for iBin = 1:size(mean_frs,2)
        if isnan(mean_occ(iBin))
            pi = 0;
        else
            pi = mean_occ(iBin)/nansum(mean_occ);
        end
        ri = mean_frs(iCell,iBin);
        ri_dash = ri/ov_mean_fr;
        bin_bit = pi*ri_dash*log2(ri_dash);
        if isnan(bin_bit)
            bin_bit = 0;
        end
        bb_array(iBin) = bin_bit;
        spatial_info = spatial_info + bin_bit;
    end
    spatial_info_array(iCell,1) = spatial_info;
end