function [field_widths] = get_FieldWidths(s_mean_frs,frac_max_fr,track_length)

%  --- get_FieldWidths
%      - Function to compute number of place fields and width of each field
%  
%  Inputs:
%         - s_mean_frs: Smoothed mean firing rate array (can also use
%                       non-smoothed array)
%         - frac_max_fr: Fraction of peak firing rate firing rate must stay
%                        above to be considered part of place field
%                        (default = 1/3)
%         - track_length: Length of track or recording room (in cm) for
%                         conversion of place-field width to cm (Default =
%                         175cm)
%  Outputs:
%          - field_widths: Cell array of place field widths and peak
%                          positions (index position) and mean firing rate
%%
%Params
if nargin < 2
    frac_max_fr = 1/3;
end
if nargin < 3
    track_length = 175;
end
d_peak_th = 1/3; %Threshold size of main field peak to be considered a 
                  %secondary peak (default is 1/3 of peak field)
%%
for iCell = 1:size(s_mean_frs,1)
    
    if sum(s_mean_frs(iCell,:)) == 0
        real_peaks = [];
    else
        %Determine whether single place field or multiple
        [v_peak,i_peak] = findpeaks(s_mean_frs(iCell,:));

    
        %Find field peak values and positions and place into new array
        if size(v_peak,2) > 1
            peak_peak = max(v_peak);
            p_count = 1;
            for iPeak = 1:size(v_peak,2)
                if v_peak(iPeak) >= d_peak_th*peak_peak
                    real_peaks(p_count,1) = v_peak(iPeak);
                    real_peaks(p_count,2) = i_peak(iPeak);
                    p_count = p_count + 1;
                else
                    continue
                end
            end
        elseif size(v_peak,2) == 1
            real_peaks(1,1) = v_peak;
            real_peaks(1,2) = i_peak;
        else
            real_peaks = [];
        end
    end
    
    if isempty(real_peaks)
        field_widths{iCell,1} = [];
        field_widths{iCell,1} = [];
        field_widths{iCell,1} = [];
    else
    
        for iPeak = 1:size(real_peaks,1)

            fr_peak = real_peaks(iPeak,1);

            if fr_peak == 0
                field_width = 0;
                field_widths{iCell,1}(iPeak) = field_width;
            else
                peak_th = fr_peak * frac_max_fr;

                peak_pos = real_peaks(iPeak,2);

                num_lt_pixels = 0;
                ab_th = true;
                check_pos = 1;
                last_pos = fr_peak; %Value of last position to check value is reducing
                while ab_th
                    if (peak_pos - check_pos) < 1
                        ab_th = false;
                    elseif s_mean_frs(iCell, (peak_pos-check_pos)) ...
                            > peak_th && (peak_pos - check_pos) >= 1 ...
                            && s_mean_frs(iCell, (peak_pos-check_pos)) ...
                            <= last_pos
                        num_lt_pixels = num_lt_pixels + 1;
                        last_pos = s_mean_frs(iCell, (peak_pos-check_pos));
                        check_pos = check_pos + 1;
                    else
                        ab_th = false;
                    end
                end

                num_rt_pixels = 1;
                ab_th = true;
                check_pos = 1;
                last_pos = fr_peak;
                while ab_th
                    if (peak_pos + check_pos) > size(s_mean_frs,2)
                        ab_th = false;
                    elseif s_mean_frs(iCell, (peak_pos+check_pos)) ...
                            > peak_th && (peak_pos + check_pos) <= ...
                            size(s_mean_frs,2) && ...
                            s_mean_frs(iCell, (peak_pos+check_pos))...
                            <= last_pos
                        num_rt_pixels = num_rt_pixels + 1;
                        last_pos = s_mean_frs(iCell, (peak_pos+check_pos));
                        check_pos = check_pos + 1;
                    else
                        ab_th = false;
                    end
                end

                field_width = num_lt_pixels + num_rt_pixels;
                field_width = (field_width / size(s_mean_frs,2)) * track_length; %Converts to cm
                Mean_Field_FR = mean(s_mean_frs(iCell,peak_pos-num_lt_pixels:...
                    peak_pos+(num_rt_pixels-1))); %Mean in-field firing rate
                field_widths{iCell,1}(iPeak,1) = field_width;
                field_widths{iCell,1}(iPeak,2) = real_peaks(iPeak,2); %Field peak position
                field_widths{iCell,1}(iPeak,3) = Mean_Field_FR;
                field_widths{iCell,1}(iPeak,4) = real_peaks(iPeak,1); %Field peak value
                
                

            end
        end
        clear real_peaks v_peak i_peak field_width
    end
    if ~isempty(field_widths{iCell,1})
        %Sort on basis of peak firing rate
        field_widths{iCell,1} = sortrows(field_widths{iCell,1},4,'descend');
    end
    
end

end