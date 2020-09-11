function [overlap_mat] = Place_Field_Overlap(mean_frs,norm,norm2,plot_ol)

% --- Place_Field_Overlap
%       - Function to compute overlap between place-field pairs

if nargin < 4
    plot_ol = 0;
end
if nargin < 3
    norm2 = 1; %Normalizes by area of reference curve to get overlap as proportion
end
if nargin < 2
    norm = 1;
end

%Create empty comparison matrix
overlap_mat = zeros(size(mean_frs,1),size(mean_frs,1));

for iCell = 1:size(overlap_mat,1)
    for cCell = 1:size(overlap_mat,2)
        iFR = mean_frs(iCell,:);
        cFR = mean_frs(cCell,:);
        
        if norm == 1
            iFR = normalize(iFR,'range');
            cFR = normalize(cFR,'range');
        end
        
        overlap = cumtrapz(1:30,min([iFR' cFR']'));
        overlap2 = overlap(end);
        
        if plot_ol == 1
            figure; hold on
            plot(1:30,iFR)
            plot(1:30,cFR)
            area(1:30,min([iFR' cFR']'));
        end
        
        overlap_mat(iCell,cCell) = overlap2;
        
        clear overlap overlap2
        
    end
end

if norm2
    overlap_mat = overlap_mat ./ diag(overlap_mat);
    %As overlap percentage depends on size of reference field this finds an
    %average over both place-field in the pair
    upper = triu(overlap_mat);
    lower = tril(overlap_mat);
    upper = upper.';

    clear overlap_mat
    overlap_mat = (upper + lower) / 2;
end