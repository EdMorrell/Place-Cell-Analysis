function [pfdist_mat] = Place_Field_Distance(field_widths)

% --- Place_Field_Distance
%       - Function to compute distance between place-field peaks


%Create empty comparison matrix
pfdist_mat = zeros(size(field_widths,1),size(field_widths,1));

for iCell = 1:size(pfdist_mat,1)
    for cCell = 1:size(pfdist_mat,2)
        if size(field_widths{iCell,1},2) < 2 | ...
                size(field_widths{cCell,1},2) < 2
            pfdist_mat(iCell,cCell) = NaN;
            
        else
            iPos = field_widths{iCell,1}(1,2);
            cPos = field_widths{cCell,1}(1,2);

            pfdist_mat(iCell,cCell) = abs(iPos-cPos);
        end        
    end
end