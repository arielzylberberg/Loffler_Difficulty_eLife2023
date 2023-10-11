function [S_concat] = concat_struct_fields(S)
%CONCAT_STRCT_FIELDS Summary of this function goes here
%   Detailed explanation goes here

fields = fieldnames(S);
for k = 1:numel(fields)
    aField     = fields{k}; % EDIT: changed to {}
    
    for i = 1:size(S,2)
        if i == 1
            S_concat.(aField) = [];
        end
        if size(S(i).(aField),1) == 1 % if single row --> add columns
            S_concat.(aField) = [S_concat.(aField) S(i).(aField)]; %S.(aField) = cat(dim, S.(aField), T.(aField));
        else % if single column; or matrix --> add rows
            S_concat.(aField) = [S_concat.(aField); S(i).(aField)]; %S.(aField) = cat(dim, S.(aField), T.(aField));
        end
    end
end


end

