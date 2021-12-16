% Data normalization by columns
% input: Raw Data (DataRaw)
% return: Normalized Data (DataNorm)
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [DataNorm] = fNormalization(DataRaw)
    DataNorm=[];
    for i = 1:size(DataRaw,2)
        DataNorm=[DataNorm (DataRaw(:,i)-min(DataRaw(:,i)))/(max(DataRaw(:,i))-min(DataRaw(:,i)))];
    end
end