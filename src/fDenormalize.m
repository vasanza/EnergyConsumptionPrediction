% Data Denormalize by columns
% input: Normalized Data(DataNorm)
% input: Max Normalized Data(MaxDataNorm)
% input: Min Normalized Data(MinDataNorm)
% return: Raw Data (DataRaw)
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [Data] = fDenormalize(DataNorm,MaxDataNorm,MinDataNorm)
    Data=[];
    for i = 1:size(DataNorm,2)
        Data=[Data (DataNorm(:,i)+ MinDataNorm(i))*(MaxDataNorm(i)-MinDataNorm(i))];
    end
end