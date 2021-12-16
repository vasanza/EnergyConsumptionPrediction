% Function that returns the averaged data in a time window.
% input (allData): Data in time series
% input (Win): Time Window
% output (allDataMean): Data averaged in time sales
% Example: 
%   win=60;% 1 minute
%   [allDataMean] = FindCSV(allData,Win)
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [allDataMean] = fData_MeanWin(allData,win)
    allDataMean=[];
    for i=1:win:length(allData)-win
        allDataMean=[allDataMean; mean(allData(i:i+win,:))];
    end
end
