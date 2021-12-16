% Function that returns the list of .CSV files
% input (path): address of the folder where the data is .CSV
% output (filenames): Complete list of .CSV file in folder
% Example: 
%   path = fullfile('./Data/');
%   filenames=FindMAT(path)
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [filenames] = FindCSV(path)
    filenames = dir(fullfile(path ,'*.csv'));
    [~, reindex] = sort( str2double( regexp( {filenames.name}, '\d+', 'match', 'once' )));
    filenames = filenames(reindex);
end
