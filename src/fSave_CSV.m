% Function that save .MAT file as .CSV
% input (filename): name of the file .MAT
% output (saved file): NewFile.csv as default name
% Example:
%   [] = fSave_CSV(allDataMean);
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [] = fSave_CSV(filename)
    allDataMean=array2table(filename);
    writetable(filename,'NewFile.csv');
end