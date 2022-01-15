% Reshape Array using max features by time windows
% input: Input Data (DataIn)
% input: Time window or timestep value (TimeStep)
% return: Reshape Array with rms features (DataRMS)
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [DataRMS] = fmax_features(DataIn,TimeStep)
    DataRMS=[];
        
        for i = 1:TimeStep:length(DataIn)-TimeStep
            DataRMS=[DataRMS; max(DataIn(i:i+TimeStep,:)-DataIn(i,:))];
        end
        
end