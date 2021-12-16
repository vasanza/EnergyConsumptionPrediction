% Reshape Array using rms features by time windows
% input: Input Data (DataIn)
% input: Time window or timestep value (TimeStep)
% return: Reshape Array with rms features (DataRMS)
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/
function [DataRMS] = frms_features_v2(DataIn,TimeStep)
    DataRMS=[];
        for i = 1:TimeStep:length(DataIn)-TimeStep
            DataRMS=[DataRMS; rms(DataIn(i:i+TimeStep,:))];
        end
        
end