% Function that selects the features based on the threshold correlation percentage.
% input (DataFeatures): Data with features, where each row is an example
% input (FeaturesLabels): Names of each column or feature
% input (threshold): Maximum correlation value allowed
% output (NewDataFeatures): New dataset without the features with correlation higher than threshold
% output(NewFeaturesLabels): New Labels without the features with correlation higher than threshold
% output(LabelsRemove): List of removed features
%
% Example: 
%   TimeStep =10080;
%   DataFeatures = frms_features_v2(allDataMean,TimeStep);
%   threshold = 0.750;%<-------Maximum correlation value allowed
%   Features_labels = {'WeekDay','Voltage (V)','Current (A)','Active Power (W)','Frecuency (Hz)','Active Energy (KWh)','Power Factor',...
%    'ESP32 Temperature (°C)','CPU Consumption (%)','CPU Power Consumption (%)','CPU Temperature (%)','GPU Consumption (%)'...
%    ,'GPU Power Consumption (%)','GPU Temperature (%)','RAM Consumption (%)','RAM Power Consumption (%)'};
%   [NewDataFeatures,NewFeaturesLabels,LabelsRemove] = Feature_Selection(DataFeatures,Features_labels,threshold)
% 
% More examples: https://github.com/vasanza/Matlab_Code
% Read more: https://vasanza.blogspot.com/

function [NewDataFeatures,NewFeaturesLabels,LabelsRemove] = Feature_Selection(DataFeatures,FeaturesLabels,threshold)
    corr_matr = [];
    corr_matr = corrcoef(DataFeatures);
    % Load nxn Fisher Z correlation matrix, stored as variable 'corrmat' 
    nROI = length(corr_matr); % # regions of interest (ROIs)
    sort_ind = [1:2:nROI,2:2:nROI]; % make all homotopic ROIs on 1 diagonal
    [~, h_corrmat, h_colorbar] = plot_corrmat([],... % leave timeSeries input empty, since already calculated corrmat
    'corrmat', corr_matr,... % the correlation matrix
    'title', 'Electrical Consumption Parameters',... % plot title
    'labels', FeaturesLabels,... % correlation matrix cell labels
    'sort_ind', sort_ind,... % sort ROIs for plotting
    'label_FontSize', 10,... % long labels get cutoff with larger fonts
    'outline', 1); % outline the cells (appearance improvement)
    upper_corr_matr = triu(corr_matr, 1);
    ind=[1:length(upper_corr_matr)];
    to_drop = ind(any(upper_corr_matr > threshold));
    LabelsRemove=FeaturesLabels(to_drop);%Features to be removed
    DataFeatures(:, to_drop) = [];%Eliminate features with values greater than the threshold
    FeaturesLabels(:, to_drop) = [];%Eliminate features with values greater than the threshold
    NewDataFeatures=DataFeatures;
    NewFeaturesLabels=FeaturesLabels;
end

