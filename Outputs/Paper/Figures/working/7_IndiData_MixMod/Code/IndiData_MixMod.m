clear, clc, close all

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromPreProcessedData';

%% Behind the scenes...

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

%%

if strcmp(AnalysisDepth,'fromRawData')
    % load data
    % save data (?)

    % AnalysisDepth = 'fromPreProcessedData';
end

if strcmp(AnalysisDepth,'fromPreProcessedData')
    % fit model
    % filename(participant) = 

    % AnalysisDepth = 'fromModelOutput';
end

if strcmp(AnalysisDepth,'fromModelOutput')

end

%% Plot data

for participant = 1:4
    plotSimilarityMatrix(x(participant).x,filename{participant},'../')
end
