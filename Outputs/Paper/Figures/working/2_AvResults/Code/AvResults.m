%function 

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

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')
    error('not coded yet')
    % load data
    % save data (?)
end

if strcmp(AnalysisDepth,'fromPreProcessedData')
    dirname = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\'; % Fix this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cleandata = combineData(dirname);
    save('../../../../../../Analyses/combinedData.mat','cleandata');
    warning('If combinedData.mat previously existed, it has been overwritten')
end

if strcmp(AnalysisDepth,'fromModelOutput')
    load('combinedData.mat','cleandata')
end

%% Plot data

moving_bias = fitMixtureModel(cleandata);
disp('Figures saved')

save('../../../../../../Analyses/moving_bias.mat','moving_bias')
warning('If moving_bias.mat previously existed, it has been overwritten')
