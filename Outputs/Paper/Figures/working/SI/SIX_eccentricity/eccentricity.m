% Before starting, `cd` to the location of this script

%% Options

clc, clear, close all

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromRawData';

%% Behind the scenes...

% Add path to required script

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..'];
addpath(genpath(repoHomeDir))
modelOutputDir = [repoHomeDir,filesep,'Analyses'];

rng(0)

%% Load/process data

% ecc = 'upper';
ecc = 'lower';

if strcmp(AnalysisDepth,'fromRawData')
    cleandata = combineData([repoHomeDir,filesep,'Data'],[],[],ecc);
end

if strcmp(AnalysisDepth,'fromRawData') || strcmp(AnalysisDepth,'fromPreProcessedData')

    rng(0) % the modelling might be probabilistic - TODO check this

    model = fitMixtureModel(cleandata);
end

%% Plot data

whichFigures.MixMod_polar    = true;

plotMixtureModel(model,...
    whichFigures,['F2_CombinedMMResults_ecc_',ecc,AnalysisDepth])

% withLabels = false;
% plotMixtureModel(model,...
%     whichFigures,['F2_CombinedMMResults_ecc_',ecc,AnalysisDepth],withLabels)
% 

