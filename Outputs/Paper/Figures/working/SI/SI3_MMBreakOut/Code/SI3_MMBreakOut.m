% This is a stripped down version of:
% Outputs\Paper\Figures\working\F2_CombinedMMResults\Code\F2_CombinedMMResults.m

% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromModelOutput';

%% Behind the scenes...

% Add path to required script
addpath(genpath('../../../../../../../Analyses/'))
modelOutputDir = '../../../../../../../Analyses';

rng(0)

%% Load data

if strcmp(AnalysisDepth,'fromModelOutput')
    ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
        'combined_*.mat']);
    if length(ModelFile) > 1
        warning('Multiple model files. Using most recent.')
        [~,idx] = sort([ModelFile.datenum]);
        ModelFile = ModelFile(idx);
        ModelFile = ModelFile(end);
    end
    load([modelOutputDir,'/MixtureModels/',ModelFile.name],'model')
end

%% Plot data

whichFigures.mixMod_BreakOut = true;

plotMixtureModel(model,...
    whichFigures,'../')

