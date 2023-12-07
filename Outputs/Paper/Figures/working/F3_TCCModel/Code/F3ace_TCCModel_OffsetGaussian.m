%% Generate simulated data based on a skewed gaussian model

clear, clc, close all

rng(0)

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data

AnalysisDepth = 'fromPreProcessedData';

%%

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

DataDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Data'];

%%

if strcmp(AnalysisDepth,'fromRawData')
    addpath(genpath('../../../../../../Data/'))
    cleandata = combineData(DataDir);
end

if strcmp(AnalysisDepth,'fromPreProcessedData')
    load([DataDir,filesep,'combinedData.mat']);
end

%% 

lengthOfSlidingWindow = 9; %Extra smoothing to simplify visual interpretation of instructive cartoon figures
model = fitMixtureModel(cleandata,lengthOfSlidingWindow);
moving_bias = model.moving_bias;

whichFigures.MixMod_polar = true;
plotMixtureModel(model,...
    whichFigures,['F3_TCCModel_og_Input_',AnalysisDepth])
plotMixtureModel(model,...
    whichFigures,['F3_TCCModel_og_Input_',AnalysisDepth],false)

%% Generate a set of paramters for skewed gaussians that would create this data structure

gaussianWidth = 25;

[~, OGdata] = GenerativeModel([],'offsetGaussians',moving_bias,...
    'gaussianWidth',gaussianWidth,'nTrials',size(cleandata.trialdata.cues,1));
OGdata.trialdata.chosen = OGdata.trialdata.chosen';

save('OGdata.mat')

OG_model = fitMixtureModel(OGdata,lengthOfSlidingWindow);
OG_moving_bias = OG_model.moving_bias;

whichFigures.MixMod_polar = true;
plotMixtureModel(OG_model,...
    whichFigures,['F3_TCCModel_og_Output_',AnalysisDepth])
plotMixtureModel(OG_model,...
    whichFigures,['F3_TCCModel_og_Output_',AnalysisDepth],false)

%%

plotSimilarityMatrix(OGdata.trialdata.similarityMatrix,...
    'og','../')
plotSimilarityMatrix(OGdata.trialdata.similarityMatrix,...
    'og','../',[],false)

