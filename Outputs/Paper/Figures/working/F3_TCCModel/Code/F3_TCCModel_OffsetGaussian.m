%% Generate simulated data based on a skewed gaussian model

clear, clc, close all

rng(0)

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:                          % Generate figures from the raw data (slowest)
% fromPreProcessedData_preCombined:     % Generate figures from the pre-processed data (before it has been combined across participants)
% fromPreProcessedData_postCombined:    % Generate figures from the pre-processed data (after it has been combined across participants) (fastest)

AnalysisDepth = 'fromPreProcessedData_postCombined';

%%

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

DataDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Data'];

%%

if strcmp(AnalysisDepth,'fromRawData')
    error('option not implemented yet') %!!!!!!!!!!
end

if strcmp(AnalysisDepth,'fromPreProcessedData_preCombined')
    addpath(genpath('../../../../../../Data/'))
    data = combineData(DataDir);
end

if strcmp(AnalysisDepth,'fromPreProcessedData_postCombined')
    data = readtable([DataDir,filesep,'combinedData.csv']);
end

%% 

lengthOfSlidingWindow = 9; %Extra smoothing to simplify visual interpretation of instructive cartoon figures
model = fitMixtureModel(data,[],lengthOfSlidingWindow);
moving_bias = model.moving_bias;

whichFigures.MixMod_polar = true;
plotMixtureModel(model,...
    whichFigures,['F3_TCCModel_og_Input_',AnalysisDepth])

%% Generate a set of paramters for skewed gaussians that would create this data structure

gaussianWidth = 25;

[~, OGdata] = GenerativeModel([],'offsetGaussians',moving_bias,...
    'gaussianWidth',gaussianWidth,'nTrials',size(data,1),...
    'pltSimFigs', true);
OGdata.trialdata.chosen = OGdata.trialdata.chosen';

save('OGdata.mat')

figure(3)
saveas(gcf,['../og_SimilarityFunction_',datestr(now,'yymmdd-HHMMSS'),'.svg']);

OG_model = fitMixtureModel(OGdata,[],lengthOfSlidingWindow);
OG_moving_bias = OG_model.moving_bias;

whichFigures.MixMod_polar = true;
plotMixtureModel(OG_model,...
    whichFigures,['F3_TCCModel_og_Output_',AnalysisDepth])

%%

plotSimilarityMatrix(OGdata.trialdata.similarityMatrix,...
    'og','../')

