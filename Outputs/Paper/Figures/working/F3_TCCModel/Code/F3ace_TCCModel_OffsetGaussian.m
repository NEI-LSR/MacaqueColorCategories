%% Generate simulated data based on a cognitive bias model
% (using offset gaussians)

% This code generates a simulated data set that when analyzed with a 
% mixture model captures the average pattern of results made by macaque 
% monkeys in a match-to-sample color-matching task. 
%
% The code fits a mixture model (Zhang and Luck, 2008, Bae et al, 2015)
% to behavioral data, and then simulates data to match the biases found in 
% that data.
%
% It then fits both a mixture model and a modified Target 
% Confusability Competition (TCC) model (Schurgin et al, 2020) to that
% simulated data and plots visualisations of both.
%
% The simulated data set captures the results of an agent with a true 
% cognitive bias when assessed using a paradigm that uniformly samples 
% colors from a truly uniform perceptual color space.

clear, clc, close all

rng(0)

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromPreProcessedData';

%% Add path to required scripts

addpath(genpath('../../../../../../Analyses/'))

DataDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Data'];

%% Load data

if strcmp(AnalysisDepth,'fromRawData')
    addpath(genpath('../../../../../../Data/'))
    cleandata = combineData(DataDir);
end

if strcmp(AnalysisDepth,'fromPreProcessedData')
    load([DataDir,filesep,'combinedData.mat'],'cleandata');
end

%% Fit models

if strcmp(AnalysisDepth,'fromRawData') || strcmp(AnalysisDepth,'fromPreProcessedData')

    % % Fit a mixture model to the real data (so that we can replicate the bias)

    lengthOfSlidingWindow = 9; % Extra smoothing to simplify visual interpretation of instructive cartoon figures
    model = fitMixtureModel(cleandata,lengthOfSlidingWindow);
    moving_bias = model.moving_bias;

    whichFigures.MixMod_polar = true;
    plotMixtureModel(model,...
        whichFigures,['F3_TCCModel_og_Input_',AnalysisDepth])
    plotMixtureModel(model,...
        whichFigures,['F3_TCCModel_og_Input_',AnalysisDepth],false)

    % % Generate simulated data that matches the bias in the real data

    gaussianWidth = 25;

    [~, OGdata] = GenerativeModel([],'offsetGaussians',moving_bias,...
        'gaussianWidth',gaussianWidth,'nTrials',size(cleandata.trialdata.cues,1));
    OGdata.trialdata.chosen = OGdata.trialdata.chosen';

    OG_model = fitMixtureModel(OGdata,lengthOfSlidingWindow);

    save('OGdata.mat','OGdata','OG_model')

else
    load('OGdata.mat')
end

%% Plot figures

whichFigures.MixMod_polar = true;
plotMixtureModel(OG_model,...
    whichFigures,['F3_TCCModel_og_Output_',AnalysisDepth])
plotMixtureModel(OG_model,...
    whichFigures,['F3_TCCModel_og_Output_',AnalysisDepth],false)

plotSimilarityMatrix(OGdata.trialdata.similarityMatrix,...
    'og','../')
plotSimilarityMatrix(OGdata.trialdata.similarityMatrix,...
    'og','../',[],false)

