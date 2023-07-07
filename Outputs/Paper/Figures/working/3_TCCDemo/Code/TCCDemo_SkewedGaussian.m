%% Generate simulated data based on a skewed gaussian model

clear, clc, close all
rng(0)

% !
% Before starting, `cd` to the location of this script

%%

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

load C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\combinedData.mat 
% regenerate this, and write in options for regenerating it within this file 
% !!!!!!!!!!!!!

%% 

moving_bias = fitMixtureModel(cleandata);

% figure, plot(moving_bias)

%% Generate a set of paramters for skewed gaussians that would create this data structure

% sd = 30;
m = -34.7129397313997; %((64a78a13-c347-4503-999c-1f6a5f5ff0f9))
c = 17.3564698678105;

skewedGaussians = (moving_bias - c)/m;

[~, SGdata] = GenerativeModel([],'skewedGaussians',skewedGaussians,...
    'nTrials',size(cleandata.trialdata.cues,1));

SG_moving_bias = fitMixtureModel(SGdata);

%%
figure, 
hold on
plot(moving_bias)
plot(SG_moving_bias)

plotSimilarityMatrix(SGdata.trialdata.similarityMatrix,...
    'sg','../')