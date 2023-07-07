%% Generate simulated data based on stimulus-space non-uniformity

clear, clc, close all
rng(42)

% !
% Before starting, `cd` to the location of this script

%%

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

load('../../../../../../Analyses/combined/combined_TCC-fullremap-workspace_2300706.mat')


%%

[~,simdata] = f(x); % Compute similarity matrix

plotSimilarityMatrix(simdata.trialdata.similarityMatrix,...
    'ssnu','../')


%%

fitMixtureModel(simdata);