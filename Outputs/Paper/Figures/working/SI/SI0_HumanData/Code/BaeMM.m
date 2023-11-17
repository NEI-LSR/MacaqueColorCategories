% clear, clc, close all

%%

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath(repoHomeDir));

%%
data = processBaeData();

%%

Lab = 1;
excludeCorrect = false;
lengthOfSlidingWindow = 9; % we set our sliding window to 3, and their stimulus interval is roughly 3 times as small as ours, so this means that they get roughly the same smoothing treatment

model = fitMixtureModel(data,Lab,lengthOfSlidingWindow,excludeCorrect);

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

plotMixtureModel(model,...
    whichFigures,'BaeMM')

