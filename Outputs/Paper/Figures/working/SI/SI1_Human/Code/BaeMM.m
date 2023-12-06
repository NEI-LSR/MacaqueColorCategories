clear, clc, close all

%%

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath(repoHomeDir));

%%
data = processBaeData();

%%

Lab = 1;
includeCorrect = true;
lengthOfSlidingWindow = 9; % we set our sliding window to 3, and their stimulus interval is roughly 3 times as small as ours, so this means that they get roughly the same smoothing treatment

model = fitMixtureModel(data,Lab,lengthOfSlidingWindow,includeCorrect);

%% For direct comparison to Bae+ 2015 F8b

% figure, 
% plot(0:2:358,-model.bias)
% xticks(0:60:360)
% yticks([-17,-11,-6,0,6,11,17])

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

plotMixtureModel(model,...
    whichFigures,'BaeMM')

