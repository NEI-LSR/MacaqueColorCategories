clear, clc, close all

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];

addpath(genpath([repoHomeDir,'Analyses']))

%% Load data

loadedData = readtable([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Panichello_2019',filesep,'41467_2019_11298_MOESM4_ESM.xlsx'],...
    'Sheet','Figure1C','Range',"A2:B17777");

%%

% I'm making the assumption that the numbers are hue angles, rather than
% index numbers (with with interval of 1degree it won't matter much either
% way)

for i = 1:size(loadedData,1)
    data.trialdata.cues(i,1)     = {loadedData.target(i)};
    data.trialdata.choices(i,:)  = {1:360};
    if loadedData.response(i) > 360
        data.trialdata.chosen(i,1)   = {loadedData.response(i) - 360};
    else
        data.trialdata.chosen(i,1)   = {loadedData.response(i)};
    end
end

data.nBig   = 360;
data.nSmall = 360;

%%

Lab = 1;
includeCorrect = true;
lengthOfSlidingWindow = 17; % we set our sliding window to 3, and their stimulus interval is roughly 6 times as small as ours, so this means that they get roughly the same smoothing treatment

model = fitMixtureModel(data,Lab,lengthOfSlidingWindow,includeCorrect);

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

plotMixtureModel(model,...
    whichFigures,'Panichello_')

