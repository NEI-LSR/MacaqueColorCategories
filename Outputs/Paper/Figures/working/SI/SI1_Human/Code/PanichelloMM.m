clear, clc, close all

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];

%%

load([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Panichello_2019',filesep,'Processed_data',filesep,...
    'MonkeyE',filesep,'monkeye.mat'])
% figure, histogram(cleandata.trialdata.chosen,1000)
data = cleandata;

load([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Panichello_2019',filesep,'Processed_data',filesep,...
    'MonkeyW',filesep,'monkeyw.mat'])
% figure, histogram(cleandata.trialdata.chosen,1000)

data.trialdata.cues     = [data.trialdata.cues;     cleandata.trialdata.cues];
data.trialdata.chosen   = [data.trialdata.chosen;   cleandata.trialdata.chosen];
data.trialdata.cond     = []; % removing to avoid confusion
data.trialdata.dirname  = []; % removing to avoid confusion
data.trialdata.paradigm = []; % removing to avoid confusion

data.trialdata.allchoices{1,1}  = 1:360;
data.trialdata.choices = repmat({1:360},length(data.trialdata.cues),1);

% figure, histogram(data.trialdata.chosen,1000)

%% Preprocessing

% It's odd that there are only 64 unique cues
% but thousands of unique choices.

data.trialdata.chosen = round(data.trialdata.chosen);
data.trialdata.chosen(data.trialdata.chosen == 0) = 360;


%%

Lab = 1;
excludeCorrect = false;
lengthOfSlidingWindow = 17; % we set our sliding window to 3, and their stimulus interval is roughly 6 times as small as ours, so this means that they get roughly the same smoothing treatment

model = fitMixtureModel(data,Lab,lengthOfSlidingWindow,excludeCorrect);

%%

whichFigures.MixMod_polar    = true;
% whichFigures.MixMod_linear   = true;

plotMixtureModel(model,...
    whichFigures,['bla_'])

