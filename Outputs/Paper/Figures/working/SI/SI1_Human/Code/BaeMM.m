clear, clc, close all

%%

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath(repoHomeDir));

%%

loadedData = readtable([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Bae_2015',filesep,'raw_data',filesep,'DelayedEstimation data.csv'],...
    'Range',"D1:F10801");

% loadedData = readtable([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Bae_2015',filesep,'raw_data',filesep,'UndelayedEstimation data.csv'],...
%     'Range',"D1:F10801");

cue_hueAngle = loadedData.TargetColor*2;

error = loadedData.ClickAngle - loadedData.TargetAngle;
choice_hueAngle = cue_hueAngle + error;
choice_hueAngle(choice_hueAngle<0)   = choice_hueAngle(choice_hueAngle<0)+360;
choice_hueAngle(choice_hueAngle>360) = choice_hueAngle(choice_hueAngle>360)-360;

cue_index = (cue_hueAngle/2) + 1; % to account for the fact that the first angle is 2deg whereas our analyses assume that the first angle is 0deg
cue_index(cue_index == 181) = 1;

choice_index = (choice_hueAngle/2) + 1; % to account for the fact that the first angle is 2deg whereas our analyses assume that the first angle is 0deg
choice_index(choice_index == 181) = 1;


for i = 1:size(loadedData,1)
    data.trialdata.cues(i,1)      = {cue_index(i)};
    data.trialdata.choices(i,:)   = {1:180};
    data.trialdata.chosen(i,1)    = {choice_index(i)};
end

data.nBig = 180;

%%

Lab = 1;
includeCorrect = true;
lengthOfSlidingWindow = 15; % picked by hand

model = fitMixtureModel(data,lengthOfSlidingWindow,includeCorrect);

%% For direct comparison to Bae+ 2015 F8b

% figure, 
% plot(0:2:358,-model.bias)
% xticks(0:60:360)
% yticks([-17,-11,-6,0,6,11,17])

%%

model.stimColorSpace    = 'CIELAB';
model.stimCols          = [70,38];      %L*, chroma

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

axlims = 30;

withLabels = false;
DKL = 'Bae';
plotMixtureModel(model,...
    whichFigures,'Bae_CIELAB_',withLabels,DKL,axlims)
