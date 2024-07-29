clear, clc, close all


%% Add required paths

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath(repoHomeDir));

%%

loadedData = readtable([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Bae_2015',filesep,'raw_data',filesep,'UndelayedEstimation data.csv'],...
    'Range',"D1:F8641");

cue_hueAngle = loadedData.TargetColor*2;

error = loadedData.ClickAngle - loadedData.TargetAngle;
choice_hueAngle = cue_hueAngle + error;
choice_hueAngle(choice_hueAngle<0)   = choice_hueAngle(choice_hueAngle<0)+360;
choice_hueAngle(choice_hueAngle>360) = choice_hueAngle(choice_hueAngle>360)-360;

cue_index = (cue_hueAngle/2) + 1; % to account for the fact that the first angle is 2deg whereas our analyses assume that the first angle is 0deg
cue_index(cue_index == 181) = 1;

choice_index = round(choice_hueAngle/2) + 1; % to account for the fact that the first angle is 2deg whereas our analyses assume that the first angle is 0deg
choice_index(choice_index == 181) = 1;


%% convert into "cleandata" format

for i = 1:size(loadedData,1)
    data.trialdata.cues(i,1)      = {cue_index(i)};
    data.trialdata.choices(i,:)   = {1:180};
    data.trialdata.chosen(i,1)    = {choice_index(i)};
end

data.nBig = 180;
%%

includeCorrect = true;
lengthOfSlidingWindow = 15; % picked by hand

model = fitMixtureModel(data,lengthOfSlidingWindow,includeCorrect);

%% For direct comparison to Bae+ 2015 F8b

% figure,
% plot(0:2:358,-model.bias)
% xticks(0:60:360)
% yticks([-17,-11,-6,0,6,11,17])

%%

try
    choiceProb = model.choice_probability'; % transposing to match similarity matrix (so cue on x-axis, choice on y-axis)
catch
    % choiceProb_diag = model{1,1}.choice_probability;
    error('model variable structure is nested') % TODO work out why
end

for i = 1:size(choiceProb,1)
    choiceProb_diag(i,:) = circshift(choiceProb(i,:),i-(size(choiceProb,1))/2);
end

choiceProb_diag = choiceProb_diag/max(choiceProb_diag(:));

plotSimilarityMatrix(choiceProb_diag*2) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
