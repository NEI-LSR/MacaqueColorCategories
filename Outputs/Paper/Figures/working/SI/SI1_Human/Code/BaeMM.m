clear, clc, close all

%%

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath(repoHomeDir));

%%
data = processBaeData();
% 
% % this shifts all values by 1, such that cue #1 is at 0 degrees
% for i = 1:size(data.trialdata.cues,1) 
%     if data.trialdata.cues{i} == 180
%         data.trialdata.cues(i) = {1};
%         disp('woop1')
%     else
%         data.trialdata.cues(i)     = {data.trialdata.cues{i}   + 1};
%     end
%     if data.trialdata.chosen{i} == 180
%         data.trialdata.chosen(i) = {1};
%         disp('woop2')
%     else
%         data.trialdata.chosen(i)   = {data.trialdata.chosen{i} + 1};
%     end
% end

% - % max(cell2mat(data.trialdata.chosen))  = 182!!!???

%%

Lab = 1;
includeCorrect = true;
lengthOfSlidingWindow = 9; % we set our sliding window to 3, and their stimulus interval is roughly 3 times as small as ours, so this means that they get roughly the same smoothing treatment

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

plotMixtureModel(model,...
    whichFigures,'BaeMM')

