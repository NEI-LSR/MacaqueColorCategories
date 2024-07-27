% Conway Lab MTurk Data

clc, clear, close all

%%

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
dataDir = [repoHomeDir,'Data'];

addpath(genpath([repoHomeDir,'Analysis']));
addpath(genpath(dataDir));

%%

load([dataDir,filesep,'211112--220507_MTurk.mat'])

%%

Lab = 0;
lengthOfSlidingWindow = [];
includeCorrect = false;

model = fitMixtureModel(cleandata,lengthOfSlidingWindow,includeCorrect);

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;
whichFigures.GaussianWidth   = true;

axlims = 30;

withLabels = false;
DKL = 'NIH';
plotMixtureModel(model,...
    whichFigures,'MTurk',withLabels,DKL,axlims)

%% Plot choice probability matrix
% Code copied from SI6_choiceMatrices.m

try
    choiceProb_diag = model.choice_probability'; % transposing to match similarity matrix (so cue on x-axis, choice on y-axis)
catch
    % choiceProb_diag = model{1,1}.choice_probability;
    error('model variable structure is nested') % TODO work out why
end

for i = 1:size(choiceProb_diag,1)
    choiceProb_diag(i,:) = circshift(choiceProb_diag(i,:),i-(size(choiceProb_diag,1))/2);
end

choiceProb_diag = choiceProb_diag/max(choiceProb_diag(:));

filename = 'CP_MTurk';

% plotSimilarityMatrix(model.choice_probability) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
plotSimilarityMatrix(choiceProb_diag,filename,'../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

% TODO Relabel simaility as choice probability

% h = findobj;
% h(n).Label = 'Choice Probability'; % doesn't work


%% How many trials from each participant

[C,ia,ic] = unique(cleandata.trialdata.WorkerID);

for i = 1:length(C)
    t(i) = sum(strcmp(cleandata.trialdata.WorkerID,C(i)));
end

figure, histogram(categorical(t))
% figure, plot(t)

