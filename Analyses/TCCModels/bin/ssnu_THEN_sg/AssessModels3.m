clear, clc, close all

currentValues.dPrime           = 1.4899;
currentValues.gaussianWidth    = 39.0070;

load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu\5_230803-161847.mat','x')
currentValues.stimulusRemapping = x;
clear x

load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu_THEN_sg\5_230804-174020.mat','x')
currentValues.skewedGaussians = x;
clear x

%% Load Data

DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
if ~exist('data','var')
    data = combineData(DataDir); % TODO Switch to csv
end

data.trialdata.nBig = 64;
data.trialdata.nSmall = 4;
data.trialdata.nTrials = 98104;

%% Define function

choiceInds = cell2mat(data.trialdata.choices);
cueInd = cell2mat(data.trialdata.cues);
response = cell2mat(data.trialdata.chosen);

f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
    'choiceInds',       choiceInds',...
    'cueInd',           cueInd,...
    'response',         response,...
    'nTrials',          data.trialdata.nTrials,...
    'nBig',             data.trialdata.nBig, ...
    'nSmall',           data.trialdata.nSmall,...
    'dPrime',               currentValues.dPrime,...
    'gaussianWidth',        currentValues.gaussianWidth,...
    'stimulusRemappingPol', currentValues.stimulusRemapping,...
    'skewedGaussians',      currentValues.skewedGaussians); % (add stimCols to speed up slightly)


%% Run function
clc

nll = f([]);

