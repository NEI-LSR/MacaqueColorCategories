clc, clear, close all

%% Load data

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))
DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
% data = combineData(DataDir);
data = load([DataDir,'210517--211108_Castor_data.mat']);

if any(any(isnan(cell2mat(data.trialdata.choices))))
    warning('Detected trials with NaN choices. Removing.')
    [a,~] = ind2sub([size(data.trialdata.cues,1),size(data.trialdata.choices{1},2)],...
        find(isnan(cell2mat(data.trialdata.choices))));
    data.trialdata.allchoices(unique(a)) = [];
    data.trialdata.dirname(unique(a)) = [];
    data.trialdata.paradigm(unique(a)) = [];
    data.trialdata.choices(unique(a)) = [];
    data.trialdata.chosen_idx(unique(a)) = [];
    data.trialdata.cues(unique(a)) = [];
    data.trialdata.stimCols_raw(unique(a)) = [];
    data.trialdata.chosen(unique(a)) = [];
    data.trialdata.stimCols(unique(a)) = [];
    data.trialdata.nTrials = size(data.trialdata.cues,1);
end

%% Load model

% load('iterative230815-203401.mat') % Combined
load('iterative230830-100240.mat') % Castor

x = [currentValues.dPrime,...
    currentValues.stimulusRemapping,...
    currentValues.gaussianWidth,...
    currentValues.skewedGaussians'];

optimisationMeta = double(zeros([6,2]));

optimisationMeta(:,1) = [0,1,0,1,1,1];

nBig = 64;

% How many values does each paramter have?
optimisationMeta(:,2) = [...
    nBig*nBig;...   % Free Similarity Matrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    nBig;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig,...        % skewedGaussians
    ];

%%

choiceInds = cell2mat(data.trialdata.choices);
cueInd = cell2mat(data.trialdata.cues);
response = cell2mat(data.trialdata.chosen);

f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
    'choiceInds',       choiceInds',...
    'cueInd',           cueInd,...
    'response',         response,...
    'nTrials',          size(data.trialdata.cues,1),...
    'nBig',             size(data.trialdata.stimCols{1},1), ...
    'nSmall',           size(data.trialdata.choices{1},2),...
    'optimisationMeta', optimisationMeta);

[~,data_t] = f(x);

%%

plotSimilarityMatrix(data_t.trialdata.similarityMatrix)

%%
figure, hold on
plot(x(67:end),'k')
axis tight
ylim([0,1])
yline(0.5,'k:')

%%
x_mod = x;
x_mod(67:end) = 0.5;

[~,data_t2] = f(x_mod);
plotSimilarityMatrix(data_t2.trialdata.similarityMatrix)

%%