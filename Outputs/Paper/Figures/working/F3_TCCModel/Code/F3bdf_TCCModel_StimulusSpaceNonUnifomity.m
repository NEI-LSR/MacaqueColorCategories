%% Generate simulated data based on stimulus-space non-uniformity

clear, clc, close all

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromModelOutput';

%% Behind the scenes...

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromPreProcessedData')
    error('not yet implemented')
    % run model
    % save output
end

if strcmp(AnalysisDepth,'fromModelOutput')
    load('../../../../../../Analyses/TCCModels/Combined/combined_TCC-fullremap-workspace_2300708.mat') % TODO Investigate "Warning: Could not find appropriate function on path loading function handle" (still seems to work though...)
end % TODO Remove this dependency

% figure, hold on
% plot(x)
x = movmean(x,9);
% plot(x)

% _Could use fancier moving mean from
% https://github.com/NEI-LSR/MacaqueColorCategories/blob/b08441688de7e10b5d3c325aca12b5ac442d84e5/Analyses/fitMixtureModel.m#L182-L183
% but quick testing suggests it doesn't do much, and so for cartoons, let's
% keep the code simple!_

%% Compute similarity matrix

f = @(x)GenerativeModel(x,'choiceInds',choiceInds,'cueInd',cueInd,'response',response,'nTrials',nTrials,'nBig',nBig,'nSmall',nSmall,'dprime',2.443912865562119,'optimisationMeta',optimisationMeta,...
    'pltSimFigs', true, 'gaussianWidth', 25); % this is available in the load file, but I want to add `pltSimFigs` to it

[~,simdata] = f(x); 
figure(5)
saveas(gcf,['../ssnu_SimilarityFunction_',datestr(now,'yymmdd-HHMMSS'),'.svg']);

plotSimilarityMatrix(simdata.trialdata.similarityMatrix,...
    'ssnu','../',[])

%%

lengthOfSlidingWindow = 9; %Extra smoothing to simplify visual interpretation of instructive cartoon figures

model = fitMixtureModel(simdata,lengthOfSlidingWindow);

whichFigures.MixMod_polar = true;
plotMixtureModel(model,...
    whichFigures,['TCCDemo_ssnu_',AnalysisDepth])
plotMixtureModel(model,...
    whichFigures,['TCCDemo_ssnu_',AnalysisDepth],false)
