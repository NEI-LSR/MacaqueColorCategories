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
    load('../../../../../../Analyses/combined/combined_TCC-fullremap-workspace_2300708.mat') % !!!!!!!! Investigate "Warning: Could not find appropriate function on path loading function handle" (still seems to work though...)
end 

% figure, hold on
% plot(x)
x = movmean(x,9);
% plot(x)

% _Could use fancier moving mean from
% https://github.com/NEI-LSR/MacaqueColorCategories/blob/b08441688de7e10b5d3c325aca12b5ac442d84e5/Analyses/fitMixtureModel.m#L182-L183
% but quick testing suggests it doesn't do much, and so for cartoons, let's
% keep the code simple!_

%% Compute similarity matrix

[~,simdata] = f(x); 
plotSimilarityMatrix(simdata.trialdata.similarityMatrix,...
    'ssnu','../')

%%

lengthOfSlidingWindow = 9; %Extra smoothing to simplify visual interpretation of instructive cartoon figures
fitMixtureModel(simdata,[],[],lengthOfSlidingWindow);

%% Kernel plot

% Pulled from `GenerativeModel.m`
x = -180:0.1:180;
default_SimFunc_sd = 60;
SplitGauss = @(x,sd_left,sd_right) [exp(-((x(x<=0).^2)/(2*sd_left^2))), exp(-((x(x>0).^2)/(2*sd_right^2)))];
simFunc = SplitGauss(x, 0.5*default_SimFunc_sd, 0.5*default_SimFunc_sd);

figure, hold on
pltCols(:,4) = 0.8;
plot(x,simFunc,'k--','LineWidth',3)
xline(0,'k:')
xlim([-180,180]);
xticks([-180,0,180]);
ylim([0,1]);
yticks([0,1]);
xlabel('Degrees')
ylabel('Similarity')

saveas(gcf,fullfile('../',['ssnu-kernel_', datestr(now,'yymmdd-HHMMSS'), '.svg']))
disp('Figures saved')

