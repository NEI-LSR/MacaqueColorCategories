%% Generate simulated data based on a skewed gaussian model

clear, clc, close all

rng(0)

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:                          % Generate figures from the raw data (slowest)
% fromPreProcessedData_preCombined:     % Generate figures from the pre-processed data (before it has been combined across participants)
% fromPreProcessedData_postCombined:    % Generate figures from the pre-processed data (after it has been combined across participants) (fastest)

AnalysisDepth = 'fromPreProcessedData_postCombined';

%%

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

DataDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Data'];

%%

if strcmp(AnalysisDepth,'fromRawData')
    error('option not implemented yet') %!!!!!!!!!!
end

if strcmp(AnalysisDepth,'fromPreProcessedData_preCombined')
    addpath(genpath('../../../../../../Data/'))
    data = combineData(DataDir);
end

if strcmp(AnalysisDepth,'fromPreProcessedData_postCombined')
    data = readtable([DataDir,filesep,'combinedData.csv']);
end

%% 

lengthOfSlidingWindow = 9; %Extra smoothing to simplify visual interpretation of instructive cartoon figures
model = fitMixtureModel(data,[],lengthOfSlidingWindow);
moving_bias = model.moving_bias;

whichFigures.MixMod_polar = true;
plotMixtureModel(model,...
    whichFigures,['TCCDemo_sg_Input_',AnalysisDepth])

%% Generate a set of paramters for skewed gaussians that would create this data structure

gaussianWidth = 25;
[m,c] = computeRelationshipBetweenSkewAndFittedGaussian(gaussianWidth);

skewedGaussians = (moving_bias - c)/m;

[~, SGdata] = GenerativeModel([],'skewedGaussians',skewedGaussians,...
    'gaussianWidth',gaussianWidth,'nTrials',size(data,1),...
    'pltSimFigs', true);
SGdata.trialdata.chosen = SGdata.trialdata.chosen';

figure(3)
saveas(gcf,['../sg_SimilarityFunction_',datestr(now,'yymmdd-HHMMSS'),'.svg']);

SG_model = fitMixtureModel(SGdata,[],lengthOfSlidingWindow);
SG_moving_bias = SG_model.moving_bias;

whichFigures.MixMod_polar = true;
plotMixtureModel(SG_model,...
    whichFigures,['TCCDemo_sg_Output_',AnalysisDepth])

%%

plotSimilarityMatrix(SGdata.trialdata.similarityMatrix,...
    'sg','../')

%%

function [m,c] = computeRelationshipBetweenSkewAndFittedGaussian(sd)

x = -180:0.1:180;

SkewRange = 0.1:0.05:0.9;

SplitGauss = @(x,sd_left,sd_right) ...
    [exp(-((x(x<=0).^2)/(2*sd_left^2))),...
     exp(-((x(x> 0).^2)/(2*sd_right^2)))];

for i = 1:length(SkewRange)
    sg(:,i) = SplitGauss(x,...
            SkewRange(i)     * (2 * sd),...
            (1-SkewRange(i)) * (2 * sd));
end

% figure, 
% plot(x,sg,'k')
% hold on
% xline(0)

% Fit gaussians
gaussEqn = @(a,b,c,d,x) a*exp(-(((x-b).^2)/(2*c^2)))+d;
startingPoints = [0.5 0 78 0];

for i = 1:length(SkewRange)
    f{i} = fit(x',sg(:,i),gaussEqn,'start',startingPoints,...
        'Lower',[0 -180 0 0],'Upper',[Inf 180 Inf 1]);
end

% for i = 1:length(SkewRange)
%     figure,
%     hold on
%     plot(x,sg(:,i))
%     plot(x,gaussEqn(f{1,i}.a,f{1,i}.b,f{1,i}.c,f{1,i}.d,x))
% end

% hacky extract
for i = 1:length(SkewRange)
    MixModBias(i) = f{i}.b;
end

lm = fitlm(SkewRange,MixModBias);
m = lm.Coefficients.Estimate(2);
c = lm.Coefficients.Estimate(1);

% figure,
% hold on
% plot(SkewRange,MixModBias,'k-*')
% y = predict(lm,SkewRange');
% plot(SkewRange,y,'r--o')

end