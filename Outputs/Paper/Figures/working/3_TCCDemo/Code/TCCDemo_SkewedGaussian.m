%% Generate simulated data based on a skewed gaussian model

clear, clc, close all
rng(0)

% !
% Before starting, `cd` to the location of this script

%%

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

load C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\combinedData.mat 
% regenerate this, and write in options for regenerating it within this file 
% !!!!!!!!!!!!!

%% 

moving_bias = fitMixtureModel(cleandata);

% figure, plot(moving_bias)

%% Generate a set of paramters for skewed gaussians that would create this data structure

SimFuc_sd = 60;
[m,c] = computeRelationshipBetweenSkewAndFittedGaussian(SimFuc_sd);

% % sd = 30;
% m = -34.7129397313997; %((64a78a13-c347-4503-999c-1f6a5f5ff0f9))
% c = 17.3564698678105;

skewedGaussians = (moving_bias - c)/m;

[~, SGdata] = GenerativeModel([],'skewedGaussians',skewedGaussians,...
    'SimFuc_sd',SimFuc_sd,'nTrials',size(cleandata.trialdata.cues,1));

SG_moving_bias = fitMixtureModel(SGdata);

%%

% figure, 
% hold on
% plot(moving_bias)
% plot(SG_moving_bias)

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
    sg(:,i) = SplitGauss(x, SkewRange(i)*sd, (1-SkewRange(i))*sd);
end

% figure, 
% plot(x,sg,'k')
% hold on
% xline(0)

% Fit gaussians
gaussEqn = @(a,b,c,d,x) a*exp(-(((x-b).^2)/(2*c^2)))+d';
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