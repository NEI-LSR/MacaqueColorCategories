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

%%

if strcmp(AnalysisDepth,'fromRawData')
    error('option not implemented yet') %!!!!!!!!!!
end

if strcmp(AnalysisDepth,'fromPreProcessedData_preCombined')
    dirname = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\'; % Fix this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cleandata = combineData(dirname);
    save('../../../../../../Analyses/combinedData.mat','cleandata');
    warning('If combinedData.mat previously existed, it has been overwritten')
end

if strcmp(AnalysisDepth,'fromPreProcessedData_postCombined')
    load('combinedData.mat','cleandata')
end

%% 

lengthOfSlidingWindow = 9; %Extra smoothing to simplify visual interpretation of instructive cartoon figures
moving_bias = fitMixtureModel(cleandata,[],[],lengthOfSlidingWindow);
disp('Figures saved')

% figure, plot(moving_bias)

%% Generate a set of paramters for skewed gaussians that would create this data structure

SimFunc_sd = 50;
[m,c] = computeRelationshipBetweenSkewAndFittedGaussian(SimFunc_sd);

skewedGaussians = (moving_bias - c)/m;

[~, SGdata] = GenerativeModel([],'skewedGaussians',skewedGaussians,...
    'SimFunc_sd',SimFunc_sd,'nTrials',size(cleandata.trialdata.cues,1));

SG_moving_bias = fitMixtureModel(SGdata,[],[],lengthOfSlidingWindow);
disp('Figures saved')

%%

% figure, 
% hold on
% plot(moving_bias)
% plot(SG_moving_bias)

plotSimilarityMatrix(SGdata.trialdata.similarityMatrix,...
    'sg','../')

nBig = 64;
for i = 1:nBig
    sm_cs(i,:) = circshift(SGdata.trialdata.similarityMatrix(i,:), (nBig/2)-i+1);
end

% figure, imagesc(sm_cs)
% axis square

[~,orderByBias] = sort(SG_moving_bias);

figure, hold on
pltCols = parula(nBig);
pltCols(:,4) = 0.8;
for i = 1:nBig
    plot(-180:360/nBig:180,sm_cs(orderByBias(i),[1:end,1]),'Color',pltCols(i,:))
end
plot(-180:360/nBig:180,mean(sm_cs(:,[1:end,1])),'k--','LineWidth',3)
xline(0,'k:')
xlim([-180,180]);
xticks([-180,0,180]);
ylim([0,1]);
yticks([0,1]);
xlabel('Degrees')
ylabel('Similarity')

[~,minloc_SG_moving_bias] = min(abs(SG_moving_bias)); % find the element of SG_moving_bias that is closest to 0
SG_moving_bias_norm = SG_moving_bias - min(SG_moving_bias);
SG_moving_bias_norm = SG_moving_bias_norm/max(SG_moving_bias_norm);

colorbar('Ticks',[0,SG_moving_bias_norm(minloc_SG_moving_bias),1],...
         'TickLabels',[min(SG_moving_bias),0,max(SG_moving_bias)])
%clim([-15,15]) % This would be nice to do but it would take some fiddling so I'll leave it for now...

saveas(gcf,fullfile('../',['sg-kernel_', datestr(now,'yymmdd-HHMMSS'), '.svg']))
disp('Figures saved')

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