clc, clear, close all

addpath(genpath('../../../../../../Analyses/'))

%% Load data

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',];
addpath(genpath([repoHomeDir,filesep,'Analyses']))
DataDir = [repoHomeDir,filesep,'Analyses',filesep,'TCCModels',filesep,'combined',filesep];

filename = 'combined_TCC-FreeSimilarityMatrix-workspace_230214';

load([DataDir,filename,'.mat'],'x')

%%

plotSimilarityMatrix(x, filename, '../')

plotSimilarityMatrix(x, [filename,'warm'], '../', 3)
plotSimilarityMatrix(x, [filename,'cool'], '../', 38)

%% 

% input(['Warning: you are about to re-fit all the models.',newline,...
%     'This will take roughly 45 minutes.',newline,...
%     'Press enter to continue.'])

% realOrSimData = 'real';
% 
% % rn = 0;
% for rn = 1:9
%     RecoveryTesting(realOrSimData,rn)
% end

%%

plotbar_NLL_AIC_BIC('Combined','../')

