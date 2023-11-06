clear, clc, close all
rng(0)

% !
% Before starting, `cd` to the location of this script

%%

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))

load('../../../../../../Analyses/TCCModels/Castor/211108_090705_Castor/220517--211108_Castor_TCC-FreeSimilarityMatrix-workspace_230225.mat')

% write in options for other levels of reproduction
% !!!!!!!!!!!!!

%%

% categoryCenter = 18;
% plotSimilarityMatrix(x,'18','../',categoryCenter)

plotSimilarityMatrix(x,'','../')

%%

plotbar_NLL_AIC_BIC('Castor','../')
