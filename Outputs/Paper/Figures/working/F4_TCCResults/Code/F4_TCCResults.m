clc, clear, close all

addpath(genpath('../../../../../../Analyses/'))

%% Load data

DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\';
filename = 'combined_TCC-FreeSimilarityMatrix-workspace_230214';

load([DataDir,filename,'.mat'],'x')

%%

plotSimilarityMatrix(x, filename, '../')

plotSimilarityMatrix(x, [filename,'warm'], '../', 3)
plotSimilarityMatrix(x, [filename,'cool'], '../', 38)

