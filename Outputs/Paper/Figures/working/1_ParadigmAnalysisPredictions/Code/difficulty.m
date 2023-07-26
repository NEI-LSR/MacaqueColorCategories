clear, clc, close all

%%

% TODO Replace with `loadData.m`

filename = 'combinedData';
% TODO Replace with CSVs

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses'))


load(['C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\',...
    filename,'.mat'])

difficulty_psychometric(cleandata,filename);


%%


