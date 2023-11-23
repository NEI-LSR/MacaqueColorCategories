clear, clc, close all

% TODO Add replication levels

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];

addpath(genpath([repoHomeDir,'Analyses']))

load([repoHomeDir,'Analyses',filesep,'TCCModels',filesep,'Combined',filesep,'combined_TCC-fullremap-workspace_2300708.mat'])

%%

plotColorspace(x,'_'); 