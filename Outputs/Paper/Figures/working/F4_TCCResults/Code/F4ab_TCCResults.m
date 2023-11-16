% Before starting, `cd` to the location of this script

%% Options

clc, clear, close all

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromModelOutput';
csv = 0; % csv or mat?
addpath(genpath('../../../../../../Analyses/'))

%% Behind the scenes...

% Add path to required script

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',];
addpath(genpath([repoHomeDir,filesep,'Analyses']))
DataDir = [repoHomeDir,filesep,'Analyses',filesep,'TCCModels',filesep,'combined',filesep];

%% Load/process data

if strcmp(AnalysisDepth,'fromPreProcessedData')
    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,dim1,dim2,varargin);
end


if strcmp(AnalysisDepth,'fromModelOutput')
    filename = 'combined_TCC-FreeSimilarityMatrix-workspace_230214';
    load([DataDir,filename,'.mat'],'x')
end

%% Plot model outputs

plotSimilarityMatrix(x, filename, '../')

plotSimilarityMatrix(x, [filename,'warm'], '../', 3)
plotSimilarityMatrix(x, [filename,'cool'], '../', 38)
