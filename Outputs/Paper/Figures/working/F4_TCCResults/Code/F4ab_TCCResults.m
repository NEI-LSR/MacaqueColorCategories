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
DataDir =  [repoHomeDir,filesep,'Data',filesep];
ModelDir = [repoHomeDir,filesep,'Analyses',filesep,'TCCModels',filesep,'combined',filesep];

%% Load/process data

if strcmp(AnalysisDepth,'fromPreProcessedData')

    load([DataDir,'combinedData.mat'],'data')
    data.trialdata.nBig = 64;
    data.trialdata.nSmall = 4;
    data.trialdata.nTrials = size(data.trialdata.cues,1);

    params = [1,0,0,0,0,0,0]; % telling the model that we want to fit a free similarity matrix model
    rn = 0; % random number for reproducibility
    dim1 = []; % parameters for dimensionality reduction, left empty here
    dim2 = [];

    x = ParameterEstimator(data,params,rn,dim1,dim2);

    filename = ['combined_TCC-FreeSimilarityMatrix_',datestr(now,'yymmdd-HHMMSS'), '.mat'];
    save(filename,'x')    
end

if strcmp(AnalysisDepth,'fromModelOutput')
    filename = 'combined_TCC-FreeSimilarityMatrix_231116-190823';
    load([filename,'.mat'],'x')
end

%% Plot model outputs

plotSimilarityMatrix(x, filename, '../')

plotSimilarityMatrix(x, [filename,'warm'], '../', 3)
plotSimilarityMatrix(x, [filename,'cool'], '../', 38)
