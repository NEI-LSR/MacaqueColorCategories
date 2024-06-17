% Before starting, `cd` to the location of this script

%% Options

clc, clear, close all

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromPreProcessedData';

%% Behind the scenes...

% Add path to required script

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..'];
addpath(genpath([repoHomeDir,filesep,'Analyses']))
DataDir =  [repoHomeDir,filesep,'Data',filesep];
ModelDir = [repoHomeDir,filesep,'Analyses',filesep,'TCCModels',filesep,'MTurk',filesep];

%% Load/process data


if strcmp(AnalysisDepth,'fromRawData')
    error('Not currently implemented') % TODO
end

if strcmp(AnalysisDepth,'fromPreProcessedData')

    load([DataDir,'211112--220507_MTurk.mat'])
    data = cleandata; clear cleandata
    data.trialdata.nBig = 64;
    data.trialdata.nSmall = 4;
    data.trialdata.nTrials = size(data.trialdata.cues,1);

    params = [1,0,0,0,0,0,0]; % telling the model that we want to fit a free similarity matrix model
    rn = 0; % random number for reproducibility
    dim1 = []; % parameters for dimensionality reduction, left empty here
    dim2 = [];

    x = ParameterEstimator(data,params,rn,dim1,dim2);

    filename = ['MTurk_TCC-FreeSimilarityMatrix_',datestr(now,'yymmdd-HHMMSS')];
    save([ModelDir,filename],'x')    
end

if strcmp(AnalysisDepth,'fromModelOutput')
    filename = 'MTurk_TCC-FreeSimilarityMatrix_240415-112215';
    load([ModelDir,filename,'.mat'],'x')
end

%% Plot model outputs

plotSimilarityMatrix(x, filename, '../')
