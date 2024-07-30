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
    
    rn = 0; % random number for reproducibility

%     params = [1,0,0,0,0,0,0]; % telling the model that we want to fit a free similarity matrix model
%     dim1 = []; % parameters for dimensionality reduction, left empty here
%     dim2 = [];
% 
%     x = ParameterEstimator(data,params,rn,dim1,dim2);
% 
%     filename = ['MTurk_TCC-FreeSimilarityMatrix_',datestr(now,'yymmdd-HHMMSS')];
%     save([ModelDir,filename],'x')    
    
    %-% D-prime and gaussian width
    
    params = [0,1,0,0,1,0,0];
    
    
    [x_dpgw,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);
    
    filename = ['MTurk_TCC-dpgw_',datestr(now,'yymmdd-HHMMSS')];
    save([ModelDir,filename],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)
    
    disp(x_dpgw(1))
    disp(x_dpgw(2))
    
    %-% SSNU model
    
    params = [0,0,0,1,0,0,0];
    
    [x_ssnu,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
        'dPrime',           x_dpgw(1),...
        'gaussianWidth',    x_dpgw(2));
    
    filename = ['MTurk_TCC-ssnu_',datestr(now,'yymmdd-HHMMSS')];
    save([ModelDir,filename],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)
     
    %-% Offset gaussian model
    
    params = [0,0,0,0,0,0,1];
    
    [x_og,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
        'dPrime',           x_dpgw(1),...
        'gaussianWidth',    x_dpgw(2));
    
    filename = ['MTurk_TCC-og_',datestr(now,'yymmdd-HHMMSS')];
    save([ModelDir,filename],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)
    
end

if strcmp(AnalysisDepth,'fromModelOutput')
    filename = 'MTurk_TCC-FreeSimilarityMatrix_240415-112215';
    load([ModelDir,filename,'.mat'],'x')
end

%% Plot model outputs

plotSimilarityMatrix(x, filename, '../')

closestToZero = [7,30,38,61]; % from MTurkDataMM.m

plotSimilarityMatrix(x,... 
    [filename,'_',num2str(closestToZero)],'../',...
    closestToZero,false) 

%%

nTrials = 46000;

load('MTurk_TCC-dpgw_240730-172702.mat')
% plotSimilarityMatrix(x, filename, '../')
nParam = 2;
[aic,bic] = aicbic(-nll_x,nParam,nTrials)

load('MTurk_TCC-ssnu_240730-172824.mat')
% plotSimilarityMatrix(x, filename, '../')
nParam = 182;
[aic,bic] = aicbic(-nll_x,nParam,nTrials)

load('MTurk_TCC-og_240730-172857.mat')
% plotSimilarityMatrix(x, filename, '../')
nParam = 182;
[aic,bic] = aicbic(-nll_x,nParam,nTrials)


