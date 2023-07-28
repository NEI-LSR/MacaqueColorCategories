%function 

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromModelOutput';

%% Behind the scenes...

% Add path to required script
addpath(genpath('../../../../../../Analyses/'))


modelOutputDir = '../../../../../../Analyses';

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')
    error('not coded yet')
    % load data
    % save data (?)
end

if strcmp(AnalysisDepth,'fromPreProcessedData')
    data = combineData();
    save('../../../../../../Data/combinedData.mat','data');
    warning('If combinedData.mat previously existed, it has been overwritten')

    rng(0) % the modelling might be probabilistic - TODO check this

    model = fitMixtureModel(data,0);

    save([modelOutputDir,'/MixtureModels/','combined_',...                  % TODO Add a better filename (that includes info about the data)
        datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        'model')
end

if strcmp(AnalysisDepth,'fromModelOutput')
    ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
        'combined_*.mat']);
    if length(ModelFile) > 1
        warning('Multiple model files. Using most recent.')
        [~,idx] = sort([ModelFile.datenum]);
        ModelFile = ModelFile(idx);
        ModelFile = ModelFile(end);
    end
    load([modelOutputDir,'/MixtureModels/',ModelFile.name],'model')
end

%% Plot data

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

plotMixtureModel(model,...
    whichFigures,['AvResults_',AnalysisDepth])

