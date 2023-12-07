% Before starting, `cd` to the location of this script

%% Options

clc, clear, close all

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromModelOutput';

%% Behind the scenes...

% Add path to required script

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',];
addpath(genpath(repoHomeDir))
modelOutputDir = [repoHomeDir,filesep,'Analyses'];

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')
    
    cleandata = combineData_mat([repoHomeDir,filesep,'Data'],5); % random number seed picked by hand
       
    saveDataFile = 1;
    if saveDataFile
        save([repoHomeDir,filesep,'Data',filesep,'combinedData.mat'],'cleandata')
    end
end

if strcmp(AnalysisDepth,'fromPreProcessedData')
    load([repoHomeDir,filesep,'Data',filesep,'combinedData.mat']);
end

if strcmp(AnalysisDepth,'fromRawData') || strcmp(AnalysisDepth,'fromPreProcessedData')

    rng(0) % the modelling might be probabilistic - TODO check this

    model = fitMixtureModel(cleandata);

    if ~exist([modelOutputDir,'/MixtureModels/'],"dir")
        mkdir([modelOutputDir,'/MixtureModels/'])
    end

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
    whichFigures,['F2_CombinedMMResults_',AnalysisDepth])

withLabels = false;
plotMixtureModel(model,...
    whichFigures,['F2_CombinedMMResults_',AnalysisDepth],withLabels)


