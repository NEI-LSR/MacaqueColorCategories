% Before starting, `cd` to the location of this script

%% Options

clc, clear, close all

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromPreProcessedData';
csv = 0; % csv or mat?

%% Behind the scenes...

% Add path to required script

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',];
addpath(genpath([repoHomeDir,filesep,'Analyses']))
modelOutputDir = [repoHomeDir,filesep,'Analyses'];

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')

    if csv
        data = combineData([repoHomeDir,filesep,'Data'])
    else
        data = combineData_mat([repoHomeDir,filesep,'Data']);
    end
       
    saveDataFile = 0;
    if saveDataFile
        save([repoHomeDir,filesep,'Data',filesep,'combinedData.mat'])
    end

end

if strcmp(AnalysisDepth,'fromPreProcessedData')
 
    if csv
        warning('Using csv method currently results in different output. It should not. Work in progress.')
        loadedData = readtable([repoHomeDir,filesep,'Data',filesep,'combinedData.csv']);
        data = loadedData;
    else
        load([repoHomeDir,filesep,'Data',filesep,'combinedData.mat']);
    end

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
    whichFigures,['F2_CombinedMMResults_',AnalysisDepth])

