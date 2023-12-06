clear, clc, close all

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromPreProcessedData';

%% Behind the scenes...

% Add path to required script
addpath(genpath('../../../../../../../Analyses/'))

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..'];
addpath(genpath(repoHomeDir))
modelOutputDir = [repoHomeDir,filesep,'Analyses'];

filename{1} = '210422--211012_Pollux_data';
filename{2} = '210517--211108_Castor_data';
filename{3} = '220322--220823_Morty_data';
filename{4} = '210428--210609_Buster_data';

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')
    error('not coded yet')
    % load data
    % pass on to 'fromPreProcessedData'
end

if strcmp(AnalysisDepth,'fromPreProcessedData')

    % Load data
    for participant = 1:4
        data{participant} = load([repoHomeDir,filesep,'Data',filesep,filename{participant},'.mat']);
        % data{participant} = readtable([DataDir,filesep,filename{participant},'.csv']);
    end

    % Fit model, save model data
    for participant = 1:4
        rng(0) % the modelling might be probabilistic - TODO check this
        model = fitMixtureModel(data{participant},0);
        save([modelOutputDir,filesep,'MixtureModels',filesep,filename{participant},'_',...
            datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            'model')
        allModels{participant} = model;
    end
end

if strcmp(AnalysisDepth,'fromModelOutput')
    % Load models
    for participant = 1:4
        ModelFile = dir([DataDir,filesep,'MixtureModels',filesep,...
            filename{participant},'*.mat']);
        if length(ModelFile) > 1
            warning('Multiple model files for this participant. Using most recent.')
            [~,idx] = sort([ModelFile.datenum]);
            ModelFile = ModelFile(idx);
            ModelFile = ModelFile(end);
        end
        load([DataDir,filesep,'MixtureModels',filesep,ModelFile.name],'model')
        allModels{participant} = model;
    end
end

%% Plot data

% Castor hack (to remove purple category that doesn't cross confidence interval)

allModels{1,2}.interp_crossing(4) = [];
allModels{1,2}.interp_ci(8) = [];
allModels{1,2}.interp_ci(7) = [];
allModels{1,2}.crossing_colvals(4,:) = [];
allModels{1,2}.ci(:,4) = [];


whichFigures.MixMod_polar  = true;
whichFigures.MixMod_linear = true;

for participant = 1:4
    plotMixtureModel(allModels{participant},...
        whichFigures,[filename{participant},'_',AnalysisDepth])
end


