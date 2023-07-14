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
addpath(genpath('../../../../../../Analyses/'))

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')
    error('not coded yet')
    % load data
    % pass on to 'fromPreProcessedData'
end

if strcmp(AnalysisDepth,'fromPreProcessedData')

    modelOutputDir = '../../../../../../Analyses';
    modelOutputFiles = {...
        '211012_124119_Pollux/210422--211012_Pollux_TCC-FreeSimilarityMatrix-workspace_230222.mat',...
        '211108_090705_Castor/220517--211108_Castor_TCC-FreeSimilarityMatrix-workspace_230225.mat',...
        '220823_081207_Morty/220322--220823_Morty_TCC-FreeSimilarityMatrix-workspace_230213.mat',...
        '210609_124628_Buster/210428--210609_Buster_TCC-FreeSimilarityMatrix-workspace_230213.mat'};

    % Load data
    filename = cell(1,4);
    for participant = 1:4
        data{participant} = load([modelOutputDir,filesep,modelOutputFiles{participant}], 'trialdata');
        [~,filename{participant}] = fileparts(modelOutputFiles{participant});
    end

    % Fit model, save model data
    for participant = 1:4
        rng(0) % the modelling might be probabilistic - TODO check this
        model{participant} = fitMixtureModel(data{participant},0);
        save([modelOutputDir,'/MixtureModels/',filename{participant},'_',...
            datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            'model')
    end
end

if strcmp(AnalysisDepth,'fromModelOutput')
    % Load models

end

%% Plot data

clc

whichFigures.MixMod_polar = true;

for participant = 1:4
    plotMixtureModel(model{participant},...
        whichFigures,filename{participant})
end


