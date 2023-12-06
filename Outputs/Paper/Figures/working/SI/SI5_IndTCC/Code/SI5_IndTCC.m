clear, clc, close all

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
addpath(genpath('../../../../../../../Analyses/'))

%%

if strcmp(AnalysisDepth,'fromRawData')
    % load data
    % save data (?)

    % AnalysisDepth = 'fromPreProcessedData';
end

if strcmp(AnalysisDepth,'fromPreProcessedData')
    % fit model
    % filename(participant) = 

    % AnalysisDepth = 'fromModelOutput';
end

if strcmp(AnalysisDepth,'fromModelOutput')

    modelOutputDir = '../../../../../../../Analyses/TCCModels';
    modelOutputFiles = {...
        'Pollux/211012_124119_Pollux/210422--211012_Pollux_TCC-FreeSimilarityMatrix-workspace_230222_x.mat',...
        'Castor/211108_090705_Castor/220517--211108_Castor_TCC-FreeSimilarityMatrix-workspace_230225_x.mat',...
        'Morty/220823_081207_Morty/220322--220823_Morty_TCC-FreeSimilarityMatrix-workspace_230213_x.mat',...
        'Buster/210609_124628_Buster/210428--210609_Buster_TCC-FreeSimilarityMatrix-workspace_230213_x.mat'};

    filename = cell(1,length(modelOutputFiles));
    for participant = 1:length(modelOutputFiles)
        x(participant) = load([modelOutputDir,filesep,modelOutputFiles{participant}], 'x');
        [~,filename{participant}] = fileparts(modelOutputFiles{participant});
    end
end

%% Plot data

for participant = 1:length(modelOutputFiles)
    plotSimilarityMatrix(x(participant).x,filename{participant},'../')
end
