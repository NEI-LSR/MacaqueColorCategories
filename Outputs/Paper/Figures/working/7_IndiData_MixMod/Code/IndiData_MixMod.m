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
    % save data (?)
end

if strcmp(AnalysisDepth,'fromPreProcessedData')

    % dataDir = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data'; %fix this!!!!!!!!!!!!!!!!!!
    % dataFiles = {...
    %     '210422--211012_Pollux_data.mat',...
    %     '220517--211108_Castor_data.mat',...
    %     '220322--220823_Morty_data.mat',...
    %     '210428--210609_Buster_data.mat'};
    % 
    % filename = cell(1,4);
    % for participant = 1:4
    %     data(participant) = load([dataDir,filesep,dataFiles{participant}]);
    %     [~,filename{participant}] = fileparts(dataFiles{participant});
    % end



    % 
    % dirname = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\'; % Fix this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % cleandata = combineData(dirname);
    % save('../../../../../../Analyses/combinedData.mat','cleandata');
    % warning('If combinedData.mat previously existed, it has been overwritten')




    % I should set it up so that it has access to the actual data files
    % rather than relying on these packaged model outputs, but this should
    % work for now...
    modelOutputDir = '../../../../../../Analyses';
    modelOutputFiles = {...
        '211012_124119_Pollux/210422--211012_Pollux_TCC-FreeSimilarityMatrix-workspace_230222.mat',...
        '211108_090705_Castor/220517--211108_Castor_TCC-FreeSimilarityMatrix-workspace_230225.mat',...
        '220823_081207_Morty/220322--220823_Morty_TCC-FreeSimilarityMatrix-workspace_230213.mat',...
        '210609_124628_Buster/210428--210609_Buster_TCC-FreeSimilarityMatrix-workspace_230213.mat'};

    filename = cell(1,4);
    for participant = 1:4
        data{participant} = load([modelOutputDir,filesep,modelOutputFiles{participant}], 'trialdata');
        [~,filename{participant}] = fileparts(modelOutputFiles{participant});
    end
end

if strcmp(AnalysisDepth,'fromModelOutput')
    error('not coded yet')
end

%% Plot data

for participant = 1:4
    fitMixtureModel(data{participant});
    disp('Figures saved')
end


% save('../../../../../../Analyses/moving_bias.mat','moving_bias')
% warning('If moving_bias.mat previously existed, it has been overwritten')
