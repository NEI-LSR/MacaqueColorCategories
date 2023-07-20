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

    DataDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
        'Data'];
    DataFiles = {...
        '210422--211012_Pollux_data.csv',...
        '210517--211108_Castor_data.csv',...
        '220322--220823_Morty_data.csv',...
        '210428--210609_Buster_data.csv'};

    % Load data
    filename = cell(1,4);
    for participant = 1:4
        data{participant} = readtable([DataDir,filesep,DataFiles{participant}]);
        [~,filename{participant}] = fileparts(DataFiles{participant});
    end

    % Fit model, save model data
    for participant = 1:4
        rng(0) % the modelling might be probabilistic - TODO check this
        model{participant} = fitMixtureModel(data{participant},0);
        save([DataDir,'/MixtureModels/',filename{participant},'_',...
            datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            'model')
    end
end

if strcmp(AnalysisDepth,'fromModelOutput')
    % Load models

end

%% Plot data

whichFigures.MixMod_polar = true;

for participant = 1:4
    plotMixtureModel(model{participant},...
        whichFigures,filename{participant})
end


