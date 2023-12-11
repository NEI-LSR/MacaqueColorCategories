% (Code extended/modified from SI1_IndiMM.m)

clear, clc, close all

% !
% Before starting, `cd` to the location of this script

%% Options

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromPreProcessedData';

includeCorrect = true; % include correct trials. Note, setting this to true is unusual for this dataset, and is only included here to visualise the difference between a similarity matrix and choice probability

%% Behind the scenes...

% Add path to required script
addpath(genpath('../../../../../../../Analyses/'))

DataDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Data'];
AnalysisDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Analyses'];

filename{1} = '210422--211012_Pollux_data';
filename{2} = '210517--211108_Castor_data';
filename{3} = '210428--210609_Buster_data';
filename{4} = '220322--220823_Morty_data';
filename{5} = 'combinedData';

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')
    error('not coded yet')
    % load data
    % pass on to 'fromPreProcessedData'
end

if strcmp(AnalysisDepth,'fromPreProcessedData')

    % Load data
    for participant = 1:length(filename)
        data{participant} = load([DataDir,filesep,filename{participant},'.mat']);
    end

    % Fit model, save model data
    for participant = 1:length(filename)
        rng(0) % the modelling might be probabilistic - TODO check this
        model = fitMixtureModel(data{participant}.cleandata,[],includeCorrect);
        save([AnalysisDir,filesep,'MixtureModels',filesep,filename{participant},'_',...
            datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            'model')
        allModels{participant} = model;
    end
    AnalysisDepth = 'fromModelOutput'; % pass on to next stage
end

if strcmp(AnalysisDepth,'fromModelOutput')
    % Load models
    for participant = 1:length(filename)
        ModelFile{participant} = dir([AnalysisDir,filesep,'MixtureModels',filesep,...
            filename{participant},'*.mat']);
        if length(ModelFile{participant}) > 1
            warning('Multiple model files for this participant. Using most recent.')
            [~,idx] = sort([ModelFile{participant}.datenum]);
            ModelFile{participant} = ModelFile{participant}(idx);
            ModelFile{participant} = ModelFile{participant}(end);
        end
        load([AnalysisDir,filesep,'MixtureModels',filesep,ModelFile{participant}.name],'model')
        allModels{participant} = model;
    end
end

%% Plot choice probability

for participant = 1:length(filename)
    try
        choiceProb_diag = allModels{participant}.choice_probability'; % transposing to match similarity matrix (so cue on x-axis, choice on y-axis)
    catch
        % choiceProb_diag = model{1,1}.choice_probability;
        error('model variable structure is nested') % TODO work out why
    end

    for i = 1:size(choiceProb_diag,1)
        choiceProb_diag(i,:) = circshift(choiceProb_diag(i,:),i-32);
    end

    filename = ['CP_',ModelFile{participant}.name(1:end-4)];

    % plotSimilarityMatrix(model.choice_probability) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
    plotSimilarityMatrix(choiceProb_diag,filename,'../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

    % TODO Relabel simaility as choice probability

    % h = findobj;
    % h(n).Label = 'Choice Probability'; % doesn't work
end

%%

