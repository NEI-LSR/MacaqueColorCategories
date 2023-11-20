clc, clear, close all

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath([repoHomeDir,'Analyses']))

modelOutputDir = [repoHomeDir,'Analyses'];

%% Load model

% (re-using code from F2bc_CombinedMMResults.m)

ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
    'combined_*.mat']);
% ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
%     '*Pollux_*.mat']);
% ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
%     '*Castor*.mat']); % TODO Work out why Castor's model matches the combined model currently
% ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
%     '*Buster_*.mat']);
% ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
%     '*Morty_*.mat']);

if length(ModelFile) > 1
    warning('Multiple model files. Using most recent.')
    [~,idx] = sort([ModelFile.datenum]);
    ModelFile = ModelFile(idx);
    ModelFile = ModelFile(end);
end
load([modelOutputDir,'/MixtureModels/',ModelFile.name],'model')

%% Plot choice probability

try
    choiceProb_diag = model.choice_probability;
catch
    choiceProb_diag = model{1,1}.choice_probability;
    warning('model variable structure is nested') % TODO work out why
end

for i = 1:size(choiceProb_diag,1)
    choiceProb_diag(:,i) = circshift(choiceProb_diag(:,i),i-32);
end

filename = ['CP_',ModelFile.name(1:end-4)];

% plotSimilarityMatrix(model.choice_probability) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
plotSimilarityMatrix(choiceProb_diag,filename,'../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

% TODO Relabel simaility as choice probability

% h = findobj;
% h(n).Label = 'Choice Probability'; % doesn't work

%%

