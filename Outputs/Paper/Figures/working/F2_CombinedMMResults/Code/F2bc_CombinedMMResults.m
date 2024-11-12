% Before starting, `cd` to the location of this script

%% Options

clc, clear, close all

% Determine the depth of the analysis that you would like to reproduce

% fromRawData:               % Generate figures from the raw data (slowest)
% fromPreProcessedData:      % Generate figures from the pre-processed data
% fromModelOutput:           % Generate figures from the model outputs only (fastest)

AnalysisDepth = 'fromPreProcessedData';

%% Behind the scenes...

% Add path to required script

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',];
addpath(genpath(repoHomeDir))
modelOutputDir = [repoHomeDir,filesep,'Analyses'];

rng(0)

%% Load/process data

if strcmp(AnalysisDepth,'fromRawData')
    
    cleandata = combineData([repoHomeDir,filesep,'Data']);
       
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
    load([modelOutputDir,filesep,'MixtureModels',filesep,ModelFile.name],'model')
end

%% Plot data

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;
whichFigures.GaussianWidth   = true;

plotMixtureModel(model,...
    whichFigures,['F2_CombinedMMResults_',AnalysisDepth])

withLabels = false;
plotMixtureModel(model,...
    whichFigures,['F2_CombinedMMResults_',AnalysisDepth],withLabels)

%% Pull out specific example gaussians

hueIndex = 0:360/64:360-(360/64);
[~,closestToZero] = min(abs(hueIndex - model.interp_crossing)')

stimCols = generateStimCols('nBig',64);
stimCols_sRGB = LabTosRGB([repelem(76.0693, 64); stimCols]);

figure, hold on
xlim([-90,90])
% xlim([-180,180])
ylim([0,0.5])
for i = [1,2]
    p = plot(model.gaussfits{closestToZero(i)});
    p.Color = stimCols_sRGB(closestToZero(i),:);
    p.LineWidth = 2;
    p.DisplayName = ['Gaussian for cue ',num2str(closestToZero(i))];
end

for i = [1,2]
    p = plot(model.gaussfits{closestToZero(i)+16});
    p.Color = stimCols_sRGB(closestToZero(i)+16,:);
    p.LineWidth = 2;
    p.LineStyle = '--';
    p.DisplayName = ['Gaussian for cue ',num2str(closestToZero(i)+16)];
end

xline(0,'k:','LineWidth', 2, 'HandleVisibility','off')

yticks([0,0.5])
ylabel('Choice Probability')
xticks(-90:45:90)
% xlabel() % TODO Is the x-axis interval or degree here?

saveas(gcf,fullfile('../',['F2_CombinedMMResults_','exGaussian_', datestr(now,'yymmdd-HHMMSS'), '.svg']))
