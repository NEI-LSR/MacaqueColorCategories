clear, clc, close all
rng(0)

% !
% Before starting, `cd` to the location of this script

%%

% TODO write in options for other levels of reproduction

addpath(genpath('../../../../../../Analyses/'))
load('../../../../../../Analyses/TCCModels/Castor/211108_090705_Castor/220517--211108_Castor_TCC-FreeSimilarityMatrix-workspace_230225_x.mat')

%%

categoryCenter = 18;
plotSimilarityMatrix(x,'18','../',categoryCenter)
plotSimilarityMatrix(x,'','../')

%%

plotbar_NLL_AIC_BIC('Castor','../')

%% Plot some gaussians to show the difference between the two situations (cherry pick for explanatory benefit)

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..'];
addpath(genpath(repoHomeDir))
modelOutputDir = [repoHomeDir,filesep,'Analyses'];

filename{2} = '210517--211108_Castor_data';

rng(0)

for participant = 2
    ModelFile = dir([modelOutputDir,filesep,'MixtureModels',filesep,...
        filename{participant},'*.mat']);
    if length(ModelFile) > 1
        warning('Multiple model files for this participant. Using most recent.')
        [~,idx] = sort([ModelFile.datenum]);
        ModelFile = ModelFile(idx);
        ModelFile = ModelFile(end);
    end
    load([modelOutputDir,filesep,'MixtureModels',filesep,ModelFile.name],'model')
end

%%

hueIndex = 0:360/64:360-(360/64);
[~,closestToZero] = min(abs(hueIndex - model.interp_crossing)')

stimCols = generateStimCols('nBig',64);
stimCols_sRGB = LabTosRGB([repelem(76.0693, 64); stimCols]);

ylim([0,0.6])
for j = 8%1:10
    figure, hold on

    xline(0,'k:','LineWidth', 2, 'HandleVisibility','off')
    xlim([-90,90])
    for i = [-j,j]
        p = plot(model.gaussfits{closestToZero(2)+i});
        p.Color = stimCols_sRGB(closestToZero(2)+i,:);
        p.LineWidth = 2;
        p.DisplayName = ['Gaussian for cue ',num2str(closestToZero(2)+i)];
    end
    % title(j)
end

yticks([0,0.6])
ylabel('Choice Probability')
xticks(-90:45:90)

saveas(gcf,fullfile('../',['F5_CastorCogBias_','exGaussian_', datestr(now,'yymmdd-HHMMSS'), '.svg']))
