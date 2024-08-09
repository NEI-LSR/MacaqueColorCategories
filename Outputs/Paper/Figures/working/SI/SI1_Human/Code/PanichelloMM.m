clear, clc, close all

convertToCIELUV = false;
bootstrap_      = false;
downsampleTo64  = false;

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];

addpath(genpath([repoHomeDir,'Analyses']))

%% Load data

loadedData = readtable([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Panichello_2019',filesep,'41467_2019_11298_MOESM4_ESM.xlsx'],...
    'Sheet','Figure1C','Range',"A2:B17777");

%%

% I'm making the assumption that the numbers are hue angles, rather than
% index numbers (with with interval of 1degree it won't matter much either
% way)

indexShiftedData = loadedData;

indexShiftedData.target     = loadedData.target   + 1; % this is needed because our model fitting assumes that stimulus #1 is at 0degrees, which does not appear to be the case for this dataset
indexShiftedData.response   = loadedData.response + 1;

indexShiftedData.target(indexShiftedData.target>360)        = indexShiftedData.target(indexShiftedData.target>360) - 360; % not actually neccessary, but here for completeness
indexShiftedData.response(indexShiftedData.response>360)    = indexShiftedData.response(indexShiftedData.response>360) - 360; % neccessary, because there are some responses above 360

cue_index       = [indexShiftedData.target];
choice_index    = [indexShiftedData.response];

%% Convert to CIELUV

if convertToCIELUV
    % What is the closest analogouse CIELUV value for each of the stimuli?
    % Requires PsychToolbox

    [cart,pol] = generateStimCols('nBig',360,'sat',52);
    stimcols_CIELAB = [ones(1,360)*60; cart];
    whitePoint = xyYToXYZ([0.3184, 0.3119, 48.64]'); % white point not specified in paper, so using same as Bae (2015)
    whitePoint = whitePoint ./ whitePoint(2);
    stimcols_XYZ = LabToXYZ(stimcols_CIELAB,whitePoint);
    stimcols_CIELUV = XYZToLuv(stimcols_XYZ,whitePoint);

    pltCols = double(LabTosRGB(stimcols_CIELAB));

    figure, hold on
    scatter(stimcols_CIELUV(2,:),stimcols_CIELUV(3,:),...
        [],pltCols./255,'filled')
    if downsampleTo64
        stimcols_CIELAB = [ones(1,64)*70; generateStimCols('nBig',64,'sat',52)];
        pltCols = double(LabTosRGB(stimcols_CIELAB));
    end
    scatter(stimcols_CIELAB(2,:),stimcols_CIELAB(3,:),...
        [],pltCols./255,'filled')
    axis equal
    % view(3)

    % What are the closest stimuli?
    for i = 1:360

        % denote equivalent stimuli (treats the inner ring as CIELAB)
        % plot([stimcols_CIELAB(2,i),stimcols_CIELUV(2,i)],[stimcols_CIELAB(3,i),stimcols_CIELUV(3,i)],'k')

        dists = sqrt(...
            (stimcols_CIELUV(2,i)-stimcols_CIELAB(2,:)).^2 +...
            (stimcols_CIELUV(3,i)-stimcols_CIELAB(3,:)).^2);
        [~,LUT(i)] = min(abs(dists));

        % denote closest stimuli (treats the inner ring as CIELUV)
        plot([stimcols_CIELAB(2,LUT(i)),stimcols_CIELUV(2,i)],[stimcols_CIELAB(3,LUT(i)),stimcols_CIELUV(3,i)],'r')
    end

    cue_index    = LUT(cue_index)';
    choice_index = LUT(choice_index)';

    [~,DKLpolesIndex] = computeDKL_XYZ('Panichello');
    DKLpolesIndex_CIELUV = LUT(DKLpolesIndex);
    DKLpolesDegrees_CIELUV = (DKLpolesIndex_CIELUV - 1);
end

%%


for i = 1:size(indexShiftedData,1)
    data.trialdata.cues(i,1)     = {cue_index(i)};
    if downsampleTo64
        data.trialdata.choices(i,:)   = {1:64};
    else
        data.trialdata.choices(i,:)   = {1:360};
    end
    data.trialdata.chosen(i,1)   = {choice_index(i)};
end

if downsampleTo64
    data.nBig = 64;
else
    data.nBig = 360;
end

%%

includeCorrect = true;

if downsampleTo64
    lengthOfSlidingWindow = 9; % picked by hand
else
    lengthOfSlidingWindow = 29; % picked by hand
end

model = fitMixtureModel(data,lengthOfSlidingWindow,includeCorrect);

%% Bootstrap

if bootstrap_

    rng(0);

    for bs = 1:100
        nTrials = size(data.trialdata.cues,1);
        idx = randi(nTrials,nTrials,1);
        tempdata.trialdata.cues = data.trialdata.cues(idx);
        tempdata.trialdata.choices = data.trialdata.choices(idx);
        tempdata.trialdata.chosen = data.trialdata.chosen(idx);
        tempdata.nBig = 360;
        tempdata.nSmall = 360;

        model(bs) = fitMixtureModel(tempdata,lengthOfSlidingWindow,includeCorrect);

        nCrossings(bs) = size(model(bs).interp_crossing,1);
    end

    figure,
    hist(nCrossings)

end

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;
whichFigures.GaussianWidth   = true;

withLabels = false;
axlims = 30;

if convertToCIELUV
    filename = 'Panichello_CIELUV_';

    plotMixtureModel(model,...
        whichFigures,filename,withLabels,DKLpolesDegrees_CIELUV,axlims)
else
    filename = 'Panichello_CIELAB_';
    model.stimColorSpace    = 'CIELAB';
    model.stimCols          = [60,52];      %L*, chroma
    DKL = 'Panichello';

    plotMixtureModel(model,...
        whichFigures,filename,withLabels,DKL,axlims)
end

%% Plot choice probability matrix
% Code copied from SI6_choiceMatrices.m

try
    choiceProb_diag = model.choice_probability'; % transposing to match similarity matrix (so cue on x-axis, choice on y-axis)
catch
    % choiceProb_diag = model{1,1}.choice_probability;
    error('model variable structure is nested') % TODO work out why
end

for i = 1:size(choiceProb_diag,1)
    choiceProb_diag(i,:) = circshift(choiceProb_diag(i,:),i-(size(choiceProb_diag,1))/2);
end

choiceProb_diag = choiceProb_diag/max(choiceProb_diag(:));

filename = 'CP_Panichello';

% plotSimilarityMatrix(model.choice_probability) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
plotSimilarityMatrix(choiceProb_diag*2,filename,'../',[],false,model.stimCols) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

% TODO Relabel simaility as choice probability

% h = findobj;
% h(n).Label = 'Choice Probability'; % doesn't work

hueIndex = 0:1:359;
[~,closestToZero] = min(abs(hueIndex - model.interp_crossing)')

plotSimilarityMatrix(choiceProb_diag*2,... % multiplying by 2 just to increase visibility
    [filename,'_',num2str(closestToZero)],'../',...
    closestToZero,false,model.stimCols) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

%% TCC Models

data.trialdata.nBig   = 360;
data.trialdata.nSmall = 360;
data.trialdata.nTrials = size(data.trialdata.cues,1);
rn = 0;

%% D-prime and gaussian width

params = [0,1,0,0,1,0,0];

[x_dpgw,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

save(['Panichello_TCC_dpgw_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

disp(x_dpgw(1))
disp(x_dpgw(2))

%% SSNU model

params = [0,0,0,1,0,0,0];

[x_ssnu,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Panichello_TCC_ssnu_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

%% Offset gaussian model

params = [0,0,0,0,0,0,1];

[x_og,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Panichello_TCC_og_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

%% SSNU model (reduced)

params = [0,0,0,1,0,0,0];

[x_ssnu16,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,16,[],...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Panichello_TCC_ssnu16_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

%% Offset gaussian model (reduced)

params = [0,0,0,0,0,0,1];

[x_og16,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],16,...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Panichello_TCC_og16_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


