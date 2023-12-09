clear, clc, close all

convertToCIELUV = false;

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

    figure, hold on
    pltCols = double(LabTosRGB(stimcols_CIELAB));
    scatter(stimcols_CIELAB(2,:),stimcols_CIELAB(3,:),...
        [],pltCols./255,'filled')
    scatter(stimcols_CIELUV(2,:),stimcols_CIELUV(3,:),...
        [],pltCols./255,'filled')
    axis equal
    % view(3)

    % What are the closest stimuli?
    for i = 1:360

        % denote equivalent stimuli (treats the inner ring as CIELAB)
        plot([stimcols_CIELAB(2,i),stimcols_CIELUV(2,i)],[stimcols_CIELAB(3,i),stimcols_CIELUV(3,i)],'k')
  
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
    data.trialdata.choices(i,:)  = {1:360};
    data.trialdata.chosen(i,1)   = {choice_index(i)};
end

data.nBig   = 360;
data.nSmall = 360;

%%

Lab = 1;
includeCorrect = true;
lengthOfSlidingWindow = 29; % picked by hand

model = fitMixtureModel(data,lengthOfSlidingWindow,includeCorrect);

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

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

