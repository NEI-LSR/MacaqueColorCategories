clear, clc, close all

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

for i = 1:size(indexShiftedData,1)
    data.trialdata.cues(i,1)     = {indexShiftedData.target(i)};
    data.trialdata.choices(i,:)  = {1:360};
    data.trialdata.chosen(i,1)   = {indexShiftedData.response(i)};
end

data.nBig   = 360;
data.nSmall = 360;

%%

Lab = 1;
includeCorrect = true;
lengthOfSlidingWindow = 29; % picked by hand

model = fitMixtureModel(data,lengthOfSlidingWindow,includeCorrect);

model.stimColorSpace    = 'CIELAB';
model.stimCols          = [60,52];      %L*, chroma

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

plotMixtureModel(model,...
    whichFigures,'Panichello_')

%%
% 
% forceCIELUV = true;
% whichFigures.MixMod_polar    = true;
% whichFigures.MixMod_linear   = false;
% 
% plotMixtureModel(model,...
%     whichFigures,'Panichello_',[],forceCIELUV)


