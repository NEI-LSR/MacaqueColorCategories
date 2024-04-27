clear, clc, close all

convertToCIELUV = false;

%% Add required paths

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath(repoHomeDir));

%%

loadedData = readtable([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Bae_2015',filesep,'raw_data',filesep,'DelayedEstimation data.csv'],...
    'Range',"D1:F10801");

% loadedData = readtable([repoHomeDir,'Data',filesep,'secondaryData',filesep,'Bae_2015',filesep,'raw_data',filesep,'UndelayedEstimation data.csv'],...
%     'Range',"D1:F10801");

cue_hueAngle = loadedData.TargetColor*2;

error = loadedData.ClickAngle - loadedData.TargetAngle;
choice_hueAngle = cue_hueAngle + error;
choice_hueAngle(choice_hueAngle<0)   = choice_hueAngle(choice_hueAngle<0)+360;
choice_hueAngle(choice_hueAngle>360) = choice_hueAngle(choice_hueAngle>360)-360;

cue_index = (cue_hueAngle/2) + 1; % to account for the fact that the first angle is 2deg whereas our analyses assume that the first angle is 0deg
cue_index(cue_index == 181) = 1;

choice_index = round(choice_hueAngle/2) + 1; % to account for the fact that the first angle is 2deg whereas our analyses assume that the first angle is 0deg
choice_index(choice_index == 181) = 1;

%% Convert to CIELUV

if convertToCIELUV
    % What is the closest analogous CIELUV value for each of the stimuli?
    % Requires PsychToolbox

    [cart,pol] = generateStimCols('nBig',180,'sat',38);
    stimcols_CIELAB = [ones(1,180)*70; cart];
    whitePoint = xyYToXYZ([0.3184, 0.3119, 48.64]');
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
    for i = 1:180

        % denote equivalent stimuli (treats the inner ring as CIELAB)
        plot([stimcols_CIELAB(2,i),stimcols_CIELUV(2,i)],[stimcols_CIELAB(3,i),stimcols_CIELUV(3,i)],'k')
  
        dists = sqrt(...
            (stimcols_CIELUV(2,i)-stimcols_CIELAB(2,:)).^2 +...
            (stimcols_CIELUV(3,i)-stimcols_CIELAB(3,:)).^2);
        [~,LUT(i)] = min(abs(dists));

        % denote closest stimuli (treats the inner ring as CIELUV)
        plot([stimcols_CIELAB(2,LUT(i)),stimcols_CIELUV(2,i)],[stimcols_CIELAB(3,LUT(i)),stimcols_CIELUV(3,i)],'r')
    end

    cue_index = LUT(cue_index)';
    choice_index = LUT(choice_index)';

    [~,DKLpolesIndex] = computeDKL_XYZ('Bae');
    DKLpolesIndex_CIELUV = LUT(DKLpolesIndex);
    DKLpolesDegrees_CIELUV = (DKLpolesIndex_CIELUV - 1) * 2;
end

%% convert into "cleandata" format

for i = 1:size(loadedData,1)
    data.trialdata.cues(i,1)      = {cue_index(i)};
    data.trialdata.choices(i,:)   = {1:180};
    data.trialdata.chosen(i,1)    = {choice_index(i)};
end
data.nBig = 180;

%%

includeCorrect = true;
lengthOfSlidingWindow = 15; % picked by hand

model = fitMixtureModel(data,lengthOfSlidingWindow,includeCorrect);

%% For direct comparison to Bae+ 2015 F8b

% figure, 
% plot(0:2:358,-model.bias)
% xticks(0:60:360)
% yticks([-17,-11,-6,0,6,11,17])

%% Bootstrap 

rng(0);

for bs = 1:100
    nTrials = size(data.trialdata.cues,1);
    idx = randi(nTrials,nTrials,1);
    tempdata.trialdata.cues = data.trialdata.cues(idx);
    tempdata.trialdata.choices = data.trialdata.choices(idx);
    tempdata.trialdata.chosen = data.trialdata.chosen(idx);
    tempdata.nBig = 180;

    model(bs) = fitMixtureModel(tempdata,lengthOfSlidingWindow,includeCorrect);

    nCrossings(bs) = size(model(bs).interp_crossing,1);
end

figure,
hist(nCrossings)

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

withLabels = false;
axlims = 30;

if convertToCIELUV
    filename = 'Bae_CIELUV_';

    plotMixtureModel(model,...
    whichFigures,filename,withLabels,DKLpolesDegrees_CIELUV,axlims)
else
    filename = 'Bae_CIELAB_';
    model.stimColorSpace    = 'CIELAB';
    model.stimCols          = [70,38];      %L*, chroma
    DKL = 'Bae';

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

filename = 'CP_Bae';

% plotSimilarityMatrix(model.choice_probability) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
plotSimilarityMatrix(choiceProb_diag,filename,'../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

% TODO Relabel similarity as choice probability

% h = findobj;
% h(n).Label = 'Choice Probability'; % doesn't work


