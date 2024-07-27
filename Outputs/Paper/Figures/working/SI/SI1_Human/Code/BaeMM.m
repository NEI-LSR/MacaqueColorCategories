clear, clc, close all

convertToCIELUV = false;
bootstrap_      = false;
downsampleTo64  = false;

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

    pltCols = double(LabTosRGB(stimcols_CIELAB));

    figure, hold on
    scatter(stimcols_CIELUV(2,:),stimcols_CIELUV(3,:),...
        [],pltCols./255,'filled')
    if downsampleTo64
        stimcols_CIELAB = [ones(1,64)*70; generateStimCols('nBig',64,'sat',38)];
        pltCols = double(LabTosRGB(stimcols_CIELAB));
    end
    scatter(stimcols_CIELAB(2,:),stimcols_CIELAB(3,:),...
        [],pltCols./255,'filled')

    axis equal
    % view(3)

    % What are the closest stimuli?
    for i = 1:180

        % denote equivalent stimuli (treats the inner ring as CIELAB)
        % plot([stimcols_CIELAB(2,i),stimcols_CIELUV(2,i)],[stimcols_CIELAB(3,i),stimcols_CIELUV(3,i)],'k')

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
    if downsampleTo64
        data.trialdata.choices(i,:)   = {1:64};
    else
        data.trialdata.choices(i,:)   = {1:180};
    end
    data.trialdata.chosen(i,1)    = {choice_index(i)};
end
if downsampleTo64
    data.nBig = 64;
else
    data.nBig = 180;
end

%%

includeCorrect = true;
if downsampleTo64
    lengthOfSlidingWindow = 5; % picked by hand
else
    lengthOfSlidingWindow = 15; % picked by hand
end

model = fitMixtureModel(data,lengthOfSlidingWindow,includeCorrect);

%% For direct comparison to Bae+ 2015 F8b

% figure,
% plot(0:2:358,-model.bias)
% xticks(0:60:360)
% yticks([-17,-11,-6,0,6,11,17])

%% Bootstrap

if bootstrap_

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

end

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;
whichFigures.GaussianWidth   = true;

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
    choiceProb = model.choice_probability'; % transposing to match similarity matrix (so cue on x-axis, choice on y-axis)
catch
    % choiceProb_diag = model{1,1}.choice_probability;
    error('model variable structure is nested') % TODO work out why
end

for i = 1:size(choiceProb,1)
    choiceProb_diag(i,:) = circshift(choiceProb(i,:),i-(size(choiceProb,1))/2);
end

choiceProb_diag = choiceProb_diag/max(choiceProb_diag(:));

filename = 'CP_Bae';

% plotSimilarityMatrix(model.choice_probability) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
plotSimilarityMatrix(choiceProb_diag*2,filename,'../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

% TODO Relabel similarity as choice probability

% h = findobj;
% h(n).Label = 'Choice Probability'; % doesn't work

%% Plot gaussians for specific areas

% Modified from `plotMixtureModel.m`
hueIndex = 0:2:358;
[~,closestToZero] = min(abs(hueIndex - model.interp_crossing)')

for i = 1:length(closestToZero)
    % plotSimilarityMatrix(choiceProb_diag*2,...
    %     [filename,'_',num2str(closestToZero(i))],'../',...
    %     closestToZero(i),false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

    % figure,
    % hold on
    % plot(hueIndex-180,mean(choiceProb),'k--','DisplayName','Whole space average')
    % plot(hueIndex-180,mean(choiceProb(closestToZero(i)-25:1:closestToZero(i)+25,:)),'DisplayName','Around attractor')
    % plot(hueIndex-180,mean(choiceProb(closestToZero(i)+1:1:closestToZero(i)+25,:)),'DisplayName','After attractor')
    % plot(hueIndex-180,mean(choiceProb(closestToZero(i)-1:-1:closestToZero(i)-25,:)),'DisplayName','Before attractor')
    % 
    % legend('AutoUpdate','off')
    % xline(0,'k:')
    % xlim([-60,60])
end

plotSimilarityMatrix(choiceProb_diag*2,...
    [filename,'_',num2str(closestToZero)],'../',...
    closestToZero,false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)


%%


figure,
hold on
axis tight
% ylim([0,1])

for cueIndex = closestToZero


figure,
hold on
axis tight
% ylim([0,1])

    s = scatter(model.PotentialDistances, model.choice_probability(:,cueIndex),'filled');

    % plot(f,PotentialDistances,choice_probability(:,cueIndex),'k.')
    plot(model.gaussfits{cueIndex},...
        model.PotentialDistances(~isnan(model.choice_probability(:,cueIndex))),...
        model.choice_probability(~isnan(model.choice_probability(:,cueIndex)),cueIndex),'k.');

    p = gca;
    p.Children(2).Marker = 'none'; % turn off data, so that we can replot it how we like...
    % p.Children(1).LineWidth = 3;

    p11 = predint(model.gaussfits{cueIndex},model.PotentialDistances,0.95,'functional','off');                         %TODO: Check whether this is the appropriate type of interval: https://www.mathworks.com/help/curvefit/confidence-and-prediction-bounds.html
    plot(model.PotentialDistances,p11,'k:',...
        'DisplayName','Nonsimultaneous Functional Bounds')
    % p.Children(1).LineWidth = 1;
    % p.Children(2).LineWidth = 1;

    xline(0,'k--')

    xlabel('Error distance')
    ylabel('Choice Probability')
    % text(-90,0.9,num2str(cueIndex))
    legend('off')
    drawnow

end
% saveas(gcf,fullfile([filename,num2str(cueIndex),'_mixMod_BreakOut.svg'])) % TODO Add an argument to num2str that means that each number has 2 significant figures

%% TCC Models

data.trialdata.nBig   = 180;
data.trialdata.nSmall = 180;
data.trialdata.nTrials = size(data.trialdata.cues,1);
rn = 0;

%% D-prime and gaussian width

params = [0,1,0,0,1,0,0];

[x_dpgw,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

save(['Bae_TCC_dpgw_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

disp(x_dpgw(1))
disp(x_dpgw(2))

%% SSNU model

params = [0,0,0,1,0,0,0];

[x_ssnu,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Bae_TCC_ssnu_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

%% Offset gaussian model

params = [0,0,0,0,0,0,1];

[x_og,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Bae_TCC_og_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

%% SSNU model (reduced)

params = [0,0,0,1,0,0,0];

[x_ssnu16,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,16,[],...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Bae_TCC_ssnu16_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

%% Offset gaussian model (reduced)

params = [0,0,0,0,0,0,1];

[x_og16,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],16,...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Bae_TCC_og16_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


%% Both (reduced)

params = [0,0,0,1,0,0,1];

[x_both16,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,16,16,...
    'dPrime',           x_dpgw(1),...
    'gaussianWidth',    x_dpgw(2));

save(['Bae_TCC_both_',num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).')
% save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


%% Plot model outputs

% plotSimilarityMatrix(x, filename, '../')


