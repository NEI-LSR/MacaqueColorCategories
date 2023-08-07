clear, clc, close all

addpath(genpath(['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Analyses']))

%% Load the sg model from F3_TCCModel
% TODO LATER: Swap out with something based on natural scene statistics /
% background/object color statistics?

load(['..',filesep,'..',filesep,'F3_TCCModel',filesep,'Code',filesep,'SGdata.mat'],...
    'SGdata')

%% Fit ssnu model to data

clc

params.dPrime = true;
params.gaussianWidth = true;

rn = 0;

% [x,~,~,nll_x,~] = ParameterEstimator(SGdata,params,rn);

%%

% data.trialdata.dPrime           = x(1);
% data.trialdata.gaussianWidth    = x(2);

SGdata.trialdata.dPrime           = 0.979593659581261;
SGdata.trialdata.gaussianWidth    = 24.755050722629140;

params = rmfield(params,'dPrime'); % TODO It would be nice if we could just set these to false rather than have to delete them
params = rmfield(params,'gaussianWidth');

params.stimulusRemapping = true;

[ssnuModelOfSGData,aic,bic,nll_x,x0,f] = ParameterEstimator(SGdata,params,rn);

save(['ssnuModelOfSGData',num2str(rn),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(SGdata)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


%% Plot similarity matrix of model

[~,tempdata] = f(ssnuModelOfSGData);
plotSimilarityMatrix(tempdata.trialdata.similarityMatrix,...
    'ssnuModelOfSGData','../');

%% Plot fitted ssnu model as a colorspace

plotColorspace(ssnuModelOfSGData,'../colspace_ssnuModelOfSGData');





