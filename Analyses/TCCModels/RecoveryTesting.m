function [nll,aic,bic] = RecoveryTesting(realOrSimData,rn)

% Comparing model outputs

if ~exist('realOrSimData','var')
    realOrSimData = 'sim';
    warning('No value for `realOrSimData` provided. Defaulting to simulated data')
end

if strcmp(realOrSimData,'sim')

    data = generateSimulatedData;

elseif strcmp(realOrSimData,'real')

    addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))
    DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
    if ~exist('data','var')
        data = combineData(DataDir); % TODO Switch to csv
    end

    % data = load([DataDir,'210517--211108_Castor_data.mat']);
    % data = load([DataDir,'210422--211012_Pollux_data.mat']);
    % data = load([DataDir,'210428--210609_Buster_data.mat']);
    % data = load([DataDir,'220322--220823_Morty_data.mat']);
    
    % data.trialdata.nBig = 64;
    % data.trialdata.nSmall = 4;
    % data.trialdata.nTrials = 98104;
    % warning('Assuming nBig, nSmall, and nTrials')

    data.trialdata.nBig = size(data.trialdata.stimCols{1},1);
    data.trialdata.nSmall = size(data.trialdata.choices{1},2);
    data.trialdata.nTrials = size(data.trialdata.cues,1);

end

%% 

if ~exist('rn','var')
    rn = 0;
    warning('Using default rn (0)')
end

%%

[~,~,~,nll.dPrimeAndGaussianWidth,~] = ParameterEstimator_caller(rn,data,'dPrimeAndGaussianWidth');
[aic.dPrimeAndGaussianWidth,bic.dPrimeAndGaussianWidth] = aicbic(-nll.dPrimeAndGaussianWidth,2,data.trialdata.nTrials);

%%

% [~,~,~,nll.single_ssnu,~] = ParameterEstimator_caller(rn,data,'single-ssnu');
% [aic.single_ssnu,bic.single_ssnu] = aicbic(-nll.single_ssnu,64,data.trialdata.nTrials);

%%

% [~,~,~,nll.single_og,~] = ParameterEstimator_caller(rn,data,'single-og');
% [aic.single_og,bic.single_og] = aicbic(-nll.single_og,64,data.trialdata.nTrials);

%%

% [~,~,~,nll.single_sg,~] = ParameterEstimator_caller(rn,data,'single-sg');
% [aic.single_sg,bic.single_sg] = aicbic(-nll.single_sg,64,data.trialdata.nTrials);

%%

% [~,~,~,nll.simultaneous,~] = ParameterEstimator_caller(rn,data,'simultaneous');
% [aic.simultaneous,bic.simultaneous] = aicbic(-nll.simultaneous,130,data.trialdata.nTrials);

%%

% [~,~,~,nll.iterative,~] = ParameterEstimator_caller(rn,data,'iterative');
% [aic.iterative,bic.iterative] = aicbic(-nll.iterative,130,data.trialdata.nTrials);

end