clear, clc, close all

%% Pulled from `recoveryTesting.m`

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))
DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';

data = load([DataDir,'210517--211108_Castor_data.mat']);

data.trialdata.nBig = size(data.trialdata.stimCols{1},1);
data.trialdata.nSmall = size(data.trialdata.choices{1},2);
data.trialdata.nTrials = size(data.trialdata.cues,1);

%% Pulled from `ParameterEstimator_caller.m`

rn = 0;
dims = 1:2:30;

SaveDir = '.';
fittingType = 'single-ssnu';

for i = 1:length(dims)
    dim = dims(i);

    params = [0,1,0,0,1,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);
    disp(x);

    params = [0,0,0,1,0,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,dim,...
        'dPrime',           x(1),...
        'gaussianWidth',    x(2));

    save([SaveDir, filesep,...
        fittingType,num2str(rn),'_',num2str(dim),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% Re-load data

d = dir('single-ssnu0*');

for i = 1:length(d)
    load(d(i).name,'x','nll_x','dim')
    nll_train(dim) = nll_x;
    xt(dim).stimulusRemapping = x;
end

nll_train = nll_train(1:2:end);
x = xt(1:2:end);

%% Pulled from `fitModel2` (in the Color Space Geometry repo)

figure, hold on
for i = 1:length(x)
    plot(0:360/length(x(i).stimulusRemapping):360,...
        [x(i).stimulusRemapping,x(i).stimulusRemapping(1)]/...
        mean([x(i).stimulusRemapping,x(i).stimulusRemapping(1)]),... % TODO: This normalisation feels hacky
        'DisplayName',num2str(length(x(i).stimulusRemapping)));
end
axis tight
legend
ylabel('Normalised angular offset')
% yticks([])
xticks([0,90,180,270,360])

figure, 
scatter(dims,nll_train,...
    'ok','filled','MarkerFaceAlpha',0.5)
ylabel('NLL')
xlabel('Number of Parameters')

[aic,bic] = aicbic(-nll_train,...
    dims+2,...
    data.trialdata.nTrials); % TODO Remember to swap back out with data_train later

figure, 
scatter(dims,aic,...
    'ok','filled','MarkerFaceAlpha',0.5)
ylabel('AIC')
xlabel('Number of Parameters')

figure, 
scatter(dims,bic,...
    'ok','filled','MarkerFaceAlpha',0.5)
ylabel('BIC')
xlabel('Number of Parameters')

