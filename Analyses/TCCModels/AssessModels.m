clear, clc, close all

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))

% from ParameterEstimator_caller.m
DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
data = combineData_mat(DataDir); % TODO Switch to csv

%% From PE_caller_caller

% cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu')
% cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\sg')
% cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\sg_ssnu')
cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu_THEN_sg')
d = dir('*.mat');

% https://www.mathworks.com/matlabcentral/answers/397385-how-to-sort-a-structure-array-based-on-a-specific-field#answer_317198
T = struct2table(d);
sortedT = sortrows(T, 'date');
d = table2struct(sortedT);

for i = 1:length(d)
    load(d(i).name)
    all_aic(i) = aic;
    all_nll(i) = nll_x;
    all_x(:,i) = x;
end

%%

clc

nBig = 64;

% From ParameterEstimator_caller
data.trialdata.nBig = 64;
data.trialdata.nTrials = 98104;
data.trialdata.nSmall = 4;
data.trialdata.dPrime           = 1.4899;
data.trialdata.gaussianWidth    = 39.0070;

load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu\4_230803-161235.mat',...
    'x')
data.trialdata.stimulusRemappingPol = x;

% From ParameterEstimator.m
optimisationMeta = double(zeros([6,2]));

% Which parameters?
optimisationMeta(:,1) = [...
    isfield(params,'freeSimilarityMatrix');...
    isfield(params,'dPrime');...
    isfield(params,'stimulusRemappingCart');...
    isfield(params,'stimulusRemapping');...
    isfield(params,'gaussianWidth');...
    isfield(params,'skewedGaussians'),...
    ];

% How many values does each paramter have?
optimisationMeta(:,2) = [...
    nBig*nBig;...   % free similarityMatrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    nBig;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig,...        % skewedGaussians
    ];

for i = 1:length(d)
    choiceInds = cell2mat(data.trialdata.choices);
    cueInd = cell2mat(data.trialdata.cues);
    response = cell2mat(data.trialdata.chosen);

    f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
        'choiceInds',   choiceInds',...
        'cueInd',       cueInd,...
        'response',     response,...
        'nTrials',      data.trialdata.nTrials,...
        'nBig',         data.trialdata.nBig, ...
        'nSmall',       data.trialdata.nSmall,...
        'dPrime',       data.trialdata.dPrime,...
        'gaussianWidth',data.trialdata.gaussianWidth,...
        'stimulusRemappingPol',data.trialdata.stimulusRemappingPol,...
        'optimisationMeta',optimisationMeta); % (add stimCols to speed up slightly)

    [~,tempdata] = f(all_x(:,i)');
    
    filename = d(i).name(1:end-4);
    OutPutFileDir = './';
    plotSimilarityMatrix(tempdata.trialdata.similarityMatrix,filename,OutPutFileDir);

end


