clear, clc, close all

% which params actually makes the difference?

load('simultaneous1230814-224418.mat','x','nll_x')
disp(nll_x)

x_simul = x;
clear x
figure, hold on
plot(x_simul,'DisplayName','x\_simultaneous')

load('iterative230814-220752.mat','currentValues','nll_x')
x_iter = [currentValues.dPrime,...
    currentValues.stimulusRemapping,...
    currentValues.gaussianWidth,...
    currentValues.skewedGaussians'];
disp(nll_x)
plot(x_iter,'DisplayName','x\_iterative')
legend
ylim([0,3]) % TODO Normalise stimulusRemapping values for easier comparison

%% Compute nlls for each part 

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))
DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
if ~exist('data','var')
    data = combineData(DataDir); % TODO Switch to csv
end

data.trialdata.nBig = 64;
nBig = 64;
data.trialdata.nSmall = 4;
data.trialdata.nTrials = 98104;
warning('Assuming nBig, nSmall, and nTrials')

if size(data.trialdata.chosen,1) < size(data.trialdata.chosen,2) % TODO This is a temp hack fix for a difference between real and simulated data
    data.trialdata.chosen = data.trialdata.chosen';
    warning('transposing choices');
end

choiceInds = cell2mat(data.trialdata.choices);
cueInd = cell2mat(data.trialdata.cues);
response = cell2mat(data.trialdata.chosen);

optimisationMeta = double(zeros([6,2]));

% Which parameters?
optimisationMeta(:,1) = [0,1,0,1,1,1];

% How many values does each paramter have?
optimisationMeta(:,2) = [...
    nBig*nBig;...   % Free Similarity Matrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    nBig;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig,...        % skewedGaussians
    ];

f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
    'choiceInds',       choiceInds',...
    'cueInd',           cueInd,...
    'response',         response,...
    'nTrials',          data.trialdata.nTrials,...
    'nBig',             data.trialdata.nBig, ...
    'nSmall',           data.trialdata.nSmall,...
    'optimisationMeta', optimisationMeta);

nll_simul = f(x_simul);
disp(nll_simul)
nll_iter = f(x_iter);
disp(nll_iter)

% Good. They match.

%%

% Now let's switch out bits and see what effect they have.

nll_mod(1) = f([x_iter(1), x_simul(2:end)]);
nll_mod(2) = f([x_simul(1),  x_iter(2:65),  x_simul(66:end)]);
nll_mod(3) = f([x_simul(1),  x_simul(2:65),  x_iter(66), x_simul(67:end)]);
nll_mod(4) = f([x_simul(1),  x_simul(2:65),  x_simul(66), x_iter(67:end)]);

disp(nll_mod)

% Huh, it looks like actually the big difference is in switching in the
% value for gaussian width

disp(x_simul(66))
disp(x_iter(66))

% Oh wow yeah they are quite different
% What's the likelihood surface liek for gaussian width?

for i = 1:100
    nll_gTest(i) = f([x_simul(1),  x_simul(2:65),  i, x_simul(67:end)]);
end

figure, plot(nll_gTest)

% smooth as a baby's...

xline(x_simul(66))
xline(x_iter(66))

% and clearly the x_simul value is bad (**for this set of other parameters**)

% why?

%%

% Could it just be that lsqnonlin is real bad for gaussian width?
% What happens if I use lsqnonlin for everything in the iterative? Does it fuck it all up? 

% % temporarily in ParameterEstimator:

% if 0%length(x0) < 10 % use bads for low dimensional models
%     x = bads(f,x0,lb,ub);
% else
%     x = lsqnonlin(f,x0,lb,ub,options);
% 
%     % TODO Work out how to compute confidence intervals
%     % [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqnonlin(f,x0,lb,ub,options);
%     % CI = nlparci(x,FVAL,'jacobian',JACOB);
% end

% rn = 0;
% ParameterEstimator_caller(rn,data,'iterative');

d = dir('iterative230815-*');

for i = 1:length(d)
    load(d(i).name,'nll_x')
    all_nll(i) = nll_x;
end

figure, plot(all_nll)
axis tight
xticks(1:length(d))

disp(currentValues.gaussianWidth)

%%

d = dir('sim*230815*');

for i = 1:length(d)
    load(d(i).name,'nll_x')
    all_nll(i) = nll_x;
end

figure, plot(all_nll,'ko')
axis tight
xticks(1:length(d))

