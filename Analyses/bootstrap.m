clc, clear, close all

addpath(genpath('..'))

%%

dirname = ['..',filesep,'Data'];
rng(0)
bootstrap_ = true;

for i = 1:100
    rn = i; % use i as a random number generator so that new starting values are used each time
    data = combineData_mat(dirname,rn,bootstrap_);

    data.trialdata.nBig = 64;
    data.trialdata.nSmall = 4;
    data.trialdata.nTrials = 98104;

    [~,~,~,nll(i).single_og,~]      = ParameterEstimator_caller(rn,data,'single-og');
    [~,~,~,nll(i).single_ssnu,~]    = ParameterEstimator_caller(rn,data,'single-ssnu');

end