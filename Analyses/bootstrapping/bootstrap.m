clc, clear, close all

addpath(genpath('..'))

%%

dirname = ['..',filesep,'Data'];
rng(0)
bootstrap_ = true;

for i = 1:100
    rn = i; % use i as a random number generator so that new starting values are used each
    data = combineData_mat(dirname,rn,bootstrap_);

    data.trialdata.nBig = 64;
    data.trialdata.nSmall = 4;
    data.trialdata.nTrials = 98104;

    [~,~,~,nll(i).single_og,~]      = ParameterEstimator_caller(rn,data,'single-og');
    [~,~,~,nll(i).single_ssnu,~]    = ParameterEstimator_caller(rn,data,'single-ssnu');

    close all

end

%% Either use the nll from above, or reload everything in from the saved model fits

d = dir('*single-og*');
for i = 1:length(d)
    load(d(i).name)
    nll_reloaded(i,1) = nll_x;
end

d = dir('*single-ssnu*');
for i = 1:length(d)
    load(d(i).name)
    nll_reloaded(i,2) = nll_x;
end

figure, hold on
plot(nll_reloaded(:,1),'DisplayName','Offset Gaussian')
plot(nll_reloaded(:,2),'DisplayName','SSNU')
legend('Location','best')