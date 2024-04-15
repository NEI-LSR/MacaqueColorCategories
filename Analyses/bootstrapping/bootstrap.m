clc, clear, close all

addpath(genpath('../..'))

monkey_or_humanMTurk = 'humanMTurk'; % which dataset?

%%

dirname = ['..',filesep,'..',filesep,'Data'];
rng(0)
bootstrap_ = true;

for i = 1:100
    rn = i; % use i as a random number generator so that new starting values are used each
    if strcmp(monkey_or_humanMTurk,'monkey')
        data = combineData(dirname,rn,bootstrap_);

        data.trialdata.nBig = 64;
        data.trialdata.nSmall = 4;
        data.trialdata.nTrials = 98104;

    elseif strcmp(monkey_or_humanMTurk,'humanMTurk')

        repoHomeDir = ['..',filesep,'..'];
        DataDir = [repoHomeDir,filesep,'Data',filesep];
        load([DataDir,'211112--220507_MTurk.mat'])
        data = cleandata; clear cleandata
        nTrials = size(data.trialdata.cues,1);

        idx = randi(nTrials,nTrials,1);

        data.trialdata.cues = data.trialdata.cues(idx);
        data.trialdata.choices = data.trialdata.choices(idx);
        data.trialdata.chosen = data.trialdata.chosen(idx);

        data.trialdata.nBig = 64;
        data.trialdata.nSmall = 4;
        data.trialdata.nTrials = nTrials;

    end

    [~,~,~,nll(i).single_og,~]      = ParameterEstimator_caller(rn,data,'single-og');
    [~,~,~,nll(i).single_ssnu,~]    = ParameterEstimator_caller(rn,data,'single-ssnu');

    close all

end

%% Reload everything in from the saved model fits and plot

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

