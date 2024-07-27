% Before starting, `cd` to the location of this script

clc, clear, close all

% To regenerate from scratch, see:
% MacaqueColorCategories\Analyses\bootstrapping\bootstrap.m

%% 

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..'];
addpath(genpath([repoHomeDir,filesep,'Analyses']))

% d = dir([repoHomeDir,filesep,'Analyses',filesep,'bootstrapping',filesep,'models',filesep,'monkey',filesep,'*single-og*']);
d = dir([repoHomeDir,filesep,'Analyses',filesep,'bootstrapping',filesep,'models',filesep,'mechTurk',filesep,'*single-og*']);

for i = 1:length(d)
    load(d(i).name)
    nll_reloaded(i,1) = nll_x;
end

% d = dir([repoHomeDir,filesep,'Analyses',filesep,'bootstrapping',filesep,'models',filesep,'monkey',filesep,'*single-ssnu*']);
d = dir([repoHomeDir,filesep,'Analyses',filesep,'bootstrapping',filesep,'models',filesep,'mechTurk',filesep,'*single-ssnu*']);

for i = 1:length(d)
    load(d(i).name)
    nll_reloaded(i,2) = nll_x;
end

%%

figure, hold on
scatter(nll_reloaded(:,1),nll_reloaded(:,2))
xlabel('og')
ylabel('ssnu')
title('NLL')

plot([min(xlim),max(xlim)],[min(xlim),max(xlim)])

axis equal square

%%

clear bic
% [~,bic(:,1)] = aicbic(-nll_reloaded(:,1),66,98104);
% [~,bic(:,2)] = aicbic(-nll_reloaded(:,2),66,98104);
[~,bic(:,1)] = aicbic(-nll_reloaded(:,1),66,46000);
[~,bic(:,2)] = aicbic(-nll_reloaded(:,2),66,46000);

figure, hold on
scatter(bic(:,1),bic(:,2),'k')
xlabel('og')
ylabel('ssnu')
title('BIC')

xlim_ = xlim;
% plot([1.955,2]*10^5,[1.955,2]*10^5,'k')
plot([min(xlim_),max(xlim_)],[min(xlim_),max(xlim_)],'k')

axis equal square tight

%%

saveas(gcf,['../','bootstrap_',datestr(now,'yymmdd-HHMMSS'),'.svg'])


