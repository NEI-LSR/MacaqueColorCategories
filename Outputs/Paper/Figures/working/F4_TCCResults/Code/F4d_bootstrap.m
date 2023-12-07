% Before starting, `cd` to the location of this script

clc, clear, close all

%% 

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',];
addpath(genpath([repoHomeDir,filesep,'Analyses']))

d = dir([repoHomeDir,filesep,'Analyses',filesep,'bootstrapping',filesep,'models',filesep,'*single-og*']);

for i = 1:length(d)
    load(d(i).name)
    nll_reloaded(i,1) = nll_x;
end

d = dir([repoHomeDir,filesep,'Analyses',filesep,'bootstrapping',filesep,'models',filesep,'*single-ssnu*']);

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
[~,bic(:,1)] = aicbic(-nll_reloaded(:,1),66,98104);
[~,bic(:,2)] = aicbic(-nll_reloaded(:,2),66,98104);

figure, hold on
scatter(bic(:,1),bic(:,2),'k')
xlabel('og')
ylabel('ssnu')
title('BIC')

plot([1.955,2]*10^5,[1.955,2]*10^5,'k')

axis equal square tight

%%

saveas(gcf,['../','bootstrap_',datestr(now,'yymmdd-HHMMSS'),'.svg'])


