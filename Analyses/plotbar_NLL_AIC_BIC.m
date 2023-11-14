function plotbar_NLL_AIC_BIC(p,filedir)

%% Load data

% d = ['.',filesep,'TCCModels',filesep];

d = ['Y:\PROJECTS\MacaqueColorCategories_test\MacaqueColorCategories\Analyses',filesep,'TCCModels',filesep];

% p = 'Combined'; % participant
% p = 'Castor';
% p = 'Pollux';
% p = 'Buster';
% p = 'Morty';

if strcmp(p,'Combined')
    NLLfiles = {'\NLL_230826-101803.csv','\NLL_230831-111225.csv','\NLL_230907-093101.csv'};
elseif strcmp(p,'Castor')
    NLLfiles = {'\NLL_230830-091609.csv','\NLL_230830-105121.csv'};
elseif strcmp(p,'Pollux')
    NLLfiles = {'\NLL_230901-034950.csv','\NLL_230902-203220.csv'};
elseif strcmp(p,'Buster')
    NLLfiles = {'\NLL_230901-125452.csv','\NLL_230901-143654.csv'};
elseif strcmp(p,'Morty')
    NLLfiles = {'\NLL_230903-000444.csv','\NLL_230903-023202.csv'};
end

for i = 1:4
    t = readtable([d,p,NLLfiles{1}]);
    NLL(i) = min(table2array(t(:,i)));
end
NLL(end+1) = table2array(readtable([d,p,NLLfiles{2}]));

if strcmp(p,'Combined')
    NLL = [NLL,table2array(readtable([d,p,NLLfiles{3}]))];
    [AIC,BIC] = aicbic(-NLL,[2,64,64,130,130,64^2],102738);
else
    [AIC,BIC] = aicbic(-NLL,[2,64,64,130,130],102738);
end


%% Plot data

figure('Name','NLL'),
try
    b = bar(diag(NLL([1,3,2])),'stacked');
catch
    b = bar(diag(NLL([1,3,2])),'stacked');
end
if strcmp(p,'Combined')
    b = bar(diag(NLL([1,3,2])),'stacked');
    ylim([98e3,100e3])
    yticks(98e3:1e3:100e3)
    yticklabels({'98K',[],'100K'})
    % set(b(6),'facecolor',[0.3,0.3,0.3])
elseif strcmp(p,'Castor')
    ylim([46e3,50e3])
    yticks(46e3:2e3:50e3)
    yticklabels({'46K',[],'50K'})
elseif strcmp(p,'Pollux')
    ylim([70e3,74e3])
    yticks(70e3:2e3:74e3)
    yticklabels({'70K',[],'74K'})
elseif strcmp(p,'Buster')
    ylim([25e3,27e3])
    yticks(25e3:2e3:27e3)
    yticklabels({'25K','27K'})
elseif strcmp(p,'Morty')
    ylim([59e3,61e3])
    yticks(59e3:2e3:61e3)
    yticklabels({'59K','61K'})
end
legend({'Null Hypothesis','Cognitive Bias','Non-uniform Space','Simultaneous MMM','Iterative MMM','Free Similarity Matrix'},...
    'Location','northeastoutside')
set(b(1),'facecolor',[1,1,1])
set(b(2),'facecolor',[0.7,0.7,0.7])
set(b(3),'facecolor',[0,0,0])
% set(b(2),'EdgeColor','none')
xticks([])
grid()
ylabel('NLL')

try
    text(1:length(NLL([1,3,2])), NLL([1,3,2])', ...
        num2str(NLL([1,3,2])'/1000,'%0.1fK'),'HorizontalAlignment','center','VerticalAlignment','bottom')
catch
    text(1:length(NLL([1,3,2])), NLL([1,3,2])', ...
        num2str(NLL([1,3,2])'/1000,'%0.1fK'),'HorizontalAlignment','center','VerticalAlignment','bottom')
end

figure('Name','AIC'),
try
    b = bar(diag(AIC([1,3,2])),'stacked');
catch
    b = bar(diag(AIC([1,3,2])),'stacked');
end
if strcmp(p,'Combined')
    b = bar(diag(AIC([1,3,2])),'stacked');
    ylim([196e3,200e3])
    yticks(196e3:1e3:200e3)
    yticklabels({'196K',[],[],[],'200K'})
    % set(b(6),'facecolor',[0.3,0.3,0.3])
elseif strcmp(p,'Castor')
    ylim([94e3,100e3])
    yticks(94e3:1e3:100e3)
    yticklabels({'94K',[],[],[],[],[],'100K'})
elseif strcmp(p,'Pollux')
    ylim([140e3,149e3])
    yticks(140e3:1e3:149e3)
    yticklabels({'140K',[],[],[],[],[],[],[],[],'149K'})
elseif strcmp(p,'Buster')
    ylim([51e3,52.2e3])
    yticks(51e3:1e3:52.2e3)
    yticklabels({'51K','52K'})
elseif strcmp(p,'Morty')
    ylim([119e3,122e3])
    yticks(119e3:1e3:122e3)
    yticklabels({'119K','122K'})
end
legend({'Null Hypothesis','Cognitive Bias','Non-uniform Space','Simultaneous MMM','Iterative MMM','Free Similarity Matrix'},...
    'Location','northeastoutside')
set(b(1),'facecolor',[1,1,1])
set(b(2),'facecolor',[0.7,0.7,0.7])
set(b(3),'facecolor',[0,0,0])
% set(b(2),'EdgeColor','none')
xticks([])
grid()
ylabel('AIC')
% text(1:length(AIC([1,3,2])), AIC([1,3,2])', ...
% num2str(AIC([1,3,2])'/1000,'%0.1fK'),'HorizontalAlignment','center','VerticalAlignment','bottom')
try
    text(1:length(AIC([1,3,2])), AIC([1,3,2])', ...
        num2str(AIC([1,3,2])'/1000,'%0.1fK'),'HorizontalAlignment','center','VerticalAlignment','bottom')
catch
    text(1:length(AIC([1,3,2])), AIC([1,3,2])', ...
        num2str(AIC([1,3,2])'/1000,'%0.1fK'),'HorizontalAlignment','center','VerticalAlignment','bottom')
end

figure('Name','BIC'),
try
    b = bar(diag(BIC([1,3,2])),'stacked');
catch
    b = bar(diag(BIC([1,3,2])),'stacked');
end
if strcmp(p,'Combined')
    ylim([194e3,200e3])
    yticks(194e3:1e3:200e3)
    yticklabels({'194K',[],[],[],[],[],'200K'})
    % set(b(6),'facecolor',[0.3,0.3,0.3])
elseif strcmp(p,'Castor')
    ylim([94e3,100e3])
    yticks(94e3:1e3:100e3)
    yticklabels({'94K',[],[],[],[],[],'100K'})
elseif strcmp(p,'Pollux')
    ylim([142e3,148e3])
    yticks(142e3:1e3:148e3)
    yticklabels({'142K',[],[],[],[],[],'148K'})
elseif strcmp(p,'Buster')
    ylim([52e3,53e3])
    yticks(52e3:1e3:53e3)
    yticklabels({'52K','53K'})
elseif strcmp(p,'Morty')
    ylim([119e3,122e3])
    yticks(119e3:1e3:122e3)
    yticklabels({'119K','122K'})
end
legend({'Null Hypothesis','Cognitive Bias','Non-uniform Space','Simultaneous MMM','Iterative MMM','Free Similarity Matrix'},...
    'Location','northeastoutside')
set(b(1),'facecolor',[1,1,1])
set(b(2),'facecolor',[0.7,0.7,0.7])
set(b(3),'facecolor',[0,0,0])
% set(b(2),'EdgeColor','none')
xticks([])
grid()
ylabel('BIC')
try
    text(1:length(BIC([1,3,2])), BIC([1,3,2])', ...
        num2str(BIC([1,3,2])'/1000,'%0.1fK'),'HorizontalAlignment','center','VerticalAlignment','bottom')
catch
    text(1:length(BIC([1,3,2])), BIC([1,3,2])', ...
        num2str(BIC([1,3,2])'/1000,'%0.1fK'),'HorizontalAlignment','center','VerticalAlignment','bottom')
end

%% Save

saveas(findobj( 'Type', 'Figure', 'Name', 'NLL'),[filedir,'NLL_',p,'_',datestr(now,'yymmdd-HHMMSS'),'.svg'])
saveas(findobj( 'Type', 'Figure', 'Name', 'AIC'),[filedir,'AIC_',p,'_',datestr(now,'yymmdd-HHMMSS'),'.svg'])
saveas(findobj( 'Type', 'Figure', 'Name', 'BIC'),[filedir,'BIC_',p,'_',datestr(now,'yymmdd-HHMMSS'),'.svg'])


