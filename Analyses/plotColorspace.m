clear, clc, close all


%%
load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\combined_TCC-0att_fullremap-workspace_230510.mat')


%%

figure,
hold on
axis equal
axis off

nBig = 64;
[stimCols,pol] = generateStimCols('nBig',nBig);
stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols(1,:);stimCols(2,:)]);
colvals = im2double(stimCols_sRGB);

x_norm  =  x(1:nBig)*(360 /sum(x(1:nBig)));
[stimCols_x(1,:),stimCols_x(2,:)] = pol2cart(deg2rad(cumsum([0, x_norm(1:end-1)])), pol(2,:));

%scatter(stimCols_x(1,:),stimCols_x(2,:),100,colvals,"filled")
scatter(stimCols_x(1,1:2:end),stimCols_x(2,1:2:end),200,colvals(1:2:end,:),"filled")

%saveas(gcf,fullfile(['behaviorally-derived-colorspace_', datestr(now,'yymmdd'), '.svg']))
saveas(gcf,fullfile(['behaviorally-derived-colorspace_everySecond', datestr(now,'yymmdd'), '.svg']))

%%

cla
% scatter(stimCols(1,:),stimCols(2,:),100,colvals,"filled")
scatter(stimCols(1,1:2:end),stimCols(2,1:2:end),200,colvals(1:2:end,:),"filled")

saveas(gcf,fullfile(['colorspace_everySecond', datestr(now,'yymmdd'), '.svg']))



