clear, clc, close all


%%
load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\combined_TCC-0att_fullremap-workspace_230510.mat')

%%

nBig = 64;
[stimCols,pol] = generateStimCols('nBig',nBig);
stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols(1,:);stimCols(2,:)]);
colvals = im2double(stimCols_sRGB);

x_norm  =  x(1:nBig)*(360 /sum(x(1:nBig)));
stimCols_x_polar = cumsum([0, x_norm(1:end-1)]);

rotVals = 1:10000; % arbitrary value, determines precision
for rotVal = rotVals
    stimCols_x_polar_temp = stimCols_x_polar + rotVals(rotVal)/max(rotVals)*360;
    sq_error(rotVal) = sum(angdiff(deg2rad(pol(1,:)),deg2rad(stimCols_x_polar_temp)).^2); %angdiff from submodule: spatialmath-matlab (https://github.com/petercorke/spatialmath-matlab)
end

figure, 
plot(rotVals,sq_error,'-o')
axis tight

[~,minLoc] = min(sq_error);
rotationVal = rotVals(minLoc)/max(rotVals)*360;
stimCols_x_polar_rotated = stimCols_x_polar + rotationVal;

[stimCols_x(1,:),stimCols_x(2,:)]                   = pol2cart(deg2rad(stimCols_x_polar),           pol(2,:));
[stimCols_x_rotated(1,:),stimCols_x_rotated(2,:)]   = pol2cart(deg2rad(stimCols_x_polar_rotated),   pol(2,:));


%%
figure,
hold on
axis equal
axis off

scatter(stimCols_x(1,:),stimCols_x(2,:),100,colvals,"filled")

%saveas(gcf,fullfile(['behaviorally-derived-colorspace_unrotated', datestr(now,'yymmdd'), '.svg']))

%%

cla

scatter(stimCols_x_rotated(1,:),stimCols_x_rotated(2,:),100,colvals,"filled")
%scatter(stimCols_x(1,1:2:end),stimCols_x(2,1:2:end),200,colvals(1:2:end,:),"filled")

saveas(gcf,fullfile(['behaviorally-derived-colorspace_', datestr(now,'yymmdd'), '.svg']))
%saveas(gcf,fullfile(['behaviorally-derived-colorspace_everySecond', datestr(now,'yymmdd'), '.svg']))

%%

cla
scatter(stimCols(1,:),stimCols(2,:),100,colvals,"filled")
%scatter(stimCols(1,1:2:end),stimCols(2,1:2:end),200,colvals(1:2:end,:),"filled")

%saveas(gcf,fullfile(['colorspace_everySecond', datestr(now,'yymmdd'), '.svg']))

