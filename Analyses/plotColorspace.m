function plotColorspace(x,filename)

%% Load default data

if nargin == 0
    error('filename input required')
end

%%

nBig = length(x);
[stimCols,pol] = generateStimCols('nBig',nBig);
stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols(1,:);stimCols(2,:)]);
colvals = im2double(stimCols_sRGB);

x_norm  =  x(1:nBig)*(360 /sum(x(1:nBig)));
stimCols_x_polar = cumsum([0, x_norm(1:end-1)]);

rotVals = 1:10000; % arbitrary value, determines precision
for rotVal = rotVals
    stimCols_x_polar_temp = stimCols_x_polar + rotVals(rotVal)/max(rotVals)*360;
    sq_error(rotVal) = sum(rad2deg(angdiff(deg2rad(pol(1,:)),deg2rad(stimCols_x_polar_temp))).^2); %angdiff from submodule: spatialmath-matlab (https://github.com/petercorke/spatialmath-matlab)
end

% figure, 
% plot(rotVals/max(rotVals)*360,sq_error,'-o')
% xlabel('Rotation (degrees)')
% ylabel('Sum squared error')
% axis tight

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

%saveas(gcf,fullfile(['behaviorally-derived-colorspace_unrotated', datestr(now,'yymmdd-HHMMSS'), '.svg']))

%%

cla

% scatter(stimCols_x_rotated(1,:),stimCols_x_rotated(2,:),100,colvals,"filled")
scatter(stimCols_x(1,1:2:end),stimCols_x(2,1:2:end),200,colvals(1:2:end,:),"filled")

% saveas(gcf,fullfile([filename,'_', datestr(now,'yymmdd-HHMMSS'), '.svg']))
saveas(gcf,['../',filename,'_behaviorally-derived-colorspace_everySecond', datestr(now,'yymmdd-HHMMSS'), '.svg'])

%%

cla
% scatter(stimCols(1,:),stimCols(2,:),100,colvals,"filled")
scatter(stimCols(1,1:2:end),stimCols(2,1:2:end),200,colvals(1:2:end,:),"filled")

saveas(gcf,fullfile(['../','colorspace_everySecond', datestr(now,'yymmdd-HHMMSS'), '.svg']))

%%

% TODO Inversion - specify in behavioral space, plot in CIELUV?