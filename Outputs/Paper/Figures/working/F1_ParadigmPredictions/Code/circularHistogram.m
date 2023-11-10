% requires:
% circstats (https://github.com/circstat/circstat-matlab)
% circhist (https://github.com/zifredder/CircHist)

% to add circhist to the path I had to delete the '@' from the start of the folder name
%
% to tweak the size and location of the points and bars you can change:
% `baseLineOffset` - higher makes the bars start further out
% `pointSize` - makes the colored circled larger or smaller

clear, clc, close all

rng(0)

addpath(genpath(['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Analyses']))

nBig = 64;

PointCue = 20; % selected cue
PointChoice = 20; % peak of response histogram

stimCols = generateStimCols('nBig',nBig,'sat',37);
stimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols]);


figure,

%%

%%%%%%% FAKE DATA %%%%%%%%%%%%%%%%
sDist1 = mod(rad2deg(circ_vmrnd(pi/3, 2, 1000)), 360); % generate random sample, convert to deg
sDist2 = histcounts(sDist1,nBig);
sDist2(PointCue+47) = 0;
%%%%%%% FAKE DATA END %%%%%%%%%%%%%%%%

% Plot the circular histogram:
baseLineOffset = 450;
obj1 = CircHist(sDist2,'dataType','histogram',... % if you already have the histcounts you'll need to add 'dataType','histogram'
    'drawAvgAng',false,'drawAvgAngCi',false,'drawR',false,... % get rid of stuff we don't want... (drawR might be useful)
    'baseLineOffset',baseLineOffset,'colorBar','k'); % add stuff we _do_ want

% get rid of more stuff
axis off
title('')
colorbar off
obj1.baseLineH.Visible = 'off'; % this requires changing line 412 of CircHist.m from `properties (SetAccess = private)` to `properties (SetAccess = public)`
obj1.polarAxs.Position = [0.1300 0.1100 0.7400 0.7800];

%%

axes('Position',[0.2 0.2 0.6 0.6])

th = (360/nBig)/2;
stimCols_rotatedToMatchHistogram = stimCols' * [cosd(th),-sind(th);sind(th),cosd(th)];% THIS IS A CLUDGE THAT WE SHOULD NOT USE FOR PLOTTING REAL DATA (...without thinking more about it)

pointSize = 70;
scatter(stimCols_rotatedToMatchHistogram(:,1),stimCols_rotatedToMatchHistogram(:,2),...
    pointSize,double(stimCols_sRGB)/255,'filled')
%set(gcf,'color',[0.7,0.7,0.7]);
axis equal
axis off

hold on

% add arc % h/t to https://www.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle#answer_8752

startvar = rad2deg(cart2pol(stimCols_rotatedToMatchHistogram(PointCue,1),stimCols_rotatedToMatchHistogram(PointCue,2)));
endvar = rad2deg(cart2pol(stimCols_rotatedToMatchHistogram(PointChoice,1),stimCols_rotatedToMatchHistogram(PointChoice,2)));
radius = 20;

t = linspace(startvar,endvar);
x = radius*cosd(t);
y = radius*sind(t);
x = [x 0 x(1)];
y = [y 0 y(1)];
P = fill(x,y,[0.7,0.7,0.7],'LineStyle','none');

% add outline to point
% scatter(stimCols_rotatedToMatchHistogram(whichPoint,1),stimCols_rotatedToMatchHistogram(whichPoint,2),...
%     70,double(stimCols_sRGB(whichPoint,:))/255,'filled',...
%     'MarkerEdgeColor','k')

% add arrow
quiver(0,0,stimCols_rotatedToMatchHistogram(PointCue,1),stimCols_rotatedToMatchHistogram(PointCue,2),...
    'k','LineWidth',2)

quiver(0,0,stimCols_rotatedToMatchHistogram(PointChoice,1),stimCols_rotatedToMatchHistogram(PointChoice,2),...
    'k','LineWidth',2)

scatter(0,0,4,'k','filled') % rounds off the arrows at the centre



%% save

saveas(gcf,'circHist.svg')
%saveas(gcf,'circHist.png') % turned off so that we don't add new versions of the png to the repo every time we run this and then commit