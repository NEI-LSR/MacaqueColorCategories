% Generated panels a and b of figure 1 for the color categories paper
% A) stimuli in CIE u*v*
% B) stimuli in CIE 1931 xy

clc, clear, close all

uv = generateStimCols('nBig',64); % requires causal globs code

Lstar = 70;
LUV = [ones(1,size(uv,2))*Lstar;uv];

RGB = LuvTosRGB(LUV); % requires causal globs code

XYZ_D65 = [0.95047;1.0;1.08883];
XYZ = LuvToXYZ(LUV,XYZ_D65); % requires PsychToolbox
xyY = XYZToxyY(XYZ); % requires PsychToolbox

%% plot CIE u*v*

f = figure;

s = tiledlayout(1,2); % requires MATLAB R2019b. Use `subplot` otherwise (will break the use of multiple axes below, but should be fixable)

ax1 = axes(s);
scatter(LUV(2,:),LUV(3,:),...
    [],double(RGB)/255,'filled')
axis equal
xlabel('u*')
ylabel('v*')

xlim([-40,40])
ylim([-40,40])


%% second axes on u*v* to show u'v'

% This code would add a second set of axes which could be used to denote
% u'v' (in addition to u*v*, since for our conditions one is a linear
% scaling of the other)
% However, getting those second axes to show the right values would require
% some work, so I'm going to not do that until we're sure we actually want
% it.
% If we DO decide we want it:
% - calculate the scaling factors for u*v* -> u'v'
% - put those values on the second set of axes. Something like:
%       ax2.xlim = ax1.xlim * scalingFactor
%       ax2.ylim = ax1.ylim * scalingFactor 
%       (2 different scaling factors? not sure)
% - work out whether the fact that the two sets of axes aren't locked (and
% thus change whenever the figure is resized is an issue, and if so:
% whether that issue is solvable.

% ax2 = axes(s);
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.Color = 'none';
% ax1.Box = 'off';
% ax2.Box = 'off';
% 
% xlabel('u''')
% ylabel('v''')


%% plot 1931 xy

nexttile
DrawChromaticity('1931') % requires PsychToolbox
scatter(xyY(1,:),xyY(2,:),...
    [],double(RGB)/255,'filled')

%% save

saveas(gcf,'F1_AB_stimuliChromaticities_1.pdf')
saveas(gcf,'F1_AB_stimuliChromaticities_1.fig')


