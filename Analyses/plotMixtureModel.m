clc, clear, close all

%%

load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\combined_TCC-FreeSimilarityMatrix-workspace_230509.mat',...
    'cleandata')

%%

[moving_bias, lower_95, upper_95] = categorybias_analysis(cleandata);

warning('files saved - go delete them')

%% Plot

stimCols = generateStimCols('nBig',64);

f2 = figure('Position',[360 123 582 495]);
ax1 = axes();
hold on

interval = 360/64;

plot((moving_bias/interval)'+(0:63),1:64,'k','LineWidth',1)

% plot((lower_95/interval)'+(0:63),1:64,'k:');
% plot((upper_95/interval)'+(0:63),1:64,'k:');
h = fill([(lower_95/interval)'+(0:63), fliplr((upper_95/interval)'+(0:63))],[1:64, fliplr(1:64)], 'k');
set(h, 'facealpha', .1, 'LineStyle', 'none');


plot([1,64],[1,64],'k--')

axis square tight
ax1.YDir = 'reverse';
xticks([])
yticks([])
xticklabels([]);
yticklabels([]);

%
ax2 = axes('Visible','off');
axis equal tight
axis off
ax2.Box = 'off';

cb2 = colorbar;

lambda = 1:64;

stimCols_sRGB = LuvTosRGB([ones(1,64)*76.0693;stimCols]);
ax = gca;
ax.Colormap = stimCols_sRGB;
ax.CLim = [min(lambda) max(lambda)];

cb2.Location = 'southoutside';

cb2.Ticks = [];
%cb2.Label.String = "Choice";
cb2.Color = 'none';
cb2.Box = 'off';
% cb2.TickDirection = "out";
ax.XTickLabels = [];
ax.XLabel = [];
ax.YTickLabels = [];
ax.YLabel = [];

%
% Add colorbar
% h/t: https://blogs.mathworks.com/steve/2020/08/18/making-color-spectrum-plots-part-3/

ax3 = axes('Visible','off');
axis equal tight
axis off
ax3.Box = 'off';

cb3 = colorbar;

lambda = 1:64;

ax = gca;
ax.Colormap = stimCols_sRGB(end:-1:1,:);
ax.CLim = [min(lambda) max(lambda)];

cb3.Location = 'westoutside';

cb3.Ticks = [];
%cb3.Label.String = "Cue";
cb3.Color = 'none';
cb3.Box = 'off';
ax.XTickLabels = [];
ax.XLabel = [];
ax.YTickLabels = [];
ax.YLabel = [];

%

hlink = linkprop([ax1,ax2,ax3],{'Position','DataAspectRatio'});

cb2_Position = cb2.Position;
cb3_Position = cb3.Position;

%ax_Position = ax.Position;

cb2.Position = [cb2_Position(1:3),cb3_Position(2)-cb2_Position(2)];
cb3.Position = [cb3_Position(1:2),cb2_Position(1)-cb3_Position(1),cb3_Position(4)];


%%

saveas(f2,['combined_TCC-MixtureModel',datestr(now,'yymmdd'),'.svg'])





