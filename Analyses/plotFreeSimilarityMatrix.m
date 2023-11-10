%% 

clear, clc, close all


%%

%load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\combined_TCC-FreeSimilarityMatrix-workspace_230509.mat')
load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\combined_TCC-0att_fullremap-workspace_230510.mat')

[~,tempdata] = f(x);
sm = tempdata.trialdata.similarityMatrix;
stimCols = generateStimCols('nBig',64);

%sm = reshape(x,[nBig, nBig]);

%%
f2 = figure('Position',[360 123 582 495]);
ax1 = axes();
imagesc(sm)
axis equal tight
colormap('gray')
axis off
ax1.Box = 'off';

hold on
plot([1,64],[1,64],'k--')

cb1 = colorbar;
cb1.Ticks = [];
cb1.Label.String = "Similarity";
caxis([0 1])
%xlabel('Choice')
%ylabel('Cue')
%rectangle('Position',[0.5,19.5,nBig,1],'EdgeColor','w')

xticklabels([]);
yticklabels([]);

%%
% Add colorbar
% h/t: https://blogs.mathworks.com/steve/2020/08/18/making-color-spectrum-plots-part-3/

ax2 = axes('Visible','off');
%ax2 = axes;
%im2 = imagesc(sm);
%im2.AlphaData = 1;
axis equal tight
axis off
ax2.Box = 'off';

cb2 = colorbar;

lambda = 1:nBig;

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

%%
% Add colorbar
% h/t: https://blogs.mathworks.com/steve/2020/08/18/making-color-spectrum-plots-part-3/

ax3 = axes('Visible','off');
%ax3 = axes;
%im3 = imagesc(sm);
%im3.AlphaData = 1;
axis equal tight
axis off
ax3.Box = 'off';

cb3 = colorbar;

lambda = 1:nBig;

ax = gca;
ax.Colormap = stimCols_sRGB(end:-1:1,:);
ax.CLim = [min(lambda) max(lambda)];

cb3.Location = 'westoutside';

cb3.Ticks = [];
%cb3.Label.String = "Cue";
cb3.Color = 'none';
cb3.Box = 'off';
%cb_position = cb3.Position
%cb3.Position = [cb_position(1),cb_position(2),cb_position(3)+3,cb_position(4)]
%cb3.Ruler.Color = 'k';
% cb3.TickDirection = "out";
ax.XTickLabels = [];
ax.XLabel = [];
ax.YTickLabels = [];
ax.YLabel = [];

%%

hlink = linkprop([ax1,ax2,ax3],{'Position','DataAspectRatio'});

cb1_Position = cb1.Position;
cb2_Position = cb2.Position;
cb3_Position = cb3.Position;

%ax_Position = ax.Position;

cb1.Position = [cb1_Position(1),cb1_Position(2),cb1_Position(3)/2,cb1_Position(4)-0.6];
cb2.Position = [cb2_Position(1:3),cb3_Position(2)-cb2_Position(2)];
cb3.Position = [cb3_Position(1:2),cb2_Position(1)-cb3_Position(1),cb3_Position(4)];

%%

%saveas(f2,'combined_TCC-FreeSimilarityMatrix_230509.svg')
saveas(f2,'combined_TCC-0att_fullremap-similaityMatrix_230510.svg')



