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
cb = colorbar;
cb.Ticks = [0,1];
cb.Label.String = "Similarity";
caxis([0 1])
%xlabel('Choice')
%ylabel('Cue')
%rectangle('Position',[0.5,19.5,nBig,1],'EdgeColor','w')

xticklabels([]);
yticklabels([]);

% Add colorbar
% h/t: https://blogs.mathworks.com/steve/2020/08/18/making-color-spectrum-plots-part-3/

ax2 = axes('Visible','off');
%ax2 = axes;
%im2 = imagesc(sm);
%im2.AlphaData = 1;
axis image tight
%ax2 = axes;

cb = colorbar;

lambda = 1:nBig;

stimCols_sRGB = LuvTosRGB([ones(1,64)*76.0693;stimCols]);
ax = gca;
ax.Colormap = stimCols_sRGB;
ax.CLim = [min(lambda) max(lambda)];

cb.Location = 'southoutside';

cb.Ticks = [];
cb.Label.String = "Choice";
% cb.TickDirection = "out";
ax.XTickLabels = [];
ax.XLabel = [];
ax.YTickLabels = [];
ax.YLabel = [];

% Add colorbar
% h/t: https://blogs.mathworks.com/steve/2020/08/18/making-color-spectrum-plots-part-3/

ax3 = axes('Visible','off');
%ax3 = axes;
%im3 = imagesc(sm);
%im3.AlphaData = 1;
axis image tight
%ax3 = axes;

cb = colorbar;

lambda = 1:nBig;

ax = gca;
ax.Colormap = stimCols_sRGB(end:-1:1,:);
ax.CLim = [min(lambda) max(lambda)];

cb.Location = 'westoutside';

cb.Ticks = [];
cb.Label.String = "Cue";
% cb.TickDirection = "out";
ax.XTickLabels = [];
ax.XLabel = [];
ax.YTickLabels = [];
ax.YLabel = [];

hlink = linkprop([ax1,ax2,ax3],{'Position','DataAspectRatio'});

% saveas(f2,'combined_TCC-FreeSimilarityMatrix_230509.svg')
saveas(f2,'combined_TCC-0att_fullremap-similaityMatrix_230510.svg')



