%% What if it is *just* bias?

clear, clc, close all

rng(0)

nBig = 64;
interval = 360/nBig;
offset = sin(deg2rad(0:interval:360-interval))*30;

%bias_vector_index = round(offset/interval); % vector of bias values, in units of "interval" away from direct match
% at some point we'll want to change this into degrees or radians, rather
% than these massive units, but for now this should suffice as a proof of
% concept

% % for a default similarity matrix:
%bias_vector_index = zeros(size(bias_vector_index));

load moving_bias.mat
bias_vector_index = round(moving_bias/interval)';

%% Grab default similarity matrix

[~, tempdata] = GenerativeModel([],...
    'nBig',nBig);
sm = tempdata.trialdata.similarityMatrix;
stimCols = tempdata.trialdata.stimCols;

%% Offset similarity matrix

for i = 1:nBig
    a_sm(i,:) = circshift(sm(i,:),bias_vector_index(i));
end

f2 = figure('Position',[360 123 582 495]);
ax1 = axes();
imagesc(a_sm)
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

%% Pull out and plot a single row

subset = a_sm(20,:);

ap = 501; % approximation, higher the value - the less of an approximation this is (but the slower it will run)
%     This assumes that the points are equally distributed across the bins
%     There *are* more complex ways of doing this if we're happy to assume a certain distribution of responses within the bin, but this probably suffices for our purpose
%

% points = -180+interval:interval:180;
points = 0:interval:360-interval;

bin_starts = points - interval/2;
bin_ends = points + interval/2;

t = [];
for j = 1:nBig
    t = [t,deg2rad(linspace(bin_starts(j),bin_ends(j),round(subset(j)*ap)))];
end
bias = rad2deg(circ_median(t,2));

f3 = figure('Position',[360 123 582 495]);
ax4 = axes;
hold on

xline(20,'k:','LineWidth',1,'DisplayName','Cue / Direct Match')
plot(1:nBig,subset,'k','DisplayName','Similarity Function');
[~,minloc] = min(abs(points-bias));
xline(minloc,'k--','LineWidth',1,'DisplayName','Function Median')
xlim([1,nBig])
axis square tight
xticks([])
yticks([0,1])

legend boxoff

ylabel('Similarity')

%% Add colorbar
% h/t: https://blogs.mathworks.com/steve/2020/08/18/making-color-spectrum-plots-part-3/

cb = colorbar;

lambda = 1:nBig;

colors = LuvTosRGB([ones(1,64)*76.0693;stimCols]);
ax = gca;
ax.Colormap = colors;
ax.CLim = [min(lambda) max(lambda)];

cb.Location = 'southoutside';

cb.Ticks = ax.XTick;
cb.Label.String = "Choice";
% cb.TickDirection = "out";
% ax.XTickLabels = [];
% ax.XLabel = [];

%%

hlink = linkprop([ax1,ax2,ax3,ax4],{'Position'});

%%

saveas(f2,'justBias.svg')
saveas(f3,'justBias_subset.svg')


%%

[~, tempdata] = GenerativeModel([],...
    'nBig',nBig,'similarityMatrix',a_sm,'nTrials',200000);
moving_bias = fitMixtureModel(tempdata);

