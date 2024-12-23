function plotSimilarityMatrix(x,filename,OutPutFileDir,categoryCenter,withLabels,Lab)

%clear, clc, close all

rotateMatrix = false;

if ~exist('withLabels','var')
    withLabels = true;
end

if min(size(x)) == 1 % if we pass a vector rather than a matrix, assume it needs reshaping
    nBig = sqrt(length(x));
    sm = reshape(x,[nBig, nBig]);
else
    nBig = length(x);
    sm = x; % Similarity Matrix
end

if rotateMatrix
    if exist('categoryCenter','var') && ~isempty(categoryCenter)
        sm = circshift(sm,[nBig/2-categoryCenter,nBig/2-categoryCenter]);
    end
end

%%

if exist('Lab','var')
    stimCols = generateStimCols('nBig',nBig,'sat',Lab(2));
    stimCols_sRGB = LabTosRGB([repelem(Lab(1), nBig); stimCols]);
else
    stimCols = generateStimCols('nBig',nBig);
    stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols]);
end

if rotateMatrix
    if exist('categoryCenter','var')
        stimCols = circshift(stimCols,[0,nBig/2-categoryCenter]);
    end
end

%%
f = figure;
ax1 = axes();
imagesc(sm')
axis equal tight
colormap('gray')
axis off
ax1.Box = 'off';

if withLabels
    xlabel('Cue')
    ylabel('Choice')
end

hold on
plot([1,nBig],[1,nBig],'k--')

if exist('categoryCenter','var') && ~isempty(categoryCenter) % Add center lines
    if rotateMatrix
        if exist('Lab','var')
            xline(nBig/2,'Color', stimCols_sRGB(nBig/2,:),'LineWidth',3)
            yline(nBig/2,'Color', stimCols_sRGB(nBig/2,:),'LineWidth',3)
        else
            xline(nBig/2,'Color', stimCols_sRGB(nBig/2,:),'LineWidth',3)
            yline(nBig/2,'Color', stimCols_sRGB(nBig/2,:),'LineWidth',3)
        end
    else
        for i = 1:length(categoryCenter)
            xline(categoryCenter(i),...
                'Color', stimCols_sRGB(categoryCenter(i),:),...
                'LineWidth',3)
        end
    end
end

if false
    xline(81,...
        'Color', stimCols_sRGB(81,:),...
        'LineWidth',3,'LineStyle','--')
end

cb1 = colorbar;
cb1.Ticks = [];
if withLabels
    cb1.Label.String = "Similarity";
end
caxis([0 1])

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

linkprop([ax1,ax2,ax3],{'Position','DataAspectRatio'});

ax1_position = ax1.Position;
set(ax1, 'Position', [0.175,0.15,ax1_position(3:4)]) 

cb1_Position = cb1.Position;
cb2_Position = cb2.Position;
cb3_Position = cb3.Position;

%ax_Position = ax.Position;

cb1.Position = [cb1_Position(1),cb1_Position(2),cb1_Position(3)/2,cb1_Position(4)-0.6];
cb2.Position = [cb2_Position(1:3),cb3_Position(2)-cb2_Position(2)];
cb3.Position = [cb3_Position(1:2),cb2_Position(1)-cb3_Position(1),cb3_Position(4)];

%%

%saveas(f2,'combined_TCC-FreeSimilarityMatrix_230509.svg')
%saveas(f2,'combined_TCC-0att_fullremap-similaityMatrix_230510.svg')

if exist('filename','var')
    if ~exist('OutPutFileDir','var')
        OutPutFileDir = './'; % if `OutPutFileDir` not passed, save here
    end
    saveas(f,[OutPutFileDir,'sm_',filename,'_',datestr(now,'yymmdd-HHMMSS'),'.svg'])
end




