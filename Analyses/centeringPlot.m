clear, clc, close all

rng(42)

%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\combined_TCC-FreeSimilarityMatrix-workspace_230214.mat';
%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\211012_124119_Pollux\210422--211012_Pollux_TCC-FreeSimilarityMatrix-workspace_230222.mat';
filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\211108_090705_Castor\220517--211108_Castor_TCC-FreeSimilarityMatrix-workspace_230225';
%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\210609_124628_Buster\210428--210609_Buster_TCC-FreeSimilarityMatrix-workspace_230213.mat';
%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\220823_081207_Morty\220322--220823_Morty_TCC-FreeSimilarityMatrix-workspace_230213.mat';

load(filepath)

nBig = 64;

%%

sm = reshape(x, [64,64]);

% sm = zeros(nBig); % simple simulated cognitive bias
% sm(logical(eye(nBig))) = 1;
% sm(28:32,28:32) = triu(ones(5));
% sm(33:37,33:37) = tril(ones(5));

figure,
imagesc(sm)
axis square
colormap('gray')
colorbar
set(gca,'YDir','normal')

%%

%categoryCenters = [3,38]; %combined
%categoryCenters = [9,37]; %Pollux
categoryCenters = [4,18,37,50]; %Castor
%categoryCenters = [1,43]; %Buster
%categoryCenters = [37,63]; %Morty

% categoryCenters = [32]; %simple simulated cognitive bias


for categoryCenter = categoryCenters(2)

    sm_cs = circshift(sm,[nBig/2-categoryCenter,nBig/2-categoryCenter]);
    % note that that this doesn't put it full on "in the center", but in the
    % 32nd of 64 positions (there are 31 before and 32 after)

    f2 = figure('Position',[360 123 582 495]);
    ax1 = axes();
    imagesc(sm_cs)
    axis equal tight
    colormap('gray')
    cb = colorbar; 
    axis off
    ax1.Box = 'off';

    xticklabels(str2num(cell2mat(xticklabels)) - (nBig/2-categoryCenter))
    yticklabels(str2num(cell2mat(yticklabels)) - (nBig/2-categoryCenter))

    hold on
    xline(nBig/2,'r')
    yline(nBig/2,'r')
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

    stimCols = generateStimCols('nBig',64);
    stimCols_sRGB = LuvTosRGB([ones(1,64)*76.0693;stimCols]);
    stimCols_sRGB_shifted = circshift(stimCols_sRGB,nBig/2-categoryCenter);

    ax = gca;
    ax.Colormap = stimCols_sRGB_shifted;
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
    ax.Colormap = stimCols_sRGB_shifted(end:-1:1,:);
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

    hlink = linkprop([ax1,ax2,ax3],{'Position','DataAspectRatio'});

    cb1_Position = cb1.Position;
    cb2_Position = cb2.Position;
    cb3_Position = cb3.Position;

    %ax_Position = ax.Position;

    cb1.Position = [cb1_Position(1),cb1_Position(2),cb1_Position(3)/2,cb1_Position(4)-0.6];
    cb2.Position = [cb2_Position(1:3),cb3_Position(2)-cb2_Position(2)];
    cb3.Position = [cb3_Position(1:2),cb2_Position(1)-cb3_Position(1),cb3_Position(4)];

    %     figure,
    %     plot(sm_cs(nBig/2,:))
    %     xline(nBig/2)

    saveas(gcf,fullfile([num2str(categoryCenter),'_', datestr(now,'yymmdd'), '.svg']))



end

% %% Compute symmetry metric
%
% range = ceil(nBig/2)+1:nBig+floor(nBig/2)+1;
% t2 = NaN(nBig,nBig); % for debugging/understanding selection
% for n = 1:length(range)
%     t = false(nBig,nBig);
%     for i = 1:nBig
%         for j = 1:nBig
%             if i+j == range(n) && abs(i-j) < (floor(nBig/2)+1)
%                 t(i,j) = 1;
%                 t2(i,j) = 1;
%             end
%         end
%     end
%     halfLength = floor(sum(t(:))/2);
%     sm_t = sm_cs(t); % switch out for desired sm (centred on category)
%     sym(n) = mean(sm_t(1:halfLength)) - mean(sm_t(end:-1:end-halfLength+1));
%
%     % bootstrap
%     for boot = 1:1000
%         switch_ = logical(randi([0, 1], [1, halfLength]));
%         if mod(sum(t(:)),2)
%             switch_ = [switch_,false,flip(switch_)];
%         else
%             switch_ = [switch_,flip(switch_)];
%         end
%         sm_t_bs(switch_) = flip(sm_t(switch_)); % similarity matrix, temporary, bootstrap
%         sm_t_bs(~switch_) = sm_t(~switch_);
%         sym_bs(boot) = mean(sm_t_bs(1:halfLength)) - mean(sm_t(end:-1:end-halfLength+1)); %symmetry, bootstrap
%     end
%     lowerCI(n) = prctile(sym_bs,2.5);
%     upperCI(n) = prctile(sym_bs,97.5);
% end
%
% figure, hold on
% imagesc(sm_cs,'AlphaData',t2)
% ax = gca();
% ax.YDir = 'reverse';
% plot([0,nBig],[0,nBig],'r')
% axis square tight
% colormap('gray')
% colorbar
%
% figure, hold on
% plot(sym,'k','DisplayName','Symmetry')
% %yline(0,'HandleVisibility','off')
% axis tight
% ylabel('Symmetry (one side - other side)')
%
% plot(lowerCI,'DisplayName','lowerCI')
% plot(upperCI,'DisplayName','upperCI')
%
% legend
