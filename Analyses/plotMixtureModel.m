function plotMixtureModel(model, whichFigures, filename, withLabels)

if ~exist('withLabels','var')
    withLabels = true;
end

axlims = 40;

PotentialDistances  = model.PotentialDistances;
interp_crossing     = model.interp_crossing;
interp_ci           = model.interp_ci;
presentation_counts = model.presentation_counts;
choice_probability  = model.choice_probability;
upper_95_w          = model.upper_95_w;
lower_95_w          = model.lower_95_w;
bias                = model.bias;
be_w                = model.be_w;
ci                  = model.ci;
moving_bias         = model.moving_bias;

% request *all* figures
if isfield(whichFigures,'all') && whichFigures.all == true
    whichFigures.mixMod_BreakOut = true;
    whichFigures.MixMod_linear = true;
    whichFigures.MixMod_polar = true;
end

%%
nBig = size(model.gaussfits,2);
interval = 360/nBig;

%% Colors

hue_angle = 0:interval:360; %includes wraparound
stimCols = generateStimCols('nBig',nBig,'sat',37);

rotVal = interval/2;
rotationMatrix = [cosd(rotVal), -sind(rotVal); sind(rotVal), cosd(rotVal)]; % h/t: https://www.mathworks.com/matlabcentral/answers/323483-how-to-rotate-points-on-2d-coordinate-systems#answer_253463
stimCols_rotated = rotationMatrix * stimCols;

if exist('Lab','var') %CIELAB
    stimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols_rotated]);
else % CIELUV
    stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols_rotated]);
end

colvals = im2double(stimCols_sRGB);

for i = 1:length(interp_crossing)
    [~,closestHueAngleToCrossings(i)] = min(abs(interp_crossing(i) - hue_angle));
end
interp_crossing_colvals = colvals(closestHueAngleToCrossings,:);

%% Break out figures

if isfield(whichFigures,'mixMod_BreakOut') && whichFigures.mixMod_BreakOut == true
    for cueIndex = 1:nBig

        figure('visible','off')
        hold on
        axis tight
        ylim([0,1])

        pltCol = 1 - repmat(...
            (presentation_counts(:,cueIndex)-min(presentation_counts(:,cueIndex)))...
            /(max(presentation_counts(:,cueIndex))-min(presentation_counts(:,cueIndex))),...
            1,3);

        s = scatter(PotentialDistances, choice_probability(:,cueIndex),'filled');
        s.CData = pltCol;

        % plot(f,PotentialDistances,choice_probability(:,cueIndex),'k.')
        plot(model.gaussfits{cueIndex},...
            PotentialDistances(~isnan(choice_probability(:,cueIndex))),...
            choice_probability(~isnan(choice_probability(:,cueIndex)),cueIndex),'k.');

        p = gca;
        p.Children(2).Marker = 'none'; % turn off data, so that we can replot it how we like...
        % p.Children(1).LineWidth = 3;

        p11 = predint(model.gaussfits{cueIndex},PotentialDistances,0.95,'functional','off');                         %TODO: Check whether this is the appropriate type of interval: https://www.mathworks.com/help/curvefit/confidence-and-prediction-bounds.html
        plot(PotentialDistances,p11,'k:',...
            'DisplayName','Nonsimultaneous Functional Bounds')
        % p.Children(1).LineWidth = 1;
        % p.Children(2).LineWidth = 1;

        xline(0,'k--')

        xlabel('Error distance')
        ylabel('Choice Probability')
        text(-90,0.9,num2str(cueIndex))
        legend('off')
        drawnow

        saveas(gcf,fullfile([filename,num2str(cueIndex),'_mixMod_BreakOut.svg']))
        close all
    end

end

%% Figures: Bias by Cue with Category Crossings and Confidence Intervals

if isfield(whichFigures,'MixMod_linear') && whichFigures.MixMod_linear == true

    figure

    % % % display fit type
    % axes0 = axes('Position',[0.01 0.95 0.15 0.05]);
    % text('Parent', axes0, 'Interpreter', 'none', 'String', ['Fit Type: ' fit_type], 'FontSize', 15, 'Position', [0.05 0.5 0]);
    % set(axes0,'XColor','none','YColor','none')
    %
    % % % display file name
    % axes1 = axes('Position',[0.81 0.95 0.24 0.05]);
    % text('Parent',axes1,'Interpreter','none','String',filename,'FontSize', 15, 'Position',[0.05 0.5 0]);
    % set(axes1,'XColor','none','YColor','none')

    axes('PositionConstraint','innerposition',...
        'Position',[0.13 0.19 0.82 0.75])

    hold on

    x = 0;
    for i = 2:2:length(interp_ci)
        x = x+1;
        if interp_ci(i) > interp_ci(i-1)
            fill([interp_ci(i-1) interp_ci(i) interp_ci(i) interp_ci(i-1)],...
                [-50 -50 50 50],interp_crossing_colvals(x,:),'EdgeColor','none','FaceAlpha',0.5);%opacity(x))
        elseif interp_ci(i) < interp_ci(i-1)
            fill([interp_ci(i-1) 360 360 interp_ci(i-1)],...
                [-50 -50 50 50],interp_crossing_colvals(x,:),'EdgeColor','none','FaceAlpha',0.5);
            fill([0 interp_ci(i) interp_ci(i) 0],...
                [-50 -50 50 50],interp_crossing_colvals(x,:),'EdgeColor','none','FaceAlpha',0.5);
        end
    end

    h = fill([hue_angle, fliplr(hue_angle)], [lower_95_w', fliplr(upper_95_w')], 'k');
    set(h, 'facealpha', .1, 'LineStyle', 'none');

    scatter(hue_angle,bias([1:end 1]),100,colvals([1:end 1],:),'filled');
    plot(hue_angle,be_w,'k');
    % xline(interp_ci,'--');

    for i = 1:length(interp_crossing)
        plot([interp_crossing(i), interp_crossing(i)], ylim(), 'Color', interp_crossing_colvals(i,:));
    end
    plot(interp_crossing,0,'ko','MarkerFaceColor','k');

    grid on
    yticks([-axlims,-20:20:20,axlims])
    yticklabels({'-40','-20','0','+20','+40',''})
    xticks(0:45:360);
    ax = gca;
    set(gca,'TickDir','out');
    % ax.TickLength = [0.025 0.025];
    % ax.FontSize = 10;
    xlim([0 360]);
    ylim([-axlims axlims]);
    yline(0,'LineStyle','--','Color',[0.3,0.3,0.3]);
    xlabel('Hue Angle');
    ylabel('Bias');

    saveas(gcf,fullfile('../',[filename,'_MixMod_linear_', datestr(now,'yymmdd-HHMMSS'), '.svg']))

end

%% Pie chart figure

if isfield(whichFigures,'MixMod_polar') && whichFigures.MixMod_polar == true

    axisoffset = axlims;

    axPositions = [0.02,0.02,0.96,0.96];

    figure('Position',[360, 97+(2/3), 420, 420])
    % figure('color', 'white')

    hold on
    pie(repelem(interval,nBig));% pie chart w/ nBig equally sized slices
    colormap(rstimCols_sRGB);
    ax = gca;
    ax.Position = axPositions;
    axis equal
    delete(ax.Children([1, 1:2:nBig*2])) % stop displaying % for each slice
    for i = 1:nBig
        ax.Children(i).EdgeAlpha = 0; % get rid of lines between slices
    end
    for i = 1:nBig
        ax.Children(i).FaceAlpha = 0.5;
    end
    ax.View = [90 90];

    % Mark out anything outside of confidence intervals
    if ~isempty(ci) && ~isempty(interp_ci) % && isempty(change_range) == 0

        x0=0;
        y0=0;

        CI_range = unique(interp_ci','stable');

        CI_range = CI_range([1:end 1]);
        for i = 2:2:length(CI_range)
            if i < length(CI_range)
                pie_direction = (.25*CI_range(1,i))/90 + 0.25; % starting value
                if CI_range(1,i) > CI_range(1,i+1)
                    theta = deg2rad((360 - CI_range(1,i)) + CI_range(1,i+1));
                else
                    theta = deg2rad(CI_range(1,i+1) - CI_range(1,i));
                end
                a1 = 2*pi*pie_direction; % Starting direction
                r = 1; % radius
                a2 = a1 + theta; % Ending direction
                t = linspace(a1,a2);
                x = x0 + r*cos(t);
                y = y0 + r*sin(t);
                fill([x0,x,x0],[y0,y,y0],'w','FaceAlpha',1,'EdgeAlpha',0); %0.7
            end
        end

        % % Transparency for confidence intervals

        %     unique_interpci = unique(interp_ci,'stable');
        %     %unique_interpci = interp_ci;
        %
        %     x = 0;
        %     for i = 2:2:length(unique_interpci)
        %         x = x+1;
        %         change(x,1) = abs(unique_interpci(i) - unique_interpci(i-1));
        %     end
        %
        %     opacity = (change ./min(change)) *0.15;
        %     opacity_number = 0;
        %     %CI_range = unique(CI_range,'stable');
        %
        %     unique_opacity = unique(opacity,'stable');

        % for i = 1:2:length(CI_range)-1
        %     %opacity_number = opacity_number+1;
        %     pie_direction = (.25*CI_range(1,i))/90 + 0.25;
        %     if CI_range(1,i+1) > CI_range(1,i)
        %         theta = deg2rad(CI_range(1,i+1) - CI_range(1,i));
        %     elseif CI_range(1,i+1) < CI_range(1,i)
        %         theta = deg2rad((360 - CI_range(1,i)) + CI_range(1,i+1));
        %     end
        %     a1 = 2*pi*pie_direction; % Starting direction
        %     a2 = a1 + theta; % Ending direction
        %     t = linspace(a1,a2);
        %     x = x0 + r*cos(t);
        %     y = y0 + r*sin(t);
        %     %fill([x0,x,x0],[y0,y,y0],'w','FaceAlpha',unique_opacity(opacity_number),'EdgeAlpha',0);
        %     fill([x0,x,x0],[y0,y,y0],'w','FaceAlpha',0.3,'EdgeAlpha',0);
        % end
    else
        for i = 1:nBig
            ax.Children(i).FaceAlpha = 0;
        end
    end

    ax.Color = 'none';

    shift_colvals = [colvals(nBig*(3/4):end,:); colvals(1:nBig*(3/4)-1,:)];
    [cart,~] = generateStimCols('nBig',nBig); % generate values to plot cues

    scatter(cart(1,:)./36,cart(2,:)./36,150,shift_colvals,'filled'); % plot all cues around pie chart

    set(gca,'visible','off')
    axis equal tight

    rad_angle = deg2rad(hue_angle); %hue angles in radians
    axes3 = axes('Position',[0.02,0.02,0.96,0.96]);
    polarplot(rad_angle, be_w + axisoffset,'b'); % dummy holder

    rlim([0 axlims*2]);
    hold on
    polarplot(rad_angle, zeros(length(rad_angle),1) + axisoffset,'LineStyle','--','Color',[0.3,0.3,0.3]);
    polarplot(rad_angle, be_w + axisoffset,'k');
    thetaticks(0:45:360)

    % if isempty(ci) == 0
    %     polarplot(rad_angle, lower_95_w + axisoffset,':k');
    %     polarplot(rad_angle, upper_95_w + axisoffset,':k');
    % end

    % add lines
    for k = 1:length(interp_crossing)
        polarplot([deg2rad(interp_crossing(k)) deg2rad(interp_crossing(k))],[0 axlims*2],'Color',rstimCols_sRGB(round(interp_crossing(k)/interval),:));
        polarplot(deg2rad(interp_crossing(k)), 0 + axisoffset, 'ko','MarkerFaceColor','k')
    end

    plotDemoLine = 0;
    if plotDemoLine
        polarplot([deg2rad(180) deg2rad(180)],[0 axlims*2],...
            'Color',rstimCols_sRGB(round(180/interval),:),...
            'LineStyle','--','LineWidth',2);
    end

    ax = gca;
    ax.Color = 'none';
    % ax.ThetaTickLabel = {'0','','','90','','','180','','','270','',''};
    ax.ThetaTickLabel = {};
    ax.RTick = [0:20:axlims*2];
    if ~withLabels
        rticklabels({'','','','',''});
    else
        rticklabels({'','-20','0','+20',''}); % Could use `axisoffset`?
    end

    ax.RAxisLocation = 45;
    % ax.FontSize = 8;
    % ax.Units = 'normalized';
    %ax.Position = [0.3734375,0.212952799121844,0.2890625,0.609220636663008];

    theta = rad_angle;

    ax_polar = gca;
    ax_cart = axes();
    axis equal
    ax_cart.Position = ax_polar.Position;

    if isempty(ci) == 0
        rlow = lower_95_w' + axisoffset;
        rhigh = upper_95_w' + axisoffset;

        [x1,y1] = pol2cart(theta,rlow);
        [x2,y2] = pol2cart(theta,rhigh);
        patch([x1 fliplr(x2)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
    end

    xlim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]);
    ylim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]);
    %axis equal;
    set(ax_cart,'visible','off');

    % % display fit type
    % axes0 = axes('Position',[0.01 0.95 0.15 0.05]);
    % text('Parent', axes0, 'Interpreter', 'none', 'String', ['Fit Type: ' fit_type], 'FontSize', 15, 'Position', [0.05 0.5 0]);
    % set(axes0,'XColor','none','YColor','none');

    % % display file name
    % axes1 = axes('Position',[0.81 0.95 0.24 0.05]);
    % text('Parent',axes1,'Interpreter','none','String',filename,'FontSize', 15, 'Position',[0.05 0.5 0]);
    % set(axes1,'XColor','none','YColor','none');

    saveas(gcf,fullfile('../',[filename,'_MixMod_polar_', datestr(now,'yymmdd-HHMMSS'), '.svg']))

end


end