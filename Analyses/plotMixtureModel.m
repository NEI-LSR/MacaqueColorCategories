function plotMixtureModel()

%% Colors

stimCols = generateStimCols('nBig',nBig,'sat',37);

rotVal = interval/2;
rotationMatrix = [cosd(rotVal), -sind(rotVal); sind(rotVal), cosd(rotVal)]; % h/t: https://www.mathworks.com/matlabcentral/answers/323483-how-to-rotate-points-on-2d-coordinate-systems#answer_253463
stimCols_rotated = rotationMatrix * stimCols;

if Lab %CIELAB
    stimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols_rotated]);
else % CIELUV
    stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols_rotated]);
end

colvals = im2double(stimCols_sRGB);


%% Figures: Bias by Cue with Category Crossings and Confidence Intervals

figure

% % create filename
if isfield(cleandata.trialdata,'dirname')
    if numel(dirname) == 1
        filename = dirname{1};
        file_dir = dirname{1};
    elseif numel(dirname) > 1
        dates = str2double(extractBefore(dirname, 7));
        [dates, date_order] = sort(dates);
        filename = [num2str(dates(1)), '--', num2str(dates(end)), extractAfter(dirname{1}, 13)];
        file_dir = dirname{date_order(end)};
    end
else
    filename = dirname;
    file_dir = dirname;
end

% % % display fit type
% axes0 = axes('Position',[0.01 0.95 0.15 0.05]);
% text('Parent', axes0, 'Interpreter', 'none', 'String', ['Fit Type: ' fit_type], 'FontSize', 15, 'Position', [0.05 0.5 0]);
% set(axes0,'XColor','none','YColor','none')
%
% % % display file name
% axes1 = axes('Position',[0.81 0.95 0.24 0.05]);
% text('Parent',axes1,'Interpreter','none','String',filename,'FontSize', 15, 'Position',[0.05 0.5 0]);
% set(axes1,'XColor','none','YColor','none')

axes2 = axes;
hold on
if isempty(ci) == 0
    plot(hue_angle, lower_95_w, 'k:');
    plot(hue_angle, upper_95_w,'k:');
    h = fill([hue_angle, fliplr(hue_angle)], [lower_95_w', fliplr(upper_95_w')], 'k');
    set(h, 'facealpha', .1, 'LineStyle', ':');
end
scatter(hue_angle,bias([1:end 1]),100,colvals([1:end 1],:),'filled');
plot(hue_angle,be_w,'k');
if isempty(interp_ci) == 0
    xline(interp_ci,'--');
end

if isempty(interp_crossing) == 0
    for i = 1:length(interp_crossing)
        xline(interp_crossing(i),'Color',crossing_colvals(i,:),'LineWidth',2.5)
    end
    scatter(interp_crossing,0,'filled','k');
end
yticks([-50 -25 0 25 50])
xticks([0 60 120 180 240 300 360]);
ax = gca;
set(gca,'TickDir','out');
ax.TickLength = [0.025 0.025];
ax.FontSize = 10;
xlim([0 360]);
ylim([-50 50]);
yline(0);
xlabel('Hue Angle ({\circ})');
ylabel('Bias ({\circ})');

saveas(gcf,fullfile('../',[filename,'_categorybias1_', datestr(now,'yymmdd-HHMMSS'), '.svg']))

%% Pie chart figure

figure
% figure('color', 'white')

hold on
pie(repelem(interval,nBig));% pie chart w/ nBig equally sized slices
colormap(rstimCols_sRGB);
ax = gca;
axis equal
delete(ax.Children([1, 1:2:nBig*2])) % stop displaying % for each slice
for i = 1:nBig
    ax.Children(i).EdgeAlpha = 0; % get rid of lines between slices
end
for i = 1:nBig
    ax.Children(i).FaceAlpha = 1;
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

rotated_colvals = im2double(rstimCols_sRGB);

shift_colvals = [colvals(nBig*(3/4):end,:); colvals(1:nBig*(3/4)-1,:)];
[cart,~] = generateStimCols('nBig',nBig); % generate values to plot cues

stimscatter = scatter(cart(1,:)./36,cart(2,:)./36,90,shift_colvals,'filled'); % plot all cues around pie chart

set(gca,'visible','off')
axis equal tight

rad_angle = deg2rad(hue_angle); %hue angles in radians
axes3 = axes;
p = polarplot(rad_angle, be_w+40,'b'); % dummy holder



rlim([0 80]);
hold on
polarplot(rad_angle, zeros(length(rad_angle),1)+40,'LineStyle',':','Color','k');
polarplot(rad_angle, be_w+40,'k');
thetaticks(0:45:315)

% if isempty(ci) == 0
%     polarplot(rad_angle, lower_95_w+40,':k');
%     polarplot(rad_angle, upper_95_w+40,':k');
% end

% add lines
% for k = 1:length(interp_crossing)
%     %     polarplot([deg2rad(interp_crossing(k)) deg2rad(interp_crossing(k))],[0 60],'Color',[rotated_colvals(crossings(k),:) (1+min(opacity))-opacity(k)],'LineWidth',1.5);
%     polarplot([deg2rad(interp_crossing(k)) deg2rad(interp_crossing(k))],[0 80],'Color',rotated_colvals(crossings(k),:),'LineWidth',1.5);
% end

ax = gca;
ax.Color = 'none';
% ax.ThetaTickLabel = {'0','','','90','','','180','','','270','',''};
ax.ThetaTickLabel = {};
ax.RTick = 0:20:80;
rticklabels({'','-20{\circ}','0{\circ}','+20{\circ}',''});


% ax.RAxisLocation = 235;
% ax.FontSize = 8;
% ax.Units = 'normalized';
%ax.Position = [0.3734375,0.212952799121844,0.2890625,0.609220636663008];

theta = rad_angle;

ax_polar = gca;
ax_cart = axes();
axis equal
ax_cart.Position = ax_polar.Position;

if isempty(ci) == 0
    rlow = lower_95_w' + 40;
    rhigh = upper_95_w' + 40;

    [x1,y1] = pol2cart(theta,rlow);
    [x2,y2] = pol2cart(theta,rhigh);
    patch([x1 fliplr(x2)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeAlpha', 0);
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

saveas(gcf,fullfile('../',[filename,'_categorybias2_', datestr(now,'yymmdd-HHMMSS'), '.svg']))

%% Figures: Bias by Cue with Category Crossings and Confidence Intervals

figure

% % % display fit type
% axes0 = axes('Position',[0.01 0.95 0.15 0.05]);
% text('Parent', axes0, 'Interpreter', 'none', 'String', ['Fit Type: ' fit_type], 'FontSize', 15, 'Position', [0.05 0.5 0]);
% set(axes0,'XColor','none','YColor','none');
%
% % % display file name
% axes1 = axes('Position',[0.81 0.95 0.24 0.05]);
% text('Parent',axes1,'Interpreter','none','String',filename,'FontSize', 15, 'Position',[0.05 0.5 0]);
% set(axes1,'XColor','none','YColor','none');

axes2 = axes('Position', [0.13 0.13 0.8 0.78]);
hold on

hue_angle = 0:interval:360;
hue_angle = hue_angle - 180;

for i = 1:length(interp_crossing)
    if interp_crossing(i) > 180
        interp_crossing(i) = interp_crossing(i) - 360;
    end
end

if isempty(ci) == 0
    for i = 1:length(interp_ci)
        if interp_ci(i) > 180
            interp_ci(i) = interp_ci(i) - 360;
        end
    end

    %      unique_interpci = unique(interp_ci,'stable');
    %      %unique_interpci = interp_ci
    %      x = 0;
    %      for i = 2:2:length(unique_interpci)
    %          x = x+1;
    %          change(x,1) = abs(unique_interpci(i) - unique_interpci(i-1));
    %      end
    %
    %      opacity = change ./min(change) *0.15;
    %      opacity = 1 - (1+min(opacity) - opacity);
    %      %opacity = flip(opacity);
    %      opacity = 1 - opacity;

    x = 0;
    for i = 2:2:length(interp_ci)
        x = x+1;
        fill([interp_ci(i-1) interp_ci(i) interp_ci(i) interp_ci(i-1)],...
            [-50 -50 50 50],crossing_colvals(x,:),'EdgeColor','none','FaceAlpha',0.75);%opacity(x))
    end

    lower_95_w_rotated = lower_95_w([nBig/2+1:nBig 1:nBig/2]);
    upper_95_w_rotated = upper_95_w([nBig/2+1:nBig 1:nBig/2]);

    plot(hue_angle, lower_95_w_rotated([1:end 1]), 'k:');
    plot(hue_angle, upper_95_w_rotated([1:end 1]), 'k:');
    h = fill([hue_angle, fliplr(hue_angle)], [lower_95_w_rotated([1:end 1])', fliplr(upper_95_w_rotated([1:end 1])')], 'k');
    set(h, 'facealpha', .1, 'LineStyle', ':');
end

moving_bias_rotated = moving_bias([nBig/2+1:nBig 1:nBig/2]);
rotated_bias = bias([nBig/2+1:nBig 1:nBig/2]);
rotated_colvals = colvals([nBig/2+1:nBig 1:nBig/2],:);
scatter(hue_angle,rotated_bias([1:end 1]),100,rotated_colvals([1:end 1],:),'filled');
plot(hue_angle,moving_bias_rotated([1:end 1]),'k');
if isempty(interp_ci) == 0
    xline(interp_ci,'--');
end
if isempty(interp_crossing) == 0
    scatter(interp_crossing,0,'filled','k');
end
yticks([-50 -25 0 25 50])
xticks([-180 -120 -60 0 60 120 180]);
hold off
ax = gca;
set(gca,'TickDir','out');
ax.TickLength = [0.025 0.025];
ax.FontSize = 10;
xlim([-180 180]);
ylim([-50 50]);
yline(0);
xlabel('Hue Angle ({\circ})');
ylabel('Bias ({\circ})');

saveas(gcf,fullfile('../',[filename,'_categorybias3_', datestr(now,'yymmdd-HHMMSS'), '.svg']))