function [moving_bias, lower_95, upper_95] = fitMixtureModel(cleandata,Lab,lengthOfSlidingWindow)


if isfield(cleandata.trialdata,'dirname') % for real data
    dirname = unique(cleandata.trialdata.dirname);
    paradigm = unique(cleandata.trialdata.paradigm);
    if length(paradigm) > 1
        warning('More than one paradigm in data file');
    end
    try
        nBig = size(cleandata.trialdata.allchoices{1,1},2);
    catch
        nBig = size(cleandata.trialdata.stimCols{1,1},1);
    end
else
    dirname = 'simdata'; % for simdata
    nBig = size(cleandata.trialdata.stimCols,2);
end

nSmall = sum(~isnan(cleandata.trialdata.choices{end,1}));

interval = 360/nBig;

% Colors
stimCols = generateStimCols('nBig',nBig,'sat',37);
rotVal = 360/(nBig*2);
rotationMatrix = [cosd(rotVal),-sind(rotVal);sind(rotVal),cosd(rotVal)];
stimCols_biased = rotationMatrix * stimCols;
if Lab %CIELAB
    stimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols_biased]);
else % CIELUV
    stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols_biased]);
end
colvals = im2double(stimCols_sRGB);

%% Calculate bias 

cues = cell2mat(cleandata.trialdata.cues);
choices = cell2mat(cleandata.trialdata.choices);
chosen = cell2mat(cleandata.trialdata.chosen);

try
    colresp = [cues, chosen, choices];
catch
    colresp = [cues, chosen', choices];
end

% Get completed trials + incorrect trials
colresp(any(isnan(colresp),2),:) = []; % gets rid of all aborts
incorrect = colresp(colresp(:,1) ~= colresp(:,2),1:2); % cue + response on incorrect trials


% Calculate Angular Error

% Distance values
distances = (-180+interval:interval:180)';

% Calculate angular error (distance) between incorrect choice and cue
distance = zeros(size(incorrect(:,1)));
for i = 1:length(incorrect) % for each trial
    if abs(incorrect(i,2) - incorrect(i,1)) < nBig/2
        distance(i,1) = (incorrect(i,2) - incorrect(i,1)) * interval;
    elseif abs(incorrect(i,2) - incorrect(i,1)) > nBig/2 && incorrect(i,2) > incorrect(i,1)
        distance(i,1) = (-(nBig - abs(incorrect(i,2) - incorrect(i,1)))) * interval;
    else
        distance(i,1) = (nBig - abs(incorrect(i,2) - incorrect(i,1))) * interval;
    end
end

% Cue + angular error for each incorrect trial
cue_dist = [incorrect(:,1), distance];

% Count of number of times each color (by distance from cue) was chosen
distance_counts = zeros(length(distances),nBig);
for i = 1:nBig
    dist_count = cue_dist(cue_dist(:,1)==i, :); % All angular error values for cue i
    distance_counts(:,i) = histc(dist_count(:,2),distances);
end

% Count of number of times each color (by distance from cue) was an option
% (completed trials only)
for i = 1:nBig
    comp_trial = colresp(colresp(:,1)==i,3:end); % choices for completed trials
    for j = 1:nSmall %
        for k = 1:size(comp_trial,1)
            if abs(comp_trial(k,j) - i) < nBig/2
                comp_trial(k,j) = (comp_trial(k,j) - i)*interval;
            elseif abs(comp_trial(k,j) - i) > nBig/2 && comp_trial(k,j) > i
                comp_trial(k,j) = (-(nBig - abs(comp_trial(k,j) - i)))*interval;
            else
                comp_trial(k,j) = (nBig - abs(comp_trial(k,j) - i))*interval;
            end
        end
    end
    for m = 1:length(distances)
        compchoice_distances(m,i)= sum(comp_trial(:)==distances(m,1));
    end
end

% Normalize counts for each color by number of times color was an option
dist_prop = distance_counts./compchoice_distances;

% Replace value at 0 (correct choice) to exclude from curve fit
dist_prop(nBig/2,:) = NaN;
compchoice_distances(nBig/2,:) = NaN; % Replace number at distance = 0

%Fit Gaussian to error distribution for each cue to get bias value
for i = 1:nBig
    gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d';
    startingPoints = [0.5 0 78 0];
    dist_idx = find(~isnan(dist_prop(:,i))); % Index of distance counts for only colors shown
    dists = distances(dist_idx,:); % Excludes distance values for colors never shown
    weights = compchoice_distances(dist_idx,i); % Weights fit by number of times each color was an option
    f = fit(dists,dist_prop(~isnan(dist_prop(:,i)),i),gaussEqn,'start',startingPoints, 'Weights',weights,'Lower',[0 -180 0 0],'Upper',[Inf 180 Inf 1]);
    ci = confint(f,0.95);
    bias(i,1) = f.b;
    ci_lower_95(i,1) = ci(1,2);
    ci_upper_95(i,1) = ci(2,2);
    if any(isnan(ci),'all')
        disp(i)
        disp(f)
        disp(ci)
        warning('NaN in CI')
    end
end

%% Category center locations

if ~exist('lengthOfSlidingWindow','var')
    lengthOfSlidingWindow = 3; % moving average input
end

if mod(lengthOfSlidingWindow,2) == 0
    error('lengthOfSlidingWindow should be odd')
end

moving_bias = movmean(padarray(bias,lengthOfSlidingWindow,"circular"),lengthOfSlidingWindow,'Endpoints','discard');
moving_bias = moving_bias(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

be_w = moving_bias([1:end,1]); % bias estimates including wraparound

for i = 1:nBig
    crossesZero(i) = and(be_w(i)>=0, be_w(i+1)<=0);
end

crossings = find(crossesZero);

hue_angle = 0:interval:360; %includes wraparound

% Category Center Location Interpolation
if ~isempty(crossings)

    for i = 1:length(crossings)
        x = [hue_angle(crossings(i)) hue_angle(crossings(i)+1)];
        y = [be_w(crossings(i)) be_w(crossings(i)+1)];
        interp_crossing(i,1) = interp1(y,x,0);
    end

    % Category Center Colors
    polarAngs = interp_crossing'; %Polar Angles
    [a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
    crossingCols = [a;b];
    if Lab % CIELAB
        crossingCols_sRGB = LabTosRGB([repelem(76.0693, length(polarAngs)); crossingCols]);
    else % CIELUV
        crossingCols_sRGB = LuvTosRGB([repelem(76.0693, length(polarAngs)); crossingCols]);
    end
    crossing_colvals = im2double(crossingCols_sRGB);
else
    interp_crossing = [];
end

%% Confidence intervals of category centers
if isempty(ci) == 0

    lower_95 = movmean(padarray(ci_lower_95,lengthOfSlidingWindow,"circular"),lengthOfSlidingWindow,'Endpoints','discard');
    lower_95 = lower_95(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

    upper_95 = movmean(padarray(ci_upper_95,lengthOfSlidingWindow,"circular"),lengthOfSlidingWindow,'Endpoints','discard');
    upper_95 = upper_95(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

    for i = 1:nBig % check whether the confidence interval for each cue includes 0
        withinCI(i) = and(lower_95(i,1)<=0, upper_95(i,1)>=0);
    end

    with_CI = withinCI([1:end 1]);


    for i = 1:nBig
        changes(i) = or(and(with_CI(i)==0, with_CI(i+1)==1), and(with_CI(i)==1, with_CI(i+1)==0));
    end

    change_range = find(changes);

    upper_95_w = upper_95([1:end 1]);
    lower_95_w = lower_95([1:end 1]);

    if isempty(change_range) == 0 % if confidence interval spans zero for all values

        change_range = change_range([1:end 1]);

        for i = 1:length(crossings)
            for j = 1:sum(changes)
                if change_range(j) < change_range(j+1)
                    if (crossings(i) >= change_range(j) && crossings(i) <= change_range(j+1)) == 1
                        CI_range(i,:) = [change_range(j) change_range(j+1)];
                    end
                elseif change_range(j) > change_range(j+1)
                    if (crossings(i) >= change_range(j) && crossings(i) >= change_range(j+1)) == 1 ...
                            || (crossings(i) <= change_range(j) && crossings(i) <= change_range(j+1) == 1)
                        CI_range(i,:) = [change_range(j) change_range(j+1)];
                    end
                end
            end
        end

        CI_range = CI_range';
        CI_range = CI_range(:);

        % Confidence Interval Interpolation
        for i = 1:length(CI_range)
            if be_w(CI_range(i)) > 0
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
                y = [lower_95_w(CI_range(i)) lower_95_w(CI_range(i)+1)];
            elseif  be_w(CI_range(i)) < 0
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
                y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
            end

            interp_ci(i,1) = interp1(y,x,0);

            if isnan(interp_ci(i,1)) && be_w(CI_range(i)) > 0
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
                y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
                interp_ci(i,1) = interp1(y,x,0);
            elseif isnan(interp_ci(i,1)) && be_w(CI_range(i)) < 0
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
                y = [lower_95_w(CI_range(i)) lower_95_w(CI_range(i)+1)];
                interp_ci(i,1) = interp1(y,x,0);
            end
        end
%     else 
%         %interp_ci = []; %needed?
    end
else
    interp_ci = [];
end

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

end