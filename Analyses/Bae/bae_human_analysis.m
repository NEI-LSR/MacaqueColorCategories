% Bae Data Analysis 

bigN = 64;
degrees = 360/bigN;

hue_angle = 0:degrees:360;
hue_angle = hue_angle(1:bigN);

colresp = [cue_angle, choice_angle];

all_cues = length(unique(target_color));
interval = 360/all_cues;
cues = target_color*interval;

% values = 0:2:360-2;
% 
% for i = 1:length(values)
%     target_color(target_color == i) = values(i);
% end

cues(cues==360) = 0; % ????????????????????????????????????????
% How to convert values???

cues_cielab = cues;


% %% Convert CIELAB to CIELUV
% 
% sat = 38;
% 
% % Convert cues to CIELUV angles
% bae_hueangle = unique(cues);
% bae_lstar = 70;
% [bae_a, bae_b] = pol2cart(deg2rad(bae_hueangle),sat);
% bae_cielab = [ones(1,all_cues)*bae_lstar;bae_a'; bae_b'];
% %XYZ_D65 = [0.95047;1.0;1.08883];
% white_point = xyYToXYZ([0.3118, 0.3119, 48.64]');
% white_point = white_point ./ white_point(2);
% bae_xyz = LabToXYZ(bae_cielab,white_point);
% bae_cieluv = XYZToLuv(bae_xyz,white_point);
% 
% stimCols_sRGB = LuvTosRGB(bae_cieluv);
% colvals = im2double(stimCols_sRGB);
% 
% rotVal = 360/(all_cues*2);
% rotationMatrix = [cosd(rotVal),-sind(rotVal);sind(rotVal),cosd(rotVal)];
% stimCols_biased = rotationMatrix * bae_cieluv(2:3,:);
% rstimCols_sRGB = LuvTosRGB([bae_cieluv(1,:); stimCols_biased]);
% 
% 
% degree = rad2deg(cart2pol(bae_cieluv(2,:),bae_cieluv(3,:)));
% degree = degree';
% degree(degree < 0) = degree(degree < 0) + 360;
% 
% cues_cieluv = cues;
% 
% cue_vals = unique(cues_cielab);
% for i = 1:length(cue_vals)
%     val = cue_vals(i);
%     deg = degree(i);
% cues_cieluv(cues_cieluv == val) = deg;
% end
% 
% cues = cues_cieluv;

%%
error = zeros(size(cue_angle));
for i = 1:length(colresp)% for trial
    if abs(colresp(i,2) - colresp(i,1)) < 180
        error(i,1) = (colresp(i,2) - colresp(i,1));
    elseif abs(colresp(i,2) - colresp(i,1)) > 180 && colresp(i,2) < colresp(i,1)
        error(i,1) = ((360 - colresp(i,1) + colresp(i,2)));
    else
        error(i,1) = (abs(colresp(i,2) - colresp(i,1)) - 360);
    end
end

incorrect_idx = error ~= 0;
incorrect = [cues(incorrect_idx), error(incorrect_idx)];

% What should bin edges be?
edges = -degrees/2:degrees:360+degrees/2 %hue_angle; %center bins at hue values?
%edges = 0:degree:360;
%cuebin_idx = discretize(cues(:,1), edges, 'IncludedEdge', 'right');
cuebin_idx = discretize(cues(:,1),edges,'IncludedEdge','right');
cuebin_idx(cuebin_idx == 65) = 1;

cue_values = hue_angle';

%choice_edges = -180:degree:180;
choice_edges = -180-degrees/2:degrees:180+degrees/2;
%x = 0;
for i = 1:64% max(cuebin_idx) % = 1:bigN
    % x = x + 1;
    idx = find(cuebin_idx(:,1) == i);
    trial_number(1:bigN,i) = length(idx);
    okay = discretize(error(idx,:),choice_edges,'IncludedEdge','left');
    okay(okay == 65) = 1;
    for k = 1:bigN
        counts(k,i) = sum(okay==k);
    end
    %        counts(:,x) = histcounts(error(idx,:),choice_edges);
    counts(:,i) = counts(:,i) ./ trial_number(1,i);
end

bin_width = abs(choice_edges(1,1) - choice_edges(1,2));
bin_centers = (choice_edges(1:end-1) + bin_width/2)';
%
%% Calculate curve fits - no plotting
for i = 1:bigN
    gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d'; %a = amp, b = mu, c = stdv, d = vert. offset
    startingPoints = [0.5 0 78 0];
    %f = fit(bin_centers, [counts(:,i); counts(1,i)],gaussEqn,'start',startingPoints);% repeat first element at end
     f = fit(bin_centers(1:end-1,:), counts(:,i),gaussEqn,'start',startingPoints,'Lower',[-Inf -180 0 0],'Upper',[Inf 180 Inf 1]);% repeat first element at end
    ci = confint(f,0.95);
    bias(i,1) = f.b;
    stdv(i,1) = f.c;
    ci_lower_95(i,1) = ci(1,2);
    ci_upper_95(i,1) = ci(2,2);
    stdv_lower95(i,1) = ci(1,3);
    stdv_upper95(i,1) = ci(2,3);
end
%%
moving_bias  = movmean(bias([63:end 1:end 1:2]),5,'Endpoints','discard');
lower_95 = movmean(ci_lower_95([63:end 1:end 1:2]),5,'Endpoints','discard');
upper_95 = movmean(ci_upper_95([63:end 1:end 1:2]),5,'Endpoints','discard');


be_w = moving_bias([1:end,1]); % bias estimates including wraparound

for i = 1:bigN
    crossesZero(i) = and(be_w(i)>=0, be_w(i+1)<=0);
end

crossings = find(crossesZero);

% % Category Center Location Interpolation
hue_angle = 0:360/bigN:360;
%moving_bias = moving_bias([1:end 1]);

for i = 1:length(crossings)
    x = [hue_angle(crossings(i)) hue_angle(crossings(i)+1)];
    y = [be_w(crossings(i)) be_w(crossings(i)+1)];
    interp_crossing(i,1) = interp1(y,x,0);
end

% confidence interval stuff

for i = 1:bigN % check whether the confidence interval for each cue includes 0
    withinCI(i) = and(lower_95(i,1)<=0, upper_95(i,1)>=0);
end

with_CI = withinCI([1:end 1]);

for i = 1:bigN
    changes(i) = or(and(with_CI(i)==0, with_CI(i+1)==1), and(with_CI(i)==1, with_CI(i+1)==0));
end

range = find(changes);

range = range([1:end 1]);

for i = 1:length(crossings)
    for j = 1:sum(changes) 
        if range(j) < range(j+1)
            if (crossings(i) >= range(j) && crossings(i) <= range(j+1)) == 1 
                CI_range(i,:) = [range(j) range(j+1)];
            end
        elseif range(j) > range(j+1)
            if (crossings(i) >= range(j) && crossings(i) >= range(j+1)) == 1 ... 
                    || (crossings(i) <= range(j) && crossings(i) <= range(j+1) == 1)
                CI_range(i,:) = [range(j) range(j+1)];
            end
        end
    end
end


CI_range = CI_range';
CI_range = CI_range(:);

 upper_95 = upper_95([1:end 1]);
 lower_95 = lower_95([1:end 1]);
 
% Confidence Interval Interpolation 
for i = 1:length(CI_range)
    if moving_bias(CI_range(i)) > 0
        %if moving_bias(CI_range(i)) > moving_bias(CI_range(i)+1)
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [lower_95(CI_range(i)) lower_95(CI_range(i)+1)];
    elseif  moving_bias(CI_range(i)) < 0
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [upper_95(CI_range(i)) upper_95(CI_range(i)+1)];
    end
    interp_ci(i,1) = interp1(y,x,0);
    if isnan(interp_ci(i,1)) && moving_bias(CI_range(i)) > 0
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [upper_95(CI_range(i)) upper_95(CI_range(i)+1)];
        interp_ci(i,1) = interp1(y,x,0);
    elseif isnan(interp_ci(i,1)) && moving_bias(CI_range(i)) < 0
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [lower_95(CI_range(i)) lower_95(CI_range(i)+1)];
        interp_ci(i,1) = interp1(y,x,0);
    end
end


%CI_range = hue_angle(CI_range) + degrees/2;
 
%CI_range = CI_range([1:end 1]);



%% Figures: Bias by Cue with Category Crossings and Confidence Intervals

figure('WindowState', 'maximized');
% 
% % create filename
% if isfield(cleandata.trialdata, 'setSize') == 1 
%     filename = [dirname,'_SetSize',num2str(cleandata.trialdata.setSize)];
% else
%     filename = [dirname,'  - ', 'All Set Sizes'];
% end
% 
% % % display file name
% axes1 = axes('Position',[0.78 0.95 0.24 0.05]);
% text('Parent',axes1,'Interpreter','none','String',filename,'FontSize', 15, 'Position',[0.05 0.5 0]);
% set(axes1,'XColor','none','YColor','none')

moving_bias  = movmean(bias([63:end 1:end 1:2]),5,'Endpoints','discard');
axes2 = axes('Position', [0.13 0.13 0.8 0.78]);
hold on
hue_angle = 0:360/bigN:360';

lower_95 = movmean(ci_lower_95([63:end 1:end 1:2]),5,'Endpoints','discard');
upper_95 = movmean(ci_upper_95([63:end 1:end 1:2]),5,'Endpoints','discard');
plot(hue_angle, lower_95([1:end 1]), 'k:');
plot(hue_angle, upper_95([1:end 1]), 'k:');
h = fill([hue_angle, fliplr(hue_angle)], [lower_95([1:end 1])', fliplr(upper_95([1:end 1])')], 'k');
set(h, 'facealpha', .1, 'LineStyle', ':');
%scatter(hue_angle,bias([1:end 1]),100,colvals([1:end 1],:),'filled');
% scatter(unique(cues_cieluv), zeros(1,all_cues),12,colvals,'filled');%, colvals,'filled');
scatter(hue_angle, bias([1:end 1]),'filled','k');
plot(hue_angle,moving_bias([1:end 1]),'k');
xline(interp_ci(~isnan(interp_ci)),'--');
%xline(range*degrees,'--');
%xline(hue_angle(range));
%scatter(interp_crossing,0,'filled','k'); %add back in
%yticks([-50 -25 0 25 50])
%yticklabels({'-50','-25','0','25','50'})
xticks([0 60 120 180 240 300 360]);
hold off
ax = gca;
set(gca,'TickDir','out');
ax.TickLength = [0.025 0.025];
ax.FontSize = 10;
xlim([0 360]);
%ylim([-50 50]);
yline(0);
xlabel('Hue Angle ({\circ})');
ylabel('Bias ({\circ})');

%% Pie chart figure
hue_angle = 0:360/bigN:360;
hue_angle = hue_angle(1:bigN);
%%
figure('WindowState', 'maximized');

hold on
pie(repelem(degrees,bigN));% pie chart w/ bigN equally sized slices
% colormap(rstimCols_sRGB);
ax = gca;
delete(ax.Children([1, 1:2:bigN*2])) % stop displaying % for each slice
for i = 1:bigN
    ax.Children(i).EdgeAlpha = 0; % get rid of lines between slices
end
for i = 1:bigN
    ax.Children(i).FaceAlpha = 1;
end
ax.View = [90 90];


% Mark out anything outside of confidence intervals
% theta = angle in radians
x0=0;
y0=0;
CI_range = interp_ci';

%CI_range = interp_ci([1:end 1])';

CI_range = unique(CI_range,'stable');

CI_range = interp_ci([1:end 1])';

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

% Transparency for confidence intervals

unique_interpci = unique(interp_ci,'stable');
% 
unique_interpci = interp_ci
x = 0;
for i = 2:2:length(unique_interpci)
    
    x = x+1;
    change(x,1) = abs(unique_interpci(i) - unique_interpci(i-1));
end

opacity = (change ./min(change)) *0.15;

opacity_number = 0;
%CI_range = unique(CI_range,'stable');

unique_opacity = unique(opacity,'stable');
for i = 1:2:length(CI_range)-1
    opacity_number = opacity_number + 1;
    pie_direction = (.25*CI_range(1,i))/90 + 0.25;
    if CI_range(1,i+1) > CI_range(1,i)
        theta = deg2rad(CI_range(1,i+1) - CI_range(1,i));
    elseif CI_range(1,i+1) < CI_range(1,i)
        theta = deg2rad((360 - CI_range(1,i)) + CI_range(1,i+1));
    end
    a1 = 2*pi*pie_direction; % Starting direction
    a2 = a1 + theta; % Ending direction
    t = linspace(a1,a2);
    x = x0 + r*cos(t);
    y = y0 + r*sin(t);
   % fill([x0,x,x0],[y0,y,y0],'w','FaceAlpha',unique_opacity(opacity_number),'EdgeAlpha',0); %0.3
       fill([x0,x,x0],[y0,y,y0],'w','FaceAlpha',0.3,'EdgeAlpha',0); %0.3
end

% rotated_colvals = im2double(rstimCols_sRGB);

% shift_colvals = [colvals(135:end,:); colvals(1:134,:)];

[cart,~] = generateStimCols('nBig',180);
scatter(cart(1,:)./36,cart(2,:)./36,'filled'); % plot all cues around pie chart

set(gca,'visible','off')
axis equal

rad_angle = deg2rad(hue_angle); %hue angles in radians
moving_bias  = movmean(bias([63:end 1:end 1:2]),5,'Endpoints','discard');
axes3 = axes('Position', [0.318 0.1 0.4 0.835]);
polarplot(rad_angle([1:end 1]), moving_bias([1:end 1])+30,'b');
rlim([0 60]);
hold on
polarplot(rad_angle([1:end 1]), zeros(length(rad_angle)+1,1)+30,'LineStyle','--','Color','k');
polarplot(rad_angle([1:end 1]), moving_bias([1:end 1])+30,'k');
polarplot(rad_angle([1:end 1]), lower_95([1:end 1])+30,':k');
polarplot(rad_angle([1:end 1]), upper_95([1:end 1])+30,':k');
for k = 1:length(interp_crossing)
    polarplot([deg2rad(interp_crossing(k)) deg2rad(interp_crossing(k))],[0 60],'LineWidth',1.5);
end
ax = gca;
ax.Color = 'none';
 ax.ThetaTickLabel = {'0','','','90','','','180','','','270','',''};
% ax.RTick =[0 10 30 50 70 90];
% rticklabels({'','-40{\circ}','-20{\circ}','0{\circ}','20{\circ}','40{\circ}'});
ax.RTick = 0:10:60;
rticklabels({'','-20{\circ}','-10{\circ}','0{\circ}','10{\circ}','20{\circ}',''});
ax.RAxisLocation = 235;
ax.FontSize = 8;
ax.Units = 'normalized';
ax.Position = [0.3734375,0.212952799121844,0.2890625,0.609220636663008];

theta = rad_angle([1:end 1]);

ax_polar = gca;

rlow = lower_95([1:end 1])';
rhigh = upper_95([1:end 1])';

rlow = rlow + 30;
rhigh = rhigh + 30;

ax_cart = axes();
ax_cart.Position = ax_polar.Position;
[x1,y1] = pol2cart(theta,rlow);
[x2,y2] = pol2cart(theta,rhigh);
patch([x1 fliplr(x2)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
xlim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]);
ylim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]);
axis square;
set(ax_cart,'visible','off');

hold off
