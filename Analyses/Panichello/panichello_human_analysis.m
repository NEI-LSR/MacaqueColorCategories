function panichellobias_analysis(cleandata)
%% Load in data

cues = cleandata.trialdata.cues;
chosen = cleandata.trialdata.chosen;
dirname = cleandata.trialdata.dirname;
%%
cues(cues(:,1)==360) = 0; 

% for human data 
chosen(chosen > 360) = chosen(chosen > 360) - 360; % for human data only
chosen(chosen(:,1)==360) = 0; % for human data only
% 
% cues_cielab = cues;
% 
% %% Convert CIELAB to CIELUV
% 
% % Convert cues to CIELUV angles
all_cues = length(unique(cues)); % Number of unique cue colors  
% 
% if all_cues == 64 % for monkey data
%     sat = 57;
% elseif all_cues == 360 % for human data
%     sat = 52;
% end
% 
% panichello_hueangle = unique(cues);
% panichello_lstar = 60;
% [panichello_a, panichello_b] = pol2cart(deg2rad(panichello_hueangle),sat);
% panichello_cielab = [ones(1,all_cues)*panichello_lstar;panichello_a'; panichello_b'];
% XYZ_D65 = [0.95047;1.0;1.08883];
% panichello_xyz = LabToXYZ(panichello_cielab,XYZ_D65);
% panichello_cieluv = XYZToLuv(panichello_xyz,XYZ_D65);
% 
% stimCols_sRGB = LuvTosRGB(panichello_cieluv);
% colvals = im2double(stimCols_sRGB);
% 
% rotVal = 360/(all_cues*2);
% rotationMatrix = [cosd(rotVal),-sind(rotVal);sind(rotVal),cosd(rotVal)];
% stimCols_biased = rotationMatrix * panichello_cieluv(2:3,:);
% rstimCols_sRGB = LuvTosRGB([panichello_cieluv(1,:); stimCols_biased]);
% 
% degrees = rad2deg(cart2pol(panichello_cieluv(2,:),panichello_cieluv(3,:)));
% degrees = degrees';
% degrees(degrees < 0) = degrees(degrees < 0) + 360;
% 
% cues_cieluv = cues;
% 
% cue_vals = unique(cues_cielab);
% for i = 1:length(cue_vals)
%     val = cue_vals(i);
%     deg = degrees(i);
% cues_cieluv(cues_cieluv == val) = deg;
% end
% 
% cues = cues_cieluv;
% 
% %%
% % Convert chosen values to CIELUV angles 
% all_chosen = length(unique(chosen));
% 
% % Do I need to any shifting with these values?
% panichello_cielab = generateStimCols('nBig',all_chosen,'sat',sat);
% panichello_cielab =[ones(1,all_chosen)*60;panichello_cielab];
% XYZ_D65 = [0.95047;1.0;1.08883];
% panichello_xyz = LabToXYZ(panichello_cielab,XYZ_D65);
% panichello_cieluv = XYZToLuv(panichello_xyz,XYZ_D65);
% 
% panichello_chosenangle = unique(chosen);
% panichello_lstar = 60;
% [chosen_a, chosen_b] = pol2cart(deg2rad(panichello_chosenangle),sat);
% chosen_cielab = [ones(1,all_chosen)*panichello_lstar;chosen_a'; chosen_b'];
% XYZ_D65 = [0.95047;1.0;1.08883];
% chosen_xyz = LabToXYZ(chosen_cielab,XYZ_D65);
% chosen_cieluv = XYZToLuv(chosen_xyz,XYZ_D65);
% 
% chosen_degrees = rad2deg(cart2pol(chosen_cieluv(2,:),chosen_cieluv(3,:)));
% chosen_degrees = chosen_degrees';
% chosen_degrees(chosen_degrees < 0) = chosen_degrees(chosen_degrees < 0) + 360;
% 
% % Convert choices to cieluv angles 
% chosen_luv = chosen;
% chosen_lab = chosen;
% 
% chosen_vals = unique(chosen);
% for i = 1:length(chosen_vals)
%     val = chosen_vals(i);
%     deg = chosen_degrees(i);
% chosen_luv(chosen_luv == val) = deg;
% end
% 
% chosen = chosen_luv;
% 
nBig = 64;
colresp = [cues, chosen];
degree = 360/nBig;
interval = degree;

hue_angle = 0:degree:360;
hue_angle = hue_angle(1:nBig);
%% Colors

% panichello_lab = generateStimCols('nBig',64,'sat',sat);
% panichello_lab =[ones(1,64)*60;panichello_lab];
% XYZ_D65 = [0.95047;1.0;1.08883];
% panichello_xyz = LabToXYZ(panichello_cielab,XYZ_D65);
% panichello_cieluv = XYZToLuv(panichello_xyz,XYZ_D65);
% 
% stimCols_sRGB = LuvTosRGB(panichello_cieluv);
% colvals = im2double(stimCols_sRGB);
% 
% stimCols = generateStimCols('nBig',nBig,'sat',52);
% stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols]);
% colvals = im2double(stimCols_sRGB);

% rotVal = 360/(nBig*2);
% rotationMatrix = [cosd(rotVal),-sind(rotVal);sind(rotVal),cosd(rotVal)];
%  stimCols_biased = rotationMatrix * panichello_cieluv(2:3,:);
%  rstimCols_sRGB = LuvTosRGB([panichello_cieluv(1,:); stimCols_biased]);

%% Calculate Error
% Human Data (Cues = 360)
if all_cues == 360
    error = zeros(size(chosen));
    for i = 1:length(cues)% for trial
        if abs(colresp(i,2) - colresp(i,1)) < 180
            error(i,1) = (colresp(i,2) - colresp(i,1));
        elseif abs(colresp(i,2) - colresp(i,1)) > 180 && colresp(i,2) < colresp(i,1)
            error(i,1) = ((360 - colresp(i,1) + colresp(i,2)));
        else
            error(i,1) = (abs(colresp(i,2) - colresp(i,1)) - 360);
        end
    end
   
    edges = -2.8125:degree:360+degree/2; %hue_angle; %center bins at hue values?
    cuebin_idx = discretize(cues(:,1),edges,'IncludedEdge','right');
    cuebin_idx(cuebin_idx == 65) = 1;
    
    cue_values = hue_angle';
 
    choice_edges = -180-degree/2:degree:180+degree/2;
    for i = 1:max(cuebin_idx)
        idx = find(cuebin_idx(:,1) == i);
        trial_number(1:64,i) = length(idx);
        binned = discretize(error(idx,:),choice_edges,'IncludedEdge','left');
        binned(binned == 65) = 1;
        for k = 1:nBig
            counts(k,i) = sum(binned==k);
        end
        counts(:,i) = counts(:,i) ./ trial_number(1,i);
    end
    
    bin_width = abs(choice_edges(1,1) - choice_edges(1,2));
    bin_centers = (choice_edges(1:end-1) + bin_width/2)';
end

%%
% Monkey Data (Cues = 64)
if all_cues == 64
    error = zeros(size(chosen)); % distance of incorrect choice from cue, i.e. delta theta
    for i = 1:length(cues) % for each incorrect trial
        if abs(colresp(i,2) - colresp(i,1)) < 180
            error(i,1) = (colresp(i,2) - colresp(i,1));  
        elseif abs(colresp(i,2) - colresp(i,1)) > 180 && colresp(i,2) < colresp(i,1)
            error(i,1) = ((360 - colresp(i,1) + colresp(i,2)));
        else
            error(i,1) = (abs(colresp(i,2) - colresp(i,1)) - 360);  
        end
    end
    
    cue_values = unique(cues); % all cue values
    choice_edges = -180-degree/2:degree:180+degree/2; % choice edges such that centers are at our angles
    bin_width = abs(choice_edges(1,1) - choice_edges(1,2)); % bin width should be = degrees
    bin_centers = (choice_edges(1:end-1) + bin_width/2)'; % bin centers are at our values
    
    for i = 1:nBig
       idx = find(colresp(:,1) == cue_values(i)); % find trials for cue i 
       trial_number(1:nBig,i) = length(idx); % number of trials for given cue
       binned = discretize(error(idx,:),choice_edges,'IncludedEdge','left'); % sort into bins
       binned(binned == 65) = 1; % bin 65 is same as 1 
         for k = 1:nBig
             counts(k,i) = sum(binned==k); % total number in each bin
         end
    end

end

%% Calculate curve fits - no plotting
for i = 1:nBig
    gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d'; %a = amp, b = mu, c = stdv, d = vert. offset
    startingPoints = [0.5 0 78 0];
    f = fit(bin_centers(1:end-1,:), counts(:,i),gaussEqn,'start',startingPoints,'Lower',[0 -180 0 0],'Upper',[Inf 180 Inf 1]);% repeat first element at end
    ci = confint(f,0.95);
    bias(i,1) = f.b;
    ci_lower_95(i,1) = ci(1,2);
    ci_upper_95(i,1) = ci(2,2);
end

%%
moving_bias  = movmean(bias([63:end 1:end 1:2]),5,'Endpoints','discard');
lower_95 = movmean(ci_lower_95([63:end 1:end 1:2]),5,'Endpoints','discard');
upper_95 = movmean(ci_upper_95([63:end 1:end 1:2]),5,'Endpoints','discard');

be_w = moving_bias([1:end,1]); % bias estimates including wraparound

% Category Center Location Interpolation
for i = 1:nBig
    crossesZero(i) = and(be_w(i)>=0, be_w(i+1)<=0);
end

crossings = find(crossesZero);


for i = 1:length(crossings)
    x = [hue_angle(crossings(i)) hue_angle(crossings(i)+1)];
    y = [be_w(crossings(i)) be_w(crossings(i)+1)];
    interp_crossing(i,1) = interp1(y,x,0);
end


for i = 1:nBig % check whether the confidence interval for each cue includes 0
    withinCI(i) = and(lower_95(i,1)<=0, upper_95(i,1)>=0);
end

with_CI = withinCI([1:end 1]);


for i = 1:nBig
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

 upper_95_w = upper_95([1:end 1]);
 lower_95_w = lower_95([1:end 1]);
 
% Confidence Interval Interpolation 
for i = 1:length(CI_range)
    if moving_bias(CI_range(i)) > 0
        %if moving_bias(CI_range(i)) > moving_bias(CI_range(i)+1)
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [lower_95_w(CI_range(i)) lower_95_w(CI_range(i)+1)];
    elseif  moving_bias(CI_range(i)) < 0
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
    end
    interp_ci(i,1) = interp1(y,x,0);
    if isnan(interp_ci(i,1)) && moving_bias(CI_range(i)) > 0
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
        interp_ci(i,1) = interp1(y,x,0);
    elseif isnan(interp_ci(i,1)) && moving_bias(CI_range(i)) < 0
        x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
        y = [lower_95_w(CI_range(i)) lower_95_w(CI_range(i)+1)];
        interp_ci(i,1) = interp1(y,x,0);
    end
end



%% Figures: Bias by Cue with Category Crossings and Confidence Intervals

figure('WindowState', 'maximized');

% create filename
if isfield(cleandata.trialdata, 'setSize') == 1 
    filename = [dirname,'_SetSize',num2str(cleandata.trialdata.setSize)];
else
    filename = [dirname,'  - ', 'All Set Sizes'];
end

% % display file name
axes1 = axes('Position',[0.78 0.95 0.24 0.05]);
text('Parent',axes1,'Interpreter','none','String',filename,'FontSize', 15, 'Position',[0.05 0.5 0]);
set(axes1,'XColor','none','YColor','none')

moving_bias  = movmean(bias([63:end 1:end 1:2]),5,'Endpoints','discard');
axes2 = axes('Position', [0.13 0.13 0.8 0.78]);
hold on
hue_angle = 0:360/nBig:360';

lower_95 = movmean(ci_lower_95([63:end 1:end 1:2]),5,'Endpoints','discard');
upper_95 = movmean(ci_upper_95([63:end 1:end 1:2]),5,'Endpoints','discard');
plot(hue_angle, lower_95([1:end 1]), 'k:');
plot(hue_angle, upper_95([1:end 1]), 'k:');
h = fill([hue_angle, fliplr(hue_angle)], [lower_95([1:end 1])', fliplr(upper_95([1:end 1])')], 'k');
set(h, 'facealpha', .1, 'LineStyle', ':');
if all_cues == 64
    scatter(hue_angle,bias([1:end 1]),100,colvals([1:end 1],:),'filled');
elseif all_cues ==360
    scatter(1:all_cues, zeros(1,all_cues),8,'filled');
end
plot(hue_angle,moving_bias([1:end 1]),'k');
xline(interp_ci,'--');
%xline(hue_angle(range));
scatter(interp_crossing,0,'filled','k');
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

% create filename
if isfield(cleandata.trialdata, 'setSize') == 1 
    filename = [dirname,'_SetSize',num2str(cleandata.trialdata.setSize)];
else
    filename = [dirname,'  - ', 'All Set Sizes'];
end

interval = degree;
figure('WindowState', 'maximized');

hold on
pie(repelem(interval,nBig));% pie chart w/ nBig equally sized slices
% colormap(rstimCols_sRGB);
ax = gca;
delete(ax.Children([1, 1:2:nBig*2])) % stop displaying % for each slice
for i = 1:nBig
    ax.Children(i).EdgeAlpha = 0; % get rid of lines between slices
end
for i = 1:nBig
    ax.Children(i).FaceAlpha = 0.5;
end
ax.View = [90 90];

% Mark out anything outside of confidence intervals
if isempty(ci) == 0 && isempty(interp_ci) == 0 % && isempty(change_range) == 0
    
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
   
    
    for i = 1:2:length(CI_range)-1
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
        fill([x0,x,x0],[y0,y,y0],'w','FaceAlpha',0.3,'EdgeAlpha',0);
    end
else
    for i = 1:nBig
        ax.Children(i).FaceAlpha = 0;
    end
end

ax.Color = 'none';

% rotated_colvals = im2double(rstimCols_sRGB);

% shift_colvals = [colvals(all_cues*(3/4):end,:); colvals(1:all_cues*(3/4)-1,:)];
[cart,~] = generateStimCols('nBig',all_cues); % generate values to plot cues

scatter(cart(1,:)./36,cart(2,:)./36,185,'filled'); % plot all cues around pie chart

set(gca,'visible','off')
axis equal tight

rad_angle = deg2rad(hue_angle); %hue angles in radians
axes3 = axes('Position', [0.318 0.1 0.4 0.835]);
polarplot(rad_angle, be_w+40,'b');
rlim([0 80]);
hold on
polarplot(rad_angle, zeros(length(rad_angle),1)+40,'LineStyle','--','Color','k');
polarplot(rad_angle, be_w+40,'k','Linewidth',1.75);

% if isempty(ci) == 0
%     polarplot(rad_angle, lower_95_w+40,':k');
%     polarplot(rad_angle, upper_95_w+40,':k');
% end

% add lines
for k = 1:length(interp_crossing)
    polarplot([deg2rad(interp_crossing(k)) deg2rad(interp_crossing(k))],[0 80],'LineWidth',1.5);
end


ax = gca;
ax.Color = 'none';
ax.ThetaTickLabel = {'0','','','90','','','180','','','270','',''};
ax.RTick = 0:20:80;
rticklabels({'-40{\circ}','-20{\circ}','0{\circ}','20{\circ}','40{\circ}'});

ax.RAxisLocation = 235;
% ax.FontSize = 8;
% ax.Units = 'normalized';
% ax.Position = [0.3734375,0.212952799121844,0.2890625,0.609220636663008];

theta = rad_angle;

ax_polar = gca;
ax_cart = axes();
ax_cart.Position = ax_polar.Position;

if isempty(ci) == 0
    rlow = lower_95_w' + 40;
    rhigh = upper_95_w' + 40;
    
    [x1,y1] = pol2cart(theta,rlow);
    [x2,y2] = pol2cart(theta,rhigh);
    patch([x1 fliplr(x2)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
end


xlim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]);
ylim(ax_cart,[-max(get(ax_polar,'RLim')),max(get(ax_polar,'RLim'))]);
axis square;
set(ax_cart,'visible','off');

% % display file name
axes1 = axes('Position',[0.81 0.95 0.24 0.05]);
text('Parent',axes1,'Interpreter','none','String',filename,'FontSize', 15, 'Position',[0.05 0.5 0]);
set(axes1,'XColor','none','YColor','none');
end