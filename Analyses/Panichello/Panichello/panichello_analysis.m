function panichello_analysis(cleandata)
%% Load in data
cues = cleandata.trialdata.cues;
chosen = cleandata.trialdata.chosen;
paradigm = cleandata.trialdata.paradigm;
dirname = cleandata.trialdata.dirname;

if isfield(cleandata.trialdata, 'choices') == 1 % if our data
    choices = cleandata.trialdata.choices;
    colresp = [cues, chosen, choices];
    all_paradigms = unique(cleandata.trialdata.paradigm);
    all_dirnames = unique(cleandata.trialdata.dirname);
    load(length(cues),1) == 1;
    
else colresp = [cues, chosen]; % if panichello data 
    load = cleandata.trialdata.cond.load;
    delay = cleandata.trialdata.cond.delay;
    panichello_error = cleandata.trialdata.cond.error;
end

all_cues = length(unique(cues)); % number of unique cues

%% Subset data for specific paradigm
% assumes that last paradigm in list is the correct one

%{
if recursive == 0
    subset_id = paradigm(end);
    subset_idx = strcmp(subset_id, all_paradigms);
    
    cues = cues(subset_idx);
    choices = choices(subset_idx,:);
    chosen = chosen(subset_idx);
end
%}

%% Distribution of Angular Error (Fig 1B)

if isfield(cleandata.trialdata, 'load') == 1
% Cond 1 = Load 1, Short Delay
% Cond 2 = Load 1, Long Delay
% Cond 3 = Load 2/3, Short Delay
% Cond 4 = Load2/3, Long Delay
    load_values = unique(load);
    delay_values = unique(delay);

    cond1_idx = find(load == min(load) & delay == min(delay));
    cond2_idx = find(load == min(load) & delay == max(delay));
    cond3_idx = find(load == max(load) & delay == min(delay));
    cond4_idx = find(load == max(load) & delay == max(delay));

    cond1 = [load(cond1_idx,1), delay(cond1_idx,1), error(cond1_idx,1)];
    cond2 = [load(cond2_idx,1), delay(cond2_idx,1), error(cond2_idx,1)];
    cond3 = [load(cond3_idx,1), delay(cond3_idx,1), error(cond3_idx,1)];
    cond4 = [load(cond4_idx,1), delay(cond4_idx,1), error(cond4_idx,1)];

    cond1_unique = unique(cond1(:,3));
    cond2_unique = unique(cond2(:,3));
    cond3_unique = unique(cond3(:,3));
    cond4_unique = unique(cond4(:,3));

    for j = 1:length(cond1_unique);
        i = cond1_unique(j,:);
        error_count = cond1(cond1(:,3)==i, :); % sort all incorrect trials for cue i
        errorcount_cond1(j,:) = histc(error_count(:,3),i);
    end

    errorcount_cond1 = errorcount_cond1 ./sum(errorcount_cond1);

    for j = 1:length(cond2_unique)
        i = cond2_unique(j,:);
        error_count = cond2(cond2(:,3)==i, :); % sort all incorrect trials for cue i
        errorcount_cond2(j,:) = histc(error_count(:,3),i);
    end

    errorcount_cond2 = errorcount_cond2 ./sum(errorcount_cond2);

    for j = 1:length(cond3_unique)
        i = cond3_unique(j,:);
        error_count = cond3(cond3(:,3)==i, :); % sort all incorrect trials for cue i
        errorcount_cond3(j,:) = histc(error_count(:,3),i);
    end

    errorcount_cond3 = errorcount_cond3 ./sum(errorcount_cond3);

    for j = 1:length(cond4_unique)
        i = cond4_unique(j,:);
        error_count = cond4(cond4(:,3)==i, :); % sort all incorrect trials for cue i
        errorcount_cond4(j,:) = histc(error_count(:,3),i);
    end

    errorcount_cond4 = errorcount_cond4 ./sum(errorcount_cond4);
    
    % Figure 
    figure;
    hold on
    %bar(val_1short, error_1short);
    gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d'; %a = amp, b = mu, c = stdv, d = vert. offset
    startingPoints = [100 0 15 3];
    [f,gof]= fit(cond1_unique, errorcount_cond1, gaussEqn,'start',startingPoints);
    short1 = plot(f,cond1_unique, errorcount_cond1);
    set(short1(2),'LineWidth',2);
    set(short1(1),'visible','off');

    [f,gof]= fit(cond2_unique, errorcount_cond2, gaussEqn, 'start',startingPoints);
    long1 = plot(f,cond2_unique, errorcount_cond2);
    set(long1(2),'Linestyle','--','LineWidth',2);
    set(long1(1),'visible','off');

    [f,gof]= fit(cond3_unique, errorcount_cond3, gaussEqn, 'start',startingPoints);
    short3 = plot(f,cond3_unique, errorcount_cond3);
    set(short3(2),'Color','b','LineWidth',2);
    set(short3(1),'visible','off');

    [f,gof]= fit(cond4_unique, errorcount_cond4, gaussEqn, 'start',startingPoints);
    long3 = plot(f,cond4_unique, errorcount_cond4);
    set(long3(2),'Linestyle','--','Color','b','LineWidth',2);
    set(long3(1),'visible','off');
    legend('hide');
    xlabel('Angular error');
    ylabel('Density');
    title('Human');
    xlim([-180 180]);
    hold off
end

%% Figure 2C
for i = 1:length(all_cues);
    cue_count = cues(cues==i,:);
    cue_total(i,1) = length(cue_count);
end

cue_total = cue_total ./ length(cues);

for i = 1:length(chosen)
    if chosen(i,1) > 360
        chosen(i,1) = (chosen(i,1) - 360);
    end
end

for i = 1:length(all_cues);
    resp_count = chosen(chosen==i,:);
    resp_total(i,1) = length(resp_count);
end
resp_total = resp_total ./ length(cues);

figure;
hold on
plot(all_cues, cue_total);
bar(all_cues, resp_total);
hold off

%% Calculate Error 
% Human Data (Cues = 360)
if all_cues == 360
    error = zeros(size(chosen));
    for i = 1:length(cues)% for trial
        if abs(colresp(i,1) - colresp(i,2)) < 180 
            error(i,1) = (colresp(i,1) - colresp(i,2));
        elseif abs(colresp(i,1) - colresp(i,2)) > 180 && colresp(i,1) < colresp(i,2)
            error(i,1) = ((colresp(i,1) - colresp(i,2)) + 360);
        elseif abs(colresp(i,1) - colresp(i,2)) > 180 && colresp(i,1) > colresp(i,2)
            error(i,1) = (360 - colresp(i,1)) + colresp(i,2);
        end
    end
end

% Monkey Data (Cues = 64)
if all_cues == 64
    error = zeros(size(chosen)); % distance of incorrect choice from cue, i.e. delta theta
    for i = 1:length(cues) % for each incorrect trial
        if abs(colresp(i,1) - colresp(i,2)) < 180
            error(i,1) = (colresp(i,1) - colresp(i,2));
        elseif abs(colresp(i,1) - colresp(i,2)) > 180 && colresp(i,1) > colresp(i,2)
        error(i,1) = (-(180 - abs(colresp(i,1) - colresp(i,2))));
        else
        error(i,1) = (180 - abs(colresp(i,1) - colresp(i,2)));  
        end
    end
end

%%

%monkeyw = table2array(monkeyw);
all_cues = unique(cues);
cue_dist = [cues, chosen];
%% for Humans
data = cues;
edges = [0:4:360];
bin_idx = discretize(data, edges);
%rad_chosen = deg2rad(circ_chosen);
rad_chosen = deg2rad(chosen);
for i = 1:max(bin_idx);
%for i = 1;
    idx = find(bin_idx == i);
    cue_dist = [cues(idx,:),rad_chosen(idx,:)];
    %cue_dist = [cues(idx,:), chosen(idx,:)];
    cue_avg(i,1) = mean(cue_dist(:,1));
    resp_mean(i,1) = circ_mean(cue_dist(:,2));
    %resp_mean(i,1) = mean(cue_dist(:,2));
end

for i = 1:length(resp_mean)
    resp_mean(i,1) = rad2deg(resp_mean(i,1));
    resp_mean(i,1) = resp_mean(i,1) - cue_avg(i,1);
   % resp_mean(i,1) = resp_mean(i,1) - cue_avg(i,1);
end

for j = 46:length(resp_mean)
       resp_mean(j,1) = resp_mean(j,1) + 360;
   end

%%
figure;
plot(1:90, resp_mean);

%fsadsdjlfsd


%% 

%sdfjlsjdkldfs
%% for monkey data
data = cues;
edges = [0:6:360];
bin_idx = discretize(data, edges);
%rad_chosen = deg2rad(circ_chosen);
rad_chosen = deg2rad(chosen);
for i = 1:max(bin_idx);
%for i = 1;
    idx = find(bin_idx == i);
    cue_dist = [cues(idx,:),rad_chosen(idx,:)];
    %cue_dist = [cues(idx,:), chosen(idx,:)];
    cue_avg(i,1) = mean(cue_dist(:,1));
    resp_mean(i,1) = circ_mean(cue_dist(:,2));
    %resp_mean(i,1) = mean(cue_dist(:,2));
end

for i = 1:length(resp_mean)
    resp_mean(i,1) = rad2deg(resp_mean(i,1));
    resp_mean(i,1) = resp_mean(i,1) - cue_avg(i,1);
   % resp_mean(i,1) = resp_mean(i,1) - cue_avg(i,1);
end

for j = 28:length(resp_mean)
       resp_mean(j,1) = resp_mean(j,1) + 360;
   end


%%
for i = 1:length(all_cues)
    cue_value = all_cues(i,1);
    dist_count = cue_dist(cue_dist(:,1)==cue_value, 2); % sort all incorrect trials for cue i
    error_mean(i,:) = mean(dist_count);
    %error_counts(:,i) = histc(dist_count(:,2),total_distances);
end


 
%% Plot biases (Figure 4d)
figure('WindowState','maximized');
% human 4d
subplot(1,3,1)
hold on
plot(1:length(human4d),human4d(:,1));
plot(1:length(human4d), human4d(:,2));
hold off
legend('Empirical', 'Model');
ylim([-20 20]);
xlim([0 90]);
ylabel('Bias (mean response - target)');
xlabel('Target color');
title('Human');

% monkey W 4d
subplot(1,3,2)
hold on
plot(1:length(monkeyw4d),monkeyw4d(:,1));
plot(1:length(monkeyw4d), monkeyw4d(:,2));
hold off
legend('Empirical', 'Model');
ylim([-50 50]);
ylabel('Bias (mean response - target)');
xlabel('Target color');
title('Monkey W');

subplot(1,3,3)
hold on
plot(1:length(monkeye4d),monkeye4d(:,1));
plot(1:length(monkeye4d), monkeye4d(:,2));
hold off
legend('Empirical', 'Model');
ylim([-50 50]);
ylabel('Bias (mean response - target)');
xlabel('Target color');
title('Monkey E');


end
