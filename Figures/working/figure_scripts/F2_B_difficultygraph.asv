% Generated figure 2B (combined difficulty panel) for color categories
% paper

clc, clear, close all

smallN = 4; % number of choices on screen
bigN = 64; % total number of stimuli
monkeyN = 3; % total number of monkeys

%% load data
M1_cleandata = load(fullfile(pwd, 'AllData', 'pollux_cleandata.mat'));
M2_data = load(fullfile(pwd, 'AllData', 'Buster_colresp.mat'));
M3_cleandata = load(fullfile(pwd, 'AllData', 'Castor_cleandata.mat'));

% a little bit of messy processing to get everything in the same format
M1_data(:,1) = cell2mat(M1_cleandata.cleandata.trialdata.cues);
M1_data(:,2) = cell2mat(M1_cleandata.cleandata.trialdata.chosen);
M1_data(:,3:smallN+2) = cell2mat(M1_cleandata.cleandata.trialdata.choices);
M2_data = M2_data.colresp;
M3_data(:,1) = cell2mat(M3_cleandata.cleandata.trialdata.cues);
M3_data(:,2) = cell2mat(M3_cleandata.cleandata.trialdata.chosen);
M3_data(:,3:smallN+2) = cell2mat(M3_cleandata.cleandata.trialdata.choices);


%% Calculate difficulty values & model fit for each monkey

pct_corr = zeros(30,monkeyN);
b_glmfit = zeros(2,monkeyN);

alldata = {M1_data, M2_data, M3_data};
difficulties = 1:30;

for m = 1:monkeyN
    
    colresp = alldata{m};
    colresp(any(isnan(colresp),2),:) = []; % gets rid of aborted trials
    compltrial_choices = colresp(:,3:smallN+2); % get choices for only completed trials
    
    % For all completed trials
    difficulty = zeros(length(colresp),1);
    for i = 1:length(compltrial_choices)
        for j = 1:smallN
            compl_distance(1,j) = abs(colresp(i,1) - compltrial_choices(i,j)); % abs. distance from cue to choice j
            if compl_distance(1,j) > bigN/2
                compl_distance(1,j) = bigN - abs(colresp(i,1) - compltrial_choices(i,j));
            end
        end
        compl_distance(compl_distance<=0) = [];
        difficulty(i,1) = min(compl_distance);
    end
    
    diff_count = zeros(length(difficulties), 1);
    for i = 1:max(unique(difficulty))
        diff_count(i,1) = sum(difficulty==i);
    end
    
    % For correct trials only
    correct_data = colresp(colresp(:,1) == colresp(:,2),:);
    correct_choices = correct_data(:,3:smallN+2);
    correct_difficulty = zeros(length(colresp),1);
    
    for i = 1:length(correct_choices)
        for j = 1:smallN
            correct_distance(1,j) = abs(correct_data(i,1) - correct_choices(i,j)); % abs. distance from cue for incorrect trial i
            if correct_distance(1,j) > bigN/2
                correct_distance(1,j) = bigN - abs(correct_data(i,1) - correct_choices(i,j));
            end
        end
        correct_distance(correct_distance<=0) = [];
        correct_difficulty(i,1) = min(correct_distance);
    end
    
    cordiff_count = zeros(length(difficulties), 1);
    for i = 1:max(unique(difficulty))
        corrdiff_count(i,1) = sum(correct_difficulty==i);
    end
    
    % Percent correct by distance of closest distractor
    pct_corr(:,m) = corrdiff_count(:)./diff_count(:);
    
    
    % Fit model
    %y = K./(1+exp(-G*(x-Dm)));
    %f = fit(x,y,'1./(1+exp(b1+b2*x))','start',[ -1 1]);
    
    b_glmfit(:,m) = glmfit(difficulties,[corrdiff_count diff_count],'binomial','Link','logit');
    
end

%% FIGURES

figure;
hold on

% M1 data
M1col = 66/256;
yfit1 = glmval(b_glmfit(:,1),difficulties,'logit','Size',diff_count);
plot(difficulties,pct_corr(:,1),'.',difficulties,yfit1./diff_count,'-', ... 
    'Color', [M1col M1col M1col], 'LineWidth', 2, 'MarkerSize', 10);

% M2 data
M2col = 200/256;
yfit2 = glmval(b_glmfit(:,2),difficulties,'logit','Size',diff_count);
plot(difficulties,pct_corr(:,2),'.',difficulties,yfit2./diff_count,'-', ... 
'Color', [M2col M2col M2col], 'LineWidth', 2, 'MarkerSize', 10);


% M3 data
M3col = 132/256;
yfit3 = glmval(b_glmfit(:,3),difficulties,'logit','Size',diff_count);
plot(difficulties,pct_corr(:,3),'.',difficulties,yfit3./diff_count,'-', ...
    'Color', [M2col M2col M2col], 'LineWidth', 2, 'MarkerSize', 10);


% M4 data
M4col = 0;

% chance line
yline(1/smallN, 'k--'); 

ylim([0 1]);
xlim([min(difficulties), max(difficulties)]);
ylabel('% Correct');
xticks([]);
yticks([0,1]);
yticklabels({'0', '100%'});
pbaspect([5 6 1]);

saveas(gcf,'F2_B_difficultygraph_v2.pdf')
