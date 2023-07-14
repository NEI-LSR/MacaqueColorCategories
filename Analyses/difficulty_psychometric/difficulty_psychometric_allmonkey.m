
%% Colors
% Hard-coded - need to add more colors for additional monkeys
colvals = [175 171 171; 0 0 0; 103 99 99; 132 132 132; 185 185 185];
colvals = colvals ./ 255;
%% Process Inputs
figure('WindowState', 'maximized');

axes3 = axes('Position', [0.04, 0.31, 0.2, 0.55]);
for r = 1:length(cleandata_all) %for each monkey

    cleandata = cleandata_all{r};

    cues = cell2mat(cleandata.trialdata.cues);
    choices = cell2mat(cleandata.trialdata.choices);
    chosen = cell2mat(cleandata.trialdata.chosen);

    nBig = size(cleandata.trialdata.stimCols{1,1},1);
    nSmall = sum(~isnan(cleandata.trialdata.choices{end,1}));


    completed_idx = ~isnan(chosen);
    correct_idx = cues(completed_idx) == chosen(completed_idx);

    distance = abs(cues(completed_idx) - choices(completed_idx,:));
    distance(distance > nBig/2) = nBig - distance(distance > nBig/2);

    for trial = 1:size(cues(completed_idx),1)
        try
            difficulty(trial) = min(nonzeros(distance(trial,:)));
        catch
            if ~any(distance(trial,:) ~= 0) % if there are no choices where distance ~= 0
                difficulty(trial) = 0;
            else
                error('unknown error'); % I can't predict any situation where this would be reached, but just in case...
            end
        end
    end

    % total_difficulties = 1:30;
    %
    %     difficulties{r} = total_difficulties;

    unique_difficulties = 1:max(unique(difficulty));

    for i = 1:max(unique_difficulties)
        difficulties_completed_counts(i,1) = sum(difficulty==i);
    end

    for i = 1:max(unique_difficulties)
        difficulties_correct_counts(i,1) = sum(difficulty(correct_idx)==i);
    end

    % Percent correct by distance of closest distractor
    pct_correct = difficulties_correct_counts./difficulties_completed_counts;
    pct_correct(isnan(pct_correct)) = 0;

    %pct_corr{r} = pct_correct(1:max(unique_difficulties),:);

    %pct_corr{r} = pct_correct;

    % Fit Weibull function with 4 params

    Weibull = fittype(@(slope, inflect, floor, ceil, x) (floor+(1-floor-ceil).*(1-exp(-(x./slope).^inflect)))); %Define version of Weibull with 4 args

    %Weighted fit
    [f, gof] = fit(unique_difficulties',pct_correct,...
        Weibull,...
        'weights', difficulties_completed_counts,...
        'Lower', [0 -20 0 0],...
        'Upper', [20 20 1 1],...
        'StartPoint', [10 0.01 0.5 0.5]);
%end

% FIGURES
%axes3 = gca;
hold on

plot(axes3, unique_difficulties, pct_correct, '.','MarkerEdgeColor','none');
fitplot = plot(f); set(fitplot,'color',colvals(r,:)); set(fitplot,'LineWidth',2);
clear pct_correct
 clear difficulties_correct_counts
 clear difficulties_completed_counts
clear f
clear difficulty 
end

%%
ylim([0.25 1]);
yticks([0.25 1]);
yticklabels({'25%', '100%'});
xticks([1 30]);
xticklabels({'',''});
%xlim([min(difficulties{1}), max(difficulties{1})]);
xlim([min(unique_difficulties), max(unique_difficulties)])
title('Accuracy by Trial Difficulty');
    xlabel('Distance of Closest Distractor');
    ylabel('Accuracy');
legend1 = legend('','M1','','M2','','M3','','M4','','Location','southeast');
