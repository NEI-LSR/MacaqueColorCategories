function [unique_difficulties,pct_correct,...
    difficulties_completed_counts, f]...
    = difficulty_psychometric(cleandata)
%% PURPOSE:
%Error analysis by difficulty
%Get difficulty for each trial, plot errors by difficulty
%Difficulty calculated as distance of closest distractor to cue

%% Process Inputs
cues = cell2mat(cleandata.trialdata.cues);
choices = cell2mat(cleandata.trialdata.choices);
chosen = cell2mat(cleandata.trialdata.chosen);
paradigm = unique(cleandata.trialdata.paradigm);

% If real data, get dirname from cleandata
if isfield(cleandata.trialdata, 'dirname') == 1
    dirname = unique(cleandata.trialdata.dirname);
end

nBig = size(cleandata.trialdata.stimCols{1,1},1);
nSmall = sum(~isnan(cleandata.trialdata.choices{end,1}));


%% Calculate difficulty

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


%% Fit Weibull function with 4 params

Weibull = fittype(@(slope, inflect, floor, ceil, x) (floor+(1-floor-ceil).*(1-exp(-(x./slope).^inflect)))); %Define version of Weibull with 4 args

%% Weighted fit
[f, gof] = fit(unique_difficulties',pct_correct,...
    Weibull,...
    'weights', difficulties_completed_counts,...
    'Lower', [0 -20 0 0],...
    'Upper', [20 20 1 1],...
    'StartPoint', [10 0.01 0.5 0.5]);

%% FIGURES

if isfield(cleandata.trialdata, 'dirname') == 1 % if real data

    % create filename
    if numel(dirname) == 1
        filename = dirname{1};
        file_dir = dirname{1};
    elseif numel(dirname) > 1
        dates = str2double(extractBefore(dirname, 7)); % DG comment: does this clash with the section above?
        [dates, date_order]  = sort(dates);
        filename = [num2str(dates(1)), '--', num2str(dates(end)), extractAfter(dirname{1}, 13)];
        file_dir = dirname{date_order(end)};
    end

    figure('WindowState', 'maximized');

    % difficulty plot ("psychometric function")
    axes1 = axes('Position', [0.04, 0.31, 0.2, 0.55]);
    hold on
    plot(axes1, unique_difficulties, pct_correct, '.','MarkerEdgeColor','none');
    fitplot = plot(f); set(fitplot,'color',[0 0.6 0.3]); set(fitplot,'LineWidth',2);

    ylim([0.25 1]);
    yticks([0.25 1]);
    yticklabels({'25%', '100%'});
    xlim([min(unique_difficulties), max(unique_difficulties)]);
    xticks([min(unique_difficulties) max(unique_difficulties)]);
    xticklabels({'',''});
    title('Accuracy by Trial Difficulty');
    xlabel('Distance of Closest Distractor');
    ylabel('Accuracy');
    axes1.YAxis(1).Color = [0 0 0];
    legend('off')

end

end
