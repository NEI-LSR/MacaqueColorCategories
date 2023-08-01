function [unique_difficulties,pct_correct,...
    difficulties_completed_counts, f]...
    = difficulty_psychometric(cleandata,filename)
%% PURPOSE:
%Error analysis by difficulty
%Get difficulty for each trial, plot errors by difficulty
%Difficulty calculated as distance of closest distractor to cue

%% Process Inputs
cues = cell2mat(cleandata.trialdata.cues);
choices = cell2mat(cleandata.trialdata.choices);
chosen = cell2mat(cleandata.trialdata.chosen);

nBig = size(cleandata.trialdata.stimCols{1,1},1);

%% Calculate difficulty

completed_idx = ~isnan(chosen);
correct_idx = cues(completed_idx) == chosen(completed_idx);

distance = abs(cues(completed_idx) - choices(completed_idx,:));
distance(distance > nBig/2) = nBig - distance(distance > nBig/2);

interval = 360/nBig;
distance = distance * interval;

difficulty = zeros(1,size(cues(completed_idx),1));
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

unique_difficulties = unique(difficulty);

difficulties_completed_counts = zeros(size(unique_difficulties,2),1);
for i = 1:length(unique_difficulties)
    difficulties_completed_counts(i,1) = sum(difficulty == unique_difficulties(i));
end

difficulties_correct_counts = zeros(size(unique_difficulties,2),1);
for i = 1:length(unique_difficulties)
    difficulties_correct_counts(i,1) = sum(difficulty(correct_idx) == unique_difficulties(i));
end

% Percent correct by distance of closest distractor
pct_correct = difficulties_correct_counts./difficulties_completed_counts.*100;

%% Fit Weibull function with 4 params

Weibull = fittype(@(slope, inflect, floor, ceil, x)...
    (floor+(100-floor-ceil).*(1-exp(-(x./slope).^inflect))));

%% Weighted fit
f = fit(unique_difficulties',pct_correct,...
    Weibull,...
    'weights', difficulties_completed_counts,...
    'Lower', [0 -20 0 0],...
    'Upper', [200 20 100 100],...
    'StartPoint', [50 2 30 10]);

%% FIGURES

% figure(1) % Uncomment this to get them all on the same graph
figure,

axes('PositionConstraint','innerposition',...
    'Position',[0.13 0.19 0.82 0.75])

hold on
plot(unique_difficulties, pct_correct, '.','MarkerEdgeColor','none'); % this is invisible, yet neccessary, for reasons I don't understand
xlim([min(unique_difficulties), max(unique_difficulties)]); % this determines the limits within which the function is plotted in the following line

plot(f,'k')

p22 = predint(f,unique_difficulties,0.95,'functional','on');
x_plot =[unique_difficulties, fliplr(unique_difficulties)]; % h/t https://www.mathworks.com/matlabcentral/answers/425206-plot-of-confidence-interval-with-fill
y_plot=[p22(:,1)', flipud(p22(:,2))'];
fill(x_plot, y_plot, 1,'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.2);

xlim([0,180])
xticks([0,interval,45:45:180])
xtickangle(90)

ylim([25 100])
ytickformat('percentage')
yticks([25,50,75,100])
set(gca,'TickDir','out');

% title('Accuracy by Trial Difficulty');
xlabel('Distance of Closest Distractor');
ylabel('Accuracy');
legend('off')

if 1
    startDate = datetime(str2double([num2str(20),(cleandata.trialdata.dirname{1}(1:2))]),...
        str2double(cleandata.trialdata.dirname{1}(3:4)),...
        str2double(cleandata.trialdata.dirname{1}(5:6)));

    endDate = datetime(str2double([num2str(20),(cleandata.trialdata.dirname{end}(1:2))]),...
        str2double(cleandata.trialdata.dirname{end}(3:4)),...
        str2double(cleandata.trialdata.dirname{end}(5:6)));

    ci = confint(f); % TODO Check CI type (simultaneous or not etc.)

    text(100,60,[...
        upper(filename(16:17)),newline...
        'Trials: ',num2str(size(cues,1)),newline...
        'Sessions: ', num2str(length(unique(cleandata.trialdata.dirname))),newline,...
        'Duration: ', char(between(startDate, endDate)), newline,...
        'Ceiling: ', num2str(100-f.ceil','%.1f'),'% (+/- ',num2str(f.ceil-ci(1,4),'%.1f'),')'])

end
%% 

saveas(gcf,['../','difficulty_', filename, ...
    datestr(now,'yymmdd-HHMMSS'), '.svg'])

end
