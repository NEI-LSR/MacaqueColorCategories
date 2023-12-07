%% Figures: Bias by Cue with Category Crossings and Confidence Intervals

% Need:
% 1. moving bias, lower_95 (smoothed), and upper_95 (smoothed) for each
% monkey

figure('WindowState', 'maximized');

bigN = 64;
hue_angle = (0:360/bigN:360)';

hold on

% for i = 1:length(monkey);
% [moving_bias, lower_95, upper_95] = categorybias(cleandata);
% end

plot(hue_angle,moving_bias_pollux([1:end 1]),'k');
plot(hue_angle,moving_bias_castor([1:end 1]),'Color', [106/255 102/255 102/255]);
plot(hue_angle,moving_bias_bustercomb([1:end 1]),'Color', [106/255 102/255 102/255]);

plot(hue_angle, lower_95_pollux([1:end 1]), 'k');
plot(hue_angle, upper_95_pollux([1:end 1]), 'k');
h = fill([hue_angle, fliplr(hue_angle)], [lower_95_pollux([1:end 1])', fliplr(upper_95_pollux([1:end 1])')], 'k');
set(h, 'facealpha', 1);

plot(hue_angle, lower_95_castor([1:end 1]), 'Color', [106/255 102/255 102/255]);
plot(hue_angle, upper_95_castor([1:end 1]), 'Color', [106/255 102/255 102/255]);
h = fill([hue_angle, fliplr(hue_angle)], [lower_95_castor([1:end 1])', fliplr(upper_95_castor([1:end 1])')], [106/255 102/255 102/255]);
set(h, 'facealpha', 1);

plot(hue_angle, lower_95_bustercomb([1:end 1]), 'Color',[51/255 50/255 49/255]);
plot(hue_angle, upper_95_bustercomb([1:end 1]), 'Color',[51/255 50/255 49/255]);
h = fill([hue_angle, fliplr(hue_angle)], [lower_95_bustercomb([1:end 1])', fliplr(upper_95_bustercomb([1:end 1])')], [51/255 50/255 49/255]);
set(h, 'facealpha', 1);

yticks([-50 -25 0 25 50])
yticklabels({'-50','-25','0','25','50'})
xticks([0 60 120 180 240 300 360]);
hold off
ax = gca;
set(gca,'TickDir','out');
ax.TickLength = [0.025 0.025];
ax.FontSize = 10;
xlim([0 360]);
ylim([-50 50]);
yline(0);
xlabel('Hue Angle ({\circ})');
ylabel('Bias ({\circ})');
legend('','','','','','Pollux','','','Castor','','','Buster','');
