clear, clc, close all

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))

for i = 1%:10
    ParameterEstimator_caller(i)
end

%%

% cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu')
% cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\sg')
% cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\sg_ssnu')
cd('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu_THEN_sg')
d = dir('*.mat');

% https://www.mathworks.com/matlabcentral/answers/397385-how-to-sort-a-structure-array-based-on-a-specific-field#answer_317198
T = struct2table(d);
sortedT = sortrows(T, 'date');
d = table2struct(sortedT);

for i = 1:length(d)
    load(d(i).name)
    all_aic(i) = aic;
    all_nll(i) = nll_x;
    all_x(:,i) = x;
end

%%

fmt = '%.9f';

figure, 
ax1 = axes;
plot(all_aic,'k*')
title('AIC')
ytickformat(fmt)
ax1.YAxis.Exponent = 0;

figure, 
ax2 = axes;
plot(all_nll,'k*')
title('nll')
ytickformat(fmt)
ax2.YAxis.Exponent = 0;

