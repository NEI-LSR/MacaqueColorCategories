clear, clc, close all

modelTypes = {'dPrime*','single-ssnu*','single-sg*','sim*'};

NLL = zeros(10,4);
for j = 1:length(modelTypes)
    d = dir(modelTypes{j});
    for i = 1:length(d)
        load(d(i).name,'nll_x')
        NLL(i,j) = nll_x;
    end
end

writetable(array2table(NLL,'VariableNames',modelTypes),...
    ['NLL','_',datestr(now,'yymmdd-HHMMSS'),'.csv'])
