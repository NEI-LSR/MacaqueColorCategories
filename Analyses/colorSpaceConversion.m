function indexedLocation = colorSpaceConversion(requestedHueAngles)

% It is assumed that the first value in `requestedHueAngles` will be 0

%%

%clear, clc, close all

%%

% load('.\combined\combined_TCC-0att_fullremap-workspace_230510.mat',...
%     'x');

x = readmatrix('combined_TCC-0att_fullremap-x_230510.csv');

%%

nBig = 64;
x_norm  =  x(1:nBig)*(360 /sum(x(1:nBig))); % normalised to sum to 360
x_cum = [0, cumsum(x_norm(1:end-1))]; % convert to degrees matched to index (first value is fixed at 0)

%% 

if ~exist('requestedHueAngles','var') % if this function is being run without input, use these defaults
    requestedNumberOfPoints = 100;
    [~,samplePoints] = generateStimCols('nBig',requestedNumberOfPoints);
    samplePoints = samplePoints(1,:);
else
    samplePoints = requestedHueAngles;
end

for sp = 1:length(samplePoints) % compute indices for locations
    [~,closestIndexLoc] = min(abs(samplePoints(sp) - x_cum)); 
    if samplePoints(sp) < x_cum(closestIndexLoc)
        indexedLocation(sp,1) = closestIndexLoc - 1; % this is the lower neighbour (the requested sample point falls in between this and the next point)
    else
        indexedLocation(sp,1) = closestIndexLoc;
    end
    if indexedLocation(sp,1) + 1 <= nBig
        diffBetweenNeighbours = x_cum(indexedLocation(sp,1) + 1) - x_cum(indexedLocation(sp,1));
        lowerDifference = samplePoints(sp) - x_cum(indexedLocation(sp,1));
        indexedLocation(sp,2) = lowerDifference/diffBetweenNeighbours;
    else
        diffBetweenNeighbours = 360 - x_cum(indexedLocation(sp,1));
        lowerDifference = samplePoints(sp) - x_cum(indexedLocation(sp,1));
        indexedLocation(sp,2) = lowerDifference/diffBetweenNeighbours;
    end
end

end
