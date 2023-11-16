% for extended notes see:
% Y:\PROJECTS\Causal Globs\log\2023\01\20210108_Calibration.md

clear, clc, close all
    
%% Making test measurements with the paradigm

retval = PR655init('COM3') % connect to PR655 (check device manager to find out which COM port the PR is connected to)

figure, hold on

%%

% a whole if statement so that I can run this code regardless of whether it
% is the first time or restarting after a "Bad return code -2 from meter"

while true
    tic
    
    if ~exist('XYZ','var')
        n = 1;
    else
        n = size(XYZ,2)+1;
    end
    
    XYZ(:,n) = PR655measxyz;
    SPD(:,n) = PR655measspd([380 4 101]);
    disp(n)
    
    plot(SPD(:,n))
    drawnow
    
    save StimuliMeasurements20230108
    
    pause(68-toc)
end


