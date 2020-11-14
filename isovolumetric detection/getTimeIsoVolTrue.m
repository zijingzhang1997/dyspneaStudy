
% True Iso Vol time based on script settings of detIsoVol.m, Refer to
% analysis 3 data for operating that.

% Reference time is 320s (Update this based on detIsoVol.m)
function [tIsoVolTrue,tCalib] = getTimeIsoVolTrue(caseNum)
tCalib = [1,5]; % Calibration time by default, can be overwritten.
switch(caseNum)
    case 7
        tIsoVolTrue = [22.4, 29.57;39.36,44.85];
    case 8 
        tIsoVolTrue = [24,31; 39,47.11]; 
    case 9
        tIsoVolTrue = [7.18,16.09; 22.2, 30.05; 39.13,46.34]; tCalib = [32.6,37];
    case 10
        tIsoVolTrue = [6.2,10.99; 23.35, 27.15; 39.53, 43.75];
    case 11
        tIsoVolTrue = [13.3,19.5; 22.88, 28.1; 39.89, 45.65];
    case 12
%         tIsoVolTrue = [7.06,11.87; 22.65, 27.78; 39.95, 44.38];
        tIsoVolTrue = [7.06,12.6; 22.65, 28.7; 39.2, 44.38];
    case 13
        tIsoVolTrue = [5.92,9.742; 21.87, 27.91; 38.54, 45.59];
    % Case 14 Ignored
    case 15
        tIsoVolTrue = [6.77,15.39; 22, 31.25; 39.68, 48.1];
    % Case 16 Ignored
    case 17
        % Not sure if ther eare any Isovolumetric cases at all.. some
        % random motion and some breath-hold like periods
        tIsoVolTrue = [.588,11.89; 21.7,28.51; 34.63, 37.28];
        tCalib = [41,51];
    case 18
        % There is no clear isovolumetric motion here.
       tIsoVolTrue = [7.03,12.72; 23.42, 27.25; 40.5, 44.63];
%        tIsoVolTrue = [];
       tCalib=[13.5,22.8];
    case 19
        tIsoVolTrue = [10.24,13.55; 23.29,30.6;39.89,45.17];
    case 20
        % Not sure if IV performed or just breath holds
        tIsoVolTrue = [6, 11.83;21.06,28.44;39.34,44.49];
    case 21
        % Not sure if IV performed or just breath holds
        tIsoVolTrue = [6,10.5; 22, 26.9; 39.1,44.12];
    case 22
        % Clear IV in Biopac, not in NCS
        tIsoVolTrue = [7.35,12.78; 23.17,27.87;38.42,43.29];
    %case 23 Ignored
    case 24 
        % Slight IV in Biopac, not in NCS
        tIsoVolTrue = [5.76,12.28;22.88,27.34;38.18,46.5];
    case 25
        % Clearer in NCS
        tIsoVolTrue = [12.64,16.21;22.28,26.7;37.79,47.16];
    case 26
        tIsoVolTrue = [5.69,13.15;22.88,29.65;38.35,47.75];
    case 27
        tIsoVolTrue = [6.3,13.89;25.13,30.24;39.35,47.62];
    case 28
        tIsoVolTrue = [9.408,13.47; 22.94, 27.71; 39.62,46.2];
    case 29
        tIsoVolTrue = [6.2,14;21.28,29.12;38.99,45.64];
    %case 30 Ignored
    otherwise
        fprintf('Enter valid caseNum for isovol indicator \n');
end

end