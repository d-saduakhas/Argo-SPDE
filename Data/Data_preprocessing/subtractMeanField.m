% subtractMeanField.m 
% =========================================
% Subtracts mean field using all beta coefficients at actual observation locations
% Saves estimated means for later reconstruction (CV, potential density calculations)
%
% OUTPUT VARIABLES SAVED:
%   - differenceTemp, differencePsal: residuals (observation - mean)
%   - estimatedMeanTemp, estimatedMeanPsal: mean at each observation's location/time
%   - profYearDayRatio: the time ratio used for each observation
%   - interpLat, interpLong, interpTemp, interpPsal, interpJulDay, interpFloatID, interpYear
%   - meanGridTemp, meanGridPsal: grid-center means (for backward compatibility/visualization)
%   - useActualDate: flag recording which method was used
%
% RECONSTRUCTION:
%   original_temp = estimatedMeanTemp + differenceTemp  (should equal interpTemp)
%   predicted_temp = estimatedMeanTemp + predicted_anomaly_from_model

close all; clear;

%% =====================================================
%% CONFIGURATION
%% =====================================================
presLevel = 10;           % Pressure level: 10 or 1000
windowSizeMargined = 5;    % Must match estimateMeanField.m
betaLength = 18;           % Must match estimateMeanField.m
targetMonth = [1];         % Month to process (1 = January)

% OPTION: Set to true to use actual observation date for seasonal terms
%         Set to false to use fixed 15th of month
useActualDate = false;

%% =====================================================
%% LOAD DATA
%% =====================================================
load(['./Results/interpolated_',num2str(presLevel),'.mat']);
load(['./Results/meanField',num2str(presLevel),'_','_w',num2str(windowSizeMargined),'_','month',num2str(targetMonth),'_', num2str(betaLength),'.mat']);

% Subset the data for specific month/s
dateAggr = datevec(interpJulDay);

filterIdx = ismember(dateAggr(:,2), targetMonth);
interpLat = interpLat(filterIdx);
interpLong = interpLong(filterIdx);
interpTemp = interpTemp(filterIdx);
interpPsal = interpPsal(filterIdx);
interpFloatID = interpFloatID(filterIdx);
interpJulDay = interpJulDay(filterIdx);
interpYear = interpYear(filterIdx);
profLatRounded = roundHalf(interpLat);
profLongRounded = roundHalf(interpLong);
nProf = length(profLatRounded);
disp(['Number of profiles: ', num2str(nProf)]);

% Define the grid (MUST match estimateMeanField.m)
latGrid = -89.5:1:89.5;
longGrid = 20.5:1:379.5;

differenceTemp = zeros(nProf, 1);
differencePsal = zeros(nProf, 1);
estimatedMeanTemp = zeros(nProf, 1);
estimatedMeanPsal = zeros(nProf, 1);

% Fixed seasonal ratio for 15th of target month (used if useActualDate=false)
profYearDayRatio_fixed = (15+30*(targetMonth-1))/365;
profYearDayRatio = zeros(nProf, 1);

%% =====================================================
%% MAIN LOOP - SUBTRACT MEAN FIELD
%% =====================================================
for iProf = 1:nProf
    if ~mod(iProf, floor(nProf/20))
        disp([int2str(iProf), '/', int2str(nProf)]);
    end

    % Find the closest grid point for the current profile
    [~, iLat] = min(abs(latGrid - profLatRounded(iProf)));
    [~, iLong] = min(abs(longGrid - profLongRounded(iProf)));

    % Get beta coefficients for this grid cell
    betaTemp = squeeze(betaGridTemp(iLat, iLong, :));
    betaPsal = squeeze(betaGridPsal(iLat, iLong, :));

    % Compute spatial offsets from grid center
    latDiff = interpLat(iProf) - latGrid(iLat);
    longDiff = interpLong(iProf) - longGrid(iLong);

    % Determine which seasonal ratio to use
    if useActualDate
        profYearDay = fromJulDayToYearDay(interpJulDay(iProf));
        profYearLength = yearLength(interpJulDay(iProf));
        profYearDayRatio(iProf) = profYearDay / profYearLength;
    else
        profYearDayRatio(iProf) = profYearDayRatio_fixed;
    end

    % Compute mean at ACTUAL observation location using ALL beta coefficients
    meanTempAtObs = computeMeanAtLocation(betaTemp, latDiff, longDiff, profYearDayRatio(iProf));
    meanPsalAtObs = computeMeanAtLocation(betaPsal, latDiff, longDiff, profYearDayRatio(iProf));

    estimatedMeanTemp(iProf) = meanTempAtObs;
    estimatedMeanPsal(iProf) = meanPsalAtObs;
    differenceTemp(iProf) = interpTemp(iProf) - meanTempAtObs;
    differencePsal(iProf) = interpPsal(iProf) - meanPsalAtObs;
end

nInterp = sum(~isnan(differenceTemp));
disp(['Non-NaN residuals: ', num2str(nInterp)]);

%% =====================================================
%% VERIFICATION
%% =====================================================
reconErrorTemp = max(abs(estimatedMeanTemp + differenceTemp - interpTemp));
reconErrorPsal = max(abs(estimatedMeanPsal + differencePsal - interpPsal));
disp(['Reconstruction error (Temp): ', num2str(reconErrorTemp)]);
disp(['Reconstruction error (Psal): ', num2str(reconErrorPsal)]);
if reconErrorTemp > 1e-10 || reconErrorPsal > 1e-10
    warning('Reconstruction error too large! Check computeMeanAtLocation.');
else
    disp('Reconstruction verified: estimatedMeanTemp + differenceTemp = interpTemp');
end

%% =====================================================
%% SUMMARY STATISTICS
%% =====================================================
disp(' ');
disp('--- Residual Statistics ---');
disp(['Temperature residuals: mean = ', num2str(nanmean(differenceTemp), '%.4f'), ...
      ' C, std = ', num2str(nanstd(differenceTemp), '%.4f'), ' C']);
disp(['Salinity residuals:    mean = ', num2str(nanmean(differencePsal), '%.6f'), ...
      ' psu, std = ', num2str(nanstd(differencePsal), '%.6f'), ' psu']);

%% =====================================================
%% PLOTTING (grid-center means for visualization)
%% =====================================================
plotWorldMap(latGrid, longGrid, meanGridTemp, ...
    ['meanField', num2str(presLevel), '_', num2str(targetMonth,'%02d'), '_Temp_w', num2str(windowSizeMargined)]);
plotWorldMap(latGrid, longGrid, meanGridPsal, ...
    ['meanField', num2str(presLevel), '_', num2str(targetMonth,'%02d'), '_Psal_w', num2str(windowSizeMargined)]);

%% =====================================================
%% SAVE
%% =====================================================
save(['./Results/residual_', num2str(presLevel), '_', num2str(targetMonth,'%02d'), '.mat'], ...
    'interpLat', 'interpLong', 'interpTemp', 'interpPsal', ...
    'differenceTemp', 'differencePsal', ...
    'estimatedMeanTemp', 'estimatedMeanPsal', ...
    'profYearDayRatio', ...
    'meanGridTemp', 'meanGridPsal', ...
    'interpJulDay', 'interpFloatID', 'interpYear', ...
    'useActualDate');

disp(' ');
disp(['Saved to: ./Results/residual_', num2str(presLevel), '_', num2str(targetMonth,'%02d'), '.mat']);

%% =====================================================
%% HELPER FUNCTIONS
%% =====================================================

function meanVal = computeMeanAtLocation(beta, latDiff, longDiff, profYearDayRatio)
% Computes mean field at a specific location using ALL beta coefficients
    betaLen = length(beta);

    if betaLen < 6 || all(beta == 0)
        meanVal = beta(1);
        return;
    end

    % Spatial polynomial terms (beta2-6)
    meanVal = beta(1) + ...
              beta(2) * latDiff + ...
              beta(3) * longDiff + ...
              beta(4) * latDiff * longDiff + ...
              beta(5) * latDiff^2 + ...
              beta(6) * longDiff^2;

    % Seasonal harmonic terms (beta7-18)
    if betaLen >= 18
        for k = 1:6
            meanVal = meanVal + ...
                beta(6 + 2*k - 1) * sin(2*pi*k*profYearDayRatio) + ...
                beta(6 + 2*k)     * cos(2*pi*k*profYearDayRatio);
        end
    end
end

% Round to the nearest half integer
function r = roundHalf(x)
    r = round(x-0.5)+0.5;
end

function r = fromJulDayToYearDay(julDay)
    tempDateVec = datevec(julDay);
    r = datenum(tempDateVec) - datenum([tempDateVec(:,1) repmat([1 1 0 0 0], length(julDay), 1)]);
end

function r = yearLength(julDay)
    tempDateVec = datevec(julDay);
    r = datenum(tempDateVec(:,1), 12, 31, 24, 0, 0) - datenum(tempDateVec(:,1), 1, 1, 0, 0, 0);
end

function plotWorldMap(lat, long, values, titleStr)
    figure;
    worldmap('World');
    hold on;
    load coastlines.mat;
    plotm(coastlat, coastlon, 'k');
    pcolorm(lat, long, values);
    colorbar;
    title(titleStr);

    if ~exist('Figures', 'dir')
        mkdir('Figures');
    end
    saveas(gcf, ['./Figures/', titleStr, '.png']);
    close(gcf);
end