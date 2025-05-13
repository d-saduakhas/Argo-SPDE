% Assuming betaGridTemp and betaGridPsal are of size [numLat, numLong, betaLength]
% and you have matrices 'interpLat', 'interpLong', 'interpTemp', and 'interpPsal'
%only fit local regression
clear; close all;

presLevel = 300;
idx_pres_RG = 25; %2- 10dbar; 25 -300dbar 44- 1000dbar 
load(['./Results/interpolated_',num2str(presLevel),'.mat']);

% Load land mask
mask = ncread('Data/RG_climatology/RG_ArgoClim_Temperature_2019.nc','BATHYMETRY_MASK',[1 1 idx_pres_RG],[Inf Inf 1]);
mask = [NaN*ones(360,25) mask NaN*ones(360,10)]; 
mask(mask == 0) = 1; 
% compute monthly mean data
january_data = estimateMean(interpLat, interpLong, interpTemp, interpPsal, interpJulDay, 1, 18, mask,presLevel);
february_data= estimateMean(interpLat, interpLong, interpTemp, interpPsal, interpJulDay, 2, 18, mask,presLevel);
march_data= estimateMean(interpLat, interpLong, interpTemp, interpPsal, interpJulDay, 3, 18, mask,presLevel);


function [betaGridTemp, betaGridPsal,meanGridTemp,meanGridPsal] = estimateMean(interpLat, interpLong, interpTemp, interpPsal, interpJulDay, targetMonth, betaLength, mask,presLevel)
% Define the grid and window size
% latGrid = linspace(-64.5,79.5,145);
% longGrid = linspace(20.5,379.5,360);

latGrid = -89.5:1:89.5;
longGrid = 20.5:1:379.5;

windowSizeMargined = 5; % 5-degree window
nGrid = numel(latGrid)* numel(longGrid);

%% FIT REGRESSION
% Initialize the betaGrid for storing regression coefficients
betaGridTemp = zeros(180, 360, betaLength);
betaGridPsal = zeros(180, 360, betaLength);
meanGridTemp = zeros(180, 360);
meanGridPsal = zeros(180, 360);
if betaLength ==6
    nMinObs = 10;
else
    nMinObs = 20;
end

%% Uncomment if want to fit for subset of data
% dateAggr = datevec(interpJulDay);
% filterIdx = ismember(dateAggr(:,2), targetMonth);
% interpLat = interpLat(filterIdx);
% interpLong = interpLong(filterIdx);
% interpTemp = interpTemp(filterIdx);
% interpPsal = interpPsal(filterIdx);
% interpFloatID = interpFloatID(filterIdx);
% interpJulDay = interpJulDay(filterIdx);
% nInterp = length(interpLong);

% Wrap-around for longitude
leftBoundaryIdx = find(interpLong <= 20 + windowSizeMargined);
rightBoundaryIdx = find(interpLong >= 380 - windowSizeMargined);
interpLong = [interpLong; interpLong(leftBoundaryIdx) + 360; interpLong(rightBoundaryIdx) - 360];
interpLat = [interpLat; interpLat(leftBoundaryIdx); interpLat(rightBoundaryIdx)];
interpTemp = [interpTemp; interpTemp(leftBoundaryIdx); interpTemp(rightBoundaryIdx)];
interpPsal = [interpPsal; interpPsal(leftBoundaryIdx); interpPsal(rightBoundaryIdx)];
interpJulDay = [interpJulDay; interpJulDay(leftBoundaryIdx); interpJulDay(rightBoundaryIdx)];

tic
for iGrid = 1:nGrid
    
    if ~mod(iGrid, floor(nGrid/20))
        disp([int2str(iGrid), '/', int2str(nGrid)]);
    end
    
    [iLat, iLong] = ind2sub([numel(latGrid), numel(longGrid)], iGrid);
    latSel = latGrid(iLat);
    longSel = longGrid(iLong);
    
    latMin = latSel - windowSizeMargined;
    latMax = latSel + windowSizeMargined;
    longMin = longSel - windowSizeMargined;
    longMax = longSel + windowSizeMargined;
    
    idx = find(interpLat > latMin & interpLat < latMax & interpLong > longMin & interpLong < longMax);
    
    % Need at least 10(no seasonal) or 20(with seasonal effect) data points to estimate the regression coefficients
    if ( length(idx) < nMinObs ) || isnan(mask(iLong,iLat) )
        if ( length(idx)>=10 )
            betaGridTemp(iLat, iLong, 1) = mean(interpTemp(idx));
            betaGridPsal(iLat, iLong, 1) = mean(interpPsal(idx));
            disp(['lat: ',num2str(latSel),', lon: ',num2str(longSel),' mean '])
            continue;
        else
            betaGrid(iGrid, :) = NaN;
            disp(['lat: ',num2str(latSel),', lon: ',num2str(longSel),' unreliable grid point'])
            continue;
        end
    end
    
    profJulDayAggrWindow = interpJulDay(idx)';
    profYearDayAggrWindow = fromJulDayToYearDay(profJulDayAggrWindow);
    profYearLengthAggrWindow = yearLength(profJulDayAggrWindow);
    profYearDayRatioWindow = profYearDayAggrWindow ./ profYearLengthAggrWindow;
    % Setup Design Matrix
    if betaLength == 6
        X = [ones(length(idx), 1) (interpLat(idx)-latSel) (interpLong(idx)-longSel) ...
            (interpLat(idx)-latSel).*(interpLong(idx)-longSel) ...
            (interpLat(idx)-latSel).^2 (interpLong(idx)-longSel).^2];
    elseif betaLength == 18
        X = [ones(length(idx), 1) (interpLat(idx)-latSel) (interpLong(idx)-longSel) ...
            (interpLat(idx)-latSel).*(interpLong(idx)-longSel) ...
            (interpLat(idx)-latSel).^2 (interpLong(idx)-longSel).^2 ...
            sin(2*pi*1*profYearDayRatioWindow) cos(2*pi*1*profYearDayRatioWindow) ...
            sin(2*pi*2*profYearDayRatioWindow) cos(2*pi*2*profYearDayRatioWindow) ...
            sin(2*pi*3*profYearDayRatioWindow) cos(2*pi*3*profYearDayRatioWindow) ...
            sin(2*pi*4*profYearDayRatioWindow) cos(2*pi*4*profYearDayRatioWindow) ...
            sin(2*pi*5*profYearDayRatioWindow) cos(2*pi*5*profYearDayRatioWindow) ...
            sin(2*pi*6*profYearDayRatioWindow) cos(2*pi*6*profYearDayRatioWindow)];
    end
    
    
    % Fit polynomial regression for temperature
    betaTemp = X \ interpTemp(idx);
    betaPsal = X \ interpPsal(idx);
    
    %% clipping data if no mask used
    %     maxTemp = max(interpTemp(idx))+20;
    %     minTemp = min(interpTemp(idx))-20;
    %     maxPsal = max(interpPsal(idx))+10;
    %     minPsal = max(interpPsal(idx))-10;
    %     if (  betaTemp(1)>maxTemp ||  betaTemp(1)<minTemp || betaPsal(1)>maxPsal || betaPsal(1)<minPsal  )
    %         betaTemp = zeros(betaLength,1);
    %         betaPsal = zeros(betaLength,1);
    %             betaGridTemp(iLat, iLong, 1) = mean(interpTemp(idx));
    %             betaGridPsal(iLat, iLong, 1) = mean(interpPsal(idx));
    %             disp(['lat: ',num2str(latSel),', lon: ',num2str(longSel),' clip '])
    %             continue;
    %     end
    betaGridTemp(iLat,iLong, :) = betaTemp;
    betaGridPsal(iLat,iLong, :) = betaPsal;
%% compute monthly mean
meanGridTemp(iLat,iLong) = computeMean(betaTemp,(15+30*(targetMonth-1))/365);
meanGridPsal(iLat,iLong) = computeMean(betaPsal,(15+30*(targetMonth-1))/365);
end

toc;

% Save the results
saveName = ['./Results/meanField',num2str(presLevel),'_','_w',num2str(windowSizeMargined),'_','month',num2str(targetMonth),'_', num2str(betaLength),'.mat'];
save(saveName,'betaGridTemp','latGrid','longGrid',"meanGridTemp","betaGridPsal","meanGridPsal");

% % Display the Mean Temperature Grid
% figure;
% imagesc(longGrid, latGrid, meanGridTemp);
% axis xy; % To ensure the latitude starts from the bottom
% colorbar;
% title('Mean Temperature');
% xlabel('Longitude');
% ylabel('Latitude');

% % Overlay Coastlines
% load coastlines; % Load coastline data
% hold on;
% plot(coastlon, coastlat, 'k', 'LineWidth', 2); % Plot coastlines in black
% hold off;
return;
end


function meanVal = computeMean(beta, profYearDayRatio)
meanVal = beta(1);
if size(beta,2) == 18
    meanVal = beta(1) + ...
        beta(7) * sin(2*pi*1*profYearDayRatio) + ...
        beta(8) * cos(2*pi*1*profYearDayRatio) + ...
        beta(9) * sin(2*pi*2*profYearDayRatio) + ...
        beta(10) * cos(2*pi*2*profYearDayRatio) + ...
        beta(11) * sin(2*pi*3*profYearDayRatio) + ...
        beta(12) * cos(2*pi*3*profYearDayRatio) + ...
        beta(13) * sin(2*pi*4*profYearDayRatio) + ...
        beta(14) * cos(2*pi*4*profYearDayRatio) + ...
        beta(15) * sin(2*pi*5*profYearDayRatio) + ...
        beta(16) * cos(2*pi*5*profYearDayRatio) + ...
        beta(17) * sin(2*pi*6*profYearDayRatio) + ...
        beta(18) * cos(2*pi*6*profYearDayRatio);
end
end

function r = fromJulDayToYearDay(julDay)

tempDateVec = datevec(julDay);
r = datenum(tempDateVec) - datenum([tempDateVec(:,1) repmat([1 1 0 0 0],length(julDay),1)]);
end
function daysPassed = fromMonthToYearDay(month)
startOfYear = datenum(2012, 1, 1);
middleOfMonth = datenum(2012, month, 15);
daysPassed = middleOfMonth - startOfYear + 1;
end


function r = yearLength(julDay)

tempDateVec = datevec(julDay);
r = datenum(tempDateVec(:,1),12,31,24,0,0)-datenum(tempDateVec(:,1),1,1,0,0,0);
end
