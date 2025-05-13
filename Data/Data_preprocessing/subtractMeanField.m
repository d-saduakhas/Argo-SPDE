% In this file, we substract the monthly average only
% fix the monthly average to 15th day of the month
close all; clear;
presLevel = 300;
windowSizeMargined = 5;
load(['./Results/interpolated_',num2str(presLevel),'.mat']);
betaLength = 18;
targetMonth = [1];
load(['./Results/meanField',num2str(presLevel),'_','_w',num2str(windowSizeMargined),'_','month',num2str(targetMonth),'_', num2str(betaLength),'.mat']);

% %% Uncomment below for specific year
% targetYear = 2012;
% filterIdx = interpYear == targetYear;
% interpLat = interpLat(filterIdx);
% interpLong = interpLong(filterIdx);
% interpTemp = interpTemp(filterIdx);
% interpPsal = interpPsal(filterIdx);
% interpJulDay = interpJulDay(filterIdx);
% nProf = length(interpLat);

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
disp(nProf);
% Define the grid
latGrid = -89.5:1:89.5;
longGrid = 20.5:1:379.5;
% latGrid = linspace(-64.5,79.5,145); %RG mask
% longGrid = linspace(20.5,379.5,360);

differenceTemp = zeros(nProf, 1);
differencePsal = zeros(nProf, 1);
monthlyMeanTemp = zeros(size(latGrid, 2), size(longGrid, 2));
monthlyMeanPsal = zeros(size(latGrid, 2), size(longGrid, 2));

for iProf = 1:nProf
    if ~mod(iProf, floor(nProf/20))
        disp([int2str(iProf), '/', int2str(nProf)]);
    end

    % Find the closest grid point for the current profile
    [~, iLat] = min(abs(latGrid - profLatRounded(iProf)));
    [~, iLong] = min(abs(longGrid - profLongRounded(iProf)));

    differenceTemp(iProf) = interpTemp(iProf) - meanGridTemp(iLat,iLong);
    differencePsal(iProf) = interpPsal(iProf) - meanGridPsal(iLat,iLong);
end

nInterp = sum(~isnan(differenceTemp));
disp(nInterp);

% Plotting the mean temperature for each profile
plotWorldMap(latGrid, longGrid, meanGridTemp, ['meanField',num2str(presLevel),'_',num2str(targetMonth,'%02d'),'_','Temp_w',num2str(windowSizeMargined)]);

% Plotting the mean salinity for each profile
plotWorldMap(latGrid, longGrid, meanGridPsal, ['meanField',num2str(presLevel),'_',num2str(targetMonth,'%02d'),'_','Psal_w',num2str(windowSizeMargined)]);


% Save the difference data
save(['./Results/residual_', num2str(presLevel), '_',num2str(targetMonth,'%02d'),'.mat'], 'interpLat', "interpLong", "interpTemp","interpPsal", "differencePsal", "differenceTemp","meanGridPsal","meanGridTemp", "interpJulDay", "interpFloatID", "interpYear")


%% Annual residuals
startYear = 2007;
endYear = 2020;
for iYear = startYear:endYear
    for plotType = {'Temp', 'Psal'}
        if strcmp(plotType, 'Temp')
            plotTarget = differenceTemp;
            titleSuffix = 'Temperature';
        else
            plotTarget = differencePsal;
            titleSuffix = 'Salinity';
        end
        mask = (interpYear == iYear & ~isnan(plotTarget));

        interpLatYear = interpLat(mask);
        interpLongYear = interpLong(mask);
        interpFloatIDYear = interpFloatID(mask);
        interpJulDayYear = interpJulDay(mask);
        interpResYear = plotTarget(mask);

        %         cLimit = max(abs(quantile(interpResYear,[0.01 0.99])));

        figure;
        handle = worldmap('World');
        %         setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coastlines.mat
        plotm(coastlat, coastlon,'k')
        scatterm(interpLatYear,interpLongYear,[],interpResYear, '.');
        h = colorbar;
        %title([num2str(presLevel),' db, ',num2str(targetMonth),'/',num2str(iYear)]);
        if strcmp(plotType, 'Temp')

            h.Label.String = 'Temperature anomaly (°C)';

            switch presLevel
                case 10
                    caxis([-3.7,2.7]);
                case 300
                    caxis([-1.1,1.1]);
                case 1000
                    caxis([-0.5,0.5]);
            end

        else
            h.Label.String = 'Salinity anomaly (psu)';
            switch presLevel
                case 10
                    caxis([-0.4,0.4]);
                case 300
                    caxis([-0.15,0.15]);
                case 1000
                    caxis([-0.09,0.09]);
            end
        end
        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units'))
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2','-r330',['./Figures/residuals_',titleSuffix,'_',num2str(presLevel),'_',num2str(targetMonth,'%02d'),'_',num2str(iYear),'.eps']);
        print('-dpng','-r330',['./Figures/residuals_',titleSuffix,'_',num2str(presLevel),'_',num2str(targetMonth,'%02d'),'_',num2str(iYear),'.png']);
        save(['./Results/residuals_',titleSuffix,'_',num2str(presLevel),'_',num2str(targetMonth,'%02d'),'_',num2str(iYear),'.mat'],'interpResYear','interpLatYear','interpLongYear','interpFloatIDYear','interpJulDayYear');
    end
end

% Round to the nearest half integer
function r = roundHalf(x)
    r = round(x-0.5)+0.5;
end


function plotWorldMap(lat, long, values, titleStr)
    figure;
    worldmap('World');
    hold on;
    load coastlines.mat;
    plotm(coastlat, coastlon, 'k'); % Plot coastlines
    pcolorm(lat, long, values);
    %scatterm(lat, long, [], values, 'x');
    colorbar;
    title(titleStr);
    
    % Save the figure in the Figures folder
    if ~exist('Figures', 'dir')
        mkdir('Figures');
    end
    saveas(gcf, ['./Figures/', titleStr, '.png']);
    close(gcf); % Close the figure after saving
end
