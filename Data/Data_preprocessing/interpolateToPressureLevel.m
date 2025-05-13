close all; clear;

%  presLevel = 10; minPresInterp = 5; maxPresInterp = 15; presMask=2;
% presLevel = 300; minPresInterp = 280; maxPresInterp = 320; presMask =25;
presLevel = 1000; minPresInterp = 950; maxPresInterp = 1050;presMask =44;

load('./Data/Argo_data_2007_2020_jan_dec.mat');

%% Spatial filtering
% Load the Roemmich-Gilson land mask
latGrid = linspace(-65,80,146);
longGrid = linspace(20,380,361);
mask = ncread('./Data/RG_climatology/RG_ArgoClim_Temperature_2019.nc','BATHYMETRY_MASK',[1 1 presMask],[Inf Inf 1]); % Load land mask from the RG climatology
mask(mask == 0) = 1;

% Filter using the mask
idx = zeros(nProf,1);
for iProf = 1:nProf
    if profLatAggr(iProf) < -65 || profLatAggr(iProf) >= 80
        continue;
    end
    iLat = find(profLatAggr(iProf) >= latGrid,1,'last');
    iLong = find(profLongAggr(iProf) >= longGrid,1,'last');
    idx(iProf) = ~isnan(mask(iLong,iLat));
end
idx = logical(idx);

% Filter your data using the mask
profPresAggr = profPresAggr(idx);
profTempAggr = profTempAggr(idx);
profPsalAggr = profPsalAggr(idx);
profYearAggr = profYearAggr(idx);
profLatAggr = profLatAggr(idx);
profLongAggr = profLongAggr(idx);
profJulDayAggr = profJulDayAggr(idx);
profFloatIDAggr = profFloatIDAggr(idx);
profCycleNumberAggr = profCycleNumberAggr(idx);
nProf = sum(idx);


%% Interpolation
interpTemp_linear = NaN(nProf,1);
interpPsal_linear = NaN(nProf,1);

% Counters
counter_regular = 0;
counter_deeper_extrap = 0;
counter_shallower_extrap = 0;

for iProf = 1:nProf
    maxPres = max(profPresAggr{iProf});
    minPres = min(profPresAggr{iProf});
    
    if minPres <= presLevel && maxPres >= presLevel
        % Regular linear interpolation
        interpTemp_linear(iProf) = interp1(profPresAggr{iProf}, profTempAggr{iProf}, presLevel);
        interpPsal_linear(iProf) = interp1(profPresAggr{iProf}, profPsalAggr{iProf}, presLevel);
        counter_regular = counter_regular + 1;
%     elseif maxPres >= minPresInterp && maxPres < presLevel 
%         % Linear extrapolation at the deeper end
%         % Find the two nearest points to presLevel
%         [~, idx] = min(abs(profPresAggr{iProf} - presLevel));
%         if idx == length(profPresAggr{iProf})
%             idx = idx - 1;
%         end
%         x = profPresAggr{iProf}(idx:idx+1);
%         yTemp = profTempAggr{iProf}(idx:idx+1);
%         yPsal = profPsalAggr{iProf}(idx:idx+1);
%         % Linear extrapolation formula
%         interpTemp_linear(iProf) = yTemp(1) + (presLevel - x(1)) * (yTemp(2) - yTemp(1)) / (x(2) - x(1));
%         interpPsal_linear(iProf) = yPsal(1) + (presLevel - x(1)) * (yPsal(2) - yPsal(1)) / (x(2) - x(1));
%         counter_deeper_extrap = counter_deeper_extrap + 1;
% % Uncomment for presLevel 10
%     elseif minPres <= maxPresInterp && minPres > presLevel 
%         % Linear extrapolation at the shallower end
%         % Find the two nearest points to presLevel
%         [~, idx] = min(abs(profPresAggr{iProf} - presLevel));
%         if idx == 1
%             idx = idx + 1;
%         end
%         x = profPresAggr{iProf}(idx-1:idx);
%         yTemp = profTempAggr{iProf}(idx-1:idx);
%         yPsal = profPsalAggr{iProf}(idx-1:idx);
%         % Linear extrapolation formula
%         interpTemp_linear(iProf) = yTemp(1) + (presLevel - x(1)) * (yTemp(2) - yTemp(1)) / (x(2) - x(1));
%         interpPsal_linear(iProf) = yPsal(1) + (presLevel - x(1)) * (yPsal(2) - yPsal(1)) / (x(2) - x(1));
%         counter_shallower_extrap = counter_shallower_extrap + 1;
    end
end
nInterp = sum(~isnan(interpTemp_linear));

% Summary
disp(['Total available profiles: ', num2str(nProf)]);
disp(['Regular interpolation: ', num2str(counter_regular)]);
disp(['Deeper end extrapolation: ', num2str(counter_deeper_extrap)]);
disp(['Shallower end extrapolation: ', num2str(counter_shallower_extrap)]);
disp(['Final dataset size after interpolation/extrapolation: ', num2str(nInterp)]);

mask = find(~isnan(interpTemp_linear)); % Find non-NaN interpolated values

interpYear = profYearAggr(mask)';
interpJulDay = profJulDayAggr(mask)';
interpLat = profLatAggr(mask)';
interpLong = profLongAggr(mask)';
interpFloatID = profFloatIDAggr(mask)';
interpTemp = interpTemp_linear(mask);
interpPsal = interpPsal_linear(mask);
save(['./Results/interpolated_',num2str(presLevel),'.mat'],'interpYear','interpJulDay','interpLat','interpLong','interpFloatID','interpTemp','interpPsal','nInterp','startYear','endYear');
% 
% %%------
% 
% interpTemp_pchip = NaN(nProf,1);
% interpPsal_pchip = NaN(nProf,1);
% 
% % Counters
% counter_regular = 0;
% counter_deeper_extrap = 0;
% counter_shallower_extrap = 0;
% 
% for iProf = 1:nProf
%     maxPres = max(profPresAggr{iProf});
%     minPres = min(profPresAggr{iProf});
%     
%     if minPres <= presLevel && maxPres >= presLevel
%         % Regular pchip interpolation
%         interpTemp_pchip(iProf) = interp1(profPresAggr{iProf}, profTempAggr{iProf}, presLevel, 'pchip');
%         interpPsal_pchip(iProf) = interp1(profPresAggr{iProf}, profPsalAggr{iProf}, presLevel, 'pchip');
%         counter_regular = counter_regular + 1;
%     elseif maxPres >= minPresInterp && maxPres < presLevel
%         % Pchip extrapolation at the deeper end
%         profPresAggr{iProf} = [profPresAggr{iProf}; presLevel];
%         profTempAggr{iProf} = [profTempAggr{iProf}; NaN];
%         profPsalAggr{iProf} = [profPsalAggr{iProf}; NaN];
%         interpTemp_pchip(iProf) = interp1(profPresAggr{iProf}, profTempAggr{iProf}, presLevel, 'pchip', 'extrap');
%         interpPsal_pchip(iProf) = interp1(profPresAggr{iProf}, profPsalAggr{iProf}, presLevel, 'pchip', 'extrap');
%         counter_deeper_extrap = counter_deeper_extrap + 1;
% %     elseif minPres <= maxPresInterp && minPres > presLevel
% %         % Pchip extrapolation at the shallower end
% %         profPresAggr{iProf} = [presLevel; profPresAggr{iProf}];
% %         profTempAggr{iProf} = [NaN; profTempAggr{iProf}];
% %         profPsalAggr{iProf} = [NaN; profPsalAggr{iProf}];
% %         interpTemp_pchip(iProf) = interp1(profPresAggr{iProf}, profTempAggr{iProf}, presLevel, 'pchip', 'extrap');
% %         interpPsal_pchip(iProf) = interp1(profPresAggr{iProf}, profPsalAggr{iProf}, presLevel, 'pchip', 'extrap');
% %         counter_shallower_extrap = counter_shallower_extrap + 1;
%     end
% end
% 
% nInterp = sum(~isnan(interpTemp_pchip));
% 
% % Summary
% disp(['Total available profiles: ', num2str(nProf)]);
% disp(['Regular interpolation: ', num2str(counter_regular)]);
% disp(['Deeper end extrapolation: ', num2str(counter_deeper_extrap)]);
% disp(['Shallower end extrapolation: ', num2str(counter_shallower_extrap)]);
% disp(['Final dataset size after interpolation/extrapolation: ', num2str(nInterp)]);
% 
% mask = find(~isnan(interpTemp_pchip)); % Find non-NaN interpolated values
% 
% interpYear_pchip = profYearAggr(mask)';
% interpJulDay_pchip = profJulDayAggr(mask)';
% interpLat_pchip = profLatAggr(mask)';
% interpLong_pchip = profLongAggr(mask)';
% interpFloatID_pchip = profFloatIDAggr(mask)';
% interpTemp_pchip = interpTemp_pchip(mask);
% interpPsal_pchip =interpPsal_pchip(mask);
% 
% % Statistical comparison
% diffTemp = interpTemp_linear - interpTemp_pchip;
% diffPsal = interpPsal_linear - interpPsal_pchip;
% MSE_Temp = mean(diffTemp.^2);
% MAE_Temp = mean(abs(diffTemp));
% RMSE_Temp = sqrt(MSE_Temp);
% 
% MSE_Psal = mean(diffPsal.^2);
% MAE_Psal = mean(abs(diffPsal));
% RMSE_Psal = sqrt(MSE_Psal);
% 
% disp('--- Temperature Comparison ---');
% disp(['Mean Squared Error: ', num2str(MSE_Temp)]);
% disp(['Mean Absolute Error: ', num2str(MAE_Temp)]);
% disp(['Root Mean Squared Error: ', num2str(RMSE_Temp)]);
% 
% disp('--- Salinity Comparison ---');
% disp(['Mean Squared Error: ', num2str(MSE_Psal)]);
% disp(['Mean Absolute Error: ', num2str(MAE_Psal)]);
% disp(['Root Mean Squared Error: ', num2str(RMSE_Psal)]);
% 
% % Histogram of differences
% figure;
% subplot(1, 2, 1);
% hist(diffTemp, 50);
% xlabel('Difference (linear - pchip)');
% ylabel('Number of Profiles');
% title('Histogram of Temperature Differences');
% 
% subplot(1, 2, 2);
% hist(diffPsal, 50);
% xlabel('Difference (linear - pchip)');
% ylabel('Number of Profiles');
% title('Histogram of Salinity Differences');
% 
% % Spatial distribution of differences
% mask = find(~isnan(diffTemp));
% interpLat = profLatAggr(mask);
% interpLong = profLongAggr(mask);
% 
% figure;
% worldmap('World');
% scatterm(interpLat_linear, interpLong_linear-20, [], diffTemp(mask), 'filled');
% colorbar;
% caxis(quantile(diffTemp,[0.01,0.99]));
% title('Spatial Distribution of Temperature Differences');
% 
% figure;
% worldmap('World');
% scatterm(interpLat_linear, interpLong_linear-20, [], diffPsal(mask), 'filled');
% colorbar;
% caxis(quantile(diffPsal,[0.01,0.99]));
% title('Spatial Distribution of Salinity Differences');



%%------
% %% Uncomment section below to see the comparison between the pchip and linear interpolation methods
% 
% % Interpolate raw data to a given pressure level
% % Interpolating Temperature and Salinity using both pchip and linear methods
% 
% close all;
% clear;
% 
% month = 1;
% presLevel = 1000;
% load('./Data/Argo_data_2007_2020_jan_dec.mat');
% interpTemp_linear = zeros(nProf,1);
% interpPsal_linear = zeros(nProf,1);
% interpTemp_pchip = zeros(nProf,1);
% interpPsal_pchip = zeros(nProf,1);
% 
% for iProf = 1:nProf
%     interpTemp_linear(iProf) = interp1(profPresAggr{iProf},profTempAggr{iProf},presLevel, 'linear');
%     interpPsal_linear(iProf) = interp1(profPresAggr{iProf},profPsalAggr{iProf},presLevel, 'linear');
%     interpTemp_pchip(iProf) = interp1(profPresAggr{iProf},profTempAggr{iProf},presLevel, 'pchip');
%     interpPsal_pchip(iProf) = interp1(profPresAggr{iProf},profPsalAggr{iProf},presLevel, 'pchip');
% end
% nInterp = sum(~isnan(interpTemp_linear));
% nInterp = sum(~isnan(interpTemp_pchip));
% disp(nInterp);
% 
% mask = find(~isnan(interpTemp_linear)); % Find non-NaN interpolated values
% interpLat = profLatAggr(mask);
% interpLong = profLongAggr(mask);
% 
% % Calculate the difference between the two interpolation methods
% diffTemp = interpTemp_pchip - interpTemp_linear;
% diffPsal = interpPsal_pchip - interpPsal_linear;
% 
% idx = 1:nInterp;
% % Plot the differences for Salinity
% figure;
% handle = worldmap('World');
% setm(handle, 'Origin', [0 200 0]);
% tightmap;
% load coastlines.mat
% plot(coastlat, coastlon,'k')
% mlabel('off');
% plabel('on');
% 
% scatterm(interpLat(idx), interpLong(idx), [], diffPsal(idx), 'x');
% colorbar;
% title('Difference in Salinity between pchip and linear interpolation');
% % caxis(quantile(diffPsal,[0.1,0.95]));
% 
% % Plot the differences for Temperature
% figure;
% handle = worldmap('World');
% setm(handle, 'Origin', [0 200 0]);
% tightmap;
% load coastlines.mat
% plot(coastlat, coastlon,'k')
% mlabel('off');
% plabel('on');
% 
% scatterm(interpLat(idx), interpLong(idx), [], diffTemp(idx), 'x');
% colorbar;
% title('Difference in Temperature between pchip and linear interpolation');
% % caxis(quantile(diffTemp,[0.1,0.95]));



% %% Uncomment below to analyze the specific profiles to see the differnce between pchip and linear
% % Load your data
% load('./Data/Argo_data_2007_2020_jan_dec.mat');
% 
% % Interpolation at the desired pressure level
% presLevel = 300;
% interpTemp_linear = zeros(nProf, 1);
% interpTemp_pchip = zeros(nProf, 1);
% 
% for iProf = 1:nProf
%     interpTemp_linear(iProf) = interp1(profPresAggr{iProf}, profTempAggr{iProf}, presLevel, 'linear');
%     interpTemp_pchip(iProf) = interp1(profPresAggr{iProf}, profTempAggr{iProf}, presLevel, 'pchip');
% end
% 
% % Find profiles where linear interpolation resulted in NaN but pchip did not
% idx_diff = find(isnan(interpTemp_linear) & ~isnan(interpTemp_pchip));
% 
% % Sample a subset of these profiles for visualization
% num_profiles_to_plot = 5; % Adjust this number as needed
% sample_idx = randsample(idx_diff, min(num_profiles_to_plot, length(idx_diff)));
% 
% figure;
% hold on;
% 
% for i = 1:length(sample_idx)
%     subplot(length(sample_idx), 1, i);
%     plot(profPresAggr{sample_idx(i)}, profTempAggr{sample_idx(i)}, '-o');
%     xline(presLevel, '--r', 'PresLevel');
%     title(['Profile ', num2str(sample_idx(i))]);
%     xlabel('Pressure');
%     ylabel('Temperature');
%     grid on;
% end
% 
% hold off;
