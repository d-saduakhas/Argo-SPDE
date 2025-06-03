%% PLOT_BEST_MODEL_SCORE
% This script generates spatial maps showing the **best-performing model**
% for Temperature or Salinity scores (CRPS, MAE, RMSE, or SCRPS) at a
% specified pressure level.
%
% It loads cross-validation results, identifies the best model for each grid
% cell, and plots the spatial distribution of these optimal models.
%
% **Outputs**: PNG figures saved in `Figures/<pressureLevel>/`.
%
% **Dependencies**: `colorBrewer`, `readtable`, `meshgrid`, `worldmap`,
% MATLAB's mapping toolbox functions.
% 
% Uncomment required score type in lines 47-61.

close all;
clear;

%======================================================================%
% 0. STYLE SETTINGS
%======================================================================%
fs           = 16;      % Base font size
fig_width_cm = 28.5;      % Figure width in cm
cbar_width   = 0.025;   % Colorbar width (normalized units)
include_colorbar = 0;   % 1 to include, 0 to skip
if include_colorbar
    fig_height_cm = 11;
else
    fig_height_cm = 8.5;
end
presLevel        = 10;
% Set default font size for all axes
set(groot, 'defaultAxesFontSize', fs);

%======================================================================%
% 1. PARAMETERS AND DATA
%======================================================================%


model       = ["gauss_indep", "gauss_cor", "nig_indep", "nig_cor"];
model_title = { ...
    'Gaussian (IID Nugget)', ...
    'Gaussian (Corr. Nugget)', ...
    'NIG (IID Nugget)',      ...
    'NIG (Corr. Nugget)'     ...
};

% variable_figure    = 'SCRPS';
% temp_score_type    = "cv_Temp_SCRPS";
% salinity_score_type= "cv_Psal_SCRPS";

% variable_figure = 'MAE';
% temp_score_type = "cv_Temp_MAE";
% salinity_score_type = "cv_Psal_MAE";

variable_figure = 'CRPS';
temp_score_type = "cv_Temp_CRPS";
salinity_score_type = "cv_Psal_CRPS";
% % 
% variable_figure = 'RMSE';
% temp_score_type = "cv_Temp_RMSE";
% salinity_score_type = "cv_Psal_RMSE";


score_title        = ["Temperature", "Salinity"];

[latGrid, longGrid] = meshgrid(linspace(-90,90,181), linspace(20,380,361));

% Load data
filename = fullfile(getenv('HOME'), 'Documents/Results', num2str(presLevel), 'cv_results.csv');
resultsTableFull = readtable(filename);
addpath(fullfile(getenv('HOME'),'Documents','Results','colorBrewer'));
load(fullfile(getenv('HOME'),'Documents','Results','Data','grid_equal.mat'), 'Grid');

% Add RMSE fields
resultsTableFull.cv_Temp_RMSE  = sqrt(resultsTableFull.cv_Temp_MSE);
resultsTableFull.cv_Psal_RMSE  = sqrt(resultsTableFull.cv_Psal_MSE);

%======================================================================%
% 2. SETUP FIGURE & LAYOUT
%======================================================================%
fig = figure('Units','centimeters');
fig.Position(3) = fig_width_cm;
fig.Position(4) = fig_height_cm;

t = tiledlayout(1,2, 'Padding','tight','TileSpacing','tight');

%======================================================================%
% 3. PROCESS & PLOT FUNCTION (Temperature & Salinity)
%======================================================================%
for varIdx = 1:2
    if varIdx == 1
        scoreType = temp_score_type;
    else
        scoreType = salinity_score_type;
    end
    fprintf('Processing %s with score: %s\n', score_title(varIdx), scoreType);

    % Build varGrid & bestModelGrid
    varGrid = zeros(size(latGrid));
    bestModelGrid = cell(1,404);
    for iGrid = 1:404
        rows = resultsTableFull.gridID == iGrid;
        gridData = resultsTableFull(rows, :);
        if ~isempty(gridData)
            [~, idx_best] = min(gridData.(scoreType));
            bestModelGrid{iGrid} = gridData.model{idx_best};
            val = find(strcmp(model, bestModelGrid{iGrid}));
        else
            val = NaN;
        end
        % fill region
        latMin = Grid(iGrid,1); latMax = Grid(iGrid,2);
        lonMin = Grid(iGrid,3); lonMax = Grid(iGrid,4);
        cols = find(latGrid(1,:) == latMin+1) : find(latGrid(1,:) == latMax);
        rowsIdx = (ceil(lonMin):floor(lonMax))+1;
        varGrid(rowsIdx, cols) = val;
    end

    % Next tile & map
    ax = nexttile(varIdx);
    worldmap('World'); tightmap; gridm off; mlabel('off'); plabel('off');
    colormap(ax, brewermap(4,'Paired')); clim([0.5,4.5]);
    surfm(latGrid, longGrid, varGrid);
    hold on; load coastlines; plotm(coastlat,coastlon,'k'); patchm(coastlat,coastlon,ones(1,3));

    title(score_title(varIdx), 'FontSize', fs, 'FontWeight','bold');
    ax.FontSize = fs;

    % Colorbar if requested (manual normalized position)
    if include_colorbar && (varIdx==2)
        cb = colorbar;
        cb.Orientation = 'horizontal';
        cb.Ticks = 1:4;
        cb.TickLabels = model_title;
        cb.FontSize = fs-2;
        cb.Label.String = '';
        cb.Label.FontSize = fs-4;
        cb.Label.FontWeight = 'bold';
        % Manual position in normalized units
        cb.Position = [0.1, 0.075, 0.8, 0.035];
    end
end

%======================================================================%
% 4. SAVE FIGURE
%======================================================================%
outDir = fullfile(getenv('HOME'),'Documents','Results','Figures', num2str(presLevel));
if ~exist(outDir,'dir'), mkdir(outDir); end

if include_colorbar
    fname = sprintf('%d_%s_best_model_colorbar.png', presLevel, variable_figure);
else
    fname = sprintf('%d_%s_best_model.png', presLevel, variable_figure);
end
print(fig, '-dpng', '-r330', fullfile(outDir, fname));
fprintf('Figure saved to %s\n', fullfile(outDir, fname));

