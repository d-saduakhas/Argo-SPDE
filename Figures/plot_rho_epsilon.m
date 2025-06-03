function plot_rho_epsilon(presLevel, i_order, model, model_title, outputDir)
% PLOT_RHO_EPSILON Plot spatial maps of the 'rho_epsilon' parameter.
%
%   plot_rho_epsilon(presLevel, i_order, model, model_title, outputDir)
%
%   Generates a 1x2 grid of spatial maps specifically for the 'rho_epsilon'
%   parameter, typically for the correlated Gaussian and correlated NIG models.
%
%   Inputs:
%   * presLevel   – Scalar pressure level at which to load data.
%   * i_order     – A 1x2 index vector specifying which two models from the 'model'
%                   cell array contain the 'rho_epsilon' parameter to be plotted
%                   (e.g., [2 4] for Gaussian_Correlated and NIG_Correlated).
%   * model       – Cell array of model identifiers (e.g., {'Gaussian', 'Gaussian_Correlated', ...}).
%   * model_title – Cell array of display titles for each model (same length as 'model').
%   * outputDir   – Base directory where the generated figures will be saved.
%
%   Output:
%   A PNG figure is saved to outputDir/presLevel/ named after the pressure level
%   (e.g., '300_rho_epsilon.png'). The figure includes:
%   * Two panels showing the 'rho_epsilon' parameter for the two selected models.
%   * A single common colorbar for 'rho_epsilon' values.
%
%   Dependencies:
%   * Requires the 'colorBrewer' toolbox (specifically `brewermap`) to be
%     available on the MATLAB path.
%   * Relies on helper functions: `get_parameter_grid` (which should be renamed
%     to `makeGridParam` for consistency), `drawWorld`, `addCoast`, and
%     `latexLabel` (which are typically placed in a common 'helpers' directory).
%   * Data files: Reads 'main_results.csv' from `~/Documents/Results/<presLevel>/`
%     and 'grid_equal.mat' from `~/Documents/Results/Data/`.
%--------------------------------------------------------------------------%
% 0. SETTINGS
%======================================================================%
addpath('~/Documents/Results/colorBrewer/');
lat  = linspace(-90,  90,181);
lon  = linspace( 20, 380,361);
[latGrid,lonGrid] = meshgrid(lat,lon);
fs   = 16;                                   % base font size
cmapMod  = brewermap([],'YlOrBr');           % Colormap for model panels

%======================================================================%
% 1. LOAD THE TWO MODEL GRIDS & GET LIMITS
%======================================================================%
fn = ['~/Documents/Results/' num2str(presLevel) '/main_results.csv'];
if ~isfile(fn), error('File not found: %s',fn); end
T  = readtable(fn);
load('~/Documents/Results/Data/grid_equal.mat','Grid');   % 404 boxes

% Load rho_epsilon data for the specified models (model(2) and model(4))
gauss_cor_rho_e = get_parameter_grid(T, model(i_order(1)), 'rho_e', latGrid, lonGrid, Grid); % Model 2 in i_order
nig_cor_rho_e   = get_parameter_grid(T, model(i_order(2)), 'rho_e', latGrid, lonGrid, Grid); % Model 4 in i_order

% Combine all data for global limits calculation
all_rho_e_grids = {gauss_cor_rho_e, nig_cor_rho_e};
all_vals_combined = cell2mat(cellfun(@(x) x(:), all_rho_e_grids, 'UniformOutput', false));
all_vals_combined = all_vals_combined(~isnan(all_vals_combined));

% Determine common color limits for these plots (similar to get_clims_for_parameter for models)
if isempty(all_vals_combined)
    cLimMod = [0 1]; % Default fallback
else
    cLimMod = quantile(all_vals_combined, [0.01, 0.99]); % Using general quantile for rho_e
end

% Ensure cLimMod is not empty or singular
if isempty(cLimMod) || cLimMod(1) == cLimMod(2)
    unique_vals = unique(all_vals_combined);
    if length(unique_vals) == 1
        cLimMod = [unique_vals(1) - 0.1, unique_vals(1) + 0.1];
    else
        cLimMod = [min(all_vals_combined), max(all_vals_combined)];
    end
end
if isempty(cLimMod) || cLimMod(1) == cLimMod(2), cLimMod = [0 1]; end % Final fallback


fprintf('Global Climits for rho_epsilon: [%.4f, %.4f]\n', cLimMod(1), cLimMod(2));

%======================================================================%
% 2. PLOT 1 × 2 GRID (Manual 'axes' Positioning)
%======================================================================%
% Figure dimensions adjusted to keep individual plot width consistent with 2x3 figures
fig_width_cm = 28; % Adjusted width for 2 columns
fig_height_cm = 8; % Adjusted height for 1 row

fig = figure('Units','centimeters', 'OuterPosition',[0 0 1 1]);
fig.Position(3) = fig_width_cm;
fig.Position(4) = fig_height_cm;

% --- Define Layout Parameters (adapted from plot_correlation_common_cb) ---
numRows = 1;
numCols = 2;
left_margin = 0.02;     % Adjusted margin
right_margin = 0.08;    % More space on right for CB labels
bottom_margin = 0.08;   % More space for potential x-labels if needed, and to lift plots
top_margin = 0.08;      % More space for titles

cbar_width = 0.025;     % **Consistent colorbar width**
gap_col = 0.02;         % Gap between the two columns

total_plot_area_width = 1 - left_margin - right_margin - cbar_width;
plot_width = (total_plot_area_width - gap_col) / numCols;
plot_height = 1 - top_margin - bottom_margin;

left_col1 = left_margin;
left_col2 = left_col1 + plot_width + gap_col;

bottom_row1 = bottom_margin; % Only one row

positions = {
    [left_col1, bottom_row1, plot_width, plot_height]; % ax(1)
    [left_col2, bottom_row1, plot_width, plot_height]; % ax(2)
    };

maps   = {gauss_cor_rho_e, nig_cor_rho_e};
cmaps  = {cmapMod, cmapMod};
clims  = {cLimMod, cLimMod};
titles = {model_title(i_order(1)), model_title(i_order(2))}; % Titles for Gaussian with rho_e and NIG with rho_e

ax = gobjects(1,2);
for i = 1:2
    ax(i) = axes('Position', positions{i}); % Use axes for precise position
    drawWorld(latGrid,lonGrid,maps{i},cmaps{i},clims{i},titles{i},fs);
end
drawnow;

%======================================================================%
% 3. COMMON COLOUR-BAR (scaled for single row)
%======================================================================%
% Define scaling factor for the common colorbar
scale = 0.8; % 80% of the full height for a single row

% GAP from map to bar (in figure units)
gapX = 0.005;

% Calculate the full vertical span that colorbar should cover
posAny = get(ax(1),'Position'); % Use any plot position for reference
full_span_bottom = posAny(2);
full_span_top    = posAny(2) + posAny(4);
barFullH = full_span_top - full_span_bottom;

% Calculate new height and bottom position for scaled and centered colorbar
newH = scale * barFullH;
newY = full_span_bottom + (barFullH - newH)/2;

% ---------- COMMON colour-bar -----------------------
posLastPlot = get(ax(2),'Position'); % Position of the last plot in the row
cb = colorbar(ax(2),'Location','EastOutside'); % Attach to one of the axes
cb.Position = [posLastPlot(1) + posLastPlot(3) + gapX, ... % x-position
                newY, ...                                  % y-position (scaled and centered)
                cbar_width, ...                            % width (consistent value)
                newH];                                     % height (scaled)
cb.FontSize = fs-4;
ylabel(cb,['\rho_',char(949)],'Rotation',0,'FontSize',fs, ...
               'VerticalAlignment','middle', 'Interpreter', 'tex');
yl = cb.Label; yl.Units = 'data'; yl.Position(1) = 2.6; % **Nudge label right for consistency**

% Ensure All Plots Use Common Limits
clim(ax(1), cLimMod);
clim(ax(2), cLimMod);

%======================================================================%
% 4. SAVE
%======================================================================%
outDir = fullfile(outputDir,num2str(presLevel));
if ~exist(outDir,'dir'), mkdir(outDir); end
fname = fullfile(outDir,sprintf('%d_rho_epsilon.png',presLevel));
print(fig,'-dpng','-r330',fname);
fprintf('Figure saved to %s\n',fname);
end