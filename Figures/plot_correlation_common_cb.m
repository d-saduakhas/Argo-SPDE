function plot_correlation_common_cb( ...
    presLevel, i_order, model, model_title, outputDir)
% PLOT_CORRELATION_COMMON_CB Plot spatial maps for inter-component correlation (Rho).
%
%   plot_correlation_common_cb(presLevel, i_order, model, model_title, outputDir)
%
%   Generates a 2x3 grid of spatial maps to visualize the correlation
%   parameter (often denoted as Rho) for four different models, along with
%   two difference maps showing the effect of correlation.
%
%   Inputs:
%   * presLevel   – Scalar pressure level for data loading.
%   * i_order     – A 1x4 index vector specifying the display order of the four models
%                   from the 'model' cell array.
%   * model       – Cell array of model identifiers (e.g., {'Gaussian', 'Gaussian_Correlated', ...}).
%                   Its length should be at least max(i_order).
%   * model_title – Cell array of display titles for each model (same length as 'model').
%   * outputDir   – Base directory where the generated figures will be saved.
%
%   Output:
%   A PNG figure is saved to outputDir/presLevel/ named after the pressure level
%   (e.g., '300_correlation_combined.png'). The figure includes:
%   * Four panels showing the correlation parameter for each of the four selected models.
%   * Two "difference" panels: (Gaussian Correlated - Gaussian) and (NIG Correlated - NIG).
%   * Two common colorbars: one for model correlation values and one for differences.
%   * The colorbar label for correlation is set to $\rho_{u_1,u_2}$.
%
%   Dependencies:
%   * Requires the 'colorBrewer' toolbox (specifically `brewermap`) and the
%     `redblue()` utility to be available on the MATLAB path.
%   * Relies on helper functions: `makeGrid` (which should be renamed to `makeGridParam`
%     for consistency with other plotting functions), `drawWorld`, `addCoast`,
%     and `biv_corr`. These helpers are typically placed in a common 'helpers' directory.
%   * Data files: Reads 'main_results.csv' from `~/Documents/Results/<presLevel>/`
%     and 'grid_equal.mat' from `~/Documents/Results/Data/`.
%--------------------------------------------------------------------%
% 0. SETTINGS
%--------------------------------------------------------------------%
addpath('~/Documents/Results/colorBrewer/');
lat  = linspace(-90,  90,181);
lon  = linspace( 20, 380,361);
[latGrid,lonGrid] = meshgrid(lat,lon);
fs   = 16;                                   % base font size
cmapMod  = brewermap([],'YlOrBr');           % 4 model panels
cmapDiff = redblue(17);  cmapDiff(9,:) = 1;  % 2 difference panels

%--------------------------------------------------------------------%
% 1. LOAD THE FOUR MODEL GRIDS & GET LIMITS
%--------------------------------------------------------------------%
fn = ['~/Documents/Results/' num2str(presLevel) '/main_results.csv'];
if ~isfile(fn), error('File not found: %s',fn); end
T  = readtable(fn);
load('~/Documents/Results/Data/grid_equal.mat','Grid');   % 404 boxes
varGrid = cell(1,4);
for k = 1:4
    varGrid{k} = makeGrid(T,model(i_order(k)),latGrid,lonGrid,Grid);
end
gauss     = varGrid{1}; gauss_cor = varGrid{2};
nig       = varGrid{3}; nig_cor   = varGrid{4};
diffG     = gauss_cor - gauss;
diffN     = nig_cor   - nig;

% --- Global Limits (Ensured) ---
vals    = [gauss(:);gauss_cor(:);nig(:);nig_cor(:)];
cLimMod = quantile(vals(~isnan(vals)),[.05 .95]);
if isempty(cLimMod) || cLimMod(1) == cLimMod(2), cLimMod = [min(vals(~isnan(vals))), max(vals(~isnan(vals)))]; end
if isempty(cLimMod) || cLimMod(1) == cLimMod(2), cLimMod = [0 1]; end
dMax    = quantile(abs([diffG(:);diffN(:)]),0.95);
if isnan(dMax) || dMax == 0, dMax = max(abs([diffG(:);diffN(:)])); end
if isnan(dMax) || dMax == 0, dMax = 1; end
cLimDiff= [-dMax dMax];

fprintf('Global Climits: [%.2f, %.2f]\n', cLimMod(1), cLimMod(2));
fprintf('Global Diff Climits: [%.2f, %.2f]\n', cLimDiff(1), cLimDiff(2));

%--------------------------------------------------------------------%
% 2. PLOT 2 × 3 GRID (Manual 'axes' Positioning)
%--------------------------------------------------------------------%
fig = figure('Units','centimeters','OuterPosition',[0 0 1 1]);
fig.Position(3) = 35;  fig.Position(4) = 13;

% layout parameters
numRows=2; numCols=3;
left_margin  = 0.005; right_margin = 0.07; %More space on right for CB labels
bottom_margin= 0;     top_margin  = 0;

gap_c1_c2=0; gap_c2_c3=0.06; vert_space=0;

% Calculate widths/heights based on a 3-column layout
total_width = 1 - left_margin - right_margin;
plot_width = (total_width - gap_c1_c2 - gap_c2_c3) / numCols;
plot_height = (1 - top_margin - bottom_margin - (numRows-1)*vert_space) / numRows;

% Calculate 'left' positions based on the 3 plot columns
left_col1 = left_margin;
left_col2 = left_col1 + plot_width + gap_c1_c2;
left_col3 = left_col2 + plot_width + gap_c2_c3;

% Calculate 'bottom' positions
bottom_row2 = bottom_margin;
bottom_row1 = bottom_row2 + plot_height + vert_space;

% Define positions array
positions = {
    [left_col1, bottom_row1, plot_width, plot_height]; % ax(1)
    [left_col2, bottom_row1, plot_width, plot_height]; % ax(2)
    [left_col3, bottom_row1, plot_width, plot_height]; % ax(3)
    [left_col1, bottom_row2, plot_width, plot_height]; % ax(4)
    [left_col2, bottom_row2, plot_width, plot_height]; % ax(5)
    [left_col3, bottom_row2, plot_width, plot_height]; % ax(6)
    };

maps   = {gauss,gauss_cor,diffG,nig,nig_cor,diffN};
cmaps  = {cmapMod,cmapMod,cmapDiff,cmapMod,cmapMod,cmapDiff};
clims  = {cLimMod,cLimMod,cLimDiff,cLimMod,cLimMod,cLimDiff};
titles = {model_title(1),model_title(2),'Difference Gaussian', ...
          model_title(3),model_title(4),'Difference NIG'};

ax = gobjects(1,6);
for i = 1:6
    ax(i) = axes('Position', positions{i}); % Use axes for precise position
    drawWorld(latGrid,lonGrid,maps{i},cmaps{i},clims{i},titles{i},fs);
end
drawnow;

%--------------------------------------------------------------------%
% 3. TWO COMMON COLOUR-BARS (with scaling)
%--------------------------------------------------------------------%
% Define scaling factor
scale = 0.50; % 50% of the full height

% GAP from map to bar and bar width (in figure units)
gapX     = 0.005;
cbarW    = 0.018;

% Calculate the full vertical span that colorbars should cover
posTopAny   = get(ax(2),'Position'); posBotAny   = get(ax(5),'Position');
full_span_bottom = posBotAny(2);
full_span_top    = posTopAny(2) + posTopAny(4);
barFullH = full_span_top - full_span_bottom;

% Calculate new height and bottom position for scaled and centered colorbar
newH = scale * barFullH;
newY = full_span_bottom + (barFullH - newH)/2;

% ---------- MODELS colour-bar (column-2) -----------------------
posTop2 = get(ax(2),'Position'); % Position of the top plot for this colorbar
cb1 = colorbar(ax(2),'Location','EastOutside');
cb1.Position = [posTop2(1) + posTop2(3) + gapX, ... % x-position
                newY, ...                         % y-position (scaled and centered)
                cbarW, ...                        % width
                newH];                            % height (scaled)
cb1.FontSize = fs-4;
ylabel(cb1,'\rho_{u_1,u_2}','Rotation',0,'FontSize',fs, ...
               'VerticalAlignment','middle', 'Interpreter', 'tex');
yl1 = cb1.Label; yl1.Units = 'data'; yl1.Position(1) = 3.5; % Nudge label right


% ---------- DIFFERENCE colour-bar (column-3) -------------------
posTop3 = get(ax(3),'Position'); % Position of the top plot for this colorbar
cb2 = colorbar(ax(3),'Location','EastOutside');
cb2.Position = [posTop3(1) + posTop3(3) + gapX, ... % x-position
                newY, ...                         % y-position (scaled and centered)
                cbarW, ...                        % width
                newH];                            % height (scaled)
cb2.FontSize = fs-4;
% Ensure All Plots Use Common Limits (Important!)
clim(ax(1), cLimMod); clim(ax(4), cLimMod); clim(ax(5), cLimMod);
clim(ax(6), cLimDiff);

%--------------------------------------------------------------------%
% 4. SAVE
%--------------------------------------------------------------------%
outDir = fullfile(outputDir,num2str(presLevel));
if ~exist(outDir,'dir'), mkdir(outDir); end
    fname = fullfile(outDir,sprintf('%d_correlation_combined.png',presLevel));
    print(fig,'-dpng','-r330',fname);
    fprintf('Figure saved to %s\n',fname);
end

%=====================  helper: build one grid  ==========================%
function G = makeGrid(T,mdl,latGrid,lonGrid,Grid)
    G = nan(size(latGrid));
    S = T(strcmp(T.model,mdl),:);
    for g = 1:404
        r = S(S.gridID==g,:);
        if isempty(r), continue, end
        if size(r,1)>1
            [~,j] = max(sum(~cellfun(@isempty,table2cell(r)),2)); r=r(j,:);
        end
        v = biv_corr(r); % Ensure biv_corr is defined
        latMin=Grid(g,1); latMax=Grid(g,2);
        lonMin=Grid(g,3); lonMax=Grid(g,4);
        rows = (ceil(lonMin):floor(lonMax))+1;
        c0_list   = find(latGrid(1,:)==latMin+1);
        if isempty(c0_list)
            [~, c0] = min(abs(latGrid(1, :) - (latMin + 1)));
        else
            c0 = c0_list(1);
        end
        cols = c0:(c0+latMax-latMin-1);
        rows = rows(rows > 0 & rows <= size(G, 1));
        cols = cols(cols > 0 & cols <= size(G, 2));
        if ~isempty(rows) && ~isempty(cols)
            G(rows,cols)=v;
        end
    end
end
%=====================  helper: draw map panel  ==========================%
function drawWorld(lat,lon,Z,cmap,cLimVals,titleStr,fs)
    worldmap('World'); tightmap; mlabel off; plabel off; gridm off
    setm(gca, 'MapProjection', 'robinson');
    surfm(lat,lon,Z);
    colormap(gca,cmap); clim(gca,cLimVals);
    title(titleStr,'FontSize',fs, 'Interpreter', 'tex','FontWeight', 'bold'); addCoast
end
%=====================  helper: biv_corr  ===============================%
function c = biv_corr(r)
    k1=r.Kappa1; k2=r.Kappa2; rho=r.Rho;
    if iscell(k1),k1=str2double(k1{1}); end
    if iscell(k2),k2=str2double(k2{1}); end
    if iscell(rho),rho=str2double(rho{1}); end
    if any(isnan([k1 k2 rho]))||k1==k2, c=NaN; else
        c=(2*k1*k2*rho*log(k1/k2))/((k1^2-k2^2)*sqrt(1+rho^2));
    end
end
%=====================  helper: coastlines  =============================%
function addCoast
    hold on; load coastlines
    plotm(coastlat,coastlon,'k', 'LineWidth', 0.5) 
    h = patchm(coastlat,coastlon,[1 1 1]); % Land
    set(h, 'EdgeColor', 'black');
    hold off;
end
% end