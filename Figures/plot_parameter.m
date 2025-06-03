function plot_parameter(presLevel,i_order,model,model_title,outputDir,paramName)
% PLOT_PARAMETER Plot spatial maps for a single model parameter across multiple models.
%
%   plot_parameter(presLevel, i_order, model, model_title, outputDir, paramName)
%
%   Generates a 2x3 grid of spatial maps displaying a chosen parameter for
%   four different models, along with two difference maps.
%
%   Inputs:
%   * presLevel   – Scalar pressure level at which to load data.
%   * i_order     – A 1x4 index vector specifying the display order of the four models.
%                   For example, [1 2 3 4] would use model{1}, model{2}, etc.
%   * model       – Cell array of model identifiers (e.g., {'Gaussian', 'Gaussian_Correlated', ...}).
%                   Its length should be at least max(i_order).
%   * model_title – Cell array of display titles for each model (same length as 'model').
%   * outputDir   – Base directory where the generated figures will be saved.
%   * paramName   – String, specifying the parameter to plot. Must be one of:
%                   'Rho', 'Sigma1', 'Sigma2', 'Kappa1', 'Kappa2', 'Mu1', 'Mu2',
%                   'Nu1', 'Nu2'.
%
%   Output:
%   A PNG figure is saved to outputDir/presLevel/ named after the pressure level
%   and parameter (e.g., '300_Rho_combined.png'). The figure includes:
%   * Four panels showing the chosen parameter for each of the four selected models.
%   * Two "difference" panels: (Gaussian Correlated - Gaussian) and (NIG Correlated - NIG).
%   * Two common colorbars: one for model parameters and one for differences.
%
%   Dependencies:
%   * Requires the 'colorBrewer' toolbox (specifically `brewermap`) and the
%     `redblue()` utility to be available on the MATLAB path.
%   * Relies on helper functions: `makeGridParam`, `drawWorld`, `addCoast`,
%     and `latexLabel` (which are typically appended within this file or
%     placed in a common 'helpers' directory).
%   * Data files: Reads 'main_results.csv' from `~/Documents/Results/<presLevel>/`
%     and 'grid_equal.mat' from `~/Documents/Results/Data/`.
%--------------------------------------------------------------------------%
addpath('~/Documents/Results/colorBrewer/');
lat = linspace(-90,  90,181);
lon = linspace( 20, 380,361);
[latGrid,lonGrid] = meshgrid(lat,lon);
fs  = 16;                                 % base font size
cmapMod  = brewermap([],'YlOrBr');        % model panels
cmapDiff = redblue(17); cmapDiff(9,:) = 1; % symmetric diff panels

%--------------------------------------------------------------------%
% 1. Load data & build four grids (1:Gauss,2:Gauss_cor,3:NIG,4:NIG_cor)
%--------------------------------------------------------------------%
fn = sprintf('~/Documents/Results/%d/main_results.csv',presLevel);
if ~isfile(fn); error('File not found: %s',fn); end
T  = readtable(fn);
load('~/Documents/Results/Data/grid_equal.mat','Grid');

varGrid = cell(1,4);
for k = 1:4
    varGrid{k} = makeGridParam(T,model{i_order(k)},paramName,latGrid,lonGrid,Grid);
end
% convenience handles
[gauss,gauss_cor,nig,nig_cor] = varGrid{:};
diffG = gauss_cor - gauss;  % Gaussian correction – base
diffN = nig_cor   - nig;    % NIG      correction – base

%--------------------------------------------------------------------%
% 2. Determine colour‑limits (global to ensure comparability)
%--------------------------------------------------------------------%
vals    = [gauss(:);gauss_cor(:);nig(:);nig_cor(:)];
vals    = vals(~isnan(vals));
if isempty(vals), error('All values are NaN for %s at %d hPa',paramName,presLevel); end
cLimMod = quantile(vals,[.05 .95]);
if diff(cLimMod)==0, cLimMod = [min(vals) max(vals)]; end
if diff(cLimMod)==0, cLimMod = [0 1]; end  % final fallback

dMax    = quantile(abs([diffG(:);diffN(:)]),0.9,'all');
if isnan(dMax) || dMax==0, dMax = max(abs([diffG(:);diffN(:)])); end
if dMax==0, dMax = 1; end
cLimDiff = [-dMax dMax];

%--------------------------------------------------------------------%
% 3. Build figure (manual axes positioning for full control)
%--------------------------------------------------------------------%
fig = figure('Units','centimeters','OuterPosition',[0 0 1 1]);
fig.Position(3) = 35;  fig.Position(4) = 13;

% layout parameters
numRows=2; numCols=3;
left_margin  = 0.005; right_margin = 0.07;
bottom_margin= 0;     top_margin  = 0;

gap_c1_c2=0; gap_c2_c3=0.06; vert_space=0;

total_w = 1-left_margin-right_margin;
plot_w  = (total_w - gap_c1_c2 - gap_c2_c3)/numCols;
plot_h  = (1-top_margin-bottom_margin-(numRows-1)*vert_space)/numRows;

left1 = left_margin;
left2 = left1 + plot_w + gap_c1_c2;
left3 = left2 + plot_w + gap_c2_c3;

bot2 = bottom_margin;
bot1 = bot2 + plot_h + vert_space;

positions = { [left1 bot1 plot_w plot_h],  ...
              [left2 bot1 plot_w plot_h],  ...
              [left3 bot1 plot_w plot_h],  ...
              [left1 bot2 plot_w plot_h],  ...
              [left2 bot2 plot_w plot_h],  ...
              [left3 bot2 plot_w plot_h]   };

maps   = {gauss,gauss_cor,diffG,nig,nig_cor,diffN};
cmaps  = {cmapMod,cmapMod,cmapDiff,cmapMod,cmapMod,cmapDiff};
clims  = {cLimMod,cLimMod,cLimDiff,cLimMod,cLimMod,cLimDiff};

titles = {model_title(1),model_title(2),'Difference Gaussian', ...
          model_title(3),model_title(4),'Difference NIG'};

ax = gobjects(1,6);
for i = 1:6
    ax(i) = axes('Position',positions{i}); 
    drawWorld(latGrid,lonGrid,maps{i},cmaps{i},clims{i},titles{i},fs);
end
drawnow;

%--------------------------------------------------------------------%
% 4. Two common colour‑bars (models & differences)
%--------------------------------------------------------------------%
scale  = 0.50;           % bar height = 50% of map stack
cbGap  = 0.005;          % horizontal gap between map and bar
cbW    = 0.018;          % colour‑bar width

% determine full vertical span from maps in col‑2
posTop = get(ax(2),'Position'); posBot = get(ax(5),'Position');
spanH  = posTop(2)+posTop(4) - posBot(2);
newH   = scale*spanH; newY = posBot(2) + (spanH-newH)/2;

% ----- models (attached to column‑2) -----
cb1 = colorbar(ax(2),'Location','EastOutside');
cb1.Position = [posTop(1)+posTop(3)+cbGap , newY , cbW , newH];
cb1.FontSize = fs-4; cb1.Label.Interpreter = 'tex';
cb1.Label.Rotation=0; cb1.Label.HorizontalAlignment='left';
cb1.Label.VerticalAlignment='middle'; cb1.Label.String = latexLabel(paramName);
cb1.Label.Position(1) = 2.5; % nudge right
cb1.Label.FontSize = fs;

% ----- differences (attached to column‑3) -----
posTop3 = get(ax(3),'Position');
cb2 = colorbar(ax(3),'Location','EastOutside');
cb2.Position = [posTop3(1)+posTop3(3)+cbGap , newY , cbW , newH];
cb2.FontSize = fs-4;

% enforce common lims
set(ax([1 2 4 5]),'CLim',cLimMod);
set(ax([3 6]),'CLim',cLimDiff);

%--------------------------------------------------------------------%
% 5. Save
%--------------------------------------------------------------------%
outDir = fullfile(outputDir,num2str(presLevel)); if ~exist(outDir,'dir'); mkdir(outDir); end
fname  = fullfile(outDir,sprintf('%d_%s_combined.png',presLevel,paramName));
print(fig,'-dpng','-r330',fname);
fprintf('Figure saved: %s\n',fname);
end

%============================= HELPERS ==============================%
function G = makeGridParam(T,mdl,paramName,latGrid,lonGrid,Grid)
    G = nan(size(latGrid));
    S = T(strcmp(T.model,mdl),:);
    for g = 1:404
        r = S(S.gridID==g,:);
        if isempty(r), continue, end
        if size(r,1)>1
            [~,j] = max(sum(~cellfun(@isempty,table2cell(r)),2)); r=r(j,:);
        end
        v = r.(paramName);
        if iscell(v), v = str2double(v{1}); end
        if isnan(v), continue, end
        latMin=Grid(g,1); latMax=Grid(g,2);
        lonMin=Grid(g,3); lonMax=Grid(g,4);
        rows = (ceil(lonMin):floor(lonMax))+1;
        c0_list = find(latGrid(1,:)==latMin+1);
        if isempty(c0_list)
            [~, c0] = min(abs(latGrid(1,:) - (latMin+1)));
        else
            c0 = c0_list(1);
        end
        cols = c0:(c0+latMax-latMin-1);
        rows = rows(rows>0 & rows<=size(G,1));
        cols = cols(cols>0 & cols<=size(G,2));
        if ~isempty(rows) && ~isempty(cols)
            G(rows,cols) = v;
        end
    end
end

function drawWorld(lat,lon,Z,cmap,clims,tStr,fs)
    worldmap('World'); tightmap; mlabel off; plabel off; gridm off
    setm(gca,'MapProjection','robinson');
    surfm(lat,lon,Z); colormap(gca,cmap); clim(gca,clims);
    title(tStr,'FontSize',fs,'Interpreter','tex','FontWeight','bold');
    addCoast;
end

function addCoast
    hold on; load coastlines
    plotm(coastlat,coastlon,'k','LineWidth',0.5);
    h = patchm(coastlat,coastlon,[1 1 1]); set(h,'EdgeColor','black'); hold off
end

function lbl = latexLabel(p)
    switch lower(p)
        case 'rho',     lbl = '\rho';
        case 'sigma1',  lbl = '\sigma_1';
        case 'sigma2',  lbl = '\sigma_2';
        case 'kappa1',  lbl = '\kappa_1';
        case 'kappa2',  lbl = '\kappa_2';
        case 'mu1',     lbl = '\mu_1';
        case 'mu2',     lbl = '\mu_2';
        case 'nu1',     lbl = '\eta_1'; 
        case 'nu2',     lbl = '\eta_2';
        otherwise,      lbl = p;
    end
end


