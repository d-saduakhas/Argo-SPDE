function plot_pair_nig(presLevel,i_order,model,model_title,outputDir,paramPrefix,paramSymbol)
% PLOT_PAIR_NIG Plot spatial maps for paired parameters (e.g., Mu1/Mu2 or Nu1/Nu2).
%
%   plot_pair_nig(presLevel, i_order, model, model_title, outputDir, paramPrefix, paramSymbol)
%
%   Generates a 2x3 grid of spatial maps for a pair of related parameters
%   (e.g., parameter 1 and parameter 2, like Mu1/Mu2 or Nu1/Nu2) for two
%   specified models (typically independent and correlated versions).
%
%   Inputs:
%   * presLevel   – Scalar pressure level for data loading.
%   * i_order     – A 1x2 index vector specifying the indices of the two models
%                   (e.g., [iid_idx, corr_idx]) from the 'model' cell array to plot.
%   * model       – Cell array of model identifiers.
%   * model_title – Cell array of display titles for each model (same order as 'model').
%   * outputDir   – Base output directory where figures will be written.
%   * paramPrefix – String, the column prefix in the data table (e.g., 'Mu' for Mu1/Mu2,
%                   or 'Nu' for Nu1/Nu2).
%   * paramSymbol – LaTeX symbol for the parameter, used in colorbar and row labels
%                   (e.g., '\mu' for Mu, '\eta' for Nu).
%
%   Output:
%   A PNG figure is saved to outputDir/presLevel/ named after the pressure level
%   and parameter prefix (e.g., '300_mu_combined.png'). The figure includes:
%   * Row 1: Panels for parameter 1 (e.g., Mu1) for both specified models,
%              plus a difference map (correlated - independent).
%   * Row 2: Panels for parameter 2 (e.g., Mu2) for both specified models,
%              plus a difference map (correlated - independent).
%   * Two common colorbars: one for model parameters and one for differences,
%     both scaled to 50% height.
%   * Row labels using the specified LaTeX symbol (e.g., $\mu_1$ and $\mu_2$).
%
%   Dependencies:
%   * Requires the 'colorBrewer' toolbox (specifically `brewermap`) and the
%     `redblue()` utility to be available on the MATLAB path.
%   * Relies on helper functions: `makeGridParam`, `drawWorld`, and `addCoast`.
%     These helpers are typically placed in a common 'helpers' directory.
%   * Data files: Reads 'main_results.csv' from `~/Documents/Results/<presLevel>/`
%     and 'grid_equal.mat' from `~/Documents/Results/Data/`.
%--------------------------------------------------------------------%
addpath('~/Documents/Results/colorBrewer/');
lat = linspace(-90,90,181); lon = linspace(20,380,361);
[latGrid,lonGrid] = meshgrid(lat,lon);
fs  = 16;
cmapMod  = brewermap([],'YlOrBr');
cmapDiff = redblue(17); cmapDiff(9,:) = 1;

%--------------------------------------------------------------------%
% 1. LOAD DATA -------------------------------------------------------%
%--------------------------------------------------------------------%
fn = sprintf('~/Documents/Results/%d/main_results.csv',presLevel);
if ~isfile(fn), error('File not found: %s',fn); end
T = readtable(fn);
load('~/Documents/Results/Data/grid_equal.mat','Grid');

params = {sprintf('%s1',paramPrefix), sprintf('%s2',paramPrefix)};
G = cell(2,2);
for r = 1:2
    for c = 1:2
        G{r,c} = makeGridParam(T,model{i_order(c)},params{r},latGrid,lonGrid,Grid);
    end
end
Diff{1} = G{1,2} - G{1,1}; Diff{2} = G{2,2} - G{2,1};

%--------------------------------------------------------------------%
% 2. COLOUR LIMITS ---------------------------------------------------%
%--------------------------------------------------------------------%
allMod = [G{1,1}(:);G{1,2}(:);G{2,1}(:);G{2,2}(:)]; allMod = allMod(~isnan(allMod));
cLimMod = quantile(allMod,[.025 .975]);
if diff(cLimMod)==0, cLimMod=[min(allMod) max(allMod)]; end
if diff(cLimMod)==0, cLimMod=[0 1]; end
allDiff=[Diff{1}(:);Diff{2}(:)]; allDiff = allDiff(~isnan(allDiff));
dMax = quantile(abs(allDiff),0.95);
if isnan(dMax)||dMax==0, dMax = max(abs(allDiff)); end
if dMax==0||isnan(dMax), dMax=1; end
cLimDiff = [-dMax dMax];

%--------------------------------------------------------------------%
% 3. FIGURE LAYOUT ---------------------------------------------------%
%--------------------------------------------------------------------%

fig = figure('Units','centimeters','OuterPosition',[0 0 1 1]);
fig.Position(3)=35; fig.Position(4)=13;
left=0.04; right=0.06; gap12=0; gap23=0.06; top=0; bot=0; vgap=0;
plot_w=(1-left-right-gap12-gap23)/3; plot_h=(1-top-bot-vgap)/2;
L1=left; L2=L1+plot_w+gap12; L3=L2+plot_w+gap23; B2=bot; B1=B2+plot_h+vgap;
pos={ [L1 B1 plot_w plot_h],[L2 B1 plot_w plot_h],[L3 B1 plot_w plot_h], ...
      [L1 B2 plot_w plot_h],[L2 B2 plot_w plot_h],[L3 B2 plot_w plot_h] };

maps={G{1,1},G{1,2},Diff{1}, G{2,1},G{2,2},Diff{2}};
cmaps={cmapMod,cmapMod,cmapDiff,cmapMod,cmapMod,cmapDiff};
clms={cLimMod,cLimMod,cLimDiff,cLimMod,cLimMod,cLimDiff};

titles={model_title{i_order(1)},model_title{i_order(2)},'Difference', ...
        model_title{i_order(1)},model_title{i_order(2)},'Difference'};

ax=gobjects(1,6);
for k=1:6
    ax(k)=axes('Position',pos{k});
    drawWorld(latGrid,lonGrid,maps{k},cmaps{k},clms{k},titles{k},fs);
end

%--------------------------------------------------------------------%
% 4. COLOUR‑BARS (50 % HEIGHT)                                       %
%--------------------------------------------------------------------%
cbGap=0.005; cbW=0.018; scale=0.50;
posTop=get(ax(2),'Position'); posBot=get(ax(5),'Position');
fullH=(posTop(2)+posTop(4))-posBot(2); newH=scale*fullH; newY=posBot(2)+(fullH-newH)/2;

cbMod=colorbar(ax(2),'Location','EastOutside');
cbMod.Position=[posTop(1)+posTop(3)+cbGap , newY , cbW , newH];
cbMod.FontSize=fs-4;
% set(cbMod.Label,'String',paramSymbol,'Interpreter','tex','Rotation',0, ...
%     'FontSize',fs,'HorizontalAlignment','left','VerticalAlignment','middle');
% cbMod.Label.Position(1)=2.5;

posDiff=get(ax(3),'Position');
cbDiff=colorbar(ax(3),'Location','EastOutside');
cbDiff.Position=[posDiff(1)+posDiff(3)+cbGap , newY , cbW , newH];
cbDiff.FontSize=fs-4; cbDiff.Label.String='';

set(ax([1 2 4 5]),'CLim',cLimMod); set(ax([3 6]),'CLim',cLimDiff);

%--------------------------------------------------------------------%
% 5. ROW LABELS ------------------------------------------------------%
%--------------------------------------------------------------------%
rowLbls={sprintf('$%s_{1}$',paramSymbol), sprintf('$%s_{2}$',paramSymbol)};
for r=1:2
    a=ax((r-1)*3+1); p=get(a,'Position'); lblW=0.035; lblH=0.05;
    xLbl=max(p(1)-lblW-0.005,0.002); yLbl=p(2)+(p(4)-lblH)/2;
    annotation(fig,'textbox',[xLbl yLbl lblW lblH],'String',rowLbls{r}, ...
        'Interpreter','latex','FontSize',fs+2,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none');
end

drawnow;

%--------------------------------------------------------------------%
% 6. SAVE ------------------------------------------------------------%
%--------------------------------------------------------------------%
outDir=fullfile(outputDir,num2str(presLevel)); if ~exist(outDir,'dir'), mkdir(outDir); end
fname=fullfile(outDir,sprintf('%d_%s_combined.png',presLevel,lower(paramPrefix)));
print(fig,'-dpng','-r330',fname);
fprintf('Figure saved: %s\n',fname);
end

%===================== helper functions =============================%
function G=makeGridParam(T,mdl,paramName,latGrid,lonGrid,Grid)
G=nan(size(latGrid)); S=T(strcmp(T.model,mdl),:);
for g=1:404
    r=S(S.gridID==g,:); if isempty(r), continue, end
    if size(r,1)>1
        [~,j]=max(sum(~cellfun(@isempty,table2cell(r)),2)); r=r(j,:);
    end
    v=r.(paramName); if iscell(v), v=str2double(v{1}); end
    if isnan(v), continue, end
    latMin=Grid(g,1); latMax=Grid(g,2); lonMin=Grid(g,3); lonMax=Grid(g,4);
    rows=(ceil(lonMin):floor(lonMax))+1;
    c0=find(latGrid(1,:)==latMin+1,1);
    if isempty(c0), [~,c0]=min(abs(latGrid(1,:)-(latMin+1))); end
    cols=c0:(c0+latMax-latMin-1);
    rows=rows(rows>0 & rows<=size(G,1)); cols=cols(cols>0 & cols<=size(G,2));
    if ~isempty(rows)&&~isempty(cols), G(rows,cols)=v; end
end
end

function drawWorld(lat,lon,Z,cmap,climVal,ttl,fs)
worldmap('World'); tightmap; mlabel off; plabel off; gridm off
setm(gca,'MapProjection','robinson'); surfm(lat,lon,Z);
colormap(gca,cmap); clim(gca,climVal);
title(ttl,'FontSize',fs,'Interpreter','tex','FontWeight','bold'); addCoast;
end

function addCoast
hold on; load coastlines; plotm(coastlat,coastlon,'k','LineWidth',0.5);
h=patchm(coastlat,coastlon,[1 1 1]); set(h,'EdgeColor','black'); hold off;
end
