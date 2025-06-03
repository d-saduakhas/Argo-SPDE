%% ========================================================================
%  Script : Make LaTeX table of “best‑model” grid counts
%           (Temperature, Salinity shown as “TempCount, SalCount”)
%  Author : Damilya Saduakhas
%  Updated: 13‑May‑2025
% ========================================================================
close all; clear;

% --------------------------- USER SETTINGS -------------------------------
pressureLevels  = [10 300 1000];                  % hPa
metrics         = {'RMSE','MAE','CRPS','SCRPS'};
vars            = {'Temp','Psal'};                % 1 = Temperature, 2 = Salinity

models          = {'gauss_indep','gauss_cor','nig_indep','nig_cor'};
modelLabels     = {'Gaussian','Gaussian','NIG','NIG'};
structureLabels = {'diagonal','general','diagonal','general'};

rootPath        = '~/Documents/Results';          % where CSV folders live
outFile         = '~/Documents/Results/Figures/best_models_summary_temp_psal.tex';
% -------------------------------------------------------------------------

% --------------------------- PRE‑ALLOCATE -------------------------------
for p = 1:numel(pressureLevels)
    results(p).Pressure = pressureLevels(p);
    % dimensions: model × metric × variable (Temp/Sal)
    results(p).Counts   = zeros(numel(models), numel(metrics), numel(vars));
end

% --------------------------- COUNT WINNERS -------------------------------
for p = 1:numel(pressureLevels)
    pres   = pressureLevels(p);
    fName  = fullfile(rootPath, num2str(pres), 'cv_results.csv');
    if ~isfile(fName)
        warning('File %s not found – skipping pressure level %d.', fName, pres);
        continue
    end
    T = readtable(fName);
    T.cv_Temp_RMSE  = sqrt(T.cv_Temp_MSE);
    T.cv_Psal_RMSE  = sqrt(T.cv_Psal_MSE);

    for g = unique(T.gridID)'                              % loop grids
        gData = T(T.gridID == g, :);
        for v = 1:numel(vars)                              % Temp / Sal
            vName = vars{v};
            for m = 1:numel(metrics)                       % MSE … SCRPS
                col = sprintf('cv_%s_%s', vName, metrics{m});
                if ~ismember(col, gData.Properties.VariableNames), continue; end
                [~, idx]  = min(gData.(col));              % best model index
                mdl       = gData.model{idx};
                mdlIdx    = find(strcmp(models, mdl));
                if ~isempty(mdlIdx)
                    results(p).Counts(mdlIdx, m, v) = results(p).Counts(mdlIdx, m, v) + 1;
                end
            end
        end
    end
end

% --------------------------- BUILD LaTeX ---------------------------------
NL  = newline;
txt = [ ...
    "\begin{table}"
    "\centering"
    "\caption{Counts of grids with the best model performance---Temperature, Salinity}"
    "\label{tab:CV_Argo_results_temperature_best}"
    "\begin{tabular}{@{}cllcccc@{}}"
    "\toprule"
    "\textbf{Pressure level} & \textbf{Model} & \textbf{Structure} & " + ...
    "\textbf{RMSE} & \textbf{MAE} & \textbf{CRPS} & \textbf{SCRPS}\\"
    "\midrule"
    ];

for p = 1:numel(results)
    C = results(p).Counts;                                % 4×4×2
    for r = 1:numel(models)
        % -------- first column (pressure level) ----------
        if r == 1
            presStr = string(results(p).Pressure);
        else
            presStr = "";
        end
        % -------- metric cells: "Temp, Sal" -------------
        cells = arrayfun(@(t,s) sprintf('%d, %d', t, s), ...
                         squeeze(C(r,:,1)), squeeze(C(r,:,2)), ...
                         'UniformOutput', false);

        row = sprintf('%s & %s & %s & %s & %s & %s & %s \\\\', ...
                      presStr, modelLabels{r}, structureLabels{r}, cells{:});
        txt(end+1) = row; %#ok<SAGROW>
    end
    txt(end+1) = "\midrule";   %#ok<SAGROW>
end

txt = [txt; "\bottomrule"; "\end{tabular}"; "\end{table}"];
latex = strjoin(txt, NL);

% --------------------------- OUTPUT --------------------------------------
disp('------------ LaTeX table string ------------');
disp(latex);

fid = fopen(outFile, 'w');
if fid == -1
    warning('Could not open %s for writing.', outFile);
else
    fprintf(fid, '%s', latex);
    fclose(fid);
    fprintf('LaTeX table saved to %s\n', outFile);
end
