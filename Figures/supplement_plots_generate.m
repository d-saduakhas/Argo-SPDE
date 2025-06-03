%% GENERATE_SUPPLEMENTAL_PLOTS
% This script generates all supplementary material plots 
%
% - `plot_correlation_common_cb` (process correlation)
% - `plot_rho_epsilon` (nugget correlation)
% - `plot_parameter_final` (Rho, Sigma, Kappa)
% - `plot_pair_nig_final` (Mu and Nu pairs)
% Figures are saved in `Figures/<pressureLevel>/`.
%
% **Dependencies**: Plotting functions, `colorBrewer`, `redblue()`, data files.
%
clear
close all
addpath('~/colorBrewer/');

pressureLevels = [10, 300, 1000];
models = ["gauss_indep", "gauss_cor","nig_indep", "nig_cor"];
modelTitles = {'Gaussian (IID Nugget)', ...
    'Gaussian (Corr. Nugget)', ...
    'NIG (IID Nugget)',      ...
    'NIG (Corr. Nugget)' ...
    };
iOrder = [1, 2, 3, 4];

% Directory for output
outputDir = 'Figures/';

% Loop through each pressure level
for i = 1:length(pressureLevels)

    plot_correlation_common_cb(pressureLevels(1), iOrder, models, modelTitles, outputDir);

    plot_rho_epsilon(pressureLevels(i), [2, 4], models, modelTitles, outputDir);

    plot_parameter(pressureLevels(i),[1, 2, 3, 4],models,modelTitles,outputDir,'Rho');
    plot_parameter(pressureLevels(i),[1, 2, 3, 4],models,modelTitles,outputDir,'Sigma1');
    plot_parameter(pressureLevels(i),[1, 2, 3, 4],models,modelTitles,outputDir,'Sigma2');
    plot_parameter(pressureLevels(i),[1, 2, 3, 4],models,modelTitles,outputDir,'Kappa1');
    plot_parameter(pressureLevels(i),[1, 2, 3, 4],models,modelTitles,outputDir,'Kappa2');

    % mu1 / mu2
    plot_pair_nig(pressureLevels(i),[3 4],models,modelTitles,outputDir,'Mu','\mu')

    % eta1 / eta2
    plot_pair_nig(pressureLevels(i),[3 4],models,modelTitles,outputDir,'Nu','\eta')

end

