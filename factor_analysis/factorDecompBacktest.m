%% params
% The number of days and markets will be determined after reading the data
nFactorsToCompute = 4;
modelType = 'PAF';
factorConstructionLookback = 60;
volLookback = 60;
tolerance = 1e-3;
iterations = 100;
kaiserNormalizeLoadings = true;  % true or false (use kaiser normalization
% for loadings?)
rotationType = '';  % '', varimax, quartimax, promax, equamax, 
% orthomax ('' = no rotation)
orthoGamma = 0.35;  % (0 < orthoGamma < 1) only used when 
% rotationType='orthomax'. coefficient that controls the correlation target
% between factors. 1=varimax, 0=quartimax (1=focus on orthogonality, 
% 0=reduce the number of significant factors rather than maintaining strict
% orthogonality)
builtInNormalizeLoadings = false;  % true or false (use built-in 
% normalization in rotation function for loadings?)
visualizeBeforeAfterRotation = '';  % '', before, after, both 
% ('' = none)
numVariablesToShow = 15;   % how many variables to show in visualizations
visualize = false;  % Whether to visualize eigenvalues and communalities
saveLoadings = true; % Whether to save factor loadings after each run
returnsFileName = 'asia_returns.xlsx'; % Flatfile for retreiving returns 
                                      % data

%% setup
clc;

% Read returns data from Excel file
currentPath = pwd; % Get the current directory
returnsPath = fullfile(currentPath, 'data', returnsFileName);
returnsData = readtable(returnsPath, 'ReadVariableNames', true);

% Extract dates from the first column
dates = returnsData{:, 1};

% Extract returns data
mktRtns = returnsData(:, 2:end);
mktRtns = table2array(mktRtns);

% Extract tickers
tickers = string(returnsData.Properties.VariableNames);
tickers = tickers(2:end);

% Ensure tickers are a column vector
tickers = tickers(:);

% Check if tickers correspond to rows in mktRtns
if length(tickers) == size(mktRtns, 1)
    % Transpose mktRtns so that tickers correspond to columns
    mktRtns = mktRtns';
elseif length(tickers) ~= size(mktRtns, 2)
    error( ...
        'Number of tickers (%d) does not match number of markets in returns data (%d)', ...
        length(tickers), size(mktRtns, 2));
end

[nDays, nMkts] = size(mktRtns);
disp(['Number of days: ', num2str(nDays)]);
disp(['Number of markets: ', num2str(nMkts)]);

% Remove any rows with missing data (NaNs)
nanRows = any(isnan(mktRtns), 2);
mktRtns = mktRtns(~nanRows, :);
dates = dates(~nanRows);

% Set positions (e.g., equal weights)
rollingDays = nDays - factorConstructionLookback;
h_deMean = @(x, dim) x - nanmean(x, dim);
myPositions = h_deMean(randn(rollingDays, nMkts), 2);

% Update params
params.nFactorsToCompute = nFactorsToCompute;
params.modelType = modelType;
params.nDays = nDays;
params.nMkts = nMkts;
params.factorConstructionLookback = factorConstructionLookback;
params.volLookback = volLookback;
params.tolerance = tolerance;
params.iterations = iterations;
params.kaiserNormalizeLoadings = kaiserNormalizeLoadings;
params.rotationType = rotationType;
params.builtInNormalizeLoadings = builtInNormalizeLoadings;
params.visualizeBeforeAfterRotation = visualizeBeforeAfterRotation;
params.orthoGamma = orthoGamma;
params.numVariablesToShow = numVariablesToShow;
params.visualize = visualize;

% Perform factor decomposition
ut = utils;
[estFactorRtns, portBetas, factorVols, factorLoadings] = ut.rolling( ...
    mktRtns, myPositions, params);

% normalized factor returns
estNormFactorRtns = bsxfun(@rdivide, estFactorRtns(:, :, end), ...
    factorVols(end, :));

% Plotting factor returns for the last iteration
figure();
for iii = 1:4
    subplot(2, 2, iii);
    plot(estNormFactorRtns(:, iii), 'b', 'LineWidth', 1.5);
    legend('Estimated Factor Returns')
end
