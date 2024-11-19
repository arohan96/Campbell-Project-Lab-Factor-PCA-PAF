%% params
% The number of days and markets will be determined after reading the data
nFactorsToCompute = 4;
modelType = 'PAF';
factorConstructionLookback = 252;
volLookback = 252;
tolerance = 1e-8;
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
visualize = true;  % Whether to visualize eigenvalues and communalities
% Add new parameters for visualizing top tickers
visualizeTopTickers = true;  % Set to true to visualize top tickers
numTopTickers = 50;  % Number of top tickers to display for each factor

%% setup
clc;

% Read returns data from CSV file
returnsData = readtable('russell_1000_retuns.csv', 'ReadVariableNames', false);

% Extract dates from the first column
dates = returnsData{:, 1};

% Extract returns data
mktRtns = returnsData{:, 2:end};

% Read tickers from a separate file
tickersData = readtable('tickers_russell_1000.csv', 'ReadVariableNames', false);

% Extract tickers from the first column
tickers = tickersData{:,1};

% Ensure tickers are a column vector
tickers = tickers(:);

% Check if tickers correspond to rows in mktRtns
if length(tickers) == size(mktRtns, 1)
    % Transpose mktRtns so that tickers correspond to columns
    mktRtns = mktRtns';
elseif length(tickers) ~= size(mktRtns, 2)
    error('Number of tickers (%d) does not match number of markets in returns data (%d)', length(tickers), size(mktRtns, 2));
end

[nDays, nMkts] = size(mktRtns);
disp(['Number of days: ', num2str(nDays)]);
disp(['Number of markets: ', num2str(nMkts)]);

% Remove any rows with missing data (NaNs)
nanRows = any(isnan(mktRtns), 2);
mktRtns = mktRtns(~nanRows, :);
dates = dates(~nanRows);

% Optionally, select top N variables based on variability
%variableStd = std(mktRtns);
%[~, idx] = sort(variableStd, 'descend');

% Select top N variables (e.g., 600)
%N = min(600, nMkts);
%selectedIdx = idx(1:N);
%mktRtns = mktRtns(:, selectedIdx);
%tickers = tickers(selectedIdx);  % Select corresponding tickers
%nMkts = size(mktRtns, 2);
%disp(['Number of markets after selecting top ', num2str(N), ' variables: ', num2str(nMkts)]);

% Set positions (e.g., equal weights)
myPositions = ones(1, nMkts) / nMkts;

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
params.tickers = tickers;  % Pass tickers to params
% Add new visualization parameters
params.visualizeTopTickers = visualizeTopTickers;
params.numTopTickers = numTopTickers;

% Perform factor decomposition
[estFactorRtns, portBetas, factorVols] = RDfactorDecomposition(mktRtns, myPositions, params);

% Visualize the estimated factor returns
figure();
for iii = 1:nFactorsToCompute
    subplot(2, 2, iii);
    plot(dates, estFactorRtns(:, iii));
    title(['Factor ', num2str(iii)]);
    datetick('x', 'yyyy');
end
