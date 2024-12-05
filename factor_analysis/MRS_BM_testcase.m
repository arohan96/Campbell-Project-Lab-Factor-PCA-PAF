% Parameters
numAssets = 5; % Number of assets
numDays = 300; % Number of days
window = 20; % Moving window size
k = 2; % Number of standard deviations for signal thresholds
dt = 1 / 252; % Daily time step (assuming 252 trading days per year)
rng(42); % Seed for reproducibility

% Simulate asset prices following Geometric Brownian Motion
initialPrices = 100 * ones(numAssets, 1); % Initial prices for all assets
mu = 0.1; % Annual drift
sigma = 0.2; % Annual volatility

% Generate Brownian motion
dW = sqrt(dt) * randn(numDays, numAssets); % Random walk increments
logReturns = (mu - 0.5 * sigma^2) * dt + sigma * dW; % GBM log returns
prices = zeros(numDays, numAssets);
prices(1, :) = initialPrices; % Set initial prices
for t = 2:numDays
    prices(t, :) = prices(t-1, :) .* exp(logReturns(t, :)); % GBM price evolution
end

% Compute daily returns
returns = diff(prices) ./ prices(1:end-1, :);

% Initialize matrices for rolling statistics and strategy signals
rollingDrift = zeros(size(returns));
rollingVolatility = zeros(size(returns));
upperThreshold = zeros(size(returns));
lowerThreshold = zeros(size(returns));
signals = zeros(numDays, numAssets);
strategyReturns = zeros(numDays-1, numAssets);

% Calculate rolling drift, volatility, and thresholds for each asset
for asset = 1:numAssets
    rollingDrift(:, asset) = movmean(log(1 + returns(:, asset)), window); % Drift
    rollingVolatility(:, asset) = movstd(log(1 + returns(:, asset)), window); % Volatility
    upperThreshold(:, asset) = rollingDrift(:, asset) + k * rollingVolatility(:, asset);
    lowerThreshold(:, asset) = rollingDrift(:, asset) - k * rollingVolatility(:, asset);
end

% Generate signals and calculate strategy returns
for asset = 1:numAssets
    for t = window:numDays-1
        if returns(t, asset) > exp(upperThreshold(t, asset)) - 1
            signals(t, asset) = -1; % Short signal
        elseif returns(t, asset) < exp(lowerThreshold(t, asset)) - 1
            signals(t, asset) = 1; % Long signal
        end
    end
    strategyReturns(:, asset) = signals(2:end, asset) .* returns(:, asset); % Fix dimension mismatch
end

% Calculate cumulative returns for each asset
cumulativeReturns = cumprod(1 + strategyReturns);

% Plot cumulative returns for all assets
figure;
plot(cumulativeReturns);
title('Cumulative Returns of Mean Reversion Strategy for Multiple Assets');
xlabel('Days');
ylabel('Cumulative Returns');
legend(arrayfun(@(x) sprintf('Asset %d', x), 1:numAssets, 'UniformOutput', false));
grid on;

% Plot price evolution of assets
figure;
plot(prices);
title('Simulated Prices of Assets (Brownian Motion)');
xlabel('Days');
ylabel('Price');
legend(arrayfun(@(x) sprintf('Asset %d', x), 1:numAssets, 'UniformOutput', false));
grid on;
