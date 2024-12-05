% Test the MeanReversionBM class
numAssets = 5; % Number of assets
numDays = 300; % Number of days
window = 20; % Moving average window size
k = 2; % Number of standard deviations for signal thresholds
dt = 1 / 252; % Daily time step
mu = 0.1; % Annual drift
sigma = 0.2; % Annual volatility
initialPrices = 100 * ones(numAssets, 1); % Initial prices for all assets

% Create the class object
strategy = MeanReversionBM(numAssets, numDays, window, k, dt, mu, sigma, initialPrices);

% Run the strategy
strategy = strategy.simulatePrices();
strategy = strategy.calculateRollingStats();
strategy = strategy.generateSignals();
strategy = strategy.calculateStrategyReturns();

% Plot results
strategy.plotPrices();
strategy.plotCumulativeReturns();
