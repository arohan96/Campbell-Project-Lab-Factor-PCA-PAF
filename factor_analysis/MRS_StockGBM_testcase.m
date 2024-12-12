% Generate Random Data
num_days = 252;  % Number of trading days (1 year)
initial_price = 100;  % Starting stock price
mu = 0.0005;  % Drift (daily mean return)
sigma = 0.02;  % Volatility (daily standard deviation)
dt = 1 / 252;  % Daily time step

% Generate random returns using Geometric Brownian Motion
random_returns = mu * dt + sigma * sqrt(dt) * randn(num_days, 1);
stock_prices = initial_price * cumprod(1 + random_returns);  % Simulated stock prices

% Create a table with datetime dates
dates = (datetime(2023, 1, 1) + caldays(0:num_days-1))';  
data = table(dates, stock_prices, 'VariableNames', {'Date', 'Factor_Returns'});  

% Initialize Strategy
strategy = MeanReversionStrategy(20, 2);

% Apply Strategy
[data, metrics] = strategy.apply_strategy(data, 'Factor_Returns', true);

% Display Metrics
disp('Performance Metrics:');
disp(metrics);
