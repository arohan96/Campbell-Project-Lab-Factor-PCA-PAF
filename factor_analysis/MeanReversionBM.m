classdef MeanReversionBM
    properties
        numAssets      % Number of assets
        numDays        % Number of days
        window         % Rolling window size
        k              % Standard deviation multiplier
        dt             % Time step (e.g., daily)
        mu             % Drift (annualized)
        sigma          % Volatility (annualized)
        initialPrices  % Initial prices for assets
        prices         % Simulated prices
        returns        % Daily returns
        rollingDrift   % Rolling drift
        rollingVolatility % Rolling volatility
        upperThreshold % Upper threshold
        lowerThreshold % Lower threshold
        signals        % Trading signals
        strategyReturns % Strategy returns
        cumulativeReturns % Cumulative returns
    end

    methods
        % Constructor
        function obj = MeanReversionBM(numAssets, numDays, window, k, dt, mu, sigma, initialPrices)
            obj.numAssets = numAssets;
            obj.numDays = numDays;
            obj.window = window;
            obj.k = k;
            obj.dt = dt;
            obj.mu = mu;
            obj.sigma = sigma;
            obj.initialPrices = initialPrices;
        end

        % Simulate Geometric Brownian Motion
        function obj = simulatePrices(obj)
            rng(42); % Seed for reproducibility
            dW = sqrt(obj.dt) * randn(obj.numDays, obj.numAssets); % Brownian motion
            logReturns = (obj.mu - 0.5 * obj.sigma^2) * obj.dt + obj.sigma * dW; % Log returns
            obj.prices = zeros(obj.numDays, obj.numAssets);
            obj.prices(1, :) = obj.initialPrices; % Set initial prices
            for t = 2:obj.numDays
                obj.prices(t, :) = obj.prices(t-1, :) .* exp(logReturns(t, :)); % GBM evolution
            end
            obj.returns = diff(obj.prices) ./ obj.prices(1:end-1, :); % Daily returns
        end

        % Calculate Rolling Drift and Volatility
        function obj = calculateRollingStats(obj)
            obj.rollingDrift = zeros(size(obj.returns));
            obj.rollingVolatility = zeros(size(obj.returns));
            obj.upperThreshold = zeros(size(obj.returns));
            obj.lowerThreshold = zeros(size(obj.returns));
            for asset = 1:obj.numAssets
                obj.rollingDrift(:, asset) = movmean(log(1 + obj.returns(:, asset)), obj.window); % Drift
                obj.rollingVolatility(:, asset) = movstd(log(1 + obj.returns(:, asset)), obj.window); % Volatility
                obj.upperThreshold(:, asset) = obj.rollingDrift(:, asset) + obj.k * obj.rollingVolatility(:, asset);
                obj.lowerThreshold(:, asset) = obj.rollingDrift(:, asset) - obj.k * obj.rollingVolatility(:, asset);
            end
        end

        % Generate Signals
        function obj = generateSignals(obj)
            obj.signals = zeros(obj.numDays, obj.numAssets); % Initialize signals
            for asset = 1:obj.numAssets
                for t = obj.window:obj.numDays-1
                    if obj.returns(t, asset) > exp(obj.upperThreshold(t, asset)) - 1
                        obj.signals(t, asset) = -1; % Short signal
                    elseif obj.returns(t, asset) < exp(obj.lowerThreshold(t, asset)) - 1
                        obj.signals(t, asset) = 1; % Long signal
                    end
                end
            end
        end

        % Calculate Strategy Returns
        function obj = calculateStrategyReturns(obj)
            obj.strategyReturns = zeros(size(obj.returns));
            for asset = 1:obj.numAssets
                obj.strategyReturns(:, asset) = obj.signals(2:end, asset) .* obj.returns(:, asset); % Use shifted signals
            end
            obj.cumulativeReturns = cumprod(1 + obj.strategyReturns); % Calculate cumulative returns
        end

        % Plot Cumulative Returns
        function plotCumulativeReturns(obj)
            figure;
            plot(obj.cumulativeReturns);
            title('Cumulative Returns of Mean Reversion Strategy');
            xlabel('Days');
            ylabel('Cumulative Returns');
            legend(arrayfun(@(x) sprintf('Asset %d', x), 1:obj.numAssets, 'UniformOutput', false));
            grid on;
        end

        % Plot Price Evolution
        function plotPrices(obj)
            figure;
            plot(obj.prices);
            title('Simulated Prices of Assets (Brownian Motion)');
            xlabel('Days');
            ylabel('Price');
            legend(arrayfun(@(x) sprintf('Asset %d', x), 1:obj.numAssets, 'UniformOutput', false));
            grid on;
        end
    end
end
