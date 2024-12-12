classdef MeanReversionStrategy
    properties
        window = 20;  % Rolling window size
        k = 2;        % Number of standard deviations for threshold
    end
    
    methods
        % Constructor
        function obj = MeanReversionStrategy(window, k)
            if nargin > 0
                obj.window = window;
                obj.k = k;
            end
        end

        % Preprocess Data
        % Inputs: 
        % DataFrame containing a Date column and factor return column
        % Outputs:
        % Processed DataFrame

        function data = preprocess_data(obj, data)
            % datetime
            if ~isa(data.Date, 'datetime') && isnumeric(data.Date)
                data.Date = datetime(data.Date, 'ConvertFrom', 'datenum'); 
            elseif ~isa(data.Date, 'datetime')
                error('Date column must be a numeric matrix or a datetime array.');
            end
        end

        % Indicators
        % Inputs: 
        % Preprocessed dataframe
        % factor_column: name of the column with factor returns
        % Outputs:
        % DataFrame with added columns for Mean, Std, Upper_Threshold, and Lower_Threshold
        function data = calculate_indicators(obj, data, factor_column)
            % Rolling mean and standard deviation
            factor_data = data.(factor_column);
            data.Mean = movmean(factor_data, obj.window);
            data.Std = movstd(factor_data, obj.window);
            
            % Thresholds
            data.Upper_Threshold = data.Mean + obj.k * data.Std;
            data.Lower_Threshold = data.Mean - obj.k * data.Std;
        end

        % Generate Signals
        % Inputs: 
        % Dataframe with indicators
        % factor_column: name of the column with factor returns
        % Outputs:
        % DataFrame with added Signal column
        function data = generate_signals(obj, data, factor_column)
            % Initialize Signal
            data.Signal = zeros(height(data), 1);
            
            % Generate Long and Short Signals
            factor_data = data.(factor_column);
            data.Signal(factor_data > data.Upper_Threshold) = -1; % Short signal
            data.Signal(factor_data < data.Lower_Threshold) = 1;  % Long signal
        end

        % Track Positions
        % Track positions based on the generated signals
        % Long position is held when Signal = 1
        % Position is exited when Signal = -1
        function data = track_positions(obj, data)
            % Initialize Position column
            data.Position = zeros(height(data), 1);
            
            % Update Positions
            for i = 2:height(data)
                if data.Signal(i) == 1  % Enter long position
                    data.Position(i) = 1;
                elseif data.Signal(i) == -1  % Exit position
                    data.Position(i) = 0;
                else
                    % Hold previous position
                    data.Position(i) = data.Position(i - 1);
                end
            end
        end

        % Calculate Returns
        function data = calculate_returns(obj, data, factor_column)
            % Calculate returns starting at 1 and deviating multiplicatively
            factor_data = data.(factor_column);
            data.Strategy_Returns = [1; 1 + data.Position(1:end-1) .* factor_data(2:end)];  % Incremental returns starting at 1
            data.Cumulative_Returns = cumprod(data.Strategy_Returns);  % Cumulative product of returns
        end

        % Plot Results
        function plot_results(obj, data, factor_column)
            % Plot Factor Returns, Mean, and Thresholds
            figure;
            plot(data.Date, data.(factor_column), 'DisplayName', 'Factor Returns');
            hold on;
            plot(data.Date, data.Mean, '--', 'DisplayName', 'Rolling Mean');
            plot(data.Date, data.Upper_Threshold, ':', 'DisplayName', 'Upper Threshold');
            plot(data.Date, data.Lower_Threshold, ':', 'DisplayName', 'Lower Threshold');
            title('Mean Reversion Strategy Using Factor Returns');
            xlabel('Date');
            ylabel('Factor Returns');
            legend;
            grid on;
            hold off;

            % Plot Cumulative Returns
            figure;
            plot(data.Date, data.Cumulative_Returns, 'DisplayName', 'Cumulative Returns', 'Color', 'green');
            title('Cumulative Returns of Mean Reversion Strategy');
            xlabel('Date');
            ylabel('Cumulative Returns');
            legend;
            grid on;
        end

        % Calculate Performance Metrics
        function metrics = calculate_performance_metrics(obj, data)
            % Average return, standard deviation, and Sharpe ratio
            avg_return = mean(data.Strategy_Returns - 1);  % Remove the 1 for baseline
            std_dev = std(data.Strategy_Returns - 1);  % Remove the 1 for baseline
            sharpe_ratio = avg_return / std_dev;
            
            metrics = struct(...
                'AverageReturn', avg_return, ...
                'StandardDeviation', std_dev, ...
                'SharpeRatio', sharpe_ratio);
        end

        % Apply Strategy
        function [data, performance_metrics] = apply_strategy(obj, data, factor_column, plotFlag)
            % Apply preprocessing, indicators, signals, and returns calculations
            data = obj.preprocess_data(data);
            data = obj.calculate_indicators(data, factor_column);
            data = obj.generate_signals(data, factor_column);
            data = obj.track_positions(data);
            data = obj.calculate_returns(data, factor_column);

            % Plot results
            if plotFlag
                obj.plot_results(data, factor_column);
            end

            % Calculate performance metrics
            performance_metrics = obj.calculate_performance_metrics(data);
        end
    end
end
