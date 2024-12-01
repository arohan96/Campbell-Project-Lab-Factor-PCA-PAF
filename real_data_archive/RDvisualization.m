classdef RDvisualization
    properties
        factorLoadings
        eigenValues
        communalities
        tickers
    end

    methods
        %% Plot Eigenvalues
        function plotEigenValues(obj)
            % Plotting eigenvalues
            % Scree Plot (Eigenvalue Plot)
            figure;
            plot(1:length(obj.eigenValues), obj.eigenValues, '-o');
            title('Scree Plot (Eigenvalues)');
            xlabel('Factor');
            ylabel('Eigenvalue');

            % Variance Explained by Factors
            totalVariance = sum(obj.eigenValues);
            explainedVariance = obj.eigenValues / totalVariance * 100;
            figure;
            bar(explainedVariance);
            title('Variance Explained by Each Factor');
            xlabel('Factor');
            ylabel('Percentage of Variance Explained');

        end

        %% Plot Communalities
        function plotCommunalities(obj)
            % Visualize Communalities
            figure;
            bar(obj.communalities);
            title('Estimated Communalities');
            xlabel('Market');
            ylabel('Communality');
            if ~isempty(obj.tickers)
                % Use tickers as labels
                xticks(1:length(obj.tickers));
                xticklabels(obj.tickers);
                xtickangle(90);  % Rotate labels for readability
            end

        end

        %% Visualize Factor Loadings Heatmap
        function visualizeFactorLoadingsHeat(obj, k, rotationType, numVariablesToShow, beforeAfterRotation)
            % Visualize the factor loadings as a heatmap
            % Sort variables based on maximum absolute factor loading
            [~, sortIdx] = sort(max(abs(obj.factorLoadings), [], 2), 'descend');
            % Select top variables
            selectedIdx = sortIdx(1:min(numVariablesToShow, size(obj.factorLoadings, 1)));
            subsetFactorLoadings = obj.factorLoadings(selectedIdx, 1:k);

            % Create heatmap
            figure;
            imagesc(subsetFactorLoadings);
            colorbar;
            colormap('jet');
            title(['Factor Loadings Heatmap - ' beforeAfterRotation]);
            xlabel('Factors');
            ylabel('Variables');

            % Set Y-axis labels to tickers or variable indices
            if ~isempty(obj.tickers)
                set(gca, 'YTick', 1:length(selectedIdx), 'YTickLabel', obj.tickers(selectedIdx));
            else
                set(gca, 'YTick', 1:length(selectedIdx), 'YTickLabel', selectedIdx);
            end

            % Set X-axis labels for factors
            set(gca, 'XTick', 1:k, 'XTickLabel', arrayfun(@(x) ['Factor ' num2str(x)], 1:k, 'UniformOutput', false));
            xtickangle(45);
        end

        %% Visualize Factor Loadings Bar Plot
        function visualizeFactorLoadingsBar(obj, k, rotationType, numVariablesToShow, beforeAfterRotation)
            % Visualize the factor loadings as bar plots for each factor
            % Sort variables based on maximum absolute factor loading
            [~, sortIdx] = sort(max(abs(obj.factorLoadings), [], 2), 'descend');
            % Select top variables
            selectedIdx = sortIdx(1:min(numVariablesToShow, size(obj.factorLoadings, 1)));
            subsetFactorLoadings = obj.factorLoadings(selectedIdx, 1:k);

            % Create bar plots for each factor
            for i = 1:k
                figure;
                bar(subsetFactorLoadings(:, i));
                title(['Factor ' num2str(i) ' Loadings - ' beforeAfterRotation]);
                xlabel('Variables');
                ylabel('Factor Loading');
                ylim([-1, 1]);  % Adjust y-axis limits if necessary
                grid on;

                % Set X-axis labels to tickers or variable indices
                if ~isempty(obj.tickers)
                    set(gca, 'XTick', 1:length(selectedIdx), 'XTickLabel', obj.tickers(selectedIdx));
                else
                    set(gca, 'XTick', 1:length(selectedIdx), 'XTickLabel', selectedIdx);
                end
                xtickangle(90);
            end
        end

        %% Plot Top Tickers for Each Factor
        function plotTopTickersForFactors(obj, k, numTopTickers)
            % Plot the top tickers explaining each factor
            % Inputs:
            %   k - number of factors to consider
            %   numTopTickers - number of top tickers to display for each factor

            for factorIdx = 1:k
                % Get the factor loadings for the current factor
                factorLoadings = obj.factorLoadings(:, factorIdx);
                % Get the absolute value of the loadings
                absLoadings = abs(factorLoadings);
                % Sort the loadings in descending order
                [~, sortedIdx] = sort(absLoadings, 'descend');
                % Select top tickers
                topIdx = sortedIdx(1:min(numTopTickers, length(sortedIdx)));
                topTickers = obj.tickers(topIdx);
                topLoadings = factorLoadings(topIdx);

                % Create bar plot
                figure;
                bar(topLoadings);
                title(['Top ' num2str(numTopTickers) ' Tickers for Factor ' num2str(factorIdx)]);
                xlabel('Tickers');
                ylabel('Factor Loading');
                grid on;

                % Set X-axis labels to top tickers
                set(gca, 'XTick', 1:length(topTickers), 'XTickLabel', topTickers);
                xtickangle(90);
            end
        end
    end
end
