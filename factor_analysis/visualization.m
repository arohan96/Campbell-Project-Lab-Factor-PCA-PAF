classdef visualization
    properties
        factorLoadings
        eigenValues
        communalities
    end
    
    methods
        
        function visualizeFactorLoadings(obj, k, rotationType, ...
                numVariablesToShow, beforeAfterRotation, type)
            % higher magnitude = variable has greater positive relationship
            % with that factor. Sort absolute value of factor loadings to 
            % determine most significant
            [~, sortIdx] = sort(abs(obj.factorLoadings), 'descend');
            % Select top numVariablesToShow based on absolute loadings
            selectedIdx = sortIdx(1:min(numVariablesToShow, size( ...
                obj.factorLoadings, 1)));
            subsetFactorLoadings = obj.factorLoadings(selectedIdx, :);
            
            if type == "heat"
                obj.visualizeFactorLoadingsHeat(k, rotationType, ...
                subsetFactorLoadings, selectedIdx, beforeAfterRotation)
            elseif type == "bar"
                obj.visualizeFactorLoadingsBar(k, rotationType, ...
                subsetFactorLoadings, selectedIdx, beforeAfterRotation)
            else
                ME = MException( ...
                    'Visualization type %s is not implemented', type);
                throw(ME);
            end
        end


        function visualizeFactorLoadingsHeat(obj, k, rotationType, ...
                subsetFactorLoadings, selectedIdx, beforeAfterRotation)
            % Create a heatmap for the factor loadings
            figure;
            imagesc(subsetFactorLoadings(:, 1:min(k, ...
                size(subsetFactorLoadings, 2))))  % Adjusted for k and 
            % dimensions
            colorbar;
            title(['Factor Loadings Heatmap - ' rotationType ' - ' ...
                beforeAfterRotation]);
            xlabel('Factors');
            ylabel('Variables');
            set(gca, 'YTick', 1:length(selectedIdx), 'YTickLabel', ...
                selectedIdx);
            % Optional: Set X-tick labels for factors if needed
            xticks(1:min(k, size(subsetFactorLoadings, 2)));
            xticklabels(arrayfun(@(x) ['Factor ' num2str(x)], 1:min(k, ...
                size(subsetFactorLoadings, 2)), 'UniformOutput', false));
        end

        function visualizeFactorLoadingsBar(obj, k, rotationType, ...
                subsetFactorLoadings, selectedIdx, beforeAfterRotation)
            % Bar plot for the first few factors
            figure;
            bar(subsetFactorLoadings(:, 1:min(k, size( ...
                subsetFactorLoadings, 2))), 'grouped');
            title(['Factor Loadings Bar Chart - ' rotationType ' - ' ...
                beforeAfterRotation]);
            xlabel('Variables');
            ylabel('Loadings');
            legend(arrayfun(@(x) ['Factor ' num2str(x)], 1:min(k, ...
                size(subsetFactorLoadings, 2)), 'UniformOutput', false));
        end
        
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

        function plotCommunalities(obj)
            % Visualize Communalities
            figure;
            bar(obj.communalities);
            title('Estimated Communalities');
            xlabel('Market');
            ylabel('Communality');            
        end

    end
end