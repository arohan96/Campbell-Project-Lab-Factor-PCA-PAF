function [eigenValues, eigenVectors, factorLoadings, communalities, explainedVariance, factorVolatility] = ...
    baseRolling(data, winSize, modelType, tolerance, iterations, nFactorsToCompute)

    rolls = size(data, 1) - winSize + 1;  % Calculate number of windows
    model = eigenValueDecomposition;      % Instantiate decomposition model
    
    % Initialize outputs as cell arrays
    eigenValues = cell(1, rolls);
    eigenVectors = cell(1, rolls);
    factorLoadings = cell(1, rolls);
    communalities = cell(1, rolls);
    explainedVariance = zeros(1, rolls);
    factorVolatility = cell(1, rolls);

    % Check model type and perform rolling decomposition
    if modelType == 'PCA'
        for i = 1 : rolls
            % Define the current window of data
            window = data(i : winSize + i - 1, :);

            % Compute correlation matrix for the current window
            model.corrMatrix = corr(window);

            % Perform PCA decomposition
            [eigenVals, eigenVecs] = model.PCA();
            eigenValues{i} = eigenVals;
            eigenVectors{i} = eigenVecs;

            % Additional metrics for PCA
            loadings = eigenVecs * diag(sqrt(eigenVals));  % Factor loadings
            factorLoadings{i} = loadings;
            communalities{i} = sum(loadings.^2, 2);        % Communalities
            explainedVariance(i) = sum(eigenVals(1:nFactorsToCompute)) / sum(eigenVals); % Explained variance
            factorVolatility{i} = std(loadings, 0, 2);     % Factor volatilities
        end

    elseif modelType == 'PAF'
        for i = 1 : rolls
            % Define the current window of data
            window = data(i : winSize + i - 1, :);

            % Compute correlation matrix for the current window
            model.corrMatrix = corr(window);

            % Perform PAF decomposition with tolerance and iterations
            [eigenVals, eigenVecs] = model.PAF(tolerance, iterations);
            eigenValues{i} = eigenVals;
            eigenVectors{i} = eigenVecs;

            % Additional metrics for PAF
            loadings = eigenVecs * diag(sqrt(max(eigenVals, 0)));  % Factor loadings with positive eigenvalues
            factorLoadings{i} = loadings;
            communalities{i} = sum(loadings.^2, 2);        % Communalities
            explainedVariance(i) = sum(eigenVals(1:nFactorsToCompute)) / sum(eigenVals); % Explained variance
            factorVolatility{i} = std(loadings, 0, 2);     % Factor volatilities
        end
    end
end
