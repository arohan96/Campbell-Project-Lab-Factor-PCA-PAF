function [estFactorRtns, portBetas, factorVols] = factorDecomposition(mktRtns, myPositions, params)
   
    % De-meaning the market returns
    mktRtns = mktRtns - mean(mktRtns);
    % Computing the covariance matrix
    factorCovMat = cov(mktRtns);
    
    % Check Model Type
    if params.modelType == "PCA"
        % Apply PCA to compute Factor Loadings
        factorLoadings = pca(factorCovMat);
        % Consider only the first k principal components
        k = params.nFactorsToCompute;
        factorLoadings = factorLoadings(:, 1:k);
        
        % Compute Factor returns, Portfolio Betas, and Factor Vols
        estFactorRtns = mktRtns*factorLoadings;
        portBetas = myPositions*factorLoadings;
        factorVols = std(estFactorRtns);

    elseif params.modelType == "PAF"
        ME = MException('Model Type %s is not implemented', ...
            params.modelType);
        throw(ME);

    else
        ME = MException('Model Type %s is not implemented', ...
            params.modelType);
        throw(ME);

    end
end