function [estFactorRtns, portBetas, factorVols] = factorDecomposition( ...
    mktRtns, myPositions, params)
    %% factorDecomposition: 
    % Take as inpiut a matrix of market returns and 
    % apply Common Factor Analysis using either Principal Component 
    % Analysis (PCA) or Principal Axis Factoring (PAF).
    %% Inputs:
    %   mktRtns: Txm matrix of market returns where T is the number of
    %   periods for which we have market returns and m is the total number 
    %   of market factors.
    %   myPositions: Static Portfolio positions. A 1xm vector indicating
    %   the portfolio's position in each of the m market factors. The 
    %   portfolio is assumed to be static across all time periods T.
    %   params: Model parameters. Should have the following variables:
    %       params.modelType: Factor decomposition model to use. Only
    %       supports PCA and PAF. Note: Only PCA has been implemented right
    %       now. TODO: Implement PAF
    %       params.nFactorsToCompute: Number of factors to compute. The
    %       first 'nFactorsToCompute' eigenvectors (sorted in
    %       descending order of corresponding eigenvalue for PCA) are
    %       considered as factor loadings.
    %       params.nDays: The total time period 'T'.
    %       params.factorConstructionLookback: Lookback period for
    %       constructing factor loadings.
    %       params.volLookback: Lookback period for computing factor
    %       volatilities.
    %% Outputs:
    %   estFactorRtns: a Txk matrix of factor returns where T is the total 
    %   number of time periods and k is the number of factor loadings.
    %   portBetas: A 1xk matrix of portfolio betas indicating the
    %   portfolio's exposure to each of the k factors.
    %   factorVols: A 1xk matrix of factor volatilities (historical).
    %   indicates the historical volatility of each factor through the 
    %   given vol lookback period.
   
    % Loading model parameters
    T = params.nDays;
    m = params.nMkts;
    k = params.nFactorsToCompute;
    factorConstructionlookback = params.factorConstructionLookback;
    volLookback = params.volLookback;

    % Adjusting for Lookback period
    mktRtnsLookbackAdj = mktRtns(T-factorConstructionlookback+1:T, :);
    % De-meaning returns
    mktRtnsLookbackAdj = mktRtnsLookbackAdj - mean(mktRtnsLookbackAdj);
    % Taking into account returns only till the factor construction
    % lookback period
    covMatrix = cov(mktRtnsLookbackAdj);
    
    % Check Model Type
    if params.modelType == "PCA"
        % Eigenvalue Decomposition
        [eigVecs, eigVals] = eig(covMatrix);
        % Sorting eigen vectors in descending order of eigen values
        [eigValsSorted, idx] = sort(diag(eigVals), 'descend'); %#ok<ASGLU>
        eigVecsSorted = eigVecs(:, idx);
        % Computing the first 'k' factors
        factorLoadings = eigVecsSorted(:, 1:k);
        
        % Compute Factor returns
        estFactorRtns = mktRtns*factorLoadings;
        % Compute factor vols
        portBetas = myPositions*factorLoadings;
        % Adjusting returns for volatility Lookback
        rtnsAdjVolLookback = estFactorRtns( ...
            factorConstructionlookback - volLookback + 1: ...
            factorConstructionlookback, :);
        % Computing Vol over Lookback period
        factorVols = std(rtnsAdjVolLookback);

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