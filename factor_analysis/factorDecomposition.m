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
    corrMatrix = corr(mktRtnsLookbackAdj);
    
    % Check Model Type
    if params.modelType == "PCA"
        % Eigenvalue Decomposition
        [eigVecs, eigVals] = eig(corrMatrix);
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
            T - volLookback + 1: ...
            T, :);
        % Computing Vol over Lookback period
        factorVols = std(rtnsAdjVolLookback);

    elseif params.modelType == "PAF"
        % The PAF Model adopts an iterative process for SVD with a reduced
        % correlation matrix. If we define our correlation matrix as C, we
        % define the reduced correlation matrix as the correlation matrix
        % but with diagonal elements replaced with 'communalities'.
        % Communalities are essentially the portion of each variable's
        % variance that can be explained by other common factors. We, thus,
        % reduce the correlation matrix iteratively in order to account for
        % these communalities and remove the 'unique variance' that each
        % variable has, taking into account only 'common variances'. We
        % start with an initial estimate of communalities as the square
        % multiple correlation for each of the underlying variable,
        % reduce the correlation matrix using these communalities, apply
        % SVD to get factor loadings, calculate subsequent communalities
        % as the sum of squared loadings, and iteratively repeat
        % this process until we get a stable solution for communalities.
        % This function stops the iterative process when the max of the
        % mabsolute value of the difference of communalities between
        % the current and previous iteration is < 10^-3

        % Estimate initial communalities as squared multiple correlation
        u_prev = sum(corrMatrix.^2) - 1;
        % Adjust diagonal of the correlation matrix with communalities
        reducedMatrix = corrMatrix - eye(size(corrMatrix)) + diag(u_prev);
        % Applying SVD
        [eigenVectors, eigenValues] = eig(reducedMatrix);
        eigenValues = diag(eigenValues);
        % Sorting in order of decreasing eigenvalues
        [eigenValues, sortIdx] = sort(eigenValues, 'descend'); %#ok<ASGLU>
        eigenVectors = eigenVectors(:, sortIdx);
        % Update the estimated unique variances as 1 - the sum of squared 
        % factor loadings
        u_curr = sum(eigenVectors.^2);
        u_curr = u_curr';
        comm_diff = abs(u_curr - u_prev);
        % Iteratively applying SVD to reduced correlation matrix until the
        % max of absolute difference between subsequent unique variances is 
        % less than 10^-3
        while max(comm_diff) > 10^-3
            u_prev = u_curr;
            reducedMatrix = corrMatrix - eye( ...
                size(corrMatrix)) + diag(u_prev);
            [eigenVectors, eigenValues] = eig(reducedMatrix);
            eigenValues = diag(eigenValues);
            [eigenValues, sortIdx] = sort(eigenValues, 'descend'); %#ok<ASGLU
            eigenVectors = eigenVectors(:, sortIdx);
            positiveEigenValues = max(eigenValues, 0);
            eigenVectors = eigenVectors*diag(sqrt(positiveEigenValues));
            u_curr = sum(eigenVectors.^2, 2);
            comm_diff = abs(u_curr - u_prev);
        end
        % Computing the first 'k' factors
        factorLoadings = eigenVectors(:, 1:k);
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

    else
        ME = MException('Model Type %s is not implemented', ...
            params.modelType);
        throw(ME);

    end
end