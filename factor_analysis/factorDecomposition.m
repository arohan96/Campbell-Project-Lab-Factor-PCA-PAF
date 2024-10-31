function [estFactorRtns, portBetas, factorVols] = factorDecomposition( ...
    mktRtns, myPositions, params)
    %% factorDecomposition: Use either correlation or covariance matrix.
    % Take as input a matrix of market returns and 
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
    %       supports PCA and PAF.
    %       params.nFactorsToCompute: Number of factors to compute.
    %       params.nDays: The total time period 'T'.
    %       params.factorConstructionLookback: Lookback period for
    %       constructing factor loadings.
    %       params.volLookback: Lookback period for computing factor
    %       volatilities.
    %       params.useCorrelation: Boolean flag, if true, use correlation
    %       matrix; if false, use covariance matrix.
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
    
    % Selecting the matrix type based on the flag
    if params.useCorrelation
        matrixType = corr(mktRtnsLookbackAdj);
    else
        matrixType = cov(mktRtnsLookbackAdj);
    end
    
    % Check Model Type
    if params.modelType == "PCA"
        % Eigenvalue Decomposition
        [eigVecs, eigVals] = eig(matrixType);
        [eigValsSorted, idx] = sort(diag(eigVals), 'descend');
        eigVecsSorted = eigVecs(:, idx);
        % Computing the first 'k' factors
        factorLoadings = eigVecsSorted(:, 1:k);

        % Optional normalization for covariance matrix
        if ~params.useCorrelation
            factorLoadings = factorLoadings ./ sqrt(diag(matrixType));
        end
        
        % Compute Factor returns
        estFactorRtns = mktRtns*factorLoadings;
        % Compute portfolio betas
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
        if params.useCorrelation
            u_curr = 1 - 1 ./ diag(inv(matrixType));
        else
            u_curr = diag(matrixType);
        end
        % Iteratively applying SVD to reduced correlation matrix until the
        % max of absolute difference between subsequent communalities is 
        % less than 10^-8
        comm_diff = 1;
        while max(comm_diff) > 10^-8
            u_prev = u_curr;
            % Reduced correlation/covariance matrix
            reducedMatrix = matrixType;
            reducedMatrix(1:size(matrixType, 1) + 1:end) = u_prev;
            % Eigen decomposition
            [eigenVectors, eigenValues] = eig(reducedMatrix);
            eigenValues = diag(eigenValues);
            % Sorting in order of decreasing eigenvalues
            [eigenValues, sortIdx] = sort(eigenValues, 'descend');
            eigenVectors = eigenVectors(:, sortIdx);
            % Only positive Eigenvalues
            positiveEigenValues = max(eigenValues, 0);
            % Scaling eigenvectors with standard deviation
            eigenVectors = eigenVectors*diag(sqrt(positiveEigenValues));
            % Update communalities
            u_curr = sum(eigenVectors.^2, 2);
            comm_diff = abs(u_curr - u_prev);
        end
        % Computing the first 'k' factors
        factorLoadings = eigenVectors(:, 1:k);

        % Compute Factor returns
        estFactorRtns = mktRtns*factorLoadings;
        % Compute portfolio betas
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