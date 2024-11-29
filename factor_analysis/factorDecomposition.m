function [factorLoadings, estFactorRtns, portBetas, factorVols] = factorDecomposition( ...
    mktRtns, myPositions, params)   
    %% factorDecomposition
    % Take as input a matrix of market returns and  apply dimensionality
    % reduction using either Principal Component Analysis (PCA) or 
    % Principal Axis Factoring (PAF).
    %% Inputs:
    %   mktRtns: Txm matrix of market returns where T is the number of
    %   periods for which we have market returns and m is the total number 
    %   of market factors.
    %   myPositions: Static Portfolio positions. A 1xm vector indicating
    %   the portfolio's position in each of the m market factors. The 
    %   portfolio is assumed to be static across all time periods T.
    %   params: Model parameters. Should have the following variables:
    %%      Mandatory Parameters
    %       params.modelType: Factor decomposition model to use. Only
    %       supports PCA and PAF.
    %       params.nFactorsToCompute: Number of factors to compute.
    %       params.nDays: The total time period 'T'.
    %       params.factorConstructionLookback: Lookback period for
    %       constructing factor loadings.
    %       params.volLookback: Lookback period for computing factor
    %       volatilities.
    %%      Optional Parameters
    %       params.tolerance: Tolerance level for convergence of
    %       communalities for PAF method. Default is 1e-3
    %       params.iterations: Number of iterations to run for PAF. Default
    %       is 100 iterations
    %       params.kaiserNormalizeLoadings: Boolean variable indicating
    %       wether to normalize factor loadings or not
    %       params.rotationType: What factor rotation type to apply for
    %       rotating factors. The supported factor rotation types are 
    %       varimax, quartimax, equamax, promax, orthomax. If an empty
    %       string is passed, no factor rotation is applied
    %       params.orthoGamma: (0 < orthoGamma < 1) only used when 
    %       rotationType='orthomax'. coefficient that controls the 
    %       correlation target between factors. 1=varimax, 0=quartimax 
    %       (1=focus on orthogonality, 0=reduce the number of significant 
    %       factors rather than maintaining strict orthogonality)
    %       params.builtInNormalizeLoadings: Boolean. Wether to use 
    %       built-in normalization in rotation function for loadings
    %       params.visualizeBeforeAfterRotation: '', before, after, both 
    %       ('' = none)
    %       params.numVariablesToShow: Number of factors to use for 
    %       visualization
    %       params.prevLoadings: An mxk matrix of factor loadings from the
    %       previous iteration.
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
    k = params.nFactorsToCompute;
    factorConstructionlookback = params.factorConstructionLookback;
    volLookback = params.volLookback;

    % Adjusting for Lookback period
    mktRtnsLookbackAdj = mktRtns(T-factorConstructionlookback+1:T, :);
    % De-meaning returns
    mktRtnsLookbackAdj = mktRtnsLookbackAdj - mean(mktRtnsLookbackAdj);
    
    % Computing correlation matrix
    corrMatrix = corr(mktRtnsLookbackAdj);

    % Setting up Eigen Value Decomposition object
    evd = eigenValueDecomposition;
    evd.corrMatrix = corrMatrix;
    
    % Check Model Type
    if params.modelType == "PCA"
        % Eigenvalue Decomposition in line with PCA
        [eigenValues, eigenVectors] = evd.PCA();

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
        % the current and previous iteration is < tolerance level
        % Estimate initial communalities as squared multiple correlation

        if isfield(params, 'iterations')
            % Use user provided number for running max iterations
            iterations = params.iterations;
        else
            iterations = 100;
        end

        if isfield(params, 'tolerance')
            % Use user provided tolerance for communality convergence
            tolerance = params.tolerance;
        else
            % default tolerance set at 1e-3
            tolerance = 1e-3;
        end
        % Computing Eigenvalues, Eigenvectors, and Communalities using
        % Principal Axis Factoring (PAF)
        [eigenValues, eigenVectors, communalities] = evd.PAF(tolerance, ...
            iterations);

    else
        ME = MException('Model Type %s is not implemented', ...
            params.modelType);
        throw(ME);
    end

    % Computing the first 'k' factors
    factorLoadings = eigenVectors(:, 1:k);

    % Check for sign indeterminancy in Factor Loadings
    factorUtils = utils;
    if isfield(params, 'prevLoadings')
        % Compute angles between previous loadings and current loadings
        angles = factorUtils.computeAngleMatrix(params.prevLoadings, ...
            factorLoadings);
        % Find all angles close to 180 degrees
        eigenVecIndices = factorUtils.findColumnsInRange(angles, 170, 190);
        % Multiplying corresponding factor loadings with -1
        factorLoadings(:, eigenVecIndices) = -1*factorLoadings( ...
            :, eigenVecIndices);
    end

    % Check if visualization is set to true
    if params.visualize == true
        % Setting up visualization object for plotting graphs
        plot = visualization;
        plot.factorLoadings = factorLoadings;
        plot.eigenValues = eigenValues;
    
        if params.modelType == "PAF"
            plot.communalities = communalities;
            plot.plotCommunalities();
        end
    
        % Plotting eigenvalues & Variance Explained by Factors
        plot.plotEigenValues();
    
        %% call visualization before rotation (heatmap, bar chart)
        if (params.visualizeBeforeAfterRotation == "before") || ( ...
                params.visualizeBeforeAfterRotation == "both")
            plot.visualizeFactorLoadings(params.nFactorsToCompute, ...
                params.rotationType, ...
                params.numVariablesToShow,  'before rotation', 'heat');
            plot.visualizeFactorLoadings(params.nFactorsToCompute, ...
                params.rotationType, params.numVariablesToShow, ...
                'before rotation', 'bar');
        end
    end
    
    % Setting up factor rotation object
    rotation = factorRotation;
    rotation.factorLoadings = factorLoadings;
    
    %% kaiser normalization (optional, can use by itself or to prepare 
    % factors for rotation)
    if params.kaiserNormalizeLoadings == true
        factorLoadings = rotation.kaiserNormalization(eigenValues(1:k));
    end

    %% call factor loading rotation function (optional)
    % Normalize - flag indicating whether the loadings matrix should be 
    % row-normalized for rotation. If 'on' (default), rows of matrix are 
    % normalized prior to rotation to have unit Euclidean norm, and
    % unnormalized right after the rotation occurs. If 'off', the raw 
    % loadings are rotated and returned.
    if isfield(params, 'rotationType')
        switch params.rotationType
            case "varimax"
                factorLoadings = rotation.varimaxRotation( ...
                    params.builtInNormalizeLoadings);
            case "quartimax"
                factorLoadings = rotation.quartimaxRotation( ...
                    params.builtInNormalizeLoadings);
            case "equamax"
                factorLoadings = rotation.equamaxRotation(...
                    params.builtInNormalizeLoadings);
            case "promax"
                factorLoadings = rotation.promaxRotation(...
                    params.builtInNormalizeLoadings);
            case "orthomax"
                factorLoadings = rotation.orthomaxRotation(...
                    params.builtInNormalizeLoadings, params.orthoGamma);
            case "parsimax"
                factorLoadings = rotation.parsimaxRotation(...
                    params.builtInNormalizeLoadings);
            case ""
            otherwise
                error( ...
                    ['Invalid rotationType: %s\nValid types: ' ...
                    'varimax, quartimax, equamax, promax, orthomax'], ...
                    rotationType);
        end
    end
    
    % Check if visualization is true
    if params.visualize == true
        %% call visualizations after rotation (heatmap, bar chart)
        if (params.visualizeBeforeAfterRotation == "after") || ( ...
                params.visualizeBeforeAfterRotation == "both")
            plot.visualizeFactorLoadings(params.nFactorsToCompute, ...
                params.rotationType, ...
                params.numVariablesToShow,  'after rotation', 'heat');
            plot.visualizeFactorLoadings(params.nFactorsToCompute, ...
                params.rotationType, params.numVariablesToShow, ...
                'after rotation', 'bar');
        end
    end

    % Compute Factor returns
    estFactorRtns = mktRtns*factorLoadings;
    % Compute portfolio betas
    portReturns = mktRtns*myPositions';
    estFactorRtnsWithIntercept = [ones( ...
        size(estFactorRtns)), estFactorRtns];
    portBetas = regress(portReturns, estFactorRtnsWithIntercept);
    % Adjusting returns for volatility Lookback
    rtnsAdjVolLookback = estFactorRtns(T - volLookback + 1:T, :);
    % Computing Vol over Lookback period
    factorVols = std(rtnsAdjVolLookback);

end