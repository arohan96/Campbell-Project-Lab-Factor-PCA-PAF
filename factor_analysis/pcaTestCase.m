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


%% params
nDays = 10000;
nMkts = 1000;
nTrueFactors = 4;
drift = 0.0001;
maxSecondFactorSize = 0.5;
nFactorsToCompute = 6;
idioVolScaler = 0.5;
seedVal = -1;       % -1 => choose a new seed value
modelType = 'PCA';
factorConstructionLookback = 10000;
volLookback = 10000;


%% setup
clc

% set random seed
if -1 == seedVal
    seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');
end
disp(['using random seed ' num2str(seedVal)]);
rng(seedVal);

% useful functions
h_deMean = @(x, dim) x - nanmean(x, dim);
h_makeRtns = @(nDays, nMkts, drift) h_deMean(randn(nDays, nMkts) ./ 100, 1) + drift;


%% generate random data
randVals = maxSecondFactorSize .* rand(nTrueFactors-1, 1);
myFactorStds = sort([1; randVals], 'descend');  % first factor is much larger
myPositions = h_deMean(randn(1, nMkts), 2);

disp('factor vol distribution');
disp(myFactorStds');

factorRtns = bsxfun(@times, h_makeRtns( nDays, nTrueFactors, drift ), myFactorStds');
idioRtns = idioVolScaler .* h_makeRtns( nDays, nMkts, drift );

myBetas = nan(nMkts, nTrueFactors);
myBetas(:, 1) = rand(nMkts, 1);                     % all markets have positive exposure
myBetas(:, 2:end) = randn(nMkts, nTrueFactors - 1); % some positive, some negative

% myBetas is nMkts x nFactors of betas
% factorRtns is nDays x nFactors of factor returns
% idioRtns is nDays x nMkts of idiosyncratic returns
% so mktRtns is nDays x nMkts of market returns (factor + idio)
mktRtns = factorRtns * myBetas' + idioRtns;


%% run PCA to get factors and portfolio betas
params.nFactorsToCompute = nFactorsToCompute;
params.modelType = modelType;
params.nDays = nDays;
params.nMkts = nMkts;
params.factorConstructionLookback = factorConstructionLookback;
params.volLookback = factorConstructionLookback;

[estFactorRtns, portBetas, factorVols] = factorDecomposition( mktRtns, ...
    myPositions, params );

%% checks
% portfolio betas
estPortVolBetas = abs(portBetas) .* factorVols;
truePortVolBetas = abs(myPositions * myBetas) .* nanstd(factorRtns);
disp('portfolio volBetas');
disp(estPortVolBetas(1:nTrueFactors));
disp(truePortVolBetas(1:nTrueFactors));

% normalized factor returns
estNormFactorRtns = bsxfun(@rdivide, estFactorRtns, factorVols);
trueNormFactorRtns = bsxfun(@rdivide, factorRtns, nanstd(factorRtns));

flipSign = (sum(abs(estNormFactorRtns(:, 1:nTrueFactors) - trueNormFactorRtns)) > ...
            sum(abs(estNormFactorRtns(:, 1:nTrueFactors) + trueNormFactorRtns)));

factorXcmp = cat(3, (-1).^flipSign .* estNormFactorRtns(:, 1:nTrueFactors), trueNormFactorRtns);
factorXcmp = permute( factorXcmp, [1, 3, 2] );

disp('normalized factor return diffs');
disp(rms(squeeze(factorXcmp(:, 1, :) - factorXcmp(:, 2, :))));

disp('factor cross-correlation matrix');
disp(corr(squeeze(factorXcmp(:, 1, :)), squeeze(factorXcmp(:, 2, :))));

figure();
for iii = 1:4
    subplot(2, 2, iii);
    scatter(factorXcmp(:, 1, iii), factorXcmp(:, 2, iii))
end
