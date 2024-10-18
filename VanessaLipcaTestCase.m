%% params
nDays = 10000;
nMkts = 1000;
nTrueFactors = 4;
drift = 0.0001;
maxSecondFactorSize = 0.5;
nFactorsToCompute = 6;
idioVolScaler = 0.5;
seedVal = -1;       % -1 => choose a new seed value


%% setup
clc

% set random seed
if -1 == seedVal
    seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');
end
disp(['using random seed ' num2str(seedVal)]);
rng(seedVal);

% useful functions
h_deMean = @(x, dim) x - mean(x, dim, 'omitnan');
h_makeRtns = @(nDays, nMkts, drift) h_deMean(randn(nDays, nMkts) ./ 100, 1) + drift;

function [estFactorRtns, portBetas, factorVols] = myPca(mktRtns, myPositions, params)
    % Perform Principal Component Analysis (PCA) on market returns
    
    % Step 1: De-mean the market returns (assuming returns are in rows)
    mktRtns = bsxfun(@minus, mktRtns, mean(mktRtns, 1));
    
    % Step 2: Compute the covariance matrix
    covMatrix = cov(mktRtns);
    
    % Step 3: Perform eigenvalue decomposition (PCA)
    [eigenVectors, eigenValues] = eig(covMatrix);
    eigenValues = diag(eigenValues);  % Convert to a vector
    
    % Step 4: Sort eigenvalues (and corresponding eigenvectors) in descending order
    [eigenValues, sortIdx] = sort(eigenValues, 'descend');
    eigenVectors = eigenVectors(:, sortIdx);
    
    % Step 5: Select the top 'nFactorsToCompute' factors
    nFactorsToCompute = params.nFactorsToCompute;
    estFactorRtns = mktRtns * eigenVectors(:, 1:nFactorsToCompute);
    
    % Step 6: Calculate portfolio betas and factor volatilities
    portBetas = myPositions * eigenVectors(:, 1:nFactorsToCompute);
    factorVols = sqrt(eigenValues(1:nFactorsToCompute));
    
    % Return the estimated factor returns, portfolio betas, and factor volatilities
end



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
[estFactorRtns, portBetas, factorVols] = myPca( mktRtns, myPositions, params );


%% checks
% portfolio betas
estPortVolBetas = abs(portBetas) .* factorVols;
truePortVolBetas = abs(myPositions * myBetas) .* std(factorRtns, 0, 1, 'omitnan');
disp('portfolio volBetas');
disp(estPortVolBetas(1:nTrueFactors));
disp(truePortVolBetas(1:nTrueFactors));

% normalized factor returns
estNormFactorRtns = estFactorRtns ./ factorVols';
trueNormFactorRtns = factorRtns ./ std(factorRtns, 0, 1, 'omitnan');

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
