%% params
nDays = 10000;
nMkts = 1000;
nTrueFactors = 4;
drift = 0.0001;
maxSecondFactorSize = 0.5;
nFactorsToCompute = 6;
idioVolScaler = 0.5;  % Note: you should get very good results if you set =0. This is the no-idiosyncratic vol case, so PCA should do a good job recovering factors
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


%% dylanPCA function creation
function [estFactorRtns, portBetas, factorVols] = dylanPCA(mktRtns, myPositions, params)
    X = mktRtns - mean(mktRtns); % demean the market returns

    covMatrix = cov(X);  % (X' * X) / (size(X, 1) - 1); compute the covariance matrix

    [eigenVectors, eigenValues] = eig(covMatrix);  % compute eigenvectors (principal components) and eigenvalues
    eigenValues = diag(eigenValues);  % extract eigenvalues from diag

    [eigenValues, idx] = sort(eigenValues, 'descend');  % sort the eigenvalues in descending order
    eigenVectors = eigenVectors(:, idx);

    k = params.nFactorsToCompute;  % select the top k eigenvectors
    principalComponents = eigenVectors(:, 1:k);

    estFactorRtns = X * principalComponents;  % compute factor returns

    portBetas = myPositions * principalComponents;  % compute portfolio betas

    factorVols = std(estFactorRtns);  % compute factor volatilities
end


%% run PCA to get factors and portfolio betas
params.nFactorsToCompute = nFactorsToCompute;
[estFactorRtns, portBetas, factorVols] = dylanPCA( mktRtns, myPositions, params );


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
