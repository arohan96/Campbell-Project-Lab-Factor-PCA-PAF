%% params
nDays = 10000;
nMkts = 10;
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

function [estFactorLoadings, estCommunalities, factorVols, eigenValues, reducedMatrix] = myPaf(mktRtns, params)
    % Perform Principal Axis Factoring (PAF) on market returns using covariance matrix

    % De-mean the market returns (assuming returns are in rows)
    mktRtns = bsxfun(@minus, mktRtns, mean(mktRtns, 1));
    
    % Compute the covariance matrix
    covMatrix = cov(mktRtns);
    
    u_curr = 1 - 1 ./ diag(inv(covMatrix));

    % Iteratively applying SVD to reduced covariance matrix until the
    % max of absolute difference between subsequent communalities is 
    % less than 10^-8
    comm_diff = 1;
    while max(comm_diff) > 10^-8
        u_prev = u_curr;
        % Adjust diagonal of the covariance matrix with communalities
        reducedMatrix = covMatrix - eye( ...
            size(covMatrix)) + diag(u_prev);
        % Applying SVD
        [eigenVectors, eigenValues] = eig(reducedMatrix);
        eigenValues = diag(eigenValues);
        % Sorting in order of decreasing eigenvalues
        [eigenValues, sortIdx] = sort(eigenValues, 'descend'); %#ok<ASGLU
        eigenVectors = eigenVectors(:, sortIdx);
        % Taking only positive Eigenvalues
        positiveEigenValues = max(eigenValues, 0);
        % Scaling eigenvectors with standard deviation
        eigenVectors = eigenVectors * diag(sqrt(positiveEigenValues));
        % Update the estimated communalities as the sum of squared 
        % factor loadings
        u_curr = sum(eigenVectors.^2, 2);
        comm_diff = abs(u_curr - u_prev);
    end
    
    % Select the top 'nFactorsToCompute' factors (loading matrix)
    nFactorsToCompute = params.nFactorsToCompute;
    estFactorLoadings = eigenVectors(:, 1:nFactorsToCompute);
    
    % Compute factor communalities (using factor loadings)
    estCommunalities = sum(estFactorLoadings.^2, 2);
    
    % Calculate factor volatilities (eigenvalues for top factors)
    factorVols = sqrt(eigenValues(1:nFactorsToCompute));
    
    % Return the estimated factor loadings, communalities, and factor volatilities
end



%% generate random data
randVals = maxSecondFactorSize .* rand(nTrueFactors-1, 1);
myFactorStds = sort([1; randVals], 'descend');  % first factor is much larger
myPositions = h_deMean(randn(1, nMkts), 2);

disp('factor vol distribution');
disp(myFactorStds');

factorRtns = bsxfun(@times, h_makeRtns(nDays, nTrueFactors, drift), myFactorStds');
idioRtns = idioVolScaler .* h_makeRtns(nDays, nMkts, drift);

myBetas = nan(nMkts, nTrueFactors);
myBetas(:, 1) = rand(nMkts, 1);                     % all markets have positive exposure
myBetas(:, 2:end) = randn(nMkts, nTrueFactors - 1); % some positive, some negative

% myBetas is nMkts x nFactors of betas
% factorRtns is nDays x nFactors of factor returns
% idioRtns is nDays x nMkts of idiosyncratic returns
% so mktRtns is nDays x nMkts of market returns (factor + idio)
mktRtns = factorRtns * myBetas' + idioRtns;

% Standardize mktRtns data before PAF
mktRtns = zscore(mktRtns);  % Standardize market returns

%% run PAF to get factor loadings, communalities, and factor volatilities
params.nFactorsToCompute = nFactorsToCompute;
[estFactorLoadings, estCommunalities, factorVols, eigenValues, reducedMatrix] = myPaf(mktRtns, params);


%% checks
disp('Estimated Factor Loadings');
disp(estFactorLoadings(:, 1:nTrueFactors));

disp('Estimated Communalities');
disp(estCommunalities);

disp('Factor Volatilities');
disp(factorVols(1:nTrueFactors));

% Visualize Communalities
figure;
bar(estCommunalities);
title('Estimated Communalities');
xlabel('Market');
ylabel('Communality');

% Scree Plot (Eigenvalue Plot)
figure;
plot(1:length(eigenValues), eigenValues, '-o');
title('Scree Plot (Eigenvalues)');
xlabel('Factor');
ylabel('Eigenvalue');

% Variance Explained by Factors
totalVariance = sum(eigenValues);
explainedVariance = eigenValues / totalVariance * 100;
figure;
bar(explainedVariance);
title('Variance Explained by Each Factor');
xlabel('Factor');
ylabel('Percentage of Variance Explained');

%% Factor Rotation (Varimax)
% Varimax rotation for better interpretability
rotatedLoadings = rotatefactors(estFactorLoadings(:, 1:nTrueFactors), 'method', 'varimax');
disp('Rotated Factor Loadings');
disp(rotatedLoadings);

%% KMO Test (Kaiser-Meyer-Olkin)
function kmoVal = kmoTest(covMatrix)
    partialCovMatrix = inv(covMatrix);
    partialCovMatrix = -partialCovMatrix ./ sqrt(diag(partialCovMatrix) * diag(partialCovMatrix)');
    kmoNumerator = sum(sum(covMatrix.^2)) - sum(diag(covMatrix).^2);
    kmoDenominator = kmoNumerator + sum(sum(partialCovMatrix.^2)) - sum(diag(partialCovMatrix).^2);
    kmoVal = kmoNumerator / kmoDenominator;
end

kmoValue = kmoTest(cov(mktRtns));
disp(['KMO Test Value: ', num2str(kmoValue)]);

%% Bartlett's Test of Sphericity
function [chiSquareValue, pValue] = bartlettTest(covMatrix, nSamples)
    nVars = size(covMatrix, 1);
    chiSquareValue = -(nSamples - 1 - (2 * nVars + 5) / 6) * log(det(covMatrix));
    df = nVars * (nVars - 1) / 2;
    pValue = 1 - chi2cdf(chiSquareValue, df);
end

[chiSquareValue, pValue] = bartlettTest(cov(mktRtns), nDays);
disp(['Bartlett''s Test Chi-Square Value: ', num2str(chiSquareValue)]);
disp(['Bartlett''s Test p-Value: ', num2str(pValue)]);

%% Residual Plot (Reproduced Covariance Matrix vs. Original)
reproducedCovMatrix = rotatedLoadings * rotatedLoadings';
residualMatrix = cov(mktRtns) - reproducedCovMatrix;
figure;
imagesc(residualMatrix);
colorbar;
title('Residual Plot (Original vs. Reproduced Covariance Matrix)');
