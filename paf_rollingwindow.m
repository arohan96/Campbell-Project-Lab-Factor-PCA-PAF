%% params
nDays = 10000;
nMkts = 10;
nTrueFactors = 4;
drift = 0.0001;
maxSecondFactorSize = 0.5;
nFactorsToCompute = 6;
idioVolScaler = 0.5;
seedVal = -1;       % -1 => choose a new seed value
windowSize = 100;   % Rolling window size (number of days)


%% setup
clc

% set random seed
if seedVal == -1
    seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');
end
disp(['using random seed ' num2str(seedVal)]);
rng(seedVal);

% useful functions
h_deMean = @(x, dim) x - mean(x, dim, 'omitnan');
h_makeRtns = @(nDays, nMkts, drift) h_deMean(randn(nDays, nMkts) ./ 100, 1) + drift;

function [estFactorLoadings, estCommunalities, factorVols, eigenValues, reducedMatrix] = myPaf(mktRtns, params)
    % Perform Principal Axis Factoring (PAF) on market returns

    % De-mean the market returns (assuming returns are in rows)
    mktRtns = bsxfun(@minus, mktRtns, mean(mktRtns, 1));
    
    % Compute the correlation matrix
    corrMatrix = corr(mktRtns);
    
    u_curr = 1 - 1 ./ diag(inv(corrMatrix));

    % Iteratively applying SVD to reduced correlation matrix until the
    % max of absolute difference between subsequent communalities is 
    % less than 10^-8
    comm_diff = 1;
    while max(comm_diff) > 10^-8
        u_prev = u_curr;
        % Adjust diagonal of the correlation matrix with communalities
        reducedMatrix = corrMatrix - eye(size(corrMatrix)) + diag(u_prev);
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


%% Rolling Window PAF
numWindows = nDays - windowSize + 1; % Calculate the number of windows
params.nFactorsToCompute = nFactorsToCompute;

% Preallocate arrays for storing rolling window results
factorLoadingsSeries = nan(numWindows, nMkts, nFactorsToCompute);
communalitiesSeries = nan(numWindows, nMkts);
factorVolsSeries = nan(numWindows, nFactorsToCompute);
eigenValuesSeries = nan(numWindows, nMkts);

% Run PAF for each rolling window
for i = 1:numWindows
    % Extract data for the current window
    currentWindowData = mktRtns(i:i + windowSize - 1, :);
    
    % Perform PAF on the current window
    [estFactorLoadings, estCommunalities, factorVols, eigenValues, ~] = myPaf(currentWindowData, params);
    
    % Store results for each window
    factorLoadingsSeries(i, :, :) = estFactorLoadings;  % Factor loadings for the window
    communalitiesSeries(i, :) = estCommunalities;       % Communalities for the window
    factorVolsSeries(i, :) = factorVols;                % Factor volatilities for the window
    eigenValuesSeries(i, :) = eigenValues;              % Eigenvalues for the window
end

%% Visualize Results Over Time

% Plot Factor Loadings over time for each factor
for j = 1:nFactorsToCompute
    figure;
    plot(squeeze(factorLoadingsSeries(:, :, j)));
    title(['Factor ' num2str(j) ' Loadings over Time']);
    xlabel('Window');
    ylabel('Loading Value');
    legend(arrayfun(@(x) ['Market ' num2str(x)], 1:nMkts, 'UniformOutput', false));
end

% Plot Communalities over time for each market
figure;
plot(communalitiesSeries);
title('Communalities over Time');
xlabel('Window');
ylabel('Communality');
legend(arrayfun(@(x) ['Market ' num2str(x)], 1:nMkts, 'UniformOutput', false));

% Plot Factor Volatilities over time for each factor
figure;
plot(factorVolsSeries);
title('Factor Volatilities over Time');
xlabel('Window');
ylabel('Volatility');
legend(arrayfun(@(x) ['Factor ' num2str(x)], 1:nFactorsToCompute, 'UniformOutput', false));

% Scree Plot (Eigenvalues over Time)
for i = 1:numWindows
    figure;
    plot(1:length(eigenValuesSeries(i, :)), eigenValuesSeries(i, :), '-o');
    title(['Scree Plot for Window ' num2str(i)]);
    xlabel('Factor');
    ylabel('Eigenvalue');
end
