% This function creates simulated market data.  Given a number of days,
% markets, an amount of drift, a number of true factors, and a percent of 
% idiosyncratic volatility, the function generates simulated market data.  
% The method outputs a matrix of market returns of size nDays x nMkts, a 
% myBetas matrix containing each market's betas to each true factor (the 
% matrix is nMkts x nTrueFactors), the true factor returns matrix which is 
% nDays x nTrueFactors, and the matrix of idiosyncratic returns which is 
% nDays x nMkts. 'Drift' specifies the amount of drift of the first factor 
% (essentially the mean of the returns). The idioVolPercent specifies the 
% amount of each market'svolatility that is purely idiosyncratic (meaning 
% orthogonal to allfactors and all other markets' idiosyncratic components).
% Setting thisparameter to 0.01 gives each market 1% idiosyncratic 
% volatility.

% This function returns orthogonal factors.
% Each market is a linear combination of the orthogonal factors plus one
% unique, orthogonal idiosyncratic component responsible for the specified
% percent of the markets volatility.


function [mktRtns, myBetas, factorRtns, idioRtns] = dataMakerVer3(nDays, nMkts, nTrueFactors, drift, idioVolPercent)
    format long
    maxSecondFactorSize = 0.5;
    seedVal = -1;       % -1 => choose a new seed value
    
    
    %% setup
    clc
    
    % set random seed
    if -1 == seedVal
        seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');
    end
    disp(['using random seed ' num2str(seedVal)]);
    rng(seedVal);

    % Generating random factor Vols with first factor having vol 1
    randVals = maxSecondFactorSize .* rand(nTrueFactors-1, 1);
    myFactorStds = sort([1; randVals], 'descend')  % first factor is much larger
    
    %% generate random data
    % Number of vectors and their dimension
    n = nMkts + nTrueFactors; % Number of vectors
    dim = nDays; % Dimension of each vector
    A = randn(dim, n); % Matrix of random vector columns

    % Demeaning the Columns
    for bbb = 1 : size(A, 2)
        A(:,bbb) = A(:,bbb) - mean(A(:,bbb));
    end

    % Adding a drift term to the first vector column
    % Drift term is scaled such that when adjusted for vol, 
    % the drift term will be correct
    A(:,1) = A(:,1) + drift * std(A(:,1)) / myFactorStds(1);

    % Saving a scale term for first factor
    scale1 = norm(A(:,1));

    % Apply Gram-Schmidt process to orthogonalize columns
    Q = zeros(dim, n);
    for i = 1:n
        v = A(:, i);
        for j = 1:i-1
            v = v - (Q(:, j)' * A(:, i)) * Q(:, j);
        end
        Q(:,i) = v / norm(v);
    end

    % Q now contains orthogonal vectors
    idioBasis = Q;
    % Applying saved scale term
    idioBasis(:,1) = idioBasis(:,1) * scale1;

    % Standardizing the vols of the Basis
    idioBasis = idioBasis ./ std(idioBasis);
    
    %% Partitioning the basis into factors and idioSyncratic components.
    
    % Setting the vols of the factor returns.
    factorRtns = [idioBasis(:,1:nTrueFactors)*diag(myFactorStds(1:nTrueFactors))];
    
    % Generating the 
    myBetas = nan(nMkts, nTrueFactors);
    myBetas(:, 1) = rand(nMkts, 1); % all markets have positive exposure
    myBetas(:, 2:end) = randn(nMkts, nTrueFactors - 1); % some positive, some negative

    % myBetas is nMkts x nFactors of betas
    % factorRtns is nDays x nFactors of factor returns
    % idioRtns is nDays x nMkts of idiosyncratic returns
    % so mktRtns is nDays x nMkts of market returns (factor + idio)
    mktRtns = factorRtns * myBetas';

    % Setting the vols of the idioSyncratic components such that they are
    % respnsible for the input percent of volatility
    idioVolScalers = std(mktRtns) * sqrt((1-(1-idioVolPercent)^2)/(1-idioVolPercent)^2);
    idioRtns =  idioBasis(:,nTrueFactors + 1 : nMkts + nTrueFactors) * diag(idioVolScalers);
    
    % Adding the idioreturns to the market returns
    mktRtns = mktRtns + idioRtns;
