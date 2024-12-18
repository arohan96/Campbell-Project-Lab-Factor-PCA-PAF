% This function creates simulated market data.  Given a number of days,
% markets, an amount of drift, a percent of idiosyncratic volatility, a
% list of factor volatilities, a list of mean exposures to each factor, a 
% list of exposure standard deviations for each factor, and a seed value 
% for replication, the function generates simulated market data.  The 
% method outputs a matrix of market returns of size nDays x nMkts, a 
% myBetas matrix containing each market's betas to each true factor (the 
% matrix is nMkts x nTrueFactors), the true factor returns matrix which is 
% nDays x nTrueFactors, and the matrix of idiosyncratic returns which is 
% nDays x nMkts. 'Drift' specifies the amount of drift of the first factor 
% (essentially the mean of the returns). The idioVolPercent specifies the 
% amount of each market's volatility that is purely idiosyncratic (meaning 
% orthogonal to all factors and all other markets' idiosyncratic 
% components). Setting thisparameter to 0.01 gives each market 1% 
% idiosyncratic volatility. The myFactorStds parameter is a list of the 
% volatilities of each factor. Each asset's exposures to each factor are 
% drawn from a normal distribution with means corresponding to mktBetaMeans
%  and standard deviations according to mktBetaVols (each entry of these 
% lists corresponds to the distribution of asset's exposures to that 
% particular factor. The seedVal paramter allows the user to specify the 
% seed to be used for random number replication. 

% This function returns orthogonal factors.
% Each market is a linear combination of the orthogonal factors plus one
% unique, orthogonal idiosyncratic component responsible for the specified
% percent of the markets volatility.

function [mktRtns, myBetas, factorRtns, idioRtns] = dataMakerVer3GivenStdsAndExposureDists(nDays, nMkts, drift, idioVolPercent, myFactorStds, mktBetaMeans, mktBetaVols, seedVal)
    format long 
    
    %% setup
    clc
    
    rng(seedVal);
    
    %% generate random data
    % Number of vectors and their dimension
    nTrueFactors = length(myFactorStds); % Number of true factors
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
    for jjj = 1:nTrueFactors
        myBetas(:, jjj) = normrnd(mktBetaMeans(jjj), mktBetaVols(jjj), [nMkts, 1]);
    end

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
