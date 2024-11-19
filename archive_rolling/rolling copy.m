function [factorRtns, portBetas, portVols] = rolling(mktRtns, ...
    myPositions, params)
    %% rolling
    % Takes as input a market returns matrix and a portfolio positions 
    % vector and returns a matrix of factor returns, portfolio betas, and
    % portfolio volatilities based on factors computed using either
    % Principal Component Analysis (PCA) or Principal Axis Factoring (PAF)
    %% Inputs:
    %   mktRtns: Txm matrix of market returns where T is the number of
    %   periods for which we have market returns and m is the total number 
    %   of market factors.
    %   myPositions: T'xm matrix of Portfolio positions, indicating
    %   the portfolio's position in each of the m market factors through 
    %   time. T' refers the time-frame for which we need to run a rolling
    %   window.
    %   params: Model Parameters. Refer to factorDecomposition.m for more
    %   details on params.
    %% Outputs:
    %   factorRtns: a Txk matrix of factor returns where T is the total 
    %   number of time periods and k is the number of factor loadings.
    %   portBetas: A T'xk matrix of portfolio betas indicating the
    %   portfolio's exposure to each of the k factors.
    %   factorVols: A T'xk matrix of factor volatilities (historical).
    %   indicates the historical volatility of each factor through the 
    %   given vol lookback period.
    
    T = size(mktRtns);
    T = T(1);
    rollingLookback = size(myPositions);
    rollingLookback = rollingLookback(1);
    portBetas = [];
    portVols = [];

    for iii = 1:rollingLookback
        [factorRtns, beta, vol] = factorDecomposition(mktRtns( ...
            1:T-rollingLookback+iii, :), myPositions(iii, :), params);
        portBetas = [portBetas; beta];
        portVols = [portVols; vol];
    end
    
    

end