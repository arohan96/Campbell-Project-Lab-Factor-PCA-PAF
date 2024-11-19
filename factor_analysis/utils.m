classdef utils

    properties
    end

    methods
        
        function [factorRtns, portBetas, portVols] = rolling(obj, ...
                mktRtns, myPositions, params)
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
                [factorLoadings, factorRtns, beta, vol] = factorDecomposition( ...
                    mktRtns(1:T-rollingLookback+iii, :), ...
                    myPositions(iii, :), params);
                portBetas = [portBetas; beta];
                portVols = [portVols; vol];
                params.prevLoadings = factorLoadings;
            end
        end
        
        function Theta = computeAngleMatrix(obj, A, B)
            %% computeAngleMatrix
            % Computes the angles between each pair of column vectors in matrices 
            % A and B.
            %% Inputs:
            %   A: An m x n matrix where each column represents a vector.
            %   B: An m x n matrix where each column represents a vector.
            %   A and B must be of the same size.
            %% Output:
            %   Theta: An n x n matrix where Theta(i, j) is the angle (in degrees)
            %   between A(:, i) and B(:, j).
        
            % Compute the norms of the columns of A and B
            A_norms = sqrt(sum(A.^2, 1));  % 1 x n vector
            B_norms = sqrt(sum(B.^2, 1));  % 1 x n vector
        
            % Compute the dot products between each pair of columns
            DP = A' * B;  % n x n matrix
        
            % Compute the outer product of the norms
            NormsProduct = (A_norms') * B_norms;  % n x n matrix
        
            % Avoid division by zero by replacing zeros with eps (very small number)
            NormsProduct(NormsProduct == 0) = eps;
        
            % Compute the cosine of the angles
            CosTheta = DP ./ NormsProduct;
        
            % Clamp values to the range [-1, 1] to avoid numerical errors
            CosTheta = max(min(CosTheta, 1), -1);
        
            % Compute the angles in radians
            Theta = acos(CosTheta);
            % Convert to degrees
            Theta = rad2deg(Theta);
        end
    end
end