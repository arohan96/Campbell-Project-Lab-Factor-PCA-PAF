classdef utils

    properties
    end

    methods
        
        function [factorRtns, portBetas, factorVols, factorLoadings] = rolling( ...
                obj, mktRtns, myPositions, params)
            %% rolling
            % Takes as input a market returns matrix and a portfolio 
            % positions vector and returns a matrix of factor returns, 
            % portfolio betas, and portfolio volatilities based on factors 
            % computed using either Principal Component Analysis (PCA) or 
            % Principal Axis Factoring (PAF)
            %% Inputs:
            %   mktRtns: Txm matrix of market returns where T is the
            %   rolling window size
            %   myPositions: Txm matrix of Portfolio positions, indicating
            %   the portfolio's position in each of the m market factors 
            %   through time.
            %   params: Model Parameters. Refer to factorDecomposition.m 
            %   for more details on params.
            %% Outputs:
            %   factorRtns: a TxkxR matrix of factor returns where k is the 
            %   number of factor loadings and R is the number of days in
            %   our sample. For each day, we get a Txk matrix of factor
            %   returns.
            %   portBetas: A Txk matrix of portfolio betas indicating the
            %   portfolio's exposure to each of the k factors.
            %   factorVols: A Txk matrix of factor volatilities 
            %   (historical). indicates the historical volatility of each 
            %   factor through the given vol lookback period.
            %   factorLoadings: An mxkxR matrix of factor loadings such
            %   that we have a matrix of factor loadings for each day the
            %   rolling script is run.
        
            
            % Set window size as the factorConstructionLookback
            windowSize = params.factorConstructionLookback;
            % Number of days to run the script for is total days minus
            % window size
            totalDays = params.nDays;
            rollingDays = totalDays - windowSize;
            % Set nDays as window size
            params.nDays = windowSize;
            % Initialise empty matrices for Betas, Vols, and Returns
            portBetas = [];
            factorVols = [];
            factorRtns = [];
            factorLoadings = [];
            
            % Iterate through rolling window
            for iii = 1:rollingDays
                % Trim market returns
                mktRet = mktRtns(1 + iii:windowSize + iii, :);
                % Run factor decomposition
                [loadings, factorRet, beta, vol] = ( ...
                    factorDecomposition( ...
                    mktRet, myPositions(iii, :), params) ...
                    );
                % Append results to corresponding matrices
                portBetas = [portBetas; beta'];
                factorVols = [factorVols; vol];
                factorRtns(:, :, iii) = factorRet;
                factorLoadings(:, :, iii) = loadings;
                % Update parameters to include previous factor loadings
                params.prevLoadings = loadings;
            end
        end
        
        function Theta = computeAngleMatrix(obj, A, B)
            %% computeAngleMatrix
            % Computes the angles between each pair of column vectors in 
            % matrices A and B.
            %% Inputs:
            %   A: An m x n matrix where each column represents a vector.
            %   B: An m x n matrix where each column represents a vector.
            %   A and B must be of the same size.
            %% Output:
            %   Theta: An n x n matrix where Theta(i, j) is the angle 
            %   (in degrees)bet ween A(:, i) and B(:, j).
        
            % Compute the norms of the columns of A and B
            A_norms = sqrt(sum(A.^2, 1));  % 1 x n vector
            B_norms = sqrt(sum(B.^2, 1));  % 1 x n vector
        
            % Compute the dot products between each pair of columns
            DP = A' * B;  % n x n matrix
        
            % Compute the outer product of the norms
            NormsProduct = (A_norms') * B_norms;  % n x n matrix
        
            % Avoid division by zero by replacing zeros with eps (very 
            % small number)
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

        function cols_in_range = findColumnsInRange(obj, matrix, ...
                lower_bound, upper_bound)
            %% findColumnsInRange 
            % finds columns where any value is within a specified range.
            %% Inputs:
            %   matrix: The input matrix.
            %   lower_bound: The lower bound of the range.
            %   upper_bound: The upper bound of the range.
            %% Output:
            %   cols_in_range: A vector containing the indices of the 
            %   columns that meet the condition.
            
            % Create a logical matrix where elements within the range are 
            % true
            within_range = matrix >= lower_bound & matrix <= upper_bound;
            
            % Check if any element in each column is within the range
            cols_logical = any(within_range, 1);
            
            % Get the indices of columns that satisfy the condition
            cols_in_range = find(cols_logical);
        end
    end
end
