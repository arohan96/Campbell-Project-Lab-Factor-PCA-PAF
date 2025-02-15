classdef eigenValueDecomposition

    properties
        corrMatrix
    end

    methods
        
        function [eigenValues, eigenVectors] = PCA(obj)
        %% PCA
        % Singular Value Decomposition for a given correlation matrix
        %% Inputs
        %   corrMatrix: correlation matrix of underlying factors
        %% Output
        %   eigenValues: Eigenvalues arranged in decreasing order of 
        %   magnitudeg
        %   eigenVectors: A matrix where each column represents an 
        %   eigenvector.
        %   Columns arranged in order of decreasing corresponding 
        %   eigenvalue
        
            % Eigenvalue Decomposition
            [eigVecs, eigVals] = eig(obj.corrMatrix);
            % Ordering eigenvalues and eigenvectors
            [eigenValues, idx] = sort(diag(eigVals), 'descend');
            eigenVectors = eigVecs(:, idx);
        end

        function [eigenValues, eigenVectors, communalities] = PAF(obj, ...
            tolerance, iterations)
        %% PAF
        % The PAF Model adopts an iterative process for SVD with a reduced
        % correlation matrix. If we define our correlation matrix as C, we
        % define the reduced correlation matrix as the correlation matrix
        % but with diagonal elements replaced with 'communalities'.
        % Communalities are essentially the portion of each variable's
        % variance that can be explained by other common factors. We, thus,
        % reduce the correlation matrix iteratively in order to account for
        % these communalities and remove the 'unique variance' that each
        % variable has, taking into account only 'common variances'. We
        % start with an initial estimate of communalities as the square
        % multiple correlation for each of the underlying variable,
        % reduce the correlation matrix using these communalities, apply
        % SVD to get factor loadings, calculate subsequent communalities
        % as the sum of squared loadings, and iteratively repeat
        % this process until we get a stable solution for communalities.
        % This function stops the iterative process when the max of the
        % mabsolute value of the difference of communalities between
        % the current and previous iteration is < tolerance level
        % Estimate initial communalities as squared multiple correlation
        %% Inputs:
        %   corrMatrix: A square correlation matrix of underlying factors
        %   tolerance: Tolerance level to use for convergence of 
        %   communalities
        %   iterations: Maximum number of iterations to run
        %% Outputs:
        %   eigenValues: Vector of eigenvalues ordered from highest to 
        %   lowest (positive only)
        %   eigenVectors: Matrix where each column represents an eigen 
        %   vectorvcorresponding to an eigenvalue. Columns are ordered in 
        %   order of decreasing eigenvalue
        %   communalities: Final communalities amongst factors
            
            corr_size = size(obj.corrMatrix);
            % Set initial communalities as zero
            u_curr = zeros(1, corr_size(1));
            u_curr = u_curr';
            % Iteratively applying SVD to reduced correlation matrix until 
            % the max of absolute difference between subsequent 
            % communalities is less than the tolerance level
            comm_diff = 1;
            iteration = 0;
            while max(comm_diff) > tolerance && iteration <= iterations
                iteration = iteration + 1;
                u_prev = u_curr;
                % Reduced correlation matrix
                reducedMatrix = obj.corrMatrix;
                n = size(obj.corrMatrix, 1);
                diagonal_indices = sub2ind([n, n], 1:n, 1:n);
                reducedMatrix(diagonal_indices) = u_prev;
                % Eigen value decomposition
                [eigenVectors, eigenValues] = eig(reducedMatrix);
                eigenValues = diag(eigenValues);
                % Sorting in order of decreasing eigenvalues
                [eigenValues, sortIdx] = sort(eigenValues, 'descend');
                eigenVectors = eigenVectors(:, sortIdx);
                % Only positive Eigenvalues
                positiveEigenValues = max(eigenValues, 0);
                % Scaling eigenvectors with standard deviation
                eigenVectors = eigenVectors*diag(sqrt( ...
                    positiveEigenValues));
                % Update communalities
                u_curr = sum(eigenVectors.^2, 2);
                comm_diff = abs(u_curr - u_prev);
            end
            
            communalities = u_curr;
        end
    end
end