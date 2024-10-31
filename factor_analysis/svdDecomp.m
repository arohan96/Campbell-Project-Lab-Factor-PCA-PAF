function [eigenValues, eigenVectors] = svdDecomp(corrMatrix)
%% svdDecomp
% Singular Value Decomposition for a given correlation matrix
%% Inputs
%   corrMatrix: correlation matrix of underlying factors
%% Output
%   eigenValues: Eigenvalues arranged in decreasing order of magnitude
%   eigenVectors: A matrix where each column represents an eigenvector.
%   Columns arranged in order of decreasing corresponding eigenvalue

    % Eigenvalue Decomposition
    [eigVecs, eigVals] = eig(corrMatrix);
    % Ordering eigenvalues and eigenvectors
    [eigenValues, idx] = sort(diag(eigVals), 'descend');
    eigenVectors = eigVecs(:, idx);
end