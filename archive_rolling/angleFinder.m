function Theta = computeAngleMatrix(A, B)
    %% computeAngleMatrix
    % Computes the angles between each pair of column vectors in matrices 
    % A and B.
    %
    %% Inputs:
    %   A: An m x n matrix where each column represents a vector.
    %   B: An m x n matrix where each column represents a vector.
    %       A and B must be of the same size.
    %
    % Output:
    %   Theta: An n x n matrix where Theta(i, j) is the angle (in degrees)
    %           between A(:, i) and B(:, j).

    % Check if A and B have the same dimensions
    if ~isequal(size(A), size(B))
        error('Matrices A and B must have the same dimensions.');
    end

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
    Theta = rad2deg(Theta);
end
