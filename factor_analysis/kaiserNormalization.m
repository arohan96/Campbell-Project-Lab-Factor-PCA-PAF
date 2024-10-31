function [normalizedLoadings] = kaiserNormalization(factorLoadings, eigenValues)
    % kaisernormalization: Scales the factor loadings by their associated eigenvalues to standardize the factors.
    % This is typically applied before rotating factor loadings because the rotations work best when the loadings
    % are scaled to keep variances balanced across factors.
    % This normalization method makes the resulting factors easier to interpret because it reduces the chance of
    % one factor overwhelming the others, which leads to a clearer picture of the underlying structure.

    % kaiserNormalization: Normalizes the factor loadings using Kaiser normalization.
    %   Input:
    %       factorLoadings - Matrix of factor loadings (nVariables x nFactors)
    %       eigenValues - Vector of eigenvalues corresponding to each factor
    %   Output:
    %       normalizedLoadings - Kaiser normalized factor loadings

    % Number of factors
    k = length(eigenValues);

    % Normalize each factor loading by the sq root of its eigenvalue
    for i = 1:k
        factorLoadings(:, i) = factorLoadings(:, i) / sqrt(eigenValues(i));
    end

    % Assign normalized loadings to the output variable
    normalizedLoadings = factorLoadings;
end