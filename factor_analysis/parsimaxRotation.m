function [rotatedLoadings] = parsimaxRotation(factorLoadings, builtInNormalizeLoadings)
    % parsimaxRotation: Applies built-in Parsimax rotation to factor loadings.
    %   Input:
    %       factorLoadings - Matrix of factor loadings (nVariables x nFactors)
    %       builtInNormalizeLoadings - Boolean indicating whether to normalize loadings matrix
    %   Output:
    %       rotatedLoadings - Parsimax rotated factor loadings

    % Use the built-in rotatefactors with 'parsimax' method
    maxit = 1000;
    rotatedLoadings = rotatefactors(factorLoadings, 'Method', 'parsimax', 'Maxit', maxit, 'Normalize', builtInNormalizeLoadings);
end
