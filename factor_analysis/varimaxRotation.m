function [rotatedLoadings] = varimaxRotation(factorLoadings, builtInNormalizeLoadings)
    % varimax: (gamma = 1)
    % correlation allowed: NO (orthogonal)
    % Designed to create completely uncorrelated factors. 
    % Maximizes the variance of squared loadings within each factor while minimizing the loadings across factors.
    % Minimizes the number of variables that have high loadings on each factor.
    % This method simplifies the interpretation of the factors.

    % varimaxRotation: Applies built-in Varimax rotation to factor loadings.
    %   Input:
    %       factorLoadings - Matrix of factor loadings (nVariables x nFactors)
    %   Output:
    %       rotatedLoadings - Varimax rotated factor loadings
    
    % Use the built-in rotatefactors with 'varimax' method
    rotatedLoadings = rotatefactors(factorLoadings, 'Method', 'varimax', 'Normalize', builtInNormalizeLoadings); % Normalize (default = on)
end
