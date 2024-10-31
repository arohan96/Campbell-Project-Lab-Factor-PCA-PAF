function [rotatedLoadings] = quartimaxRotation(factorLoadings, builtInNormalizeLoadings)
    % quartimax: (gamma = 0)
    % correlation allowed: can lead to uncorrelated factors (orthogonal) but its focus is more on reducing the number of significant factors
    % rather than maintaining orthogonality. It could still allow for correlated factors as it focuses on factor loadings.
    % It simplifies the loadings by maximizing the variance of the squared loadings across factors, but can lead to one factor dominating.
    % tends to load all variables onto as few factors as possible, aiming for simplicity
    % minimizes the number of factors needed to explain each variable.
    % This simplifies the interpretation of the observed variables.

    % quartimaxRotation: Applies built-in Quartimax rotation to factor loadings.
    %   Input:
    %       factorLoadings - Matrix of factor loadings (nVariables x nFactors)
    %   Output:
    %       rotatedLoadings - Quartimax rotated factor loadings
    
    % Use the built-in rotatefactors with 'quartimax' method
    rotatedLoadings = rotatefactors(factorLoadings, 'Method', 'quartimax', 'Normalize', builtInNormalizeLoadings);
end
