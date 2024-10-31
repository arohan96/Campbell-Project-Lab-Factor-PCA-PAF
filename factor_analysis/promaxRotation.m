function [rotatedLoadings] = promaxRotation(factorLoadings, builtInNormalizeLoadings)
    % promax:(gamma > 1)
    % correlation allowed: YES (oblique)
    % First applies an orthogonal rotation (varimax) and then obliquely transforms the factors to allow for correlation.
    % After orthogonal rotation, it raises the factor loadings to a power of gamma
    % (usually between 1-4) but default is set to 4 to achieve the oblique rotation.
    % This rotation is calculated quicker than the direct oblimin method.

    % promaxRotation: Applies built-in Promax rotation to factor loadings.
    %   Input:
    %       factorLoadings - Matrix of factor loadings (nVariables x nFactors)
    %   Output:
    %       rotatedLoadings - Promax rotated factor loadings
    
    % Use the built-in rotatefactors with 'promax' method
    rotatedLoadings = rotatefactors(factorLoadings, 'Method', 'promax', 'Normalize', builtInNormalizeLoadings);
end
