function [rotatedLoadings] = orthomaxRotation(factorLoadings, builtInNormalizeLoadings, orthoGamma)
    % orthomax: (0 < gamma < 1)
    % An oblique rotation that allows for correlated factors.
    % This rotation optimizes factor loadings for interpretability with
    % an adjustable parameter gamma for controlling obliqueness.
    
    % orthomaxRotation: Applies built-in Orthomax rotation to factor loadings.
    %   Input:
    %       factorLoadings - Matrix of factor loadings (nVariables x nFactors)
    %       builtInNormalizeLoadings - Flag indicating whether to normalize loadings before rotation
    %       gamma - Coefficient for the orthomax criterion (default is 1 for varimax)
    %   Output:
    %       rotatedLoadings - Orthomax rotated factor loadings
    
    % Use the built-in rotatefactors with 'orthomax' method
    rotatedLoadings = rotatefactors(factorLoadings, 'Method', 'orthomax', 'Normalize', builtInNormalizeLoadings, 'Coeff', orthoGamma);
end
