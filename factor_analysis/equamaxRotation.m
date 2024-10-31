function [rotatedLoadings] = equamaxRotation(factorLoadings, builtInNormalizeLoadings)
    % Note: typically will exceed maxIterations so set maxit = 800
    % equamax: (gamma = 0.5)
    % correlation allowed: allows for some correlation (oblique) but does not emphasize as much as promax
    % a combination of varimax and quartimax, aiming to balance the
    % simplicity of factors (varimax) and the variance of loadings across
    % both rows and columns (simplicity of variables) (quartimax)
    % The number of variables that load highly on a factor and the number of factors needed to explain a variable are minimized. 

    % equamaxRotation: Applies built-in Equamax rotation to factor loadings.
    %   Input:
    %       factorLoadings - Matrix of factor loadings (nVariables x nFactors)
    %   Output:
    %       rotatedLoadings - Equamax rotated factor loadings
    
    % Use the built-in rotatefactors with equamax method
    maxit = 800; % Option to set a higher iteration limit, default = 250, tested and seems to need 800+
    rotatedLoadings = rotatefactors(factorLoadings, 'Method', 'equamax', 'Maxit', maxit, 'Normalize', builtInNormalizeLoadings);
end
