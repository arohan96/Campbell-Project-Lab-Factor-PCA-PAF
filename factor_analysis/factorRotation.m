classdef factorRotation
    
    properties 
        factorLoadings
    end

    methods

        function [normalizedLoadings] = kaiserNormalization(obj, ...
                eigenValues)
            % kaisernormalization: Scales the factor loadings by their a
            % ssociated eigenvalues to standardize the factors.
            % This is typically applied before rotating factor loadings 
            % because the rotations work best when the loadings
            % are scaled to keep variances balanced across factors.
            % This normalization method makes the resulting factors easier
            % to interpret because it reduces the chance of
            % one factor overwhelming the others, which leads to a clearer 
            % picture of the underlying structure.
        
            % kaiserNormalization: Normalizes the factor loadings using 
            % Kaiser normalization.
            %   Input:
            %       factorLoadings - Matrix of factor loadings (
            %       nVariables x nFactors)
            %       eigenValues - Vector of eigenvalues corresponding to
            %       each factor
            %   Output:
            %       normalizedLoadings - Kaiser normalized factor loadings
        
            % Number of factors
            k = length(eigenValues);
        
            % Normalize each factor loading by the sq root of its 
            % eigenvalue
            for i = 1:k
                obj.factorLoadings(:, i) = obj.factorLoadings( ...
                    :, i) / sqrt(eigenValues(i));
            end
        
            % Assign normalized loadings to the output variable
            normalizedLoadings = obj.factorLoadings;
        end

        function [rotatedLoadings] = varimaxRotation(obj, ...
                builtInNormalizeLoadings)
            % varimax: (gamma = 1)
            % correlation allowed: NO (orthogonal)
            % Designed to create completely uncorrelated factors. 
            % Maximizes the variance of squared loadings within each factor
            % while minimizing the loadings across factors.
            % Minimizes the number of variables that have high loadings on 
            % each factor.
            % This method simplifies the interpretation of the factors.
        
            % varimaxRotation: Applies built-in Varimax rotation to factor 
            % loadings.
            %   Input:
            %       factorLoadings - Matrix of factor loadings 
            %       (nVariables x nFactors)
            %   Output:
            %       rotatedLoadings - Varimax rotated factor loadings
            
            % Use the built-in rotatefactors with 'varimax' method
            rotatedLoadings = rotatefactors(obj.factorLoadings, ...
                'Method', 'varimax', 'Normalize', ...
                builtInNormalizeLoadings); % Normalize (default = on)
        end

        function [rotatedLoadings] = parsimaxRotation(obj, ...
                builtInNormalizeLoadings)
            % parsimaxRotation: Applies built-in Parsimax rotation to 
            % factor loadings.
            %   Input:
            %       factorLoadings - Matrix of factor loadings 
            %       (nVariables x nFactors)
            %       builtInNormalizeLoadings - Boolean indicating whether 
            %       to normalize loadings matrix
            %   Output:
            %       rotatedLoadings - Parsimax rotated factor loadings
        
            % Use the built-in rotatefactors with 'parsimax' method
            maxit = 1000;
            rotatedLoadings = rotatefactors(obj.factorLoadings, ...
                'Method', 'parsimax', 'Maxit', maxit, 'Normalize', ...
                builtInNormalizeLoadings);
        end

        function [rotatedLoadings] = quartimaxRotation(obj, ...
                builtInNormalizeLoadings)
            % quartimax: (gamma = 0)
            % correlation allowed: can lead to uncorrelated factors 
            % (orthogonal) but its focus is more on reducing the number of 
            % significant factors rather than maintaining orthogonality. 
            % It could still allow for correlated factors as it focuses on 
            % factor loadings. It simplifies the loadings by maximizing the
            % variance of the squared loadings across factors, but can lead
            % to one factor dominating. tends to load all variables onto as
            % few factors as possible, aiming for simplicity minimizes the 
            % number of factors needed to explain each variable.
            % This simplifies the interpretation of the observed variables.
        
            % quartimaxRotation: Applies built-in Quartimax rotation to 
            % factor loadings.
            %   Input:
            %       factorLoadings - Matrix of factor loadings 
            %       (nVariables x nFactors)
            %   Output:
            %       rotatedLoadings - Quartimax rotated factor loadings
            
            % Use the built-in rotatefactors with 'quartimax' method
            rotatedLoadings = rotatefactors(obj.factorLoadings, 'Method', ...
                'quartimax', 'Normalize', builtInNormalizeLoadings);
        end

        function [rotatedLoadings] = promaxRotation(obj, ...
                builtInNormalizeLoadings)
            % promax:(gamma > 1)
            % correlation allowed: YES (oblique)
            % First applies an orthogonal rotation (varimax) and then 
            % obliquely transforms the factors to allow for correlation.
            % After orthogonal rotation, it raises the factor loadings to a
            % power of gamma (usually between 1-4) but default is set to 4 
            % to achieve the oblique rotation.
            % This rotation is calculated quicker than the direct oblimin 
            % method.
        
            % promaxRotation: Applies built-in Promax rotation to factor 
            % loadings.
            %   Input:
            %       factorLoadings - Matrix of factor loadings 
            %       (nVariables x nFactors)
            %   Output:
            %       rotatedLoadings - Promax rotated factor loadings
            
            % Use the built-in rotatefactors with 'promax' method
            rotatedLoadings = rotatefactors(obj.factorLoadings, ...
                'Method', 'promax', 'Normalize', builtInNormalizeLoadings);
        end

        function [rotatedLoadings] = orthomaxRotation(obj, ...
                builtInNormalizeLoadings, orthoGamma)
            % orthomax: (0 < gamma < 1)
            % An oblique rotation that allows for correlated factors.
            % This rotation optimizes factor loadings for interpretability 
            % with  an adjustable parameter gamma for controlling 
            % obliqueness.
            
            % orthomaxRotation: Applies built-in Orthomax rotation to 
            % factor loadings.
            %   Input:
            %       factorLoadings - Matrix of factor loadings 
            %       (nVariables x nFactors) builtInNormalizeLoadings - Flag
            % indicating whether to normalize loadings before rotation
            %       gamma - Coefficient for the orthomax criterion (default
            %       is 1 for varimax)
            %   Output:
            %       rotatedLoadings - Orthomax rotated factor loadings
            
            % Use the built-in rotatefactors with 'orthomax' method
            rotatedLoadings = rotatefactors(obj.factorLoadings, ...
                'Method', 'orthomax', 'Normalize', ...
                builtInNormalizeLoadings, 'Coeff', orthoGamma);
        end

        function [rotatedLoadings] = equamaxRotation(obj, ...
                builtInNormalizeLoadings)
            % Note: typically will exceed maxIterations so set maxit = 800
            % equamax: (gamma = 0.5)
            % correlation allowed: allows for some correlation (oblique) 
            % but does not emphasize as much as promax
            % a combination of varimax and quartimax, aiming to balance the
            % simplicity of factors (varimax) and the variance of loadings 
            % across both rows and columns (simplicity of variables) 
            % (quartimax) The number of variables that load highly on a 
            % factor and the number of factors needed to explain a variable
            % are minimized. 
        
            % equamaxRotation: Applies built-in Equamax rotation to factor 
            % loadings.
            %   Input:
            %       factorLoadings - Matrix of factor loadings 
            % (nVariables x nFactors)
            %   Output:
            %       rotatedLoadings - Equamax rotated factor loadings
            
            % Use the built-in rotatefactors with equamax method
            maxit = 800; % Option to set a higher iteration limit, 
            % default = 250, tested and seems to need 800+
            rotatedLoadings = rotatefactors(obj.factorLoadings, ...
                'Method', 'equamax', 'Maxit', maxit, 'Normalize', ...
                builtInNormalizeLoadings);
        end

    end
end