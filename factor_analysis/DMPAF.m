function [estFactorRtns, portBetas, factorVols] = factorDecomposition( ...
    mktRtns, myPositions, params)

    % Loading model parameters
    T = params.nDays;
    m = params.nMkts;
    k = params.nFactorsToCompute;
    factorConstructionLookback = params.factorConstructionLookback;
    volLookback = params.volLookback;

    % Adjusting for Lookback period
    mktRtnsLookbackAdj = mktRtns(T-factorConstructionLookback+1:T, :);
    % De-meaning returns
    mktRtnsLookbackAdj = mktRtnsLookbackAdj - mean(mktRtnsLookbackAdj);
    % Taking into account returns only till the factor construction
    % lookback period
    corrMatrix = corr(mktRtnsLookbackAdj);

    %% params
    nDays = 10000;
    nMkts = 1000;
    nTrueFactors = 4;
    drift = 0.0001;
    maxSecondFactorSize = 0.5;
    nFactorsToCompute = 6;
    idioVolScaler = 0.5;
    seedVal = -1;       % -1 => choose a new seed value
    factorConstructionLookback = 500;
    volLookback = 250;
    modelType = 'PAF';
    kaiserNormalizeLoadings = true;   % true or false (use kaiser normalization for loadings?)
    rotationType = 'varimax';   % '', varimax, quartimax, promax, equamax, orthomax ('' = no rotation)
    orthoGamma = 0.75;   % (0 < orthoGamma < 1) only used when rotationType='orthomax'. coefficient that controls the correlation target between factors. 1=varimax, 0=quartimax (1=focus on orthogonality, 0=reduce the number of significant factors rather than maintaining strict orthogonality)
    builtInNormalizeLoadings = false;   % true or false (use built-in normalization in rotation function for loadings?)
    visualizeBeforeAfterRotation = 'both';   % '', before, after, both ('' = none)
    numVariablesToShow = 15;   % how many variables to show in visualizations


    if params.modelType == "PAF"
        initialCommunalities = sum(corrMatrix.^2) - 1;

        communalities = initialCommunalities;
        reducedMatrix = corrMatrix - eye(size(corrMatrix)) + diag(communalities);

        convergenceThreshold = 1e-1;
        maxIterations = 1000;
        iter = 0;
        commHistory = zeros(maxIterations, 1);

        while iter < maxIterations
            iter = iter + 1;
            disp(iter);

            [eigenVectors, eigenValues] = eig(reducedMatrix);
            eigenValues = diag(eigenValues);

            [eigenValues, sortIdx] = sort(eigenValues, 'descend');
            eigenVectors = eigenVectors(:, sortIdx);

            factorLoadings = eigenVectors(:, 1:k);

            newCommunalities = sum(factorLoadings.^2, 2);
            
            %newCommunalities = newCommunalities'; % Transpose to a row vector for commHistory, tested, does not make a difference in output factors

            % Calculate the difference in communalities
            commDiff = max(abs(newCommunalities - communalities));

            if commDiff < convergenceThreshold
                disp(['Convergence reached after ' num2str(iter) ' iterations']);
                break;
            end

            communalities = newCommunalities;
            reducedMatrix = corrMatrix - eye(size(corrMatrix)) + diag(communalities);
        end

        %% kaiser normalization (optional, can use by itself or to prepare factors for rotation)
        if kaiserNormalizeLoadings == true
            factorLoadings = kaiserNormalization(factorLoadings, eigenValues(1:k));
        end

        %% call visualization before rotation (heatmap, bar chart)
        if (visualizeBeforeAfterRotation == "before") || (visualizeBeforeAfterRotation == "both")
            visualizeLoadingsHeat(factorLoadings, k, rotationType, numVariablesToShow, 'before rotation');
            visualizeLoadingsBar(factorLoadings, k, rotationType, numVariablesToShow, 'before rotation');
        end

        %% call factor loading rotation function (optional)
        % Normalize - flag indicating whether the loadings matrix should be row-normalized for rotation.
        % If 'on' (default), rows of matrix are normalized prior to rotation to have unit Euclidean norm, and
        % unnormalized right after the rotation occurs.
        % If 'off', the raw loadings are rotated and returned.
        if rotationType == "varimax"
            factorLoadings = varimaxRotation(factorLoadings, builtInNormalizeLoadings);
        elseif rotationType == "quartimax"
            factorLoadings = quartimaxRotation(factorLoadings, builtInNormalizeLoadings);
        elseif rotationType == "equamax"
            factorLoadings = equamaxRotation(factorLoadings, builtInNormalizeLoadings);
        elseif rotationType == "promax"
            factorLoadings = promaxRotation(factorLoadings, builtInNormalizeLoadings);
        elseif rotationType == "orthomax"
            factorLoadings = orthomaxRotation(factorLoadings, builtInNormalizeLoadings, orthoGamma);
        elseif rotationType == ""
        else
            error('Invalid rotationType: %s\nValid types: varimax, quartimax, equamax, promax, orthomax', rotationType);
        end

        %% call visualizations after rotation (heatmap, bar chart)
        if (visualizeBeforeAfterRotation == "after") || (visualizeBeforeAfterRotation == "both")
            visualizeLoadingsHeat(factorLoadings, k, rotationType, numVariablesToShow, 'after rotation');
            visualizeLoadingsBar(factorLoadings, k, rotationType, numVariablesToShow, 'after rotation');
        end

        %% break 
        estFactorRtns = mktRtns * factorLoadings;
        portBetas = myPositions * factorLoadings;

        rtnsVolLookbackAdj = estFactorRtns(T - volLookback + 1:T, :);
        factorVols = std(rtnsVolLookbackAdj);

    else
        ME = MException('modelType %s is not implemented', params.modelType);
        throw(ME);
    end
end

clc
if -1 == seedVal
    seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');
end
disp(['Using random seed ' num2str(seedVal)]);
rng(seedVal);

h_deMean = @(x, dim) x - nanmean(x, dim);
h_makeRts = @(nDays, nMkts, drift) h_deMean(randn(nDays, nMkts) ./ 100, 1) + drift;

randVals = maxSecondFactorSize .* rand(nTrueFactors - 1, 1);
myFactorStds = sort([1; randVals], 'descend');
myPositions = h_deMean(randn(1, nMkts), 2);

disp('Factor volatility distribution:');
disp(myFactorStds');

factorRtns = bsxfun(@times, h_makeRts(nDays, nTrueFactors, drift), myFactorStds');
idioRtns = idioVolScaler .* h_makeRts(nDays, nMkts, drift);

myBetas = nan(nMkts, nTrueFactors);
myBetas(:, 1) = rand(nMkts, 1);
myBetas(:, 2:end) = randn(nMkts, nTrueFactors - 1);

mktRtns = factorRtns * myBetas' + idioRtns;

%% run factor decomp
params.nDays = nDays;
params.nMkts = nMkts;
params.nFactorsToCompute = nFactorsToCompute;
params.factorConstructionLookback = factorConstructionLookback;
params.volLookback = volLookback;
params.modelType = modelType;
[estFactorRtns, portBetas, factorVols] = factorDecomposition(mktRtns, myPositions, params);

estPortVolBetas = abs(portBetas) .* factorVols;
truePortVolBetas = abs(myPositions * myBetas) .* nanstd(factorRtns);

disp('Portfolio volBetas:');
disp(estPortVolBetas(1:nTrueFactors));
disp(truePortVolBetas(1:nTrueFactors));

estNormFactorRtns = bsxfun(@rdivide, estFactorRtns, factorVols);
trueNormFactorRtns = bsxfun(@rdivide, factorRtns, nanstd(factorRtns));

flipSign = (sum(abs(estNormFactorRtns(:, 1:nTrueFactors) - trueNormFactorRtns)) > ...
            sum(abs(estNormFactorRtns(:, 1:nTrueFactors) + trueNormFactorRtns)));

factorXcmp = cat(3, (-1).^flipSign .* estNormFactorRtns(:, 1:nTrueFactors), trueNormFactorRtns);
factorXcmp = permute(factorXcmp, [1, 3, 2]);

disp('Normalized factor return differences:');
disp(rms(squeeze(factorXcmp(:, 1, :) - factorXcmp(:, 2, :))));

disp('Factor cross-correlation matrix:');
disp(corr(squeeze(factorXcmp(:, 1, :)), squeeze(factorXcmp(:, 2, :))));

figure();
for iii = 1:4
    subplot(2, 2, iii);
    scatter(factorXcmp(:, 1, iii), factorXcmp(:, 2, iii));
    title(['Factor ' num2str(iii) ' Comparison']);
    xlabel('Estimated Factor Returns');
    ylabel('True Factor Returns');
end


