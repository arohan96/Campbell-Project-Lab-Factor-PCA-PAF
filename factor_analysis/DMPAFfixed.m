%% params
nDays = 10000;
nMkts = 1000;
nTrueFactors = 4;
drift = 0.0001;   % 0.0001 start
maxSecondFactorSize = 0.5;
nFactorsToCompute = 6;
idioVolScaler = 0.5;
seedVal = -1;       % -1 => choose a new seed value
factorConstructionLookback = 500;
volLookback = 250;
modelType = 'PAF';
kaiserNormalizeLoadings = true;   % true or false (use kaiser normalization for loadings?)
rotationType = 'varimax';   % '', varimax, quartimax, promax, equamax, orthomax ('' = no rotation)
orthoGamma = 0.35;   % (0 < orthoGamma < 1) only used when rotationType='orthomax'. coefficient that controls the correlation target between factors. 1=varimax, 0=quartimax (1=focus on orthogonality, 0=reduce the number of significant factors rather than maintaining strict orthogonality)
builtInNormalizeLoadings = false;   % true or false (use built-in normalization in rotation function for loadings?)
visualizeBeforeAfterRotation = 'both';   % '', before, after, both ('' = none)
numVariablesToShow = 15;   % how many variables to show in visualizations

%% setup
clc
% set random seed
if -1 == seedVal
    seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');
end
disp(['Using random seed ' num2str(seedVal)]);
rng(seedVal);
% useful functions
h_deMean = @(x, dim) x - nanmean(x, dim);
h_makeRts = @(nDays, nMkts, drift) h_deMean(randn(nDays, nMkts) ./ 100, 1) + drift;

%% generate random data
randVals = maxSecondFactorSize .* rand(nTrueFactors - 1, 1);
myFactorStds = sort([1; randVals], 'descend');
myPositions = h_deMean(randn(1, nMkts), 2);

disp('Factor volatility distribution:');
disp(myFactorStds');

factorRtns = bsxfun(@times, h_makeRts( nDays, nTrueFactors, drift ), myFactorStds');
idioRtns = idioVolScaler .* h_makeRts( nDays, nMkts, drift );

myBetas = nan(nMkts, nTrueFactors);
myBetas(:, 1) = rand(nMkts, 1);
myBetas(:, 2:end) = randn(nMkts, nTrueFactors - 1);

mktRtns = factorRtns * myBetas' + idioRtns;

%% setup params struct
params.nDays = nDays;
params.nMkts = nMkts;
params.nFactorsToCompute = nFactorsToCompute;
params.factorConstructionLookback = factorConstructionLookback;
params.volLookback = volLookback;
params.modelType = modelType;
params.kaiserNormalizeLoadings = kaiserNormalizeLoadings;
params.rotationType = rotationType;
params.builtInNormalizeLoadings = builtInNormalizeLoadings;
params.visualizeBeforeAfterRotation = visualizeBeforeAfterRotation;
params.orthoGamma = orthoGamma;
params.numVariablesToShow = numVariablesToShow;

%% DM PAF function
function [estFactorRtns, portBetas, factorVols] = factorDecomposition( ...
    mktRtns, myPositions, params)

    % Adjusting for Lookback period
    mktRtnsLookbackAdj = mktRtns(params.nDays-params.factorConstructionLookback+1:params.nDays, :);
    % De-meaning returns
    mktRtnsLookbackAdj = mktRtnsLookbackAdj - mean(mktRtnsLookbackAdj);
    % Taking into account returns only till the factor construction
    % lookback period
    corrMatrix = corr(mktRtnsLookbackAdj);

    if params.modelType == "PAF"
        initialCommunalities = sum(corrMatrix.^2) - 1;

        communalities = initialCommunalities;
        reducedMatrix = corrMatrix - eye(size(corrMatrix)) + diag(communalities);

        convergenceThreshold = 1e-1;
        maxIterations = 1000;
        iter = 0;

        while iter < maxIterations
            iter = iter + 1;
            disp(iter);

            [eigenVectors, eigenValues] = eig(reducedMatrix);
            eigenValues = diag(eigenValues);

            [eigenValues, sortIdx] = sort(eigenValues, 'descend');
            eigenVectors = eigenVectors(:, sortIdx);

            factorLoadings = eigenVectors(:, 1:params.nFactorsToCompute);

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
        if params.kaiserNormalizeLoadings == true
            factorLoadings = kaiserNormalization(factorLoadings, eigenValues(1:params.nFactorsToCompute));
        end

        %% call visualization before rotation (heatmap, bar chart)
        if (params.visualizeBeforeAfterRotation == "before") || (params.visualizeBeforeAfterRotation == "both")
            visualizeLoadingsHeat(factorLoadings, params.nFactorsToCompute, params.rotationType, params.numVariablesToShow, 'before rotation');
            visualizeLoadingsBar(factorLoadings, params.nFactorsToCompute, params.rotationType, params.numVariablesToShow, 'before rotation');
        end

        %% call factor loading rotation function (optional)
        % Normalize - flag indicating whether the loadings matrix should be row-normalized for rotation.
        % If 'on' (default), rows of matrix are normalized prior to rotation to have unit Euclidean norm, and
        % unnormalized right after the rotation occurs.
        % If 'off', the raw loadings are rotated and returned.
        if params.rotationType == "varimax"
            factorLoadings = varimaxRotation(factorLoadings, params.builtInNormalizeLoadings);
        elseif params.rotationType == "quartimax"
            factorLoadings = quartimaxRotation(factorLoadings, params.builtInNormalizeLoadings);
        elseif params.rotationType == "equamax"
            factorLoadings = equamaxRotation(factorLoadings, params.builtInNormalizeLoadings);
        elseif params.rotationType == "promax"
            factorLoadings = promaxRotation(factorLoadings, params.builtInNormalizeLoadings);
        elseif params.rotationType == "orthomax"
            factorLoadings = orthomaxRotation(factorLoadings, params.builtInNormalizeLoadings, params.orthoGamma);
        elseif params.rotationType == "parsimax"
            factorLoadings = parsimaxRotation(factorLoadings, params.builtInNormalizeLoadings);
        elseif params.rotationType == ""
        else
            error('Invalid rotationType: %s\nValid types: varimax, quartimax, equamax, promax, orthomax', rotationType);
        end

        %% call visualizations after rotation (heatmap, bar chart)
        if (params.visualizeBeforeAfterRotation == "after") || (params.visualizeBeforeAfterRotation == "both")
            visualizeLoadingsHeat(factorLoadings, params.nFactorsToCompute, params.rotationType, params.numVariablesToShow, 'after rotation');
            visualizeLoadingsBar(factorLoadings, params.nFactorsToCompute, params.rotationType, params.numVariablesToShow, 'after rotation');
        end

        %% break 
        estFactorRtns = mktRtns * factorLoadings;
        portBetas = myPositions * factorLoadings;

        rtnsVolLookbackAdj = estFactorRtns(params.nDays - params.volLookback + 1:params.nDays, :);
        factorVols = std(rtnsVolLookbackAdj);

    else
        ME = MException('modelType %s is not implemented', params.modelType);
        throw(ME);
    end
end

%% run factor decomp
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


