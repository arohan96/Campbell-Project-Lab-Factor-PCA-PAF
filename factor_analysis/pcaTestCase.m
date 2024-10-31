%% params
nDays = 10000;
nMkts = 1000;
nTrueFactors = 4;
drift = 0.0001;
maxSecondFactorSize = 10;
nFactorsToCompute = 6;
idioVolScaler = 0.5;
seedVal = -1;       % -1 => choose a new seed value
modelType = 'PAF';
factorConstructionLookback = 10000;
volLookback = 10000;
tolerance=1e-3;
iterations=100;
kaiserNormalizeLoadings = false;   % true or false (use kaiser normalization
% for loadings?)
rotationType = '';   % '', varimax, quartimax, promax, equamax, 
% orthomax ('' = no rotation)
orthoGamma = 0.35;   % (0 < orthoGamma < 1) only used when 
% rotationType='orthomax'. coefficient that controls the correlation target
% between factors. 1=varimax, 0=quartimax (1=focus on orthogonality, 
% 0=reduce the number of significant factors rather than maintaining strict
% orthogonality)
builtInNormalizeLoadings = false;   % true or false (use built-in 
% normalization in rotation function for loadings?)
visualizeBeforeAfterRotation = 'before';   % '', before, after, both 
% ('' = none)
numVariablesToShow = 15;   % how many variables to show in visualizations


%% setup
clc

% set random seed
if -1 == seedVal
    seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', ...
        '2022-01-01');
end
disp(['using random seed ' num2str(seedVal)]);
rng(seedVal);

% useful functions
h_deMean = @(x, dim) x - nanmean(x, dim);
h_makeRtns = @(nDays, nMkts, drift) h_deMean( ...
    randn(nDays, nMkts) ./ 100, 1) + drift;


%% generate random data
randVals = maxSecondFactorSize .* rand(nTrueFactors-1, 1);
myFactorStds = sort([1; randVals], 'descend');  % first factor is much 
% larger
myPositions = h_deMean(randn(1, nMkts), 2);

disp('factor vol distribution');
disp(myFactorStds');

factorRtns = bsxfun(@times, h_makeRtns( nDays, nTrueFactors, drift ), ...
    myFactorStds');
idioRtns = idioVolScaler .* h_makeRtns( nDays, nMkts, drift );

myBetas = nan(nMkts, nTrueFactors);
myBetas(:, 1) = rand(nMkts, 1);                     % all markets have 
% positive exposure
myBetas(:, 2:end) = randn(nMkts, nTrueFactors - 1); % some positive, some 
% negative

% myBetas is nMkts x nFactors of betas
% factorRtns is nDays x nFactors of factor returns
% idioRtns is nDays x nMkts of idiosyncratic returns
% so mktRtns is nDays x nMkts of market returns (factor + idio)
mktRtns = factorRtns * myBetas' + idioRtns;


%% run PCA to get factors and portfolio betas
params.nFactorsToCompute = nFactorsToCompute;
params.modelType = modelType;
params.nDays = nDays;
params.nMkts = nMkts;
params.factorConstructionLookback = factorConstructionLookback;
params.volLookback = factorConstructionLookback;
params.tolerance = tolerance;
params.iterations = iterations;
params.kaiserNormalizeLoadings = kaiserNormalizeLoadings;
params.rotationType = rotationType;
params.builtInNormalizeLoadings = builtInNormalizeLoadings;
params.visualizeBeforeAfterRotation = visualizeBeforeAfterRotation;
params.orthoGamma = orthoGamma;
params.numVariablesToShow = numVariablesToShow;

[estFactorRtns, portBetas, factorVols] = factorDecomposition( mktRtns, ...
    myPositions, params );
%% checks
% portfolio betas
estPortVolBetas = abs(portBetas) .* factorVols;
truePortVolBetas = abs(myPositions * myBetas) .* nanstd(factorRtns);
disp('portfolio volBetas');
disp(estPortVolBetas(1:nTrueFactors));
disp(truePortVolBetas(1:nTrueFactors));

% normalized factor returns
estNormFactorRtns = bsxfun(@rdivide, estFactorRtns, factorVols);
trueNormFactorRtns = bsxfun(@rdivide, factorRtns, nanstd(factorRtns));

flipSign = (sum(abs(estNormFactorRtns(:, 1:nTrueFactors) - ...
    trueNormFactorRtns)) > ...
            sum(abs(estNormFactorRtns(:, 1:nTrueFactors) + ...
            trueNormFactorRtns)));

factorXcmp = cat(3, (-1).^flipSign .* estNormFactorRtns(:, ...
    1:nTrueFactors), trueNormFactorRtns);
factorXcmp = permute( factorXcmp, [1, 3, 2] );

disp('normalized factor return diffs');
disp(rms(squeeze(factorXcmp(:, 1, :) - factorXcmp(:, 2, :))));

disp('factor cross-correlation matrix');
disp(corr(squeeze(factorXcmp(:, 1, :)), squeeze(factorXcmp(:, 2, :))));

figure();
for iii = 1:4
    subplot(2, 2, iii);
    scatter(factorXcmp(:, 1, iii), factorXcmp(:, 2, iii))
end