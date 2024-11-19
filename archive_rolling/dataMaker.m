function [mktRtns, myBetas, factorRtns, idioRtns] = dataMaker(nDays, nMkts, nTrueFactors, drift)
    maxSecondFactorSize = 0.5;
    idioVolScaler = 0.5;
    seedVal = -1;       % -1 => choose a new seed value
    
    
    %% setup
    clc
    
    % set random seed
    if -1 == seedVal
        seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');
    end
    disp(['using random seed ' num2str(seedVal)]);
    rng(seedVal);
    
    % useful functions
    h_deMean = @(x, dim) x - nanmean(x, dim);
    h_makeRtns = @(nDays, nMkts, drift) h_deMean(randn(nDays, nMkts) ./ 100, 1) + drift;
    
    %% generate random data
    randVals = maxSecondFactorSize .* rand(nTrueFactors-1, 1);
    myFactorStds = sort([1; randVals], 'descend');  % first factor is much larger
    myPositions = h_deMean(randn(1, nMkts), 2);
    
    disp('factor vol distribution');
    disp(myFactorStds');
    
    factorRtns = bsxfun(@times, h_makeRtns( nDays, nTrueFactors, drift ), myFactorStds');
    idioRtns = idioVolScaler .* h_makeRtns( nDays, nMkts, drift );
    
    myBetas = nan(nMkts, nTrueFactors);
    myBetas(:, 1) = rand(nMkts, 1);                     % all markets have positive exposure
    myBetas(:, 2:end) = randn(nMkts, nTrueFactors - 1); % some positive, some negative

    % myBetas is nMkts x nFactors of betas
    % factorRtns is nDays x nFactors of factor returns
    % idioRtns is nDays x nMkts of idiosyncratic returns
    % so mktRtns is nDays x nMkts of market returns (factor + idio)
    mktRtns = factorRtns * myBetas' + idioRtns;
    for i = 1:size(myBetas,1)
        myBetas(i,:) = myBetas(i,:) / norm(myBetas(i,:));
    end
