% This function runs PCA on the real asset returns from the specified excel
% file. The function then regresses each of the asset returns onto the
% space of factor returns. The function returns the volatilities of the
% factors, assets' average exposure to each factor, and the standard
% deviation of these exposures in a matrix. The first column of the matrix
% represents the factor volatilities, the second column represents the
% average exposures, and the third column of the matrix represents the
% standard deviation of the exposures.

function [outputMat] = dataMakerParamGetter(nFactorsToCompute)

    modelType = 'PCA';
    returnsPath = "C:\Users\MRaid\Downloads\Campbell\NewRollingImplementation\russell_1000_returns.xlsx"; % Flatfile for retreiving returns 
                                          % data
    
    clc;
    
    % Read returns data from Excel file
    returnsData = readtable(returnsPath, 'ReadVariableNames', true);
    
    % Extract dates from the first column
    dates = returnsData{:, 1};
    
    % Extract returns data
    mktRtns = returnsData(:, 2:end);
    mktRtns = table2array(mktRtns);
    
    % Extract tickers
    tickers = string(returnsData.Properties.VariableNames);
    tickers = tickers(2:end);
    
    % Ensure tickers are a column vector
    tickers = tickers(:);
    
    % Check if tickers correspond to rows in mktRtns
    if length(tickers) == size(mktRtns, 1)
        % Transpose mktRtns so that tickers correspond to columns
        mktRtns = mktRtns';
    elseif length(tickers) ~= size(mktRtns, 2)
        error( ...
            'Number of tickers (%d) does not match number of markets in returns data (%d)', ...
            length(tickers), size(mktRtns, 2));
    end
    
    [nDays, nMkts] = size(mktRtns);
    disp(['Number of days: ', num2str(nDays)]);
    disp(['Number of markets: ', num2str(nMkts)]);
    
    % Remove any rows with missing data (NaNs)
    nanRows = any(isnan(mktRtns), 2);
    mktRtns = mktRtns(~nanRows, :);
    dates = dates(~nanRows);
    
    factorConstructionLookback = nDays;
    volLookback = nDays;
    params.factorConstructionLookback = factorConstructionLookback;
    params.volLookback = volLookback;
    params.nFactorsToCompute = nFactorsToCompute;
    params.modelType = modelType;
    params.nDays = nDays;
    params.nMkts = nMkts;
    params.visualizeBeforeAfterRotation = '';
    params.kaiserNormalizeLoadings = false;
    
    
    [estFactorRtns, ~, ~] = factorDecomposition(mktRtns, zeros(1,nMkts), params);
    
    FRI = [estFactorRtns, ones(nDays,1)];
    betas = inv(FRI' * FRI) * FRI' * mktRtns;
    
    outputMat = zeros(nFactorsToCompute, 3);
    outputMat(:,1) = std(estFactorRtns)';
    outputMat(:,2) = mean(betas(1:nFactorsToCompute,:),2);
    outputMat(:,3) = std(betas(1:nFactorsToCompute,:), 0, 2);
end
% Compute
