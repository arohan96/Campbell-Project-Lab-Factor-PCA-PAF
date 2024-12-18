% Tests PCA vs PAF using the expanding inaccuracy metric and data generated
% from dataMakerVer3. Running this file will generate a confidence interval
% which estimates the difference in the accuracies of PCA and PAF on
% simulated data. Each confidence interval represents the SSE resulting
% from using PCA minus the SSE resulting from using PAF for a specific
% factor. The SSE for a factor is the sum of squared errors between the
% true factor loadings and the factor loadings computed by the
% decomposition model. The true factor loadings are the betas found by
% regressing the factor returns onto the market returns.
% 
% The confidence interval should be interpreted as follows: If 0 is in the
% interval, there is no evidence that PCA or PAF identified that factor
% better over the specified window size (window sizes on the x axis). If
% the confidence interval is above 0, we are 95% confident that PAF
% identified the factor more accurately than PCA.  If the interval is below
% 0, we are 95% confident that PCA identified the factor more accurately
% than PAF.
%
% The parameters for the data generation are identified by running PCA on
% russel 1000 data.  dataMakerParamGetter takes in the number of true
% factors to be used, runs PCA on the russell 1000 data which is input
% inside of dataMakerParamGetter, and outputs the factor volatilities and
% the factor exposure means and standard deviations in a matrix.
% dataMakerVer3GivenStdsAndExposureDists then takes the parameters from
% this matrix (as well as the other params specified at the beginning) and
% generates data accordingly for the confidence interval generation.
%
% This function can be altered simply to show the individual accuracies of
% PCA and PAF rather than the difference between them by altering the code
% in line 111.
%
% The function takes in several parameters which are to be specified in the
% window below. 
% 
% nMkts:
%   The number of markets to be synthesized and used by the model.
%
% nTrueFactors:
%   The number of true factors to be embedded in the data and retrieved
%   using the decomposition models.
%
% idioVol:
%   This specifies the fraction of the markets' volatilities that are due
%   to each market's idiosyncratic component
%
% drift:
%   This specifies the mean of the first factor's returns.
%
% dataDays:
%   This specifies the maximum window size to be tested.
%
% initWindow:
%   This specifies the minimum window size to be tested.
%
% interval:
%   This specifies the granularity of the test, with the minimum being a
%   one day difference between each window size.
%
% nTrials:
%   This specifies the number of trials to be used to generate each
%   confidence interval. Increasing this number will increase the inference
%   power.
%

%__________________________________________________________________________

nMkts = 10
nTrueFactors = 10;
idioVol = 0.6;
drift = 0.0015;

dataDays = 400;
initWindow = 1;
interval = 10;

nTrials = 1000

%__________________________________________________________________________

[para] = dataMakerParamGetter(nTrueFactors)
factorVols = para(:,1);
mktBetaMeans = para(:,2);
mktBetaVols = para(:,3);
warning('off')
format long
nTrueFactors = length(factorVols);
outputs = cell(nTrueFactors);
rolls = (dataDays - initWindow + 1)/ interval;
for ppp = 1:nTrueFactors
    outputs{ppp} = zeros(nTrials, rolls);
end

for rrr = 1:nTrials
    % set random seed
    rng('shuffle');
    seedVal = rng;
    [data, betas, factor_rtns, idiortns] = dataMakerVer3GivenStdsAndExposureDists(dataDays, nMkts, drift, idioVol, factorVols, mktBetaMeans, mktBetaVols, seedVal);
    
    realEigenVecsTrue = inv([data,ones(dataDays,1)]' * [data,ones(dataDays,1)]) * [data,ones(dataDays,1)]' * factor_rtns;
    realEigenVecsTrue = normc(realEigenVecsTrue(1:nMkts,:));
    realEigenVecsTrue = trueEvecSorter2(realEigenVecsTrue, data);

    [inaccuraciesPCA, inaccuraciesPAF] = expandingInaccuracyByFactor(data, ...
       realEigenVecsTrue, initWindow, interval, 0.0001, 20);
    for ppp = 1:nTrueFactors

% The next line can be modified to show the pure accuracy of one model
% rather than the difference of the models accuracies by simply changing
% PCA - PAF to just PCA or PAF.

        outputs{ppp}(rrr,:) = inaccuraciesPCA(ppp, :) - inaccuraciesPAF(ppp, :);
    
    
    end
end

CM = jet(nTrueFactors);
confidence = zeros(2,rolls);
for ppp = 1:nTrueFactors
    output = outputs{ppp};
    for rrr = 1:rolls
        point = output(:,rrr);
        se = std(point)/sqrt(nTrials);
        mu = mean(point);
        confidence(1,rrr) = mu+1.96*se;
        confidence(2,rrr) = mu-1.96*se;
    end
    x = linspace(1,dataDays, rolls);
    plot(x, confidence(1, :), LineStyle = '--', color = CM(ppp,:), DisplayName = sprintf('F%d',ppp))
    hold on
    plot(x, confidence(2, :), LineStyle='--', color = CM(ppp,:), HandleVisibility = 'off')
end

plot(x, zeros(rolls,1), color = 'black', HandleVisibility = 'off')
legend
xlabel('Window Size(Days)')
ylabel('SSE_{PCA} - SSE_{PAF}')
title('Confidence That PAF is More Accurate Than PCA')
annotation('textbox', [0.1, 0.005, 0.8, 0], 'String', strcat('nMkts: ', sprintf('%d',nMkts),', factorVols: ', sprintf(' %.4f ', factorVols), ', mktBetaMeans: ', sprintf(' %.4f ',mktBetaMeans), ', mktBetaVols: ', sprintf(' %.4f ', mktBetaVols), ', idioVol: ', sprintf(' %.2f ', idioVol), ', drift: ', sprintf(' %.4f ', drift)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'EdgeColor', 'none');
