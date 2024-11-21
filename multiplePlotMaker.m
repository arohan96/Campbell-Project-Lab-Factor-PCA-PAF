% Tests PCA vs PAF using the expanding inaccuracy metric and data generated
% from dataMakerVer3. Running this file generates 9 plots, showing the
% performance of each model on the same data sythesized at 9 specified
% levels of idioVolPercentage. The plots show the average SSE per true
% factor - for the red and blue lines, lower is better.  The black line
% simply shows the difference between these two centered about y=0.5. If
% the black line is above 0.5, PAF is outperforming PCA. If the black line
% is below 0.5, PCA is outperforming PAF.

format long
nMkts = 10;
factorVols = [1,.25,.2];
dataDays = 100;
initWindow = 1;
interval = 1;
idioVol_range = linspace(.3, .03, 9);
%drift = 0
drift = exp(log(1.1) / 365) -1

% set random seed
seedVal = convertTo(datetime("now"), 'epochtime', 'Epoch', '2022-01-01');

for rrr = 1:length(idioVol_range)
    [data, betas, factor_rtns, idiortns] = dataMakerVer3GivenStds(dataDays, nMkts, drift, idioVol_range(rrr), factorVols, seedVal);
    %vols_of_factors_and_idios = ...
    %    ([diag([factor_rtns, idiortns]' * [factor_rtns,idiortns])/(dataDays-1)].^.5)';

    realEigenVecsTrue = inv(data' * data) * data' * factor_rtns;
    realEigenVecsTrue = normc(realEigenVecsTrue);
    
    [inaccuraciesPCA, inaccuraciesPAF] = expandingInaccuracy(data, ...
       realEigenVecsTrue, initWindow, interval, 0.0001, 20);
    
    x = linspace(initWindow,dataDays, floor((dataDays - initWindow + 1) / interval));
    subplot(3,3,rrr)
    plot(x, (inaccuraciesPCA)/nTrueFactors, color = 'blue', DisplayName = 'PCA')
    hold on
    plot(x, (inaccuraciesPAF)/nTrueFactors, color = 'red', DisplayName = 'PAF')
    hold on
    plot(x, ones(1, length(x))*.5, color = 'magenta')
    plot(x, (inaccuraciesPCA - inaccuraciesPAF)/nTrueFactors + .5, color = 'black', DisplayName = 'PAF>PCA')
    ylim([0 1]);
    title(sprintf('Average SSE: idioVol = %.4f', idioVol_range(rrr)))
    xlabel('Window Size; Number of Data Points')
    legend show
    ylabel('SSE')
end
annotation('textbox', [0.1, 0.005, 0.8, 0], 'String', strcat('Factor Vols: ', sprintf(' %.4f ', factorVols)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'EdgeColor', 'none');

