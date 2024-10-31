function visualizeFactorLoadingsBar(factorLoadings, k, rotationType, numVariablesToShow, beforeAfterRotation)
    % higher magnitude = variable has greater positive relationship with that factor 
    % Sort absolute value of factor loadings to determine most significant
    [~, sortIdx] = sort(abs(factorLoadings), 'descend');
    % Select top numVariablesToShow based on absolute loadings
    selectedIdx = sortIdx(1:min(numVariablesToShow, size(factorLoadings, 1)));
    subsetFactorLoadings = factorLoadings(selectedIdx, :);
    % Bar plot for the first few factors
    figure;
    bar(subsetFactorLoadings(:, 1:min(k, size(subsetFactorLoadings, 2))), 'grouped');
    title(['Factor Loadings Bar Chart - ' rotationType ' - ' beforeAfterRotation]);
    xlabel('Variables');
    ylabel('Loadings');
    legend(arrayfun(@(x) ['Factor ' num2str(x)], 1:min(k, size(subsetFactorLoadings, 2)), 'UniformOutput', false));
end
