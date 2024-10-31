function visualizeFactorLoadingsHeat(factorLoadings, k, rotationType, numVariablesToShow, beforeAfterRotation)
    % higher magnitude = variable has greater positive relationship with that factor 
    % Sort absolute value of factor loadings to determine most significant
    [~, sortIdx] = sort(abs(factorLoadings), 'descend');
    % Select top numVariablesToShow based on absolute loadings
    selectedIdx = sortIdx(1:min(numVariablesToShow, size(factorLoadings, 1)));
    subsetFactorLoadings = factorLoadings(selectedIdx, :);
    % Create a heatmap for the factor loadings
    figure;
    imagesc(subsetFactorLoadings(:, 1:min(k, size(subsetFactorLoadings, 2))))  % Adjusted for k and dimensions
    colorbar;
    title(['Factor Loadings Heatmap - ' rotationType ' - ' beforeAfterRotation]);
    xlabel('Factors');
    ylabel('Variables');
    set(gca, 'YTick', 1:length(selectedIdx), 'YTickLabel', selectedIdx);
    % Optional: Set X-tick labels for factors if needed
    xticks(1:min(k, size(subsetFactorLoadings, 2)));
    xticklabels(arrayfun(@(x) ['Factor ' num2str(x)], 1:min(k, size(subsetFactorLoadings, 2)), 'UniformOutput', false));
end
