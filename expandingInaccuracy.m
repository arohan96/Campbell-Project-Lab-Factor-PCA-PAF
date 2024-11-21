% This function takes a matrix of market returns and runs PCA and PAF on 
% the data in an expanding window fashion, saving the inaccuracy of each
% model for each window size iteratively. The function measures inaccuracy
% by the SSE between each eigenvector identified by the model and the true
% eigenvectors output by the data synthesizer. 

% The model takes in a market returns matrix with size nDays x nMkts, a
% matrix of real eigenVectors of size nMkts x nTrueFactors, an initial
% window size which dictates the minimum window size, an interval, which
% specifies the number of days to be skipped between each test, a
% tolerance, which specifies the tolerance of the iterative PAF algorithm,
% and an iterations parameter, which specifies the maximum number of
% iterations for the PAF algorithm.

% The function outputs two arrays, one containing each SSE for PCA and the
% other containing each SSE for PAF. Each element of these lists represents
% a specific windowsize which starts at initWinSize at the first index and
% increases by interval each element.

function [inaccuraciesPCA, inaccuraciesPAF] = expandingInaccuracy(data, ...
    realEigenVecs, initWinSize, interval, tolerance, iterations)
   
    rolls = floor((size(data,1) - initWinSize + 1) / interval);  %Finding number of iterations
    model = eigenValueDecomposition;   %Instantiating decomp model
    inaccuraciesPCA = zeros(1, rolls);    %Instantiating outputs
    inaccuraciesPAF = zeros(1, rolls);   %Instantiating outputs

    for xxx = 1 : rolls  %Iterating through the data
        ppp = xxx * interval;
        inaccuracyPCA = 0;
        window = data(1 : initWinSize + ppp - 1,:);   %Identifying window 
        model.corrMatrix = corr(window);    %Assigning corrMatrix
%PCA
        [eigenVals, eigenVecs] = model.PCA(); %Finding PCA decomp
        eigenVecs = normc(eigenVecs); %Normalize EigenVecs
        for yyy = 1:size(realEigenVecs, 2)
            if dot(eigenVecs(:, yyy), realEigenVecs(:, yyy)) < 0
                eigenVecs(:, yyy) = eigenVecs(:, yyy) * (-1);
            end
            incr = (norm(eigenVecs(:, yyy) - realEigenVecs(:, yyy)))^2;
            inaccuracyPCA = inaccuracyPCA + incr;
        end
        inaccuraciesPCA(xxx) = inaccuracyPCA;
%PAF
        inaccuracyPAF = 0;
        %Finding PAF decomposition
        [eigenVals, eigenVecs] = model.PAF(tolerance, iterations);
        eigenVecs = normc(eigenVecs); %Normalize EigenVecs
        for yyy = 1:size(realEigenVecs, 2)
            if dot(eigenVecs(:, yyy), realEigenVecs(:, yyy)) < 0
                eigenVecs(:, yyy) = eigenVecs(:, yyy) * (-1);
            end
            incr = (norm(eigenVecs(:, yyy) - realEigenVecs(:, yyy)))^2;
            inaccuracyPAF = inaccuracyPAF + incr;
        end
        inaccuraciesPAF(xxx) = inaccuracyPAF;
    end
end