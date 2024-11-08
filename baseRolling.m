function [eigenValues, eigenVectors] = baseRolling(data, winSize, modelType, ...
    tolerance, iterations)
   
    rolls = size(data,1) - winSize + 1;  %Finding number of iterations
    model = eigenValueDecomposition;   %Instantiating decomp model
    eigenValues = cell(1, rolls);    %Instantiating outputs
    eigenVectors = cell(1, rolls);   %Instantiating outputs
 
    if modelType == 'PCA'   %Checking model type
        for i = 1 : rolls;  %Iterating through the data
            window = data(i : winSize + i - 1,:);   %Identifying window 
            model.corrMatrix = corr(window);    %Assigning corrMatrix
            [eigenVals, eigenVecs] = model.PCA(); %Finding PCA decomp
            eigenValues{i} = eigenVals;   %Storing outputs
            eigenVectors{i} = eigenVecs;    %Storing outputs
        end
    else if modelType == 'PAF'
        for i = 1 : rolls;  %Iterating through the data
            window = data(i : winSize + i - 1,:);   %Identifying window 
            model.corrMatrix = corr(window);    %Assigning corrMatrix
            %Finding PAF decomposition
            [eigenVals, eigenVecs] = model.PAF(tolerance, iterations);
            eigenValues{i} = eigenVals;   %Storing outputs
            eigenVectors{i} = eigenVecs;    %Storing outputs
        end
    end
end