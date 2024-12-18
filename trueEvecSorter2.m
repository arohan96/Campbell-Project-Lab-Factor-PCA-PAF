% This function sorts eigenvectors. Since for the benchmark factors (the
% ones that are already imbedded in the simulated data), there is no
% inherent ordering given, we identify their order of importance by
% comparing them to PCA output using this function. Essentially, the sorter
% takes PCA of the data and then compares the resulting eigenvectors to the
% input eigenvectors, finding the true ordering which best fits the results
% of PCA.

function [sortedVecs] = trueEvecSorter2(trueVecs, data)
    model = eigenValueDecomposition;
    model.corrMatrix = corr(data); 
    [~, eigenVecs] = model.PCA();
    lis = 1 : size(trueVecs,2);
    sortedVecs = zeros(size(trueVecs));
    for iii = 1:size(trueVecs,2)
        minErr = 1000;
        minArg = -1;
        minArgSpot = false;
        for jjj = 1:length(lis)
            err = min(norm(eigenVecs(:,iii) - trueVecs(:,lis(jjj))),norm(eigenVecs(:,iii) + trueVecs(:,lis(jjj)))) ;
            if err < minErr
                minErr = err;
                minArg = lis(jjj);
                minArgSpot = jjj;
            end
        end
        sortedVecs(:,iii) = trueVecs(:,minArg);
        lis(minArgSpot) = [];
    end
end