function [anglesMats, anglesVecs] = angleFinder(vecs)

    % This code calculates the angles between each eigenvector and the previous
    % set of eigenvectors for the rolling decomposition.  This allows us to 
    % identify when eigenvectors have "flipped" arguments or when they are
    % interfering with other eigenvectors in indeterminant situations
    
    % 'angles' shows the angle in degrees between each eigenvector and the
    % previous corresponding eigenvector.
    
    % anglesMat shows the angle in degrees between each current eigenvector and
    % all other previous eigenvectors.  
    %               This matrix is interpreted as follows:
    % THE NTH ROW of this matrix corresponds to the nth eigenvector from the 
    % previous window.
    % THE MTH COLUMN of this matrix corresponds to the mth eigenvector from the
    % current window.
    % Thus, the (n,m) entry represents the angle between the previous nth 
    % eigenvector and the current mth eigenvector.

    anglesMats = cell(1, length(vecs));
    anglesVecs = cell(1, length(vecs));
    for i = 2 : length(vecs)
        angles = zeros(1,size(vecs{i},2));
        anglesMat = (vecs{i-1}' * vecs{i});
        for j = 1 : size(vecs{i},2)
            scale = norm(vecs{i-1}(:,j)) * norm(vecs{i}(:,j));
            angleRad = acos(vecs{i-1}(:,j)' * vecs{i}(:,j) / scale);
            angles(j) = angleRad * 180 / pi;
            anglesMat(j,:) = anglesMat(j,:) / norm(vecs{i-1}(:,j));
            anglesMat(:,j) = anglesMat(:,j) / norm(vecs{i}(:,j));
        end
        anglesMats{i} = acos(anglesMat) * 180 / pi;
        anglesVecs{i} = angles;
    end
end

