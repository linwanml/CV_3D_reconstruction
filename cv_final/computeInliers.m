% function inliers = computeInliers(F,match1,match2,threshold)
% Find inliers by computing perpendicular errors between the points and the epipolar lines in each image
% To be brief, we compute the Sampson distance mentioned in the lab file.
% Input: 
%   -matrix F, matched points from image1 and image 2, and a threshold (e.g. threshold=50)
% Output: 
%   -inliers: indices of inliers
function inliers = computeInliers(F,match1,match2,threshold)

    % Calculate Sampson distance for each point
    % Compute numerator and denominator at first
    numer = diag(match2' * F * match1).^2;
    
    fm1 = F * match1;
    fm2 = F' * match2;
    denom = sum( fm1(1:2,:) .^2) + sum( fm2(1:2,:) .^2);
    sd    = numer'./denom;

    % Return inliers for which sd is smaller than threshold
    inliers = find(sd<threshold);

end
