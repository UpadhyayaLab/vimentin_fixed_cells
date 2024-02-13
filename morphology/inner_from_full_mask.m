function [finalInnerMask, finalOuterMask] = inner_from_full_mask(mask, radial_fraction)

% Given a binary image mask and a radial fraction for erosion, this
% function generates a new, scaled-down mask. The number of connected
% components of the output will be the same as of the input. 
% 
% The centroids of the connected components are fixed upon 
% the transformation, and each connected component will be scaled by a
% factor radial_fraction.
%
% Note: this function is for downscaling only, i.e., radial_fraction
% should lie in [0, 1].

% image with multiple connected components for testing
% mask = zeros(100,100);
% mask(10:30,10:30) = 1;
% mask(50:80,50:80) = 1;
% mask(60:90,10:40) = 1;

% find connected components and label matrix
cc = bwconncomp(mask);
numComponents = cc.NumObjects;
labelMat = labelmatrix(cc);

% initialize output mask
finalInnerMask = zeros(size(mask));
centroids_saved = zeros(numComponents, 2);

for i = 1:numComponents
    % extract binary mask for current connected component
    currMask = labelMat == i;
    
    % find centroid of current connected component
    centroid = regionprops(currMask, 'centroid');
    centroid = round(centroid.Centroid);
    centroids_saved(i, :) = centroid;

    % extract boundary points of current connected component
    [rows,cols] = find(bwperim(currMask));
    boundaryPts = [cols rows];
    nBoundaryPts = numel(rows);
    
    % find the boundary points of the inner mask
    centroid_rep = repmat(centroid, nBoundaryPts, 1);
    innerBoundaryPts = round(centroid_rep + radial_fraction * (boundaryPts - centroid_rep));

    % create an initial mask with zeros everywhere
    innerMask = zeros(size(mask, 1), size(mask,2));
  
    % convert the boundary coordinates to linear indices
    innerBoundaryInd = sub2ind([size(mask, 1), size(mask,2)], innerBoundaryPts(:,2), innerBoundaryPts(:,1));

    % set the initial mask to 1 at the boundary coordinates
    innerMask(innerBoundaryInd) = 1;
    
    % fill the holes inside the mask
    innerMask = imfill(innerMask, 'holes');

    % update the final mask
    finalInnerMask(innerMask == 1) = 1;
end

% compute the outer mask
finalOuterMask = mask - finalInnerMask;

% check if things are working properly
% figure;
% imshow(mask);
% hold on;
% B = bwboundaries(finalInnerMask);
% for i = 1:length(B)
%     boundary = B{i};
%     plot(boundary(:,2), boundary(:,1), 'r.');
%     plot(centroids_saved(i,1), centroids_saved(i,2), 'g*');
% end
% 
% figure;
% imshow(finalInnerMask);
end