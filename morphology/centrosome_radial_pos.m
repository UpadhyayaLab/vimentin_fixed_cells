function centrosome_r_norm = centrosome_radial_pos(actin_MIP_mask, centrosome_center, psize)

% convert centrosome_center to pixels
centrosome_center(:, 1:2) = centrosome_center(:, 1:2)/psize;

actin_centroid = regionprops(actin_MIP_mask, 'centroid');
actin_centroid = round(actin_centroid.Centroid);


%% create a mask which is 1 only at the boundary
boundary_mask = zeros(size(actin_MIP_mask));

% convert the boundary coordinates to linear indices
B = bwboundaries(actin_MIP_mask);
boundary_ind = sub2ind(size(actin_MIP_mask), B{1}(:,1), B{1}(:,2));

% set the initial mask to 1 at the boundary coordinates
boundary_mask(boundary_ind) = 1;

% dilate boundary mask just to make sure that we hit the
% boundary...if this isn't done, sometimes the ray from the actin centroid
% to the projected centrosome center won't hit the boundary.
SE = strel('disk', 2);
boundary_mask = imdilate(boundary_mask, SE);

%% draw line from actin centroid through centrosome center. Determine where it
%% intersects the boundary
lineseg = [centrosome_center(1) centrosome_center(2)] - flip(actin_centroid);
lineseg_norm = norm(lineseg);

count = 0;
while true
    point_to_check = round(flip(actin_centroid) + lineseg + count*lineseg/lineseg_norm);
    % check if this is on the boundary
    if boundary_mask(point_to_check(1), point_to_check(2))
        % assign a radius for this point
        centrosome_r_norm = norm(lineseg)/(norm(lineseg)+count);
        break
    else % keep going along this line
        count = count + 1;
    end
end

end