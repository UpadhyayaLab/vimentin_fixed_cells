function [COF, zCOF, FDD_3D, FDD_3D_norm, z_FDD] = COF_3D(I, mask, bottom_slice_num, top_slice_num, psize, zstep, mask_type, thresh)

bottom_slice_num = max([bottom_slice_num 1]);
if top_slice_num == 0
    top_slice_num = size(I, 3);
end    

if strcmp(mask_type, 'actin')
    nslices = size(I, 3);
    mask_3D = repmat(mask, [1, 1, nslices]);

    % don't consider slices below and above the cell
    if bottom_slice_num > 1
        mask_3D(:,:,1:bottom_slice_num-1) = 0;
    end
    if bottom_slice_num < nslices
        mask_3D(:,:,top_slice_num+1:end) = 0;
    end
elseif strcmp(mask_type, 'nucleus')  
    mask_3D = mask;
end

% gaussian filter to reduce noise
for i = 1:nslices
    I(:, :, i) = imgaussfilt(I(:, :, i), 2);
end    

I = I.*mask_3D;

% removing all signal below a certain threshold is important, as noise can
% have a major impact here if not accounted for
I(I < thresh) = 0;

[X, Y, Z] = ndgrid(1:size(I, 1), 1:size(I, 2), 1:size(I, 3));
% convert to microns
X = X*psize; Y = Y*psize; Z = Z*zstep;

xCOF = sum(X.*I, "all") / sum(I(:));
yCOF = sum(Y.*I, "all") / sum(I(:));
zCOF = sum(Z.*I, "all") / sum(I(:));

COF = [xCOF yCOF zCOF];

distances = sqrt((X - xCOF).^2 + (Y - yCOF).^2 + (Z - zCOF).^2);
FDD_3D = sum(distances.*I, "all")/sum(I(:));

z_distances = abs(Z-zCOF);
z_FDD = sum(z_distances.*I, "all")/sum(I(:));

% normalize FDD using mask. not very meaningful if mask is not a 3D segmentation
xCOF_for_norm = sum(X.*mask, "all") / sum(mask(:)); % x-centroid of mask
yCOF_for_norm = sum(Y.*mask, "all") / sum(mask(:)); % y-centroid of mask
zCOF_for_norm = sum(Z.*mask, "all") / sum(mask(:)); % z-centroid of mask
distances_for_norm = sqrt((X - xCOF_for_norm).^2 + (Y - yCOF_for_norm).^2 + (Z - zCOF_for_norm).^2);
FDD_3D_norm = FDD_3D/(sum(distances_for_norm.*mask, "all")/sum(mask(:)));
end