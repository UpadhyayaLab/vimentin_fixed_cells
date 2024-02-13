function [MFI_around_cent_1um, frac_around_cent_1um, MFI_around_cent_2um, frac_around_cent_2um, MFI_around_cent_3um, frac_around_cent_3um] = cent_assoc_3D(I, centrosome_center, mask, psize, zstep, thresh)

% I is a background-subtracted 3D image of another channel (e.g., actin or vimentin)
nslices = size(I, 3);
mask_3D = repmat(mask, [1, 1, nslices]);

% gaussian filter to reduce noise
for i = 1:nslices
    I(:, :, i) = imgaussfilt(I(:, :, i), 2);
end

I = I.*mask_3D;

% removing all signal below a certain threshold is important, as noise can
% have a major impact here if not accounted for
I(I < thresh) = 0;
I_mask_3D = I(mask_3D == 1);

[X, Y, Z] = meshgrid(1:size(I,2),1:size(I,1),1:nslices);
x = X(mask_3D == 1); y = Y(mask_3D == 1); z = Z(mask_3D == 1);
binary_coords = [y x z]; % be very careful about keeping track of what's x and what's y!

% convert to units of microns before computing distance matrix
binary_coords(:, 1:2) = psize * binary_coords(:, 1:2);
binary_coords(:, 3) = zstep * binary_coords(:, 3);

% when in a parfor loop, pdist2 gives a somewhat cryptic warning. It
% could be difficult to resolve but does not create any clear
% issues, so we will just ignore and suppress it for now.
warning('off', 'all');
distances = pdist2(centrosome_center, binary_coords);
warning('on', 'all');
idx_1um = find(distances <= 1);
idx_2um = find(distances <= 2);
idx_3um = find(distances <= 3);

% record MFI within sphere around centrosome and fraction of
% signal within this sphere
MFI_around_cent_1um = mean(I_mask_3D(idx_1um));
frac_around_cent_1um = sum(I_mask_3D(idx_1um))/sum(I_mask_3D(:));
MFI_around_cent_2um = mean(I_mask_3D(idx_2um));
frac_around_cent_2um = sum(I_mask_3D(idx_2um))/sum(I(:));
MFI_around_cent_3um = mean(I_mask_3D(idx_3um));
frac_around_cent_3um = sum(I_mask_3D(idx_3um))/sum(I_mask_3D(:));
end