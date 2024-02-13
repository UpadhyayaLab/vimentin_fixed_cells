function I = subtract_bg(I)
% estimate bg as the 25% of all voxels in the image (quite arbitrary)
bg = quantile(I(:), .25);
I = I - bg;
I(I < 0) = 0;
end