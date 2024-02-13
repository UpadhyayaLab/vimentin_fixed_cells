function [COF, FDD, FDD_norm] = COF_2D(I, mask, psize)

% The FDD (Fluorescence Dispersion Distance) provides a measure of
% clustering. Here, I is the 2D slice under consideration (vimentin), and mask is a 2D
% mask used for normalization (here actin).
% Intended to be used on vimentin and actin slices a plane above the
% synapse.

I = I.*mask;
thresh = 10; % arbitrary
I(I < thresh) = 0;

[X, Y] = ndgrid(1:size(I, 1), 1:size(I, 2));
% convert to microns
X = X*psize; Y = Y*psize;

xCOF = sum(X.*I, "all") / sum(I(:));
yCOF = sum(Y.*I, "all") / sum(I(:));

COF = [xCOF yCOF];

distances = sqrt((X - xCOF).^2 + (Y - yCOF).^2);

FDD = sum(distances.*I, "all")/sum(I(:));

% normalize FDD using mask
xCOF_for_norm = sum(X.*mask, "all") / sum(mask(:)); % x-centroid of mask
yCOF_for_norm = sum(Y.*mask, "all") / sum(mask(:)); % y-centroid of mask
distances_for_norm = sqrt((X - xCOF_for_norm).^2 + (Y - yCOF_for_norm).^2);
FDD_norm = FDD/(sum(distances_for_norm.*mask, "all")/sum(mask(:)));
end