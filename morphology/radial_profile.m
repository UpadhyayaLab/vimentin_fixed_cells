function [r_eff, radial_bins_means_norm] = radial_profile(I, I_mask, nbins)
I = I.*I_mask; % don't consider signal from nearby cells

%% full radial profile
radial_bins = linspace(0, 1, nbins + 1);
radial_bin_centers = (radial_bins(1:end-1) + radial_bins(2:end))/2;
nbins = numel(radial_bins);
inner_masks = cell(nbins, 1);
for i = 1:nbins
    inner_masks{i} = inner_from_full_mask(I_mask, radial_bins(i));
end    

radial_bins_means_norm = zeros(nbins-1, 1);
for i = 1:nbins-1
    mask_this_radius = inner_masks{i+1}-inner_masks{i};
    radial_bins_means_norm(i) = mean(I(mask_this_radius == 1), 'all') / mean(I(I_mask == 1), 'all');
end

r_eff = sum(radial_bins_means_norm.*radial_bin_centers')/sum(radial_bins_means_norm);

end