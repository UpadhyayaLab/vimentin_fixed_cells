function [inner_outer_ratio, inner_mask_MFI, outer_mask_MFI, total_MFI] = rad_dist(I, mask, rad_fraction)
    % I should be a bg subtracted single plane image
    I = I.*mask; % don't consider signal from adjacent cells

    [inner_mask, outer_mask] = inner_from_full_mask(mask, rad_fraction);
    
    inner_outer_ratio = mean(I(inner_mask == 1), 'all') / mean(I(outer_mask == 1), 'all');

    inner_mask_MFI = mean(I(inner_mask == 1), 'all');
    outer_mask_MFI = mean(I(outer_mask == 1), 'all');
    total_MFI = mean(I(mask == 1), 'all');
end