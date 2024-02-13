function sig_frac = axial_signal_binning(I, mask, actin_bottom_slice_num, actin_top_slice_num, nbins)
    % I_vim should be a bg subtracted 3D image
    mask_rep = repmat(mask, 1, 1, size(I, 3));
    I = I.*mask_rep;
    
    edges = linspace(actin_bottom_slice_num, actin_top_slice_num, nbins+1);
    [N, ~, bin] = histcounts(actin_bottom_slice_num:actin_top_slice_num, edges);
    
    nslices_cell = numel(actin_bottom_slice_num:actin_top_slice_num);
    total_sig_by_slice = zeros(nslices_cell, 1);
    for i = 1:nslices_cell
        total_sig_by_slice(i) = sum(I(:, :, i+actin_bottom_slice_num-1), "all");
    end    
    total_sig_by_bin = accumarray(bin', total_sig_by_slice);
    mean_sig_by_bin = total_sig_by_bin./N; % not all bins will necessarily have the same # of slices
    sig_frac = mean_sig_by_bin/sum(mean_sig_by_bin);
end