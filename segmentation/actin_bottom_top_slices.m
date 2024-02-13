function [bottom_slice_num, top_slice_num, height, actin_above_FOV] = actin_bottom_top_slices(I_actin, actin_MIP_mask, zstep, min_thresh, thresh_bottom_of_cell)
    % I_actin is a background-subtracted, non-deconvolved 3D stack
    % actin_MIP_mask is a binary mask of the X-Y actin MIP, acquired with
    % segment_actin_MIP.m, for example.
    % zstep is the distance between consecutive z-slices (often .3 μm)
    % min_thresh is the 1st percentile of all intensities within the mask
    % - the median intensity outside the mask, as computed in
    % segment_actin_MIP.m
    % thresh_bottom_of_cell determines what fraction of the pixels in the
    % filtered image corresponding to a given plane need to be positive for 
    % the plane to be considered below the cell. The value to set depends
    % on the topology of the cell and on whether you want to consider small
    % actin "pockets" to constitute a "plane". Recommended values for this
    % parameter are between .7 and .98. The higher the threshold, the
    % more sensitive the algorithm will be to small actin pockets. It may
    % be appropriate to use a higher threshold for the soft gels than for
    % the stiffer gels or glass. On stiff gels or glass, you should opt for
    % a threshold of more like .7 or .8 if you want results to align closely with the
    % "max MFI" criterion.
    %
    % bottom_slice_num is an estimate for the lowest plane of the actin
    % network
    % top_slice_num is an estimate for the highest plane of the actin
    % network
    % height (in μm) is the difference between these multiplied by the "zstep"
    % actin_above_FOV will be 1 if top_slice_num = the number of
    % z-slices. In this case, the actin appears to go above the
    % top plane of imaging, and the cell should perhaps not be considered
    % for further analysis.

    if thresh_bottom_of_cell > 1 % don't allow, this corresponds to a fraction
        error('thresh_bottom_of_cell cannot exceed 1!')
    end    

    highestI_fraction_to_consider = .1; % arbitrary, set by what has worked best in practice
    thresh_top_of_cell = .1; % arbitrary, set by what has worked best in practice
    min_frac_pix_well_above_bg = .005; % any slice with below this fraction of pixels well above bg will not be counted

    I_actin_gauss_filt = zeros(size(I_actin));
    nslices = size(I_actin, 3);
    for i = 1:nslices
        % apply gaussian filter to smooth the image
        I_actin_gauss_filt(:, :, i) = imgaussfilt(I_actin(:, :, i), 2);
    end

    % attempt to only consider the cell of interest
    nslices = size(I_actin, 3);
    actin_MIP_mask_3D = repmat(actin_MIP_mask, [1, 1, nslices]);
    I_actin_gauss_filt = I_actin_gauss_filt.*actin_MIP_mask_3D; 

    % define the "δI/δz" filter
    z_filter = (1/9)*cat(3, ones(3, 3), zeros(3, 3), -ones(3, 3));

    % convolve with this filter
    actin_zfilt = convn(I_actin_gauss_filt, z_filter, 'same');

    slice_num = nan(nslices, 1);
    num_pix_well_above_bg = nan(nslices, 1);
    actin_zfilt_highestI_prop_gr_0 = nan(nslices, 1);

    % loop through all slices for which "actin_zfilt" was computed
    % without zero padding
    for i = 2:nslices-1
        slice_num(i) = i;

        actin_gauss_filt_slice = I_actin_gauss_filt(:, :, i);
        actin_zfilt_slice = actin_zfilt(:, :, i);

        well_above_bg = actin_gauss_filt_slice > min_thresh;
        num_pix_well_above_bg(i) = sum(well_above_bg(:));

        idx = actin_gauss_filt_slice > quantile(actin_gauss_filt_slice(well_above_bg), 1-highestI_fraction_to_consider);
        actin_zfilt_highestI = actin_zfilt_slice(idx);
        actin_zfilt_highestI_prop_gr_0(i) = sum(actin_zfilt_highestI > 0) / numel(actin_zfilt_highestI);
    end

    % in each plane, compute the fraction of pixels well above the bg
    frac_pix_well_above_bg = num_pix_well_above_bg / sum(actin_MIP_mask(:));

    bottom_slice_num = find(actin_zfilt_highestI_prop_gr_0 <= thresh_bottom_of_cell & frac_pix_well_above_bg >= min_frac_pix_well_above_bg, 1, 'first') - 1;
    if isempty(bottom_slice_num)
        bottom_slice_num = 1;
    end    

    top_slice_num = find(actin_zfilt_highestI_prop_gr_0 >= thresh_top_of_cell & frac_pix_well_above_bg >= min_frac_pix_well_above_bg, 1, 'last') + 1;
    if isempty(top_slice_num)
        top_slice_num = find(frac_pix_well_above_bg >= min_frac_pix_well_above_bg, 1, 'last');
    end

    % the height (in μm) is the difference between the top and bottom
    % slices multiplied by the zstep
    height = zstep*(top_slice_num - bottom_slice_num);

    % actin_above_FOV will be 1 if top_slice_num = the number of
    % z-slices. In this case, the actin appears to go above the
    % top plane of imaging, and the cell should perhaps not be considered
    % for further analysis.
    if top_slice_num == nslices
        % actin goes above the top of the FOV.
        actin_above_FOV = 1;
    else
        actin_above_FOV = 0;
    end    

end