function [total_sig, percent_below_nuc, percent_above_nuc] = total_signal(I, mask, nuc_bottom_slice_num, nuc_top_slice_num)
    % I should be a bg subtracted 3D image
    mask_rep = repmat(mask, 1, 1, size(I, 3));
    I = I.*mask_rep;
    total_sig = sum(I(:));

    % compute percent of signal below and above the nucleus
    if nargin == 4
        percent_below_nuc = 100*sum(I(:,:,1:nuc_bottom_slice_num-1), "all")/total_sig;
        percent_above_nuc = 100*sum(I(:,:,nuc_top_slice_num+1:end), "all")/total_sig;
    end    
end