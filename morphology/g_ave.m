function g_ave = g_ave(vim_slice, actin_mask, params)

% g_ave provides a measure of clustering.
% Here, I is the 2D slice under consideration (vimentin), and mask is a 2D
% mask used for normalization (actin).

g_vim_lower_rbound_pix = round(params.g_vim_lower_rbound/params.psize);
g_vim_upper_rbound_pix = round(params.g_vim_upper_rbound/params.psize);
redges_vim = 0:max(g_vim_upper_rbound_pix); % bounds for computing g

[~, ~, g] = get_autocorr(vim_slice, actin_mask, max(redges_vim), 0);
g_ave = mean(g(g_vim_lower_rbound_pix+1:end));
end