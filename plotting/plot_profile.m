function plot_profile(results_cell, field_name, profile_type, goodcells_saved, colors, condition_names, save_dir, savename, x_label, y_label)
% define bins
if strcmp(profile_type, 'radial')
    nbins = results_cell{1}.params.num_radial_bins;
    edges = linspace(0, 1, nbins+1);
elseif strcmp(profile_type, 'axial')
    nbins = results_cell{1}.params.num_axial_bins;
    edges = linspace(0, 1, nbins+1);
elseif strcmp(profile_type, 'chrom_compaction')
    nbins = results_cell{1}.params.chrom_intensity_8bit_num_bins;
    edges = linspace(0, 256, nbins+1);
elseif strcmp(profile_type, 'dist_from_nuc')
    nbins = results_cell{1}.params.num_radial_bins_around_nuc;
    max_dist = results_cell{1}.params.max_dist_around_nuc;
    edges = linspace(0, max_dist, nbins+1);  
elseif strcmp(profile_type, 'g_3D')
    nbins = floor(results_cell{1}.params.rmax_vim_3D/results_cell{1}.params.psize);
    edges = linspace(0, results_cell{1}.params.psize*(nbins+1), nbins+1);    
end
mid = (edges(1:end-1)+edges(2:end))/2;

% determine mean and SE for each bin and condition
nconditions = numel(results_cell);
sig_frac_cell = cell(nconditions, 1);
sig_frac_mean = zeros(nconditions, nbins);
sig_frac_SE = zeros(nconditions, nbins);

for i = 1:nconditions
    sig_frac_cell{i} = zeros(numel(goodcells_saved{i}), nbins);
    sig_frac_this_condition = getfield(results_cell{i}, field_name);
    sig_frac_cell{i} = sig_frac_this_condition(goodcells_saved{i}, :);
    sig_frac_mean(i, :) = mean(sig_frac_cell{i}, 1);
    sig_frac_SE(i, :) = std(sig_frac_cell{i}, 1)/sqrt(numel(goodcells_saved{i}));
end    

% do plotting

figure('Position', [1 1 0.5 .85].*get(0, 'Screensize')); 
for i = 1:nconditions
    errorbar(mid, sig_frac_mean(i, :), sig_frac_SE(i, :), ".-", "MarkerSize", 20, 'Color',  colors(i, :), 'LineWidth', 2)
    hold on;
end

xlabel(x_label)
ylabel(y_label)
% xlim([0 1])
set(gca,'linewidth',2,'fontweight','bold','fontsize',40);
axis square; 
legend(condition_names, 'Location', 'best')

saveas(gca, [save_dir, field_name, savename], 'fig')
saveas(gca, [save_dir, field_name, savename], 'tif')
close
end