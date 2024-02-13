clear; close all

%% set output folder and folder in which the results reside
save_dir = 'G:\FF\Nucleus_Data\3D Nucleus\Fixed\Blebbistatin\20240123_Fixed_E6-1_blebbistatin_Vimentin\vim_figures\';
base_dirs = {
    'G:\FF\Nucleus_Data\3D Nucleus\Fixed\Blebbistatin\20240123_Fixed_E6-1_blebbistatin_Vimentin\Well2_aCD3_DMSO_640Vim_535Actin_488Centrin_405Hoechst\cells\channels\', 
    'G:\FF\Nucleus_Data\3D Nucleus\Fixed\Blebbistatin\20240123_Fixed_E6-1_blebbistatin_Vimentin\Well1_aCD3_bleb_640Vim_535Actin_488Centrin_405Hoechst\cells\channels\'
};

base_dirs = ensure_path_separator(base_dirs);

% create the directory for saving if it doesn't exist already
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
cd(save_dir)

savename = []; % this can be used to mark which conditions or fields you're using 
% if you'll have multiple in the same folder

nconditions = numel(base_dirs);
results_cell = cell(nconditions, 1);

condition_names = ["Ctrl", "Bleb"];

colors(1, :) = [0 0 1];
colors(2, :) = [0.9290 0.6940 0.1250];

% load in results
for i = 1:nconditions
    results_cell{i} = load([base_dirs{i}, 'vimentin_fixed_cell_results.mat']);
end

%% determine fields to extract
field_names_all = fieldnames(results_cell{1}); % should be same for all datasets.
% if data have been analyzed with different versions of the code, this may lead
% to issues. 
extract_field = ones(numel(field_names_all), 1);

for i = 1:nconditions
    for j = 1:numel(field_names_all)
        field = getfield(results_cell{i}, field_names_all{j});
        if ~isnumeric(field) || ~isvector(field) % don't extract if field isn't a vector
            extract_field(j) = 0;
        elseif ~any(field ~= 0) % don't extract if field is all 0s
            extract_field(j) = 0;
        end    
    end    
end
field_names = field_names_all(extract_field == 1);

%% determine cells to eliminate. This procedure will have to be varied depending on
% the application. If the cell number did not exist in the folder, it
% should definitely not be considered.

goodcells_saved = cell(nconditions, 1);
n_goodcells = zeros(nconditions, 1);
fraction_goodcells = zeros(nconditions, 1);
for i = 1:nconditions
    ncells = numel(results_cell{i}.(field_names{1}));
    goodcells = [];
    for j = 1:ncells
        % the specific criteria of inclusion for cells can be varied!
        if results_cell{i}.cell_found(j) && ~(results_cell{i}.actin_above_FOV(j)) && ~(results_cell{i}.actin_MIP_touching_border(j))
        % if results_cell{i}.cell_found(j) && ~(results_cell{i}.actin_above_FOV(j)) && ~(results_cell{i}.actin_MIP_touching_border(j)) && (results_cell{i}.centrosome_found(j))
            goodcells = [goodcells j];
        else
            continue
        end
    end 
    goodcells_saved{i} = goodcells;
    n_goodcells(i) = numel(goodcells);
    fraction_goodcells(i) = numel(goodcells)/ncells;
end

%% distill data for each condition with the fields and cells of interest
data = cell(nconditions, 1);
for i = 1:nconditions
    for j = 1:numel(field_names)
        field = getfield(results_cell{i}, field_names{j});
        try
            data{i}.(field_names{j}) = field(goodcells_saved{i});
        end    
    end
end   

%% make a violin plot for each field
p_all = cell(numel(field_names), 1);
parfor i = 1:numel(field_names)
    % violin_plot will give errors if the data don't have enough distinct
    % values. try/catch is just a way of avoiding the issue...
    try
        p_all{i} = violinplot_different_conditions(condition_names, data, field_names{i}, colors);
        % ylabel(field_names{i})
        ylabel(field_names{i}, 'interpreter', 'none')

        saveas(gca, [save_dir, field_names{i}, savename], 'fig')
        saveas(gca, [save_dir, field_names{i}, savename], 'tif')
        close
    catch
        close
    end
end

%% radial and axial profiles of actin and vimentin
plot_profile(results_cell, 'actin_bottom_slice_rad_profile', 'radial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Radial Position', 'Normalized Signal')
plot_profile_all_cells(results_cell, 'actin_bottom_slice_rad_profile', 'radial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Radial Position', 'Normalized Signal')

plot_profile(results_cell, 'actin_axial_profile', 'axial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Axial Position', 'Normalized Signal')
plot_profile_all_cells(results_cell, 'actin_axial_profile', 'axial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Axial Position', 'Normalized Signal')

if ismember('vimentin', results_cell{1}.params.channels.names)
    plot_profile(results_cell, 'vim_above_bottom_rad_profile', 'radial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Radial Position', 'Normalized Signal')
    plot_profile_all_cells(results_cell, 'vim_above_bottom_rad_profile', 'radial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Radial Position', 'Normalized Signal')

    plot_profile(results_cell, 'vim_axial_profile', 'axial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Axial Position', 'Normalized Signal')
    plot_profile_all_cells(results_cell, 'vim_axial_profile', 'axial', goodcells_saved, colors, condition_names, save_dir, savename, 'Normalized Axial Position', 'Normalized Signal')
end        