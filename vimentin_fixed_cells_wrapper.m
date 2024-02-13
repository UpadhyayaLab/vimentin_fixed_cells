clear; close all; 

% base directories for each condition
base_dirs = {
    'G:\FF\Nucleus_Data\3D Nucleus\Fixed\Blebbistatin\20240123_Fixed_E6-1_blebbistatin_Vimentin\Well2_aCD3_DMSO_640Vim_535Actin_488Centrin_405Hoechst\cells\channels\', 
    'G:\FF\Nucleus_Data\3D Nucleus\Fixed\Blebbistatin\20240123_Fixed_E6-1_blebbistatin_Vimentin\Well1_aCD3_bleb_640Vim_535Actin_488Centrin_405Hoechst\cells\channels\'
};

base_dirs = ensure_path_separator(base_dirs);

%% Set channels, including names and numbers
params.channels.names = {'vimentin', 'actin', 'centrosome'}; % the script won't work without an actin channel
params.channels.prefixes = {'C1-cell', 'C2-cell', 'C3-cell'};
params.channels.subtract_bg = [1 1 1]; % bg should be subtracted if it hasn't been already through another method
params.fname_suffix = '.tif';

params.psize = .065; % in µm
params.zstep = .3; % in µm

params.save_prog = 1; % decide whether you want to save figures. Some figures will be saved in any case.

params.resize_to_isotropic = 0; % set to 1 if you wish to characterize vimentin around the centrosome. greatly increases computational time

params.thresh_bottom_of_cell = .9; % for locating the bottom actin plane. may require adjustment; see actin_bottom_top_slices.m for more info.

params.vim_plane_nslices_above_bottom = 1; % the vimentin network this number of planes above the bottom of the cell will be considered

params.num_radial_bins = 20; % for constructing radial profiles
params.num_axial_bins = 10; % for constructing axial profiles
params.rad_fraction = .5; % for extracting the "inner-outer ratio"

params.g_vim_lower_rbound = .5; % in μm (highly arbitrary!)
params.g_vim_upper_rbound = 1; % in μm 

params.vim_COF_3D_thresh = 5; % to reduce effects of noise on the COF calculations. If this were set to 0, bg voxels might have a large influence.

params.cent_thresh = 50; % cells for which the centrosome signal is too low should perhaps be discarded later
params.cent_intensity_quantile_thresh = .999975; % used to locate the centrosome

for i = 1:numel(base_dirs)
    vimentin_fixed_cells(base_dirs{i}, params);
end
