function vimentin_fixed_cells(base_dir, params)

% move to the base directory
cd(base_dir)

% record which channels are present
channel_flags = create_channel_flags(params.channels.names);

% create mappings from channel names to prefixes
channel_to_prefix_map = containers.Map(params.channels.names, params.channels.prefixes);
channel_to_bg_map = containers.Map(params.channels.names, params.channels.subtract_bg);

%% Create a folder in which intermediate steps will be saved
progress_dir = fullfile(base_dir, 'prog\');
if ~exist(progress_dir, 'dir')
    mkdir(progress_dir);
end

%% Locate files
nchannels = numel(params.channels.names);
sorted_files = cell(nchannels, 1);
numeric_identifiers = cell(nchannels, 1);
% process each channel
for c = 1:nchannels
    % construct the file pattern for the current channel
    pattern = [params.channels.prefixes{c} '*' params.fname_suffix];
    
    % get sorted files for the current channel
    [sorted_files{c}, numeric_identifiers{c}] = get_sorted_files(base_dir, pattern, params.channels.prefixes{c}, params.fname_suffix);
end    

% check if channels are missing for any particular cells
check_cell_identifier_consistency(numeric_identifiers, params.channels.names)
ncells = max(numeric_identifiers{1});
results_all_cells = cell(ncells, 1);

% loop through all cells
% install the Parallel Computing Toolbox to use parfor loops - it will
% speed things up a lot
parfor i = 1:ncells
    results_this_cell = struct(); % temporary struct for this iteration

    results_this_cell.cell_found = isfile([base_dir, params.channels.prefixes{1}, num2str(i), params.fname_suffix]);
    if ~results_this_cell.cell_found % continue to next iteration
        results_all_cells{i} = results_this_cell; % record results for this cell
        continue
    end

    % initialize temporary variables
    I_actin = []; I_vim = []; I_cent = []; 
    I_actin_resized = []; I_vim_resized = []; I_cent_resized = []; 
    actin_MIP_mask = []; actin_above_bottom_mask = [];

    if ~isfield(channel_flags, 'actin')
        error('An actin channel is required for this version of the code!')
    end
    [I_actin, I_actin_resized] = read_image_subtract_bg(channel_to_prefix_map, 'actin', i, channel_to_bg_map, params.resize_to_isotropic, params.psize, params.zstep);

    if isfield(channel_flags, 'vimentin')
        [I_vim, I_vim_resized] = read_image_subtract_bg(channel_to_prefix_map, 'vimentin', i, channel_to_bg_map, params.resize_to_isotropic, params.psize, params.zstep);
    end

    if isfield(channel_flags, 'centrosome')
        [I_cent, I_cent_resized] = read_image_subtract_bg(channel_to_prefix_map, 'centrosome', i, channel_to_bg_map, params.resize_to_isotropic, params.psize, params.zstep);
    end

    %% segment actin MIP and determine bottom and top slices of cell
    actin_prog = fullfile(progress_dir, 'actin\');
    if ~exist(actin_prog, 'dir')
        mkdir(actin_prog)
    end

    % segment actin MIP
    [actin_MIP_mask, results_this_cell.actin_MIP_area, results_this_cell.actin_MIP_major_axis_length, results_this_cell.actin_MIP_minor_axis_length, actin_min_thresh, results_this_cell.actin_MIP_touching_border] = segment_actin_MIP(I_actin, params.psize, 1, 1, actin_prog, i);

    % estimate total actin signal within the cell. The raw
    % values are not meaningful, but relative values may be.
    results_this_cell.actin_total_sig  = total_signal(I_actin, actin_MIP_mask);

    % determine bottom and top slices of the actin network
    [results_this_cell.actin_bottom_slice_num, results_this_cell.actin_top_slice_num, results_this_cell.actin_height, results_this_cell.actin_above_FOV] = actin_bottom_top_slices(I_actin, actin_MIP_mask, params.zstep, actin_min_thresh, params.thresh_bottom_of_cell);

    % segment actin at the bottom slice
    [actin_bottom_slice_mask, results_this_cell.actin_bottom_slice_mask_area, results_this_cell.actin_bottom_slice_MFI, results_this_cell.actin_bottom_slice_total_sig] = segment_actin_single_plane(I_actin, actin_MIP_mask, results_this_cell.actin_bottom_slice_num, params.psize, params.save_prog, actin_prog, i);

    % get actin "deformation ratio"
    results_this_cell.actin_deform_ratio = compute_actin_deform_ratio(results_this_cell.actin_MIP_major_axis_length, results_this_cell.actin_MIP_minor_axis_length, results_this_cell.actin_height);

    % segment slice at a single, user-selected slice relative to the bottom of the cell. The vimentin network
    % will be considered at this plane.
    [actin_above_bottom_mask, results_this_cell.actin_above_bottom_mask_area, results_this_cell.actin_above_bottom_MFI, results_this_cell.actin_above_bottom_total_sig] = segment_actin_single_plane(I_actin, actin_MIP_mask, results_this_cell.actin_bottom_slice_num+params.vim_plane_nslices_above_bottom, params.psize, params.save_prog, actin_prog, i);

    % describe actin's radial distribution at the "bottom" slice
    [results_this_cell.actin_bottom_slice_r_eff, results_this_cell.actin_bottom_slice_rad_profile] = radial_profile(I_actin(:, :, results_this_cell.actin_bottom_slice_num), actin_bottom_slice_mask, params.num_radial_bins);
    results_this_cell.actin_bottom_slice_inner_outer_ratio = rad_dist(I_actin(:, :, results_this_cell.actin_bottom_slice_num), actin_bottom_slice_mask, params.rad_fraction);

    % describe actin's axial distribution
    results_this_cell.actin_axial_profile = axial_signal_binning(I_actin, actin_MIP_mask, results_this_cell.actin_bottom_slice_num, results_this_cell.actin_top_slice_num, params.num_axial_bins);
   
    if isfield(channel_flags, 'vimentin')
        vim_prog = fullfile(progress_dir, 'vim\');
        if ~exist(vim_prog, 'dir')
            mkdir(vim_prog)
        end

        % segment MIP
        [~, results_this_cell.vim_MIP_area] = segment_vim_MIP(I_vim, params.psize, 1, 1, vim_prog, i);

        % estimate total vimentin signal within the cell. The raw
        % values are not meaningful, but relative values may be.
        results_this_cell.vim_total_sig  = total_signal(I_vim, actin_MIP_mask);

        % characterize vimentin clustering at the plane of interest
        vim_slice = I_vim(:, :, results_this_cell.actin_bottom_slice_num + params.vim_plane_nslices_above_bottom);

        % g_ave
        results_this_cell.g_ave_vim_above_bottom = g_ave(vim_slice, actin_above_bottom_mask, params);

        % FDD
        [~, results_this_cell.vim_FDD_above_bottom, results_this_cell.vim_FDD_above_bottom_actin_norm] = COF_2D(vim_slice, actin_above_bottom_mask, params.psize);

        % vimentin radial distribution at this slice
        [results_this_cell.vim_above_bottom_r_eff, results_this_cell.vim_above_bottom_rad_profile] = radial_profile(vim_slice, actin_above_bottom_mask, params.num_radial_bins);
        results_this_cell.vim_above_bottom_inner_outer_ratio = rad_dist(vim_slice, actin_above_bottom_mask, params.rad_fraction);

        %% characterize vimentin's axial distribution
        [~, results_this_cell.vim_zCOF, results_this_cell.vim_FDD_3D, ~, results_this_cell.vim_z_FDD] = COF_3D(I_vim, actin_MIP_mask, results_this_cell.actin_bottom_slice_num, results_this_cell.actin_top_slice_num, params.psize, params.zstep, 'actin', params.vim_COF_3D_thresh);
        results_this_cell.vim_zCOF_actin_scale = (results_this_cell.vim_zCOF - params.zstep*results_this_cell.actin_bottom_slice_num) / results_this_cell.actin_height;
        results_this_cell.vim_zCOF_cell_bottom_distance = results_this_cell.vim_zCOF - params.zstep*results_this_cell.actin_bottom_slice_num;

        % describe actin's axial distribution
        results_this_cell.vim_axial_profile = axial_signal_binning(I_vim, actin_MIP_mask, results_this_cell.actin_bottom_slice_num, results_this_cell.actin_top_slice_num, params.num_axial_bins);
    end

    if isfield(channel_flags, 'centrosome')
        cent_prog = fullfile(progress_dir, 'cent\');
        if ~exist(cent_prog, 'dir')
            mkdir(cent_prog)
        end

        % cells for which the centrosome signal is too low should perhaps be
        % discarded later
        results_this_cell.centrosome_found = max(I_cent.*actin_MIP_mask, [], "all") > params.cent_thresh;

        % localize the centrosome
        centrosome_center = localize_cent_3D(I_cent, actin_MIP_mask, params.cent_intensity_quantile_thresh, params.psize, params.zstep, params.save_prog, cent_prog, i);
        results_this_cell.centrosome_center_z = centrosome_center(3);
        results_this_cell.centrosome_center_z_rel_bottom_actin_plane = results_this_cell.centrosome_center_z - params.zstep*results_this_cell.actin_bottom_slice_num;
        if results_this_cell.centrosome_center_z_rel_bottom_actin_plane < 2 % quite arbitrary
            results_this_cell.cent_polarized = 1;
        else
            results_this_cell.cent_polarized = 0;
        end

        % centrosome's radial position within the cell
        results_this_cell.centrosome_r_norm = centrosome_radial_pos(actin_MIP_mask, centrosome_center, params.psize);

        % resize before localizing the centrosome
        if params.resize_to_isotropic
            centrosome_center_resized = localize_cent_3D(I_cent_resized, actin_MIP_mask, params.cent_intensity_quantile_thresh, params.psize, params.psize, 0);
            results_this_cell.centrosome_center_z_resized = centrosome_center_resized(3);

            % characterize actin around the centrosome
            [results_this_cell.actin_MFI_around_cent_1um, results_this_cell.actin_frac_around_cent_1um, results_this_cell.actin_MFI_around_cent_2um, results_this_cell.actin_frac_around_cent_2um, results_this_cell.actin_MFI_around_cent_3um, results_this_cell.actin_frac_around_cent_3um] = cent_assoc_3D(I_actin_resized, centrosome_center_resized, actin_MIP_mask, params.psize, params.psize, 5);

            if isfield(channel_flags, 'vimentin')
                % characterize vimentin around the centrosome
                [results_this_cell.vim_MFI_around_cent_1um, results_this_cell.vim_frac_around_cent_1um, results_this_cell.vim_MFI_around_cent_2um, results_this_cell.vim_frac_around_cent_2um, results_this_cell.vim_MFI_around_cent_3um, results_this_cell.vim_frac_around_cent_3um] = cent_assoc_3D(I_vim_resized, centrosome_center_resized, actin_MIP_mask, params.psize, params.psize, 5);
            end
        end    
   
    end    

    % record results for this cell
    results_all_cells{i} = results_this_cell;

    close all
    fprintf('Finished cell %i\n', i)
end

results = consolidate_results(results_all_cells, params);

save('vimentin_fixed_cell_results.mat', '-struct', 'results')
disp(['Finished folder ', base_dir])

end

%% identify channels present
function channel_flags = create_channel_flags(params_channels)
    % Initialize an empty structure for channel flags
    channel_flags = struct();
    
    % Iterate over the list of channels in params.channels
    for idx = 1:length(params_channels)
        channel = params_channels{idx};
        % Set a flag for the current channel
        channel_flags.(channel) = true;
    end
end

function [sorted_files, numeric_identifiers] = get_sorted_files(base_dir, pattern, prefix, suffix)
    % List all files matching the pattern
    files = dir(fullfile(base_dir, pattern));
    
    % Initialize an array to store the numeric identifiers
    numeric_identifiers = zeros(1, length(files));
    
    % Prepare regular expression for extracting numeric identifier
    % This assumes prefix and suffix do not contain regex special characters.
    % If they do, they need to be escaped accordingly.
    reg_ex_pattern = [prefix, '(\d+)', suffix];
    
    % Initialize counter for files that match the exact pattern
    valid_file_count = 0;
    
    % Extract the numeric identifier from each filename
    for i = 1:length(files)
        filename = files(i).name;
        
        % Use regular expression to match and extract numeric part
        matches = regexp(filename, reg_ex_pattern, 'tokens');
        
        if ~isempty(matches)
            % Increment the count of valid files
            valid_file_count = valid_file_count + 1;
            
            % Convert the matched numeric part to a number and store
            numeric_identifiers(valid_file_count) = str2double(matches{1}{1});
        end
    end

    % Check if no valid files were found
    if valid_file_count == 0
        error('No files matching the specified pattern were found.');
    end
    
    % Adjust arrays to include only the matched files and their identifiers
    numeric_identifiers = numeric_identifiers(1:valid_file_count);
    sorted_files = files(1:valid_file_count);
    
    % Sort files based on the numeric identifiers
    [numeric_identifiers, sort_idx] = sort(numeric_identifiers);
    sorted_files = sorted_files(sort_idx);
end

%% check if channels are missing for any particular cells
function check_cell_identifier_consistency(numeric_identifiers, channels)
    % Find the union of all numeric identifiers across channels
    all_identifiers = unique(horzcat(numeric_identifiers{:}));
    
    % Initialize a structure to hold missing identifiers for each channel
    missing_identifiers = struct();
    
    % Compare each channel's identifiers against the union of all identifiers
    for c = 1:length(channels)
        current_identifiers = numeric_identifiers{c};
        missing = setdiff(all_identifiers, current_identifiers); % Identifiers missing in the current channel
        
        if ~isempty(missing)
            missing_identifiers.(channels{c}) = missing;
        else
            missing_identifiers.(channels{c}) = []; % No missing identifiers
        end
    end
    
    % Output results
    for c = 1:length(channels)
        if ~isempty(missing_identifiers.(channels{c}))
            error('Channel %s is missing cells: %s\n', channels{c}, mat2str(missing_identifiers.(channels{c})));
        end    
    end
end

function [I, I_resized] = read_image_subtract_bg(channel_to_prefix_map, channel_name, cell_num, channel_to_bg_map, resize_to_isotropic, psize, zstep)
    file_prefix = channel_to_prefix_map(channel_name);
    subtract_BG = channel_to_bg_map(channel_name);
    fname = [file_prefix, num2str(cell_num), '.tif'];
    I = double(ReadTifStack(fname));  

    % resize to be isotropic in x, y, and z
    if resize_to_isotropic
        I_resized = resize_3D(I, psize, zstep, 'cubic');
        if subtract_BG
            I_resized = subtract_bg(I_resized);
        end
    else
        I_resized = [];
    end

    if subtract_BG
        I = subtract_bg(I);
    end  
end

function results = consolidate_results(results_all_cells, params)
    ncells = numel(results_all_cells);
    results = struct();
    
    if ncells > 0
        fields = fieldnames(results_all_cells{1});
        
        for f = 1:numel(fields)
            field_name = fields{f};
            sample_data = results_all_cells{1}.(field_name);
            
            % Handling scalar values
            if isscalar(sample_data)
                results.(field_name) = nan(ncells, 1);
            % Handling vector values
            elseif isvector(sample_data)
                % Determine vector length (m)
                if isrow(sample_data)
                    m = length(sample_data);
                    results.(field_name) = nan(ncells, m);
                else % column vector
                    m = length(sample_data);
                    results.(field_name) = nan(ncells, m);
                end
            % Handling cell arrays
            elseif iscell(sample_data)
                results.(field_name) = cell(ncells, 1);
            end
            
            % Consolidate data from all cells
            for i = 1:ncells
                if ~results_all_cells{i}.cell_found % this cell not in folder
                    results.cell_found(i) = 0;
                    continue
                else
                    % For scalar and vector values
                    if ~iscell(sample_data)
                        % Direct assignment for scalars and vectors
                        if isvector(sample_data) && isrow(sample_data)
                            results.(field_name)(i, :) = results_all_cells{i}.(field_name);
                        else % column vector or scalar
                            results.(field_name)(i, :) = reshape(results_all_cells{i}.(field_name), 1, []);
                        end
                        % For cell arrays
                    else
                        results.(field_name){i} = results_all_cells{i}.(field_name);
                    end
                end
            end
        end
    end

    % include params in results
    results.params = params;

    % create a datetime object representing the current date and time
    current_date_time = datetime('now');
    results.datetime = current_date_time;
end