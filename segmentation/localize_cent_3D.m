function centrosome_center = localize_cent_3D(I_cent, actin_MIP_mask, cent_intensity_quantile_thresh, psize, zstep, save_prog, progress_folder, cellnum)

    %% locate centrosome
    nslices = size(I_cent, 3);
    actin_MIP_mask_3D = repmat(actin_MIP_mask, [1, 1, nslices]);
    I_cent = I_cent.*actin_MIP_mask_3D; % don't consider centrosomes from nearby cells

    thresh = quantile(I_cent(:), cent_intensity_quantile_thresh); % using really high threshold because this is 3D
    BG = find(I_cent<=thresh);
    mask = I_cent;
    mask(mask<=thresh) = 0; mask(mask>thresh) = 1; % initial mask
    mask = bwareaopen(mask, 20);

    stats = regionprops3(mask, 'Volume', 'Centroid');
    % if no regions are found, lower area thresh to make sure that the code
    % doesn't crash...
    if numel(stats) == 0
        mask = I_cent;
        mask(mask<=thresh) = 0; mask(mask>thresh) = 1; % initial mask
        % mask = bwareaopen(mask, 5);
        stats = regionprops3(mask, 'Volume', 'Centroid');
    end

    volumes = stats.Volume;
    centroids_cent_regions = stats.Centroid;

    % determine centrosome center, weighting these regions by volumes
    weights = volumes/sum(volumes);
    centrosome_center = sum(centroids_cent_regions.*weights, 1);

    % be really careful about difference btw x and y!
    temp = centrosome_center(1);
    centrosome_center(1) = centrosome_center(2);
    centrosome_center(2) = temp;

    ind3 = centrosome_center(3);

    if save_prog
        % plot centrosome segmentation at slice where its center supposedly is
        figure; imagesc(I_cent(:, :, round(ind3)))
        colormap gray
        set(gca,'dataAspectRatio',[1 1 1])
        axis off
        title(['z = ', num2str(zstep*ind3), ' µm'])
        saveas(gca, [progress_folder, 'Cell_', num2str(cellnum), '_raw'], 'jpg');
        close;

        figure; imagesc(mask(:, :, round(ind3)))
        colormap gray
        set(gca,'dataAspectRatio',[1 1 1])
        axis off
        title(['z = ', num2str(zstep*ind3), ' µm'])
        saveas(gca, [progress_folder, 'Cell_', num2str(cellnum), '_mask'], 'jpg');
        close;
    end

    % convert to microns
    centrosome_center(:, 1:2) = psize * centrosome_center(:, 1:2);
    centrosome_center(:, 3) = zstep * centrosome_center(:, 3);
end