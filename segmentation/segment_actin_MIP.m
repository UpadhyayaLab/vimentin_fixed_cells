function [actin_MIP_mask, actin_MIP_area, major_axis_length, minor_axis_length, min_thresh, touching_border] = segment_actin_MIP(I, psize, subtract_bg, save_prog, progress_folder, cellnum)
% I is a background-subtracted 3D stack
% psize is the pixel size (often 0.065 μm)
% progress_folder is a folder in which the mask can be saved for later
% inspection
% cellnum is the cell number, used only to define the filename
% set save_prog = 1 if you want to save the mask and 0 otherwise
%
% actin_MIP_mask is a binary mask of the MIP
% actin_MIP_area is the area of the mask in μm^2
% major_axis_length is the length (in μm) of the major axis of the ellipse
% with the same second moments as the masked region
% minor_axis_length is the length (in μm) of the minor axis of the ellipse
% with the same second moments as the masked region
% min_thresh is the 1st percentile of all intensities within the mask - the
% median intensity outside the mask, utilized in actin_bottom_top_slices.m
% touching_border can be used to exclude cells at a later stage

% compute X-Y MIP
actin_MIP = max(I, [], 3);

if subtract_bg
    actin_MIP = actin_MIP - quantile(actin_MIP(:), .25);
    actin_MIP(actin_MIP<0) = 0;
end    

% filter and provide an initial k-means segmentation
gauss = imgaussfilt(actin_MIP, 2) + 1; % + 1 so we don't give 0 input to log...
kseg = imsegkmeans(im2single(log(gauss)), 2);

% identify brightest pixel to see if it was classified as 1 or 2
[~, indx] = max(gauss(:));
kseg_vec = kseg(:);
if kseg_vec(indx) == 1
    kseg = logical(2 - kseg);
elseif kseg_vec(indx) == 2
    kseg = logical(kseg - 1);
end

initial_mask = logical(kseg);

BWarea = bwareaopen(initial_mask, 100);
BWfill = imfill(BWarea, 'holes');

SE1 = strel('disk', 3);
BWdilated = imdilate(BWfill, SE1);
SE2 = strel('disk', 5);
BWeroded = imerode(BWdilated, SE2);
BWeroded_filled = imfill(BWeroded, 'holes');
actin_MIP_mask = bwareafilt(BWeroded_filled, 1);

actin_MIP_area = psize^2*bwarea(actin_MIP_mask);

stats = regionprops(actin_MIP_mask, "MajorAxisLength", "MinorAxisLength");
major_axis_length = psize*stats.MajorAxisLength;
minor_axis_length = psize*stats.MinorAxisLength;

% define a "min_thresh" as used in actin_bottom_top_slices.m. The actual
% definition is quite arbitrary...
min_thresh = quantile(actin_MIP(actin_MIP_mask), .01) - median(actin_MIP(~actin_MIP_mask));

% check if mask is connected to border. If so, we may wish to exclude this
% cell.
if isequal(actin_MIP_mask, imclearborder(actin_MIP_mask)) % not touching border
    touching_border = 0;
else
    touching_border = 1;
end    

%% overlay result onto raw image
if save_prog
    Perimeter = bwperim(actin_MIP_mask);
    [rows, columns]=find(Perimeter); % store edge coordinates of cell

    figure;
    imagesc(actin_MIP);
    hold on;
    plot(columns, rows, 'r.');
    set(gca,'dataAspectRatio',[1 1 1])
    axis off

    saveas(gca, [progress_folder, '\Cell_', num2str(cellnum), '_MIP'], 'jpg');
    close;
end
end