function [mask, area, MFI, total_sig] = segment_actin_single_plane(I, actin_MIP_mask, chosen_plane, psize, save_prog, progress_folder, cellnum)
I = I(:, :, chosen_plane); % extract plane of interest
I = I.*actin_MIP_mask; % don't consider nearby cells

gauss = imgaussfilt(I, 2) + 1; % + 1 so we don't give 0 input to log...
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
SE2 = strel('disk', 8);
BWeroded = imerode(BWdilated, SE2);
BWeroded_filled = imfill(BWeroded, 'holes');
mask = bwareafilt(BWeroded_filled, 1); % extract just one connected component

area = psize^2*bwarea(mask);
MFI = mean(I(mask));
total_sig = sum(I(mask));

%% overlay result onto raw image
if save_prog
    Perimeter = bwperim(mask);
    [rows, columns] = find(Perimeter); % store edge coordinates of cell

    figure;
    imagesc(I);
    hold on;
    plot(columns, rows, 'r.');
    set(gca,'dataAspectRatio',[1 1 1])
    axis off

    saveas(gca, [progress_folder, 'Cell_', num2str(cellnum), '_plane_', num2str(chosen_plane)], 'jpg');
end
close;

end