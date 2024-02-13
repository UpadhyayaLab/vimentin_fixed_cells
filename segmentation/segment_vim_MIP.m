function [vim_MIP_mask, vim_MIP_area, vim_MIP] = segment_vim_MIP(I, psize, subtract_bg, save_prog, progress_folder, cellnum)

% compute X-Y MIP
vim_MIP = max(I, [], 3);

if subtract_bg
    vim_MIP = vim_MIP - quantile(vim_MIP(:), .75);
    vim_MIP(vim_MIP<0) = 0;
end    

gauss = imgaussfilt(vim_MIP, 2) + 1; % + 1 so we don't give 0 input to log...
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
SE2 = strel('disk', 11);
BWeroded = imerode(BWdilated, SE2);
BWeroded_filled = imfill(BWeroded, 'holes');
vim_MIP_mask = bwareafilt(BWeroded_filled, 1);

vim_MIP_area = psize^2*bwarea(vim_MIP_mask);

%% overlay result onto raw image
if save_prog
    Perimeter = bwperim(vim_MIP_mask);
    [rows, columns] = find(Perimeter); % store edge coordinates of cell

    figure;
    imagesc(vim_MIP);
    hold on;
    plot(columns, rows, 'r.');
    set(gca,'dataAspectRatio',[1 1 1])
    axis off

    saveas(gca, [progress_folder, '\Cell_', num2str(cellnum), '_MIP'], 'jpg');
    close;
end
end