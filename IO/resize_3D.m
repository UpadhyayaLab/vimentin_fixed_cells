function I_resized = resize_3D(I, psize, zstep, interp_method) 
    I_resized = imresize3(I, 'Scale', [1 1 zstep/psize], 'Method', interp_method);
end