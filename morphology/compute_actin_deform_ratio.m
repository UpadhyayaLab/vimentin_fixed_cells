function actin_deform_ratio = compute_actin_deform_ratio(MIP_major_axis, MIP_minor_axis, height)
actin_deform_ratio = mean([MIP_major_axis, MIP_minor_axis]) / height;
end