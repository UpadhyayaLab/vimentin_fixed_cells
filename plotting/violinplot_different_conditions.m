
function p_all = violinplot_different_conditions(condition_names, data_cell, field_name, colors)
nconditions = numel(condition_names);
tables = cell(nconditions, 1);
T = [];
x = [];
for i = 1:nconditions
    data = getfield(data_cell{i}, field_name);
    ncells = numel(data);
    condition = cell(ncells, 1);
    for j = 1:ncells
        condition{j} = condition_names{i};
    end
    tables{i} = table(condition, data);
    T = vertcat(T, tables{i});
    
    x_var = [i*ones(ncells,1)];
    x = [x; x_var];
end

y = T.data;

y(y == 0) = NaN;
NaN_indices = isnan(y);

figure('Position', [1 1 0.5 .85].*get(0, 'Screensize')); 
h = violinplot(y(~NaN_indices), x(~NaN_indices), 'ViolinAlpha', 0.1, 'ViolinColor', colors);
xticklabels(condition_names); 
set(gca,'linewidth',2,'fontweight','bold','fontsize',36);
axis square; 
p_all = [];

%% statistics
unique_combos = nchoosek(1:nconditions,2);
for k = 1:height(unique_combos)
    s1 = getfield(data_cell{unique_combos(k,1)}, field_name);
    s2 = getfield(data_cell{unique_combos(k,2)}, field_name);
    p = ranksum(s1,s2);
    p_all(k) = p;
    h = sigstar_adj([unique_combos(k,1) unique_combos(k,2)], p);
    unique_combos(k,3) = p;
end

end