% This script was written to generate the latex table appearing in the
% paper. It is fairly hardcoded, and assumes that both non-hermitian and
% hermitian saved files exist.

file_name = 'norm_bounds.mat';
file_name_herm = 'norm_bounds_herm.mat';

results = load(file_name);
results_herm = load(file_name_herm);

% get the n-values
n_names = union(sort(fieldnames(results)),sort(fieldnames(results_herm)));
num_n = length(n_names);
n_values = zeros(num_n, 1);
for k=1:num_n
    ns = n_names{k};
    ns(1)='';
    n_values(k) = str2num(ns);
end

% get the range of m values
smallest_m = 0;
biggest_m = 0;
for k =1:num_n
    
    % populate the values of m from the two separate results files
    if isfield(results,n_names{k})
        m_names = fieldnames(results.(n_names{k}));
    else
        m_names = [];
    end
    
    if isfield(results_herm,n_names{k})
        m_names = union(m_names, fieldnames(results_herm.(n_names{k})));
    end
    
    for j = 1:length(m_names)
        m = m_names{j};
        m(1) = '';
        m = str2num(m);
        
        if (smallest_m == 0) || (m < smallest_m)
            smallest_m = m;
        end
        
        if m > biggest_m
            biggest_m = m;
        end
    end
end

% construct the table header
table = '\\begin{table}\n\\centering \n\\begin{tabular}{ c ';
for k=1:(2*num_n)
    table = [table, '| c '];
end
table = [table, '} \n '];

% top line of the table
table = [table, '\t m\\textbackslash n '];
for k = 1:num_n
    table = [table, '& ', num2str(n_values(k)), '& ', num2str(n_values(k)), '-H', '\t '];
end
table = [table, '\\\\ \\hline \n'];

% construct a row for each m
for m = smallest_m:biggest_m
    ms = ['m', num2str(m)];
    
    table = [table,' \t ', num2str(m), ' \t\t '];
    for k=1:num_n
        if isfield(results, n_names{k})
            if isfield(results.(n_names{k}), ms)
                table = [table, '& ', num2str(results.(n_names{k}).(ms).val), ' \t '];
            else
                table = [table, '& ', ' \t '];
            end
        end
        if isfield(results_herm, n_names{k})
            if isfield(results_herm.(n_names{k}), ms)
                table = [table, '& ', num2str(results_herm.(n_names{k}).(ms).val), ' \t '];
            else
                table = [table, '& ', ' \t '];
            end
        end
    end
    table = [table, '\\\\ \n']; 
end

table = [table, '\\end{tabular} \n\\end{table} \n'];
fprintf(table)