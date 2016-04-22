% This is a script written to create the plots appearing in the paper.

% what to plot
n = 6;
m_bottom = n;
m_top = n^2;

file_name = 'norm_bounds.mat';
results = load(file_name);

m_vals = zeros(1, m_top - m_bottom + 1);
for k = 1:length(m_vals)
    m_vals(k) = n + k -1;
end

ns = ['n',num2str(n)];

computed_values = zeros(1, length(m_vals));
for k=1:length(m_vals)
    ms = ['m',num2str(m_vals(k))];
    
    computed_values(k) = results.(ns).(ms).val;
end

scatter(m_vals, computed_values)