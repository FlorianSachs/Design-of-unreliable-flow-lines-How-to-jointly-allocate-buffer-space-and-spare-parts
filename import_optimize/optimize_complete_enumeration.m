function [c, s, TP, costs, counter_iterations, counter_evaluations, counter_problems] = optimize_complete_enumeration(machines, C_min, Q_min, C_max, Q_max, mu, p, gamma, TP_target, cost_buffers, cost_spares, outputdir, filename)
%optimize_complete_enumeration uses complete enumeration to reach target
%   All cominations of buffers or spare part inventories are increased to
%   reach the target

if nargin == 0
    fprintf('No parameters given... using standard values...\n\n');
    addpath('import');
    machines = 3;
    TP_target = 0.85;
    C_min = 2 * ones(1, machines-1);
    Q_min = 1 * ones(1, machines);
    C_max = 10 * ones(1, machines-1);
    Q_max = 3 * ones(1, machines);
    mu = 1 * ones(1, machines);
    p = 0.005 * ones(1, machines);
    gamma = 0.025 * ones(1, machines);
    cost_buffers = 1 * ones(1, machines-1);
    cost_spares = 1 * ones(1, machines);
    outputdir = '';
    filename = '';
end


%% init

% create instances
l_word = machines - 1;
n_letters = C_max - C_min + 1;
C = cell(1, l_word);
[C{:}] = ndgrid(0:n_letters - 1);
instances_C = reshape(cat(l_word+1, C{:}), [], l_word) + C_min;

l_word = machines;
n_letters = Q_max - Q_min + 1;
C = cell(1, l_word);
[C{:}] = ndgrid(0:n_letters - 1);
instances_Q = reshape(cat(l_word+1, C{:}), [], l_word) + Q_min;

number_instances_C = size(instances_C, 1);
number_instances_Q = size(instances_Q, 1);
number_instances = number_instances_C * number_instances_Q;

% create output table columns
table_columns = [{'TP' 'costs' 'time'} ...
                        cellfun(@(c)['C_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false) ...
                        cellfun(@(c)['Q_', num2str(c)],num2cell(1:machines),'UniformOutput',false) ...
                        {'costs_adj'}];


%% call
table_values_rows = NaN(number_instances, 3 + 2*machines - 1);
counter_iterations = 0;
counter_evaluations = number_instances;
counter_problems = 0;
%parfor k = 1:number_instances
for k = 1:number_instances
    %disp(char(strcat({' '}, '---', {' '}, num2str(k), {'/'}, num2str(number_instances), {' '}, '---')));
    
    % get general instance data
    i = mod(k-1, number_instances_C) + 1;
    j = ceil(k / number_instances_C);
    
    C = instances_C(i, :);
    Q = instances_Q(j, :);
    
    % calculate
    tic();
    if (machines >= 3)
        [TP, ~, ~, dec_terminated_normally] = Spare_Decomposition(machines, C, Q, mu, p, gamma, true);
        counter_problems = counter_problems + ~dec_terminated_normally ;
    elseif (machines == 3)
        TP = Spare_3M(C(1), C(2), Q(1), Q(2), Q(3), mu(1), mu(2), mu(3), p(1), p(2), p(3), gamma(1), gamma(2), gamma(3));
    else
        TP = Spare_2M(C, Q(1), Q(2), mu(1), mu(2), p(1), p(2), gamma(1), gamma(2));
    end
    time = toc();
    
    costs = sum([C .* cost_buffers Q .* cost_spares]);
    table_values_rows(k, :) = [TP costs time C Q];
end


%% rearrange data
%lookup_table_TP = NaN([C_max - C_min + 1 Q_max - Q_min + 1]);
%lookup_table_costs = NaN([C_max - C_min + 1 Q_max - Q_min + 1]);
%for k = 1:number_instances
%    C = table_values_rows(k, 4:(4 + machines - 2));
%    Q = table_values_rows(k, (4 + machines - 1):(4 + 2 * machines - 2));
%    lookup_table_TP([C + C_min, Q + Q_min]) = table_values_rows(k, 1);
%    lookup_table_costs([C + C_min, Q + Q_min]) = table_values_rows(k, 2);
%end


%% find minimal costs
table_values_rows(:, end+1) = table_values_rows(:, 2);
table_values_rows(table_values_rows(:, 1) < TP_target, end) = NaN;

[global_min, global_min_pos] = min(table_values_rows(:, end));
global_min_indexes = table_values_rows(:, end) == global_min;
%global_min_C_Q = table_values_rows(global_min_indexes, [4:(4 + machines - 2), (4 + machines - 1):(4 + 2 * machines - 2)]);

if (isnan(global_min))
    % no feasible solution exists
    c = NaN(1, machines-1);
    s = NaN(1, machines);
    TP = NaN;
    costs = NaN;
    table_values_rows_opt = NaN(1, 3 + 2*machines - 1);
else
    % select best TP with minimum cost
    if (length(global_min_indexes) > 1)
        TP_adjusted = table_values_rows(:, 1);
        TP_adjusted(~global_min_indexes) = NaN;
        [~, global_min_pos] = max(TP_adjusted);
    end
    
    % save optimal results
    c = table_values_rows(global_min_pos, 4:(4 + machines - 2));
    s = table_values_rows(global_min_pos, (4 + machines - 1):(4 + 2 * machines - 2));
    TP = table_values_rows(global_min_pos, 1);
    costs = global_min;
    table_values_rows_opt = table_values_rows(global_min_pos, :);
end


%% output
table_values = array2table(table_values_rows, 'VariableNames', table_columns);
if ~isempty(filename)
    writetable(table_values, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
else
    disp(table_values);
end

table_opt = array2table(table_values_rows_opt, 'VariableNames', table_columns);
if ~isempty(filename)
    writetable(table_opt, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Best', 'Range', 'A1', 'WriteMode','overwritesheet');
else
    disp(table_opt);
end


end