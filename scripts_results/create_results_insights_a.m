function create_results_insights_a(instance, TP_targets, ENV)
rng(11);


%% parameters
if nargin == 0
    instance = 'insights_a_separate_3';
    TP_targets = [0.60 0.65 0.70 0.75 0.80 0.85];
    ENV = get_environment();
end


%% settings
algorithms = {'greedy-dec', 'greedy-inc-dec'};
parts = {'both', 'onlybuffers', 'onlyspares'};
parts_nice = {'Integrated', 'Only buffers', 'Only spares'};


%% call
load([ENV.INPUT_DIR, instance, '.mat'], 'machines', 'number_instances');

data_TP = cell(length(parts), length(TP_targets));
data_costs = cell(length(parts), length(TP_targets));
data_C = cell(length(parts), length(TP_targets));
data_Q = cell(length(parts), length(TP_targets));
data_C_total = cell(length(parts), length(TP_targets));
data_Q_total = cell(length(parts), length(TP_targets));
for p = 1:length(parts)
    part = parts{p};
    for j = 1:length(TP_targets)
        TP_target = TP_targets(j);
        tmp_TP = NaN(number_instances, length(algorithms));
        tmp_costs = NaN(number_instances, length(algorithms));
        tmp_C = cell(1, length(algorithms));
        tmp_Q = cell(1, length(algorithms));
        for i = 1:length(algorithms)
            tmp_C{i} = NaN(number_instances, machines-1);
            tmp_Q{i} = NaN(number_instances, machines);
            algorithm = algorithms{i};
            load([ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_', algorithm, '_', part, '.mat'], 'table_values_rows');
            tmp_TP(:, i) = table_values_rows(:, 1);
            tmp_costs(:, i) = table_values_rows(:, 2);
            idx = tmp_TP < TP_target;
            tmp_TP(idx) = NaN;
            tmp_costs(idx) = NaN;
            tmp_C{i}(:, :) = table_values_rows(:, 4:(4+machines-2));
            tmp_Q{i}(:, :) = table_values_rows(:, (4+machines-2+1):((4+machines-2+1)+machines-1));
        end
        [~, pos] = min(tmp_costs, [], 2);
        opt_TP = NaN(number_instances, 1);
        opt_costs = NaN(number_instances, 1);
        opt_C = NaN(number_instances, machines-1);
        opt_Q = NaN(number_instances, machines);
        opt_C_total = NaN(number_instances, 1);
        opt_Q_total = NaN(number_instances, 1);
        for k = 1:number_instances
            opt_TP(k) = tmp_TP(k, pos(k));
            opt_costs(k) = tmp_costs(k, pos(k));
            if isnan(opt_TP(k))
                opt_C(k, :) = NaN(1, machines-1);
                opt_Q(k, :) = NaN(1, machines);
            else
                opt_C(k, :) = tmp_C{pos(k)}(k, :);
                opt_Q(k, :) = tmp_Q{pos(k)}(k, :);
            end
            opt_C_total(k) = sum(opt_C(k, :));
            opt_Q_total(k) = sum(opt_Q(k, :));
        end
        data_TP{p, j} = opt_TP;
        data_costs{p, j} = opt_costs;
        data_C{p, j} = opt_C;
        data_Q{p, j} = opt_Q;
        data_C_total{p, j} = opt_C_total;
        data_Q_total{p, j} = opt_Q_total;
    end
end

costs_diff = cell(length(parts), length(TP_targets));
number_solvable = NaN(length(parts), length(TP_targets));
costs_avg = NaN(length(parts), length(TP_targets));
C_avg = NaN(length(parts), length(TP_targets));
Q_avg = NaN(length(parts), length(TP_targets));
costs_diff_avg = NaN(length(parts), length(TP_targets));
costs_diff_std_error = NaN(length(parts), length(TP_targets));
costs_diff_min = NaN(length(parts), length(TP_targets));
costs_diff_max = NaN(length(parts), length(TP_targets));
for p = 1:length(parts)
    for j = 1:length(TP_targets)
        costs_diff{p, j} = (data_costs{p, j} - data_costs{1, j}) ./ data_costs{1, j};
        number_solvable(p, j) = number_instances - sum(isnan(data_costs{p, j}));
        costs_avg(p, j) = mean(data_costs{p, j}, 'omitnan');
        C_avg(p, j) = mean(data_C_total{p, j}, 'omitnan');
        Q_avg(p, j) = mean(data_Q_total{p, j}, 'omitnan');
        costs_diff_avg(p, j) = mean(costs_diff{p, j}, 'omitnan');
        costs_diff_std_error(p, j) = std(costs_diff{p, j}, 'omitnan') / sqrt(number_solvable(p, j));
        costs_diff_min(p, j) = min(costs_diff{p, j});
        costs_diff_max(p, j) = max(costs_diff{p, j});
    end
end

% table for overall results
table_results_columns = {'TP^T', 'case', '\# solvable', 'Mean TBC', 'Mean TCU', 'Mean TC', 'dev. (mean)', 'dev. (std. error)', 'dev. (min)', 'dev. (max)'};
table_results_rows = NaN(length(parts) * length(TP_targets), length(table_results_columns));
for j = 1:length(TP_targets)
    for p = 1:length(parts)
        table_results_rows((j - 1) * length(parts) + p, :) = ...
            [TP_targets(j) p number_solvable(p, j) C_avg(p, j) Q_avg(p, j) costs_avg(p, j) costs_diff_avg(p, j) costs_diff_std_error(p, j) costs_diff_min(p, j) costs_diff_max(p, j)];
    end
end

% format
vector_int = [3];
vector_two_decimals = [1 4:10];
vector_four_decimals = setdiff(1:length(table_results_columns), [vector_int vector_two_decimals]);
table_results_cell = cell(length(parts) * length(TP_targets), length(table_results_columns));
for ii = vector_int
    table_results_cell(:, ii) = compose('%i', table_results_rows(:, ii));
end
for ii = vector_two_decimals
    table_results_cell(:, ii) = compose('%.2f', table_results_rows(:, ii));
end
for ii = vector_four_decimals
    table_results_cell(:, ii) = compose('%.4f', table_results_rows(:, ii));
end
for ii = 1:length(table_results_cell(:, 1))
    table_results_cell{ii, 2} = parts_nice{str2double(table_results_cell{ii, 2})};
end

% shorten table
columns_to_delete = [];
table_results_cell(:, columns_to_delete) = [];
table_results_columns(:, columns_to_delete) = [];

% remove unnecessary TP^T entries
for i = 1:length(parts)
    for j = 1:length(TP_targets)
        if i > 1
            table_results_cell{(j-1) * length(parts) + i, 1} = '';
        end
    end
end

% remove unreasonable values in MPE columns
for j = 1:length(TP_targets)
    for c = 7:10
        table_results_cell{(j-1) * length(parts) + 1, c} = '';
    end
end

% remove NaN values
for i = 2:length(parts)
    for j = 1:length(TP_targets)
        for c = 4:10
            if strcmp(table_results_cell{(j-1) * length(parts) + i, c}, 'NaN')
                table_results_cell{(j-1) * length(parts) + i, c} = '';
            end
        end
    end
end

% output
table_results = cell2table(table_results_cell, 'VariableNames', table_results_columns);
writetable(table_results, [ENV.OUTPUT_DIR_EXCEL, instance, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
title_tex = '&&&&&& \multicolumn{4}{l}{$\Delta_\text{costs}$ compared to Integrated}';
table2latex_custom(table_results, [ENV.OUTPUT_DIR_LATEX, instance, '.tex'], title_tex);
disp(table_results);

% save
save([ENV.OUTPUT_DIR_DATA, instance, '.mat']);

end
