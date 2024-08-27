function create_results_validation(instance, TP_targets, algorithms, algorithms_nice, ENV)
rng(10);


%% parameters
if nargin == 0
    instance = 'validation_unbalanced_2';
    TP_targets = [0.75 0.80 0.85 0.90];
    algorithms = {'enumeration', 'greedy-dec', 'greedy-inc', 'greedy-inc-dec', 'simulated-annealing-dec', 'simulated-annealing-inc', 'genetic-algorithm'};
    algorithms_nice = {'ENUM', 'GD', 'GI', 'GID', 'SAD', 'SAI', 'GA'};
    ENV = get_environment();
end


%% call
data = cell(length(algorithms), length(TP_targets));
data_cleaned = cell(length(algorithms), length(TP_targets));
number_unsolvable = cell(length(algorithms), length(TP_targets));
number_solvable = cell(length(algorithms), length(TP_targets));
number_infeasible = cell(length(algorithms), length(TP_targets));
number_best = cell(length(algorithms), length(TP_targets));
TP_diff = cell(length(algorithms), length(TP_targets));
costs_diff = cell(length(algorithms), length(TP_targets));
time_diff = cell(length(algorithms), length(TP_targets));
costs = cell(1, length(TP_targets));
load([ENV.INPUT_DIR, instance, '.mat'], 'number_instances');
for j = 1:length(TP_targets)
    costs{j} = NaN(number_instances, length(algorithms));
end
for i = 1:length(algorithms)
    algorithm = algorithms{i};
    for j = 1:length(TP_targets)
        TP_target = TP_targets(j);
        load([ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_', algorithm '.mat'], 'table_values_rows');
        if any(table_values_rows(:, 1) < TP_target & table_values_rows(:, 1) ~= -1)
            idx = find(table_values_rows(:, 1) < TP_target & table_values_rows(:, 1) ~= -1);
            disp(['There are ', num2str(length(idx)), ' infeasible results for ', algorithms_nice{i}, ' (TP^T = ', num2str(TP_target, '%.2f'), ')!']);
            table_values_rows(idx, :) = -2;
        end
        idx_unsolvable = find(table_values_rows(:, 1) == -1);
        idx_infeasible = find(table_values_rows(:, 1) == -2);
        data{i, j} = table_values_rows;
        table_values_rows(idx_unsolvable, :) = NaN;
        table_values_rows(idx_infeasible, :) = NaN;
        data_cleaned{i, j} = table_values_rows;
        costs{j}(:, i) = table_values_rows(:, 2);
        tmp_data_selected_algorithm = data{1, j};
        tmp_data_selected_algorithm(idx_unsolvable, :) = NaN;
        tmp_data_selected_algorithm(idx_infeasible, :) = NaN;
        number_unsolvable{i, j} = length(idx_unsolvable);
        number_infeasible{i, j} = length(idx_infeasible);
        number_solvable{i,j} = length(data_cleaned{i, j}(:, 1)) - number_unsolvable{i,j} - number_infeasible{i,j};
        TP_diff{i, j} = (data_cleaned{i, j}(:, 1) - tmp_data_selected_algorithm(:, 1)) ./ tmp_data_selected_algorithm(:, 1);
        costs_diff{i, j} = (data_cleaned{i, j}(:, 2) - tmp_data_selected_algorithm(:, 2)) ./ tmp_data_selected_algorithm(:, 2);
        time_diff{i, j} = (data_cleaned{i, j}(:, 3) - tmp_data_selected_algorithm(:, 3)) ./ tmp_data_selected_algorithm(:, 3);
    end
end

% get best results regarding cost
costs_min = NaN(number_instances, length(TP_targets));
costs_min_realized = cell(1, length(TP_targets));
for j = 1:length(TP_targets)
    costs_min(:, j) = min(costs{j}, [], 2);
    costs_min_realized{j} = NaN(number_instances, length(algorithms));
    for i = 1:length(algorithms)
        number_best{i,j} = sum(costs_min(:, j) == costs{j}(:, i));
        costs_min_realized{j}(:, i) = costs_min(:, j) == costs{j}(:, i);
    end
end

% get information on number of best and better solutions
best_solution_devision = cell(1, length(TP_targets));
better_solution_devision = cell(1, length(TP_targets));
for j = 1:length(TP_targets)
    best_solution_devision{j} = NaN(length(algorithms), length(algorithms) + 1);
    better_solution_devision{j} = NaN(length(algorithms), length(algorithms) + 1);
    for i = 1:length(algorithms)
        idx = costs_min_realized{j}(:, i) > 0;
        for ii = 1:length(algorithms)
            best_solution_devision{j}(i, ii + 1) = sum(costs_min_realized{j}(idx, ii));
            infeasible_current = sum(isnan(costs{j}(:, i)));
            infeasible_all_A = isnan(costs{j}(:, i));
            infeasible_all_B = isnan(costs{j}(:, ii));
            tmp_A = infeasible_all_A == infeasible_all_B;
            tmp_B = infeasible_all_A == true;
            tmp_C = tmp_A == tmp_B;
            infeasible_current_and_compared = sum(tmp_C);
            infeasible = infeasible_current - infeasible_current_and_compared;
            if infeasible < 0, infeasible = 0; end
            better_solution_devision{j}(i, ii + 1) = sum(costs{j}(:, i) > costs{j}(:, ii)) + infeasible;
        end
        best_solution_devision{j}(i, 1) = best_solution_devision{j}(i, i + 1);
        best_solution_devision{j}(i, 2:end) = best_solution_devision{j}(i, 2:end) ./ best_solution_devision{j}(i, 1);
        better_solution_devision{j}(i, 1) = number_solvable{i, j};
        better_solution_devision{j}(i, 2:end) = better_solution_devision{j}(i, 2:end) ./ better_solution_devision{j}(i, 1);
    end
end

%[costs_min(:, j), costs{j}, costs_min(:, j) == costs{j}(:, 1)]
%idx = find(costs_min(:, j) ~= costs{j}(:, 1) & ~isnan(costs_min(:, j)))
%[costs_min(idx, j), costs{j}(idx, :), costs_min(idx, j) == costs{j}(idx, 1)]

% table for overall results
table_results_columns = {'TP^T', 'Algorithm', '\# solvable instances', '\# infeasible solutions', '\# best solutions', [algorithms_nice{1}, ' worse'], [algorithms_nice{1}, ' better'], 'it. (avg)', 'it. (std. error)', 'it. (min)',  'it. (max)', 'ev. (avg)', 'ev. (std. error)', 'ev. (min)', 'ev. (max)', 'time (avg)', 'time (std. error)', 'time (min)', 'time (max)', 'dev. (mean)', 'dev. (std. error)', 'dev. (min)', 'dev. (max)'};
table_results_rows = NaN(length(algorithms) * length(TP_targets), length(table_results_columns));
for i = 1:length(algorithms)
    for j = 1:length(TP_targets)
        number_worse = sum(costs_diff{i, j} < 0);
        number_better = sum(costs_diff{i, j} > 0);
        it_avg = mean(data_cleaned{i,j}(:, end-2), 'omitnan');
        it_std_error = std(data_cleaned{i,j}(:, end-2), 'omitnan') / sqrt(number_solvable{i,j});
        it_min = min(data_cleaned{i,j}(:, end-2));
        it_max = max(data_cleaned{i,j}(:, end-2));
        ev_avg = mean(data_cleaned{i,j}(:, end-1), 'omitnan');
        ev_std_error = std(data_cleaned{i,j}(:, end-1), 'omitnan') / sqrt(number_solvable{i,j});
        ev_min = min(data_cleaned{i,j}(:, end-1));
        ev_max = max(data_cleaned{i,j}(:, end-1));
        time_avg = mean(data_cleaned{i,j}(:, 3), 'omitnan');
        time_std_error = std(data_cleaned{i,j}(:, 3), 'omitnan') / sqrt(number_solvable{i,j});
        time_min = min(data_cleaned{i,j}(:, 3));
        time_max = max(data_cleaned{i,j}(:, 3));
        dev_avg = mean(costs_diff{i, j}, 'omitnan');
        dev_std_error = std(costs_diff{i, j}, 'omitnan') / sqrt(number_solvable{i,j});
        dev_min = min(costs_diff{i, j});
        dev_max = max(costs_diff{i, j});
        %table_results_rows((i - 1) * length(TP_targets) + j, :) = ...
        table_results_rows((j - 1) * length(algorithms) + i, :) = ...
            [TP_targets(j) i number_solvable{i,j} number_infeasible{i,j} number_best{i,j} number_worse number_better ...
            it_avg it_std_error it_min it_max ev_avg ev_std_error ev_min ev_max time_avg time_std_error time_min time_max dev_avg dev_std_error dev_min dev_max];
    end
end

% format
vector_int = [2 3 4 5 6 7 10 11 14 15];
vector_two_decimals = [1 8 9 12 13 16 17 18 19];
vector_four_decimals = setdiff(1:length(table_results_columns), [vector_int vector_two_decimals]);
table_results_cell = cell(length(algorithms) * length(TP_targets), length(table_results_columns));
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
    table_results_cell{ii, 2} = algorithms_nice{str2double(table_results_cell{ii, 2})};
end

% shorten table
columns_to_delete = [6:7 8:11 12:15 3];  % 3: remove solvable column
table_results_cell(:, columns_to_delete) = [];
table_results_columns(:, columns_to_delete) = [];

% remove not needed TP_targets and MPE values for reference point
for j = 1:length(TP_targets)
    table_results_cell((j-1) * length(algorithms) + 1, end) = {''};
    table_results_cell((j-1) * length(algorithms) + 1, end-1) = {''};
    table_results_cell((j-1) * length(algorithms) + 1, end-2) = {''};
    table_results_cell((j-1) * length(algorithms) + 1, end-3) = {''};
    for i = 1:length(algorithms)-1
        table_results_cell((j-1) * length(algorithms) + i + 1, 1) = {''};
    end
end

% output
table_results = cell2table(table_results_cell, 'VariableNames', table_results_columns);
writetable(table_results, [ENV.OUTPUT_DIR_EXCEL, instance, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
title_tex = ['&&&& \multicolumn{4}{l}{Computational time (s)} & \multicolumn{4}{l}{$\Delta_\text{costs}$ compared to ', algorithms_nice{1}, '}'];
table2latex_custom(table_results, [ENV.OUTPUT_DIR_LATEX, instance, '.tex'], title_tex);
disp(table_results);

% table for best results
table_best_columns = ['TP^T', 'Algorithm', '\# best solutions', algorithms_nice];
table_best_cell = cell(length(algorithms) * length(TP_targets), length(table_best_columns));
for j = 1:length(TP_targets)
    table_best_cell((j-1) * length(algorithms) + 1, 1) = compose('%.2f', TP_targets(j));
    table_best_cell((j-1) * length(algorithms) + 1 : (j) * length(algorithms), 2) = algorithms_nice;
    table_best_cell((j-1) * length(algorithms) + 1 : (j) * length(algorithms), 3) = compose('%i', best_solution_devision{j}(:, 1));
    table_best_cell((j-1) * length(algorithms) + 1 : (j) * length(algorithms), 4:end) = compose('%.4f', best_solution_devision{j}(:, 2:end));
end
for i = 1:length(algorithms)
    for j = 1:length(TP_targets)
        for ii = 1:length(algorithms)
            table_best_cell((j-1) * length(algorithms) + ii, 3+ii) = {''};
        end
    end
end
table_best = cell2table(table_best_cell, 'VariableNames', table_best_columns);
writetable(table_best, [ENV.OUTPUT_DIR_EXCEL, instance, '_best', '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
title_tex = ['&&& \multicolumn{', num2str(length(algorithms)), '}{l}{Proportion of also found best solutions}'];
table2latex_custom(table_best, [ENV.OUTPUT_DIR_LATEX, instance, '_best', '.tex'], title_tex);
disp(table_best);

% table for better results
table_better_columns = ['TP^T', 'algorithm', '\# solvable instances', algorithms_nice];
table_better_cell = cell(length(algorithms) * length(TP_targets), length(table_better_columns));
for j = 1:length(TP_targets)
    table_better_cell((j-1) * length(algorithms) + 1, 1) = compose('%.2f', TP_targets(j));
    table_better_cell((j-1) * length(algorithms) + 1 : (j) * length(algorithms), 2) = algorithms_nice;
    table_better_cell((j-1) * length(algorithms) + 1 : (j) * length(algorithms), 3) = compose('%i', better_solution_devision{j}(:, 1));
    table_better_cell((j-1) * length(algorithms) + 1 : (j) * length(algorithms), 4:end) = compose('%.4f', better_solution_devision{j}(:, 2:end));
end
for i = 1:length(algorithms)
    for j = 1:length(TP_targets)
        for ii = 1:length(algorithms)
            table_better_cell((j-1) * length(algorithms) + ii, 3+ii) = {''};
        end
    end
end
table_better = cell2table(table_better_cell, 'VariableNames', table_better_columns);
writetable(table_better, [ENV.OUTPUT_DIR_EXCEL, instance, '_better', '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
title_tex = ['&&& \multicolumn{', num2str(length(algorithms)), '}{l}{Proportion of found better solutions}'];
table2latex_custom(table_better, [ENV.OUTPUT_DIR_LATEX, instance, '_better', '.tex'], title_tex);
disp(table_better);

% save
save([ENV.OUTPUT_DIR_DATA, instance, '.mat']);

end
