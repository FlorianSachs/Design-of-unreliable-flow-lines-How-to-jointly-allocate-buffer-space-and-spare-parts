function create_results_insights_b(instance, TP_targets, create_figures, ENV)
rng(12);


%% parameters
if nargin == 0
    instance = 'insights_b_costs_all_spares_5';
    TP_targets = [0.75 0.80 0.85 0.90];
    create_figures = true;
    ENV = get_environment();
end


%% settings
algorithms = {'greedy-dec', 'greedy-inc-dec'};
parts_nice = {'Very cheap spares', 'Cheap spares', 'Equally', 'Expensive spares', 'Very expensive spares'};


%% call
load([ENV.INPUT_DIR, instance, '.mat'], 'machines', 'number_instances', 'number_instances_total');
number_parts = number_instances_total / number_instances;

data_TP = cell(number_parts, length(TP_targets));
data_costs = cell(number_parts, length(TP_targets));
data_C = cell(number_parts, length(TP_targets));
data_Q = cell(number_parts, length(TP_targets));
data_C_total = cell(number_parts, length(TP_targets));
data_Q_total = cell(number_parts, length(TP_targets));
for j = 1:length(TP_targets)
    TP_target = TP_targets(j);
    tmp_TP = NaN(number_instances_total, length(algorithms));
    tmp_costs = NaN(number_instances_total, length(algorithms));
    tmp_C = cell(1, length(algorithms));
    tmp_Q = cell(1, length(algorithms));
    for i = 1:length(algorithms)
        tmp_C{i} = NaN(number_instances_total, machines-1);
        tmp_Q{i} = NaN(number_instances_total, machines);
        algorithm = algorithms{i};
        load([ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_', algorithm, '.mat'], 'table_values_rows');
        tmp_TP(:, i) = table_values_rows(:, 1);
        tmp_costs(:, i) = table_values_rows(:, 2);
        idx = tmp_TP < TP_target;
        tmp_TP(idx) = NaN;
        tmp_costs(idx) = NaN;
        tmp_C{i}(:, :) = table_values_rows(:, 4:(4+machines-2));
        tmp_Q{i}(:, :) = table_values_rows(:, (4+machines-2+1):((4+machines-2+1)+machines-1));
    end
    [~, pos] = min(tmp_costs, [], 2);
    opt_TP = NaN(number_instances_total, 1);
    opt_costs = NaN(number_instances_total, 1);
    opt_C = NaN(number_instances_total, machines-1);
    opt_Q = NaN(number_instances_total, machines);
    opt_C_total = NaN(number_instances_total, 1);
    opt_Q_total = NaN(number_instances_total, 1);
    for k = 1:number_instances_total
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
    for p = 1:number_parts
        data_TP{p, j} = opt_TP((p-1) * number_instances + 1 : p * number_instances);
        data_costs{p, j} = opt_costs((p-1) * number_instances + 1 : p * number_instances);
        data_C{p, j} = opt_C((p-1) * number_instances + 1 : p * number_instances, :);
        data_Q{p, j} = opt_Q((p-1) * number_instances + 1 : p * number_instances, :);
        data_C_total{p, j} = opt_C_total((p-1) * number_instances + 1 : p * number_instances);
        data_Q_total{p, j} = opt_Q_total((p-1) * number_instances + 1 : p * number_instances);
    end
end

number_solvable = NaN(number_parts, length(TP_targets));
costs_avg = NaN(number_parts, length(TP_targets));
C_avg = NaN(number_parts, length(TP_targets));
C_std_error = NaN(number_parts, length(TP_targets));
C_min = NaN(number_parts, length(TP_targets));
C_max = NaN(number_parts, length(TP_targets));
Q_avg = NaN(number_parts, length(TP_targets));
Q_std_error = NaN(number_parts, length(TP_targets));
Q_min = NaN(number_parts, length(TP_targets));
Q_max = NaN(number_parts, length(TP_targets));
BI_avg = cell(number_parts, length(TP_targets));
BI_std_error = cell(number_parts, length(TP_targets));
UI_avg = cell(number_parts, length(TP_targets));
UI_std_error = cell(number_parts, length(TP_targets));
for p = 1:number_parts
    for j = 1:length(TP_targets)
        number_solvable(p, j) = number_instances - sum(isnan(data_costs{p, j}));
        costs_avg(p, j) = mean(data_costs{p, j}, 'omitnan');
        C_avg(p, j) = mean(data_C_total{p, j}, 'omitnan');
        C_std_error(p, j) = std(data_C_total{p, j}, 'omitnan') / sqrt(number_solvable(p, j));
        C_min(p, j) = min(data_C_total{p, j});
        C_max(p, j) = max(data_C_total{p, j});
        Q_avg(p, j) = mean(data_Q_total{p, j}, 'omitnan');
        Q_std_error(p, j) = std(data_Q_total{p, j}, 'omitnan') / sqrt(number_solvable(p, j));
        Q_min(p, j) = min(data_Q_total{p, j});
        Q_max(p, j) = max(data_Q_total{p, j});
        BI_avg{p, j} = mean(data_C{p, j}, 'omitnan');
        BI_std_error{p, j} = std(data_C{p, j}, 'omitnan') / sqrt(number_solvable(p, j));
        UI_avg{p, j} = mean(data_Q{p, j}, 'omitnan');
        UI_std_error{p, j} = std(data_Q{p, j}, 'omitnan') / sqrt(number_solvable(p, j));
    end
end

% table for overall results
table_results_columns = {'TP^T', 'case', '\# solvable', 'mean costs', 'TBC (mean)', 'TBC (std. error)', 'TBC (min)', 'TBC (max)', 'TCU (mean)', 'TCU (std. error)', 'TCU (min)', 'TCU (max)'};
table_results_rows = NaN(number_parts * length(TP_targets), length(table_results_columns));
for j = 1:length(TP_targets)
    for p = 1:number_parts
        table_results_rows((j - 1) * number_parts + p, :) = ...
            [TP_targets(j) p number_solvable(p, j) costs_avg(p, j) ...
            C_avg(p, j) C_std_error(p, j) C_min(p, j) C_max(p, j) ...
            Q_avg(p, j) Q_std_error(p, j) Q_min(p, j) Q_max(p, j)];
    end
end

% format
vector_int = [3 7 8 11 12];
vector_two_decimals = [1 4 5 6 9 10];
vector_four_decimals = setdiff(1:length(table_results_columns), [vector_int vector_two_decimals]);
table_results_cell = cell(number_parts * length(TP_targets), length(table_results_columns));
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
columns_to_delete = [4];
table_results_cell(:, columns_to_delete) = [];
table_results_columns(:, columns_to_delete) = [];

% remove unnecessary TP^T entries
for i = 1:number_parts
    for j = 1:length(TP_targets)
        if i > 1
            table_results_cell((j-1) * number_parts + i, 1) = {''};
        end
    end
end

% output
table_results = cell2table(table_results_cell, 'VariableNames', table_results_columns);
writetable(table_results, [ENV.OUTPUT_DIR_EXCEL, instance, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
title_tex = '&&& \multicolumn{4}{l}{Total buffer capacity (TBC)} & \multicolumn{4}{l}{Total component units (TCU)}';
table2latex_custom(table_results, [ENV.OUTPUT_DIR_LATEX, instance, '.tex'], title_tex);
disp(table_results);

% table with allocations
table_results_columns_allocation = ['TP^T', 'Case', cellfun(@(c)['C_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false), cellfun(@(c)['Q_', num2str(c)],num2cell(1:machines),'UniformOutput',false)];
table_results_cell_allocation = cell(number_parts * length(TP_targets), length(table_results_columns_allocation));
for j = 1:length(TP_targets)
    for p = 1:number_parts
        ii = (j - 1) * number_parts + p;
        table_results_cell_allocation{ii, 1} = table_results_cell{ii, 1};
        table_results_cell_allocation{ii, 2} = table_results_cell{ii, 2};
        BI = reshape([BI_avg{p, j}; BI_std_error{p, j}], 1, []);
        UI = reshape([UI_avg{p, j}; UI_std_error{p, j}], 1, []);
        BI_cell = compose('%.1f~(%.1f)', BI);
        UI_cell = compose('%.1f~(%.1f)', UI);
        for m = 1:machines-1
            table_results_cell_allocation{ii, m + 2} = BI_cell{m};
        end
        for m = 1:machines
            table_results_cell_allocation{ii, m + 2+machines-1} = UI_cell{m};
        end
    end
end
table_results_allocation = cell2table(table_results_cell_allocation, 'VariableNames', table_results_columns_allocation);
writetable(table_results_allocation, [ENV.OUTPUT_DIR_EXCEL, instance, '_allocation', '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
title_tex = '&& \multicolumn{4}{l}{Mean buffer capacities (std. error)} & \multicolumn{5}{l}{Mean number of component units (std. error)}';
table2latex_custom(table_results_allocation, [ENV.OUTPUT_DIR_LATEX, instance, '_allocation', '.tex'], title_tex);
disp(table_results_allocation);

if create_figures
    % figure: buffers
    filename_spec = [instance, '_buffers'];
    label = categorical(cellfun(@(c)['C_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false));
    f = figure();
    set(f, 'Position', [200 200 1200 length(TP_targets) * 200]);
    t = tiledlayout(length(TP_targets), number_parts);
    for j = 1:length(TP_targets)
        maxBI = 0;
        for p = 1:number_parts
            maxBI = max(maxBI, max(BI_avg{p, j}));
        end
        for p = 1:number_parts
            nexttile
            bar(label, BI_avg{p, j});
            if p == 1, ylabel(num2str(TP_targets(j), '%.2f')); end
            if j == length(TP_targets), xlabel(parts_nice(p)); end
            %ylim([0 max(BI_avg{p, j}) + 2]);
            ylim([0 maxBI + 2]);
            hold on
            er = errorbar(BI_avg{p, j}, BI_std_error{p, j});
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            hold off
        end
    end
    title(t, 'Mean buffer capacities');
    xlabel(t, 'Spare part costs');
    ylabel(t, 'Target throughput');
    saveas(f, [ENV.OUTPUT_DIR_DATA, filename_spec, '.fig']);
    exportgraphics(f, [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.png'], 'Resolution',600);
    exportgraphics(f, [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.pdf'], 'ContentType','vector');
    system(['pdfcrop "', [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.pdf'], '" "', [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.pdf'], '"']);

    % figure: spares
    filename_spec = [instance, '_spares'];
    label = categorical(cellfun(@(c)['Q_', num2str(c)],num2cell(1:machines),'UniformOutput',false));
    f = figure();
    set(f, 'Position', [200 200 1200 length(TP_targets) * 200]);
    t = tiledlayout(length(TP_targets), number_parts);
    for j = 1:length(TP_targets)
        maxUI = 0;
        for p = 1:number_parts
            maxUI = max(maxUI, max(UI_avg{p, j}));
        end
        for p = 1:number_parts
            nexttile
            bar(label, UI_avg{p, j});
            if p == 1, ylabel(num2str(TP_targets(j), '%.2f')); end
            if j == length(TP_targets), xlabel(parts_nice(p)); end
            %ylim([0 max(UI_avg{p, j}) + 2]);
            ylim([0 maxUI + 2]);
            hold on
            er = errorbar(UI_avg{p, j}, UI_std_error{p, j});
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            hold off
        end
    end
    title(t, 'Mean number of component units');
    xlabel(t, 'Spare part costs');
    ylabel(t, 'Target throughput');
    saveas(f, [ENV.OUTPUT_DIR_DATA, filename_spec, '.fig']);
    exportgraphics(f, [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.png'], 'Resolution',600);
    exportgraphics(f, [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.pdf'], 'ContentType','vector');
    system(['pdfcrop "', [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.pdf'], '" "', [ENV.OUTPUT_DIR_FIGURES, filename_spec, '.pdf'], '"']);
    close('all');
    clear f;
end

% save
save([ENV.OUTPUT_DIR_DATA, instance, '.mat']);

end
