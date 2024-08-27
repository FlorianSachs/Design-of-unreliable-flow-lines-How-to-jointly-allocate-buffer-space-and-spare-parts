%% measure time: start
tic();


%% init
run('create_environment')
ENV = get_environment();
file = 'Table2_greedy';


%% parameters
algorithm = 'greedy-inc';
instances = {'explanation_greedy-inc'};
TP_target = 0.85;
check_buffers = 1;
check_spares = 1;
allow_Q_min_adjustments = 1;


%% instance
machines = 2;
base_C_min = 1;
base_Q_min = 1;
base_C_max = 25;
base_Q_max = 5;
number_instances_total = 1;
instances_mu = 1 * ones(1, machines);
instances_p = 0.005 * ones(1, machines);
instances_gamma = 0.1 * ones(1, machines);
instances_costs_buffers = ones(1, machines-1);
instances_costs_spares = ones(1, machines);
max_TP = Spare_2M(base_C_max, base_Q_max, base_Q_max, instances_mu(1), instances_mu(2), instances_p(1), instances_p(2), instances_gamma(1), instances_gamma(2));
save([ENV.INPUT_DIR, instances{1}, '.mat'], 'machines', 'base_C_min', 'base_Q_min', 'base_C_max', 'base_Q_max', 'max_TP', 'number_instances_total', 'instances_mu', 'instances_p', 'instances_gamma', 'instances_costs_buffers', 'instances_costs_spares');


%% results from instances
call_greedy_inc(algorithm, instances, TP_target, check_buffers, check_spares, allow_Q_min_adjustments, ENV)
load([ENV.DATA_DIR, instances{1}, '/', num2str(TP_target, '%.2f'), '/', algorithm, '/', '1_C1_int.mat'], 'table_values_rows');
values = table_values_rows;
values(:, [1 2 8:13]) = round(values(:, [1 2 8:13]), 4);
table_results_rows = cell(5 * 4 + 1, 10);
table_results_rows(1, :) = {0 '' '' '' '' '' '' base_C_min, ['(', int2str(base_Q_min), ',', int2str(base_Q_min), ')'], values(1, 1)};
current_C = base_C_min;
current_Q1 = base_Q_min;
current_Q2 = base_Q_min;
for n = 1:5
    for k = 1:3
        if (k == 1), show_n = n;
        else, show_n = ''; end
        if (k == 1), next_C = current_C + 1;
        else, next_C = current_C; end
        if (k == 2), next_Q1 = current_Q1 + 1;
        else, next_Q1 = current_Q1; end
        if (k == 3), next_Q2 = current_Q2 + 1;
        else, next_Q2 = current_Q2; end
        table_results_rows((n - 1) * 4 + k + 1, :) = {show_n k next_C ['(', int2str(next_Q1), ',', int2str(next_Q2), ')'] num2str(values(n, 10 + k), '%0.4f') num2str(values(n, 7 + k), '%0.4f') '' '' '' ''};
    end
    current_C = values(n, 5);
    current_Q1 = values(n, 6);
    current_Q2 = values(n, 7);
    if isnan(values(n, 3))
        k_chosen = values(n, 4) + 1;
    else
        k_chosen = values(n, 3);
    end
    table_results_rows(n * 4 + 1, :) = {'' '' '' '' '' '' k_chosen values(n, 5), ['(', int2str(values(n, 6)), ',', int2str(values(n, 7)), ')'], num2str(values(n, 2), '%0.4f')};
end


%% copy and rename
table_results = cell2table(table_results_rows, 'VariableNames', {'n' 'k' 'Cnk' 'Qnk' 'TPnk' 'fk' 'k*' 'Cn' 'Qn' 'TPn'});
writetable(table_results, [ENV.PAPER_DIR, file, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
positions = [3, 8, 10, 14, 18];
for i = positions
    tmp = table_results_rows(i, 6);
    table_results_rows(i, 6) = {['\textbf{', tmp{1}, '}']};
end
table_results = cell2table(table_results_rows, 'VariableNames', {'n' 'k' 'Cnk' 'Qnk' 'TPnk' 'fk' 'k*' 'Cn' 'Qn' 'TPn'});
table2latex_custom(table_results, [ENV.PAPER_DIR, file, '.tex']);
disp(table_results);


%% modify LaTex headers
latex_file = readlines([ENV.PAPER_DIR, file, '.tex']);
latex_file{1} = '\begin{tabular}{LLLLLLLLLL}';
latex_file{2} = 'n & k & \mathbf{C}(n, k) & \mathbf{Q}(n, k) & TP(n, k) & f(k) & k^* & \mathbf{C}(n) & \mathbf{Q}(n) & TP(n) \\';
latex_file = [latex_file(1); '\hline'; latex_file(2:4); '\hline'; latex_file(5:8); '\hline'; latex_file(9:12); '\hline'; latex_file(13:16); '\hline'; latex_file(17:20); '\hline'; latex_file(21 :end)];
counter = 0;
for i = 1:length(latex_file)
    counter = counter + contains(latex_file{i}, '\\');
    if counter == 1 
        new = '\\';
    else
        if mod(counter - 1, 4) == 1
            new = '\\[-.05cm]';
        else
            new = '\\[-.15cm]';
        end
    end
    latex_file{i} = strrep(latex_file{i}, '\\', new);
    
end
[fid, msg] = fopen([ENV.PAPER_DIR, file, '.tex'], 'w');
if fid < 1
    error('Could not write output file: %s', msg);
end
fwrite(fid, strjoin(latex_file, '\n'));
fclose(fid);


%% measure time: end
toc();

