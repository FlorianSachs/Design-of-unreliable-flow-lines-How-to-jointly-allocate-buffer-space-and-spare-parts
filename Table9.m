%% measure time: start
tic();


%% init
run('create_environment')
ENV = get_environment();
file = 'Table9_example_instance_results';


%% results from instances
instance = 'insights_a_separate_3';
TP_target = 0.65;
load([ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_greedy-inc-dec_both', '.mat'], 'table_values_rows');
table_values_rows_integrated = table_values_rows;
load([ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_greedy-inc-dec_onlybuffers', '.mat'], 'table_values_rows');
table_values_rows_onlybuffers = table_values_rows;
load([ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_greedy-inc-dec_onlyspares', '.mat'], 'table_values_rows');
table_values_rows_onlyspares = table_values_rows;


%% copy and rename
selected_example = 18;
table_results_rows = num2cell([zeros(3,1), ...
    [round(table_values_rows_integrated(selected_example, 1), 4); round(table_values_rows_onlybuffers(selected_example, 1), 4); round(table_values_rows_onlyspares(selected_example, 1), 4)], ...
    [table_values_rows_integrated(selected_example, 2); table_values_rows_onlybuffers(selected_example, 2); table_values_rows_onlyspares(selected_example, 2)], ...
    [sum(table_values_rows_integrated(selected_example, 4:5)); sum(table_values_rows_onlybuffers(selected_example, 4:5)); sum(table_values_rows_onlyspares(selected_example, 4:5))], ...
    [sum(table_values_rows_integrated(selected_example, 6:8)); sum(table_values_rows_onlybuffers(selected_example, 6:8)); sum(table_values_rows_onlyspares(selected_example, 6:8))], ...
    [table_values_rows_integrated(selected_example, 4:5); table_values_rows_onlybuffers(selected_example, 4:5); table_values_rows_onlyspares(selected_example, 4:5)], ...
    [table_values_rows_integrated(selected_example, 6:8); table_values_rows_onlybuffers(selected_example, 6:8); table_values_rows_onlyspares(selected_example, 6:8)]]);
table_results_rows(:, 1) = {'Integrated'; 'Only buffers'; 'Only spares'};
table_results_rows{1, 2} = num2str(table_values_rows_integrated(selected_example, 1), '%0.4f');
table_results_rows{2, 2} = num2str(table_values_rows_onlybuffers(selected_example, 1), '%0.4f');
table_results_rows{3, 2} = num2str(table_values_rows_onlyspares(selected_example, 1), '%0.4f');
table_results = cell2table(table_results_rows, 'VariableNames', {'Case', 'TP', 'TC', 'TBC', 'TCU', 'C_1', 'C_2', 'Q_1', 'Q_2', 'Q_3'});
writetable(table_results, [ENV.PAPER_DIR, file, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
table2latex_custom(table_results, [ENV.PAPER_DIR, file, '.tex']);
disp(table_results);


%% modify LaTex headers
latex_file = readlines([ENV.PAPER_DIR, file, '.tex']);
latex_file{1} = '\begin{tabular}{lLRRRRRRRR}';
latex_file = [latex_file(1); '\hline'; latex_file(2:end)];
[fid, msg] = fopen([ENV.PAPER_DIR, file, '.tex'], 'w');
if fid < 1
    error('Could not write output file: %s', msg);
end
fwrite(fid, strjoin(latex_file, '\n'));
fclose(fid);


%% measure time: end
toc();

