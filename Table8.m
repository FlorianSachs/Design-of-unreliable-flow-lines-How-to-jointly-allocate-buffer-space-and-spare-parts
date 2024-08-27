%% measure time: start
tic();


%% init
run('create_environment')
ENV = get_environment();
file = 'Table8_example_instance';


%% results from instances
instance = 'insights_a_separate_3';
load([ENV.INPUT_DIR, instance, '.mat'], 'instances_mu', 'instances_p', 'instances_gamma');


%% copy and rename
selected_example = 18;
table_results_rows = cell(3, 4);
for i = 1:3
    table_results_rows{i, 1} = i;
    table_results_rows{i, 2} = num2str(instances_mu(selected_example, i), '%0.6f');
    table_results_rows{i, 3} = num2str(instances_p(selected_example, i), '%0.6f');
    table_results_rows{i, 4} = num2str(instances_gamma(selected_example, i), '%0.6f');
end

table_results = cell2table(table_results_rows, 'VariableNames', {'Machine i', 'mu_i', 'p_i', 'gamma_i'});
writetable(table_results, [ENV.PAPER_DIR, file, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
table2latex_custom(table_results, [ENV.PAPER_DIR, file, '.tex']);
disp(table_results);


%% modify LaTex headers
latex_file = readlines([ENV.PAPER_DIR, file, '.tex']);
latex_file{1} = '\begin{tabular}{lLLL}';
latex_file{2} = 'Machine $i$ & \mu_i & p_i & \gamma_i \\';
latex_file = [latex_file(1); '\hline'; latex_file(2:end)];
[fid, msg] = fopen([ENV.PAPER_DIR, file, '.tex'], 'w');
if fid < 1
    error('Could not write output file: %s', msg);
end
fwrite(fid, strjoin(latex_file, '\n'));
fclose(fid);


%% measure time: end
toc();

