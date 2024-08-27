%% measure time: start
tic();


%% init
run('create_environment')
ENV = get_environment();
create_figures = 1;
file = 'Figure4_bottlenecks';
file_table = ['Table16_for_', file];


%% results from instances
instance = 'insights_d_bottlenecks_5';
TP_targets = [0.80 0.90];
create_results_insights_d(instance, TP_targets, create_figures, ENV);
modify_latex_tables(ENV);


%% copy and rename
parts = {'buffers', 'spares'};
position = {'top', 'bottom'};
copyfile([ENV.OUTPUT_DIR_EXCEL, instance, '_allocation', '.xlsx'], [ENV.PAPER_DIR, file_table, '.xlsx']);
copyfile([ENV.OUTPUT_DIR_LATEX_MOD, instance, '_allocation', '.tex'], [ENV.PAPER_DIR, file_table, '.tex']);
for i = 1:length(parts)
    copyfile([ENV.OUTPUT_DIR_FIGURES, instance, '_', parts{i}, '.pdf'], [ENV.PAPER_DIR, file, '_', position{i}, '.pdf']);
    copyfile([ENV.OUTPUT_DIR_FIGURES, instance, '_', parts{i}, '.png'], [ENV.PAPER_DIR, file, '_', position{i}, '.png']);
end

latex_file = readlines([ENV.PAPER_DIR, file_table, '.tex']);
for i = 1:length(latex_file)
    latex_file{i} = replace(latex_file{i}, 'None', '\text{None}');
end

[fid, msg] = fopen([ENV.PAPER_DIR, file_table, '.tex'], 'w');
if fid < 1
    error('Could not write output file: %s', msg);
end
fwrite(fid, strjoin(latex_file, '\n'));
fclose(fid);


%% measure time: end
toc();

