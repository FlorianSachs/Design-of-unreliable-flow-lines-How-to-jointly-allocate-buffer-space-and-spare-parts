%% measure time: start
tic();


%% init
run('create_environment')
ENV = get_environment();
create_figures = 1;
file = 'Figure2_costs_all_spares';
file_table = ['Table14_for_', file];


%% results from instances
instance = 'insights_b_costs_all_spares_5';
TP_targets = [0.80 0.90];
create_results_insights_b(instance, TP_targets, create_figures, ENV);
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



%% measure time: end
toc();

