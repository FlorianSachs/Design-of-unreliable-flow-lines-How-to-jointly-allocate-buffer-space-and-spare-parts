%% measure time: start
tic();


%% init
run('create_environment')
ENV = get_environment();
file = 'Table7_integrated';


%% results from instances
instance = 'insights_a_separate_3';
TP_targets = [0.65 0.75 0.85];
create_results_insights_a(instance, TP_targets, ENV);
modify_latex_tables(ENV);


%% copy and rename
copyfile([ENV.OUTPUT_DIR_EXCEL, instance, '.xlsx'], [ENV.PAPER_DIR, file, '.xlsx']);
copyfile([ENV.OUTPUT_DIR_LATEX_MOD, instance, '.tex'], [ENV.PAPER_DIR, file, '.tex']);



%% measure time: end
toc();

