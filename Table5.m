%% measure time: start
tic();


%% init
run('create_environment')
ENV = get_environment();
file = 'Table5_validation_2M';
file2 = 'Table10_validation_2M_best_solutions';
file3 = 'Table11_validation_2M_better_solutions';


%% results from instances
instance = 'validation_unbalanced_2';
TP_targets = [0.80 0.90];
algorithms = {'enumeration', 'greedy-dec', 'greedy-inc', 'greedy-inc-dec', 'simulated-annealing-dec', 'simulated-annealing-inc', 'genetic-algorithm'};
algorithms_nice = {'ENUM', 'GD', 'GI', 'GID', 'SAD', 'SAI', 'GA'};
create_results_validation(instance, TP_targets, algorithms, algorithms_nice, ENV);
modify_latex_tables(ENV);


%% copy and rename
copyfile([ENV.OUTPUT_DIR_EXCEL, instance, '.xlsx'], [ENV.PAPER_DIR, file, '.xlsx']);
copyfile([ENV.OUTPUT_DIR_EXCEL, instance, '_best', '.xlsx'], [ENV.PAPER_DIR, file2, '.xlsx']);
copyfile([ENV.OUTPUT_DIR_EXCEL, instance, '_better', '.xlsx'], [ENV.PAPER_DIR, file3, '.xlsx']);
copyfile([ENV.OUTPUT_DIR_LATEX_MOD, instance, '.tex'], [ENV.PAPER_DIR, file, '.tex']);
copyfile([ENV.OUTPUT_DIR_LATEX_MOD, instance, '_best', '.tex'], [ENV.PAPER_DIR, file2, '.tex']);
copyfile([ENV.OUTPUT_DIR_LATEX_MOD, instance, '_better', '.tex'], [ENV.PAPER_DIR, file3, '.tex']);



%% measure time: end
toc();

