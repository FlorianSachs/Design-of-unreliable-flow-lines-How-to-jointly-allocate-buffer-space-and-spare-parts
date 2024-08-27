%% clean up
run('create_environment')
ENV = get_environment();
delete_if_exists('Times.mat');
save('Times.mat', 'ENV');


%% check Python
code = [
    "import numpy as np"
    "import pygad"
    "import flowline"
    "import GA_flowline"
];
pyrun(code)


%% instances
datetimeInstances_start = datetime;
save('Times.mat', 'datetimeInstances_start', '-append');
tic();
generate_instances_validation_unbalanced(ENV);
generate_instances_insights_a_separate(ENV);
generate_instances_insights_b_costs_all_spares(ENV);
generate_instances_insights_c_costs_one_spare(ENV);
generate_instances_insights_d_bottlenecks(ENV);
timeInstances = toc();
datetimeInstances_end = datetime;
load('Times.mat', 'datetimeInstances_start');
durationInstances = datetimeInstances_end - datetimeInstances_start;
save('Times.mat', 'timeInstances', 'datetimeInstances_end', 'durationInstances', '-append');


%% optimization
datetimeOptimization_start = datetime;
save('Times.mat', 'datetimeOptimization_start');
tic();

% validation big
TP_targets = [0.80 0.90];
check_buffers = true;
check_spares = true;
base_C_start = false;
base_Q_start = false;
allow_Q_min_adjustements = true;
instances = {'validation_unbalanced_2', 'validation_unbalanced_5'};

call_greedy_inc('greedy-inc', instances, TP_targets, check_buffers, check_spares, allow_Q_min_adjustements, ENV);
call_greedy_inc_dec('greedy-inc-dec', instances, TP_targets, check_buffers, check_spares, 'greedy-inc', allow_Q_min_adjustements, ENV);
call_greedy_dec('greedy-dec', instances, TP_targets, check_buffers, check_spares, base_C_start, base_Q_start, allow_Q_min_adjustements, ENV);
call_simulated_annealing_inc('simulated-annealing-inc', instances, TP_targets, ENV);
call_simulated_annealing_dec('simulated-annealing-dec', instances, TP_targets, ENV);


% validation small
TP_targets = [0.80 0.90];
instances = {'validation_unbalanced_2'};
% 
call_enumeration('enumeration', instances, TP_targets, ENV);
call_genetic_algorithm('genetic-algorithm', instances, TP_targets, ENV);


% insights: a
TP_targets = [0.65 0.75 0.85];
instances = {'insights_a_separate_3'};

types = {'onlyspares', 'onlybuffers', 'both'};
for i = 1:length(types)
    type = types{i};
    switch type
        case 'onlyspares'
            check_buffers = false;
            check_spares = true;
            base_C_start = 1;
            base_Q_start = false;
			allow_Q_min_adjustments = true;
        case 'onlybuffers'
            check_buffers = true;
            check_spares = false;
            base_C_start = false;
            base_Q_start = 1;
			allow_Q_min_adjustments = false;
        case 'both'
            check_buffers = true;
            check_spares = true;
            base_C_start = false;
            base_Q_start = false;
			allow_Q_min_adjustments = true;
    end
    call_greedy_inc(['greedy-inc_', type], instances, TP_targets, check_buffers, check_spares, allow_Q_min_adjustments, ENV);
    call_greedy_inc_dec(['greedy-inc-dec_', type], instances, TP_targets, check_buffers, check_spares, ['greedy-inc_', type], allow_Q_min_adjustments, ENV);
    call_greedy_dec(['greedy-dec_', type], instances, TP_targets, check_buffers, check_spares, base_C_start, base_Q_start, allow_Q_min_adjustments, ENV);
end

% insights: b, c, d
TP_targets = [0.80 0.90];
check_buffers = true;
check_spares = true;
base_C_start = false;
base_Q_start = false;
allow_Q_min_adjustements = true;
instances = {'insights_b_costs_all_spares_5', 'insights_c_costs_one_spare_5', ...
             'insights_d_bottlenecks_5'};

call_greedy_inc('greedy-inc', instances, TP_targets, check_buffers, check_spares, allow_Q_min_adjustements, ENV);
call_greedy_inc_dec('greedy-inc-dec', instances, TP_targets, check_buffers, check_spares, 'greedy-inc', allow_Q_min_adjustements, ENV);
call_greedy_dec('greedy-dec', instances, TP_targets, check_buffers, check_spares, base_C_start, base_Q_start, allow_Q_min_adjustements, ENV);

timeOptimization = toc();
datetimeOptimization_end = datetime;
load('Times.mat', 'datetimeOptimization_start');
durationOptimization = datetimeOptimization_end - datetimeOptimization_start;
save('Times.mat', 'timeOptimization', 'datetimeOptimization_end', 'durationOptimization', '-append');


%% output
datetimeFigure2_start = datetime;
save('Times.mat', 'datetimeFigure2_start', '-append');
tic();
Figure2;
timeFigure2 = toc();
datetimeFigure2_end = datetime;
load('Times.mat', 'datetimeFigure2_start');
durationFigure2 = datetimeFigure2_end - datetimeFigure2_start;
save('Times.mat', 'timeFigure2', 'datetimeFigure2_end', 'durationFigure2', '-append');

datetimeFigure3_start = datetime;
save('Times.mat', 'datetimeFigure3_start', '-append');
tic();
Figure3;
timeFigure3 = toc();
datetimeFigure3_end = datetime;
load('Times.mat', 'datetimeFigure3_start');
durationFigure3 = datetimeFigure3_end - datetimeFigure3_start;
save('Times.mat', 'timeFigure3', 'datetimeFigure3_end', 'durationFigure3', '-append');

datetimeFigure4_start = datetime;
save('Times.mat', 'datetimeFigure4_start', '-append');
tic();
Figure4;
timeFigure4 = toc();
datetimeFigure4_end = datetime;
load('Times.mat', 'datetimeFigure4_start');
durationFigure4 = datetimeFigure4_end - datetimeFigure4_start;
save('Times.mat', 'timeFigure4', 'datetimeFigure4_end', 'durationFigure4', '-append');

datetimeTable2_start = datetime;
save('Times.mat', 'datetimeTable2_start', '-append');
tic();
Table2;
timeTable2 = toc();
datetimeTable2_end = datetime;
load('Times.mat', 'datetimeTable2_start');
durationTable2 = datetimeTable2_end - datetimeTable2_start;
save('Times.mat', 'timeTable2', 'datetimeTable2_end', 'durationTable2', '-append');

datetimeTable5_start = datetime;
save('Times.mat', 'datetimeTable5_start', '-append');
tic();
Table5;
timeTable5 = toc();
datetimeTable5_end = datetime;
load('Times.mat', 'datetimeTable5_start');
durationTable5 = datetimeTable5_end - datetimeTable5_start;
save('Times.mat', 'timeTable5', 'datetimeTable5_end', 'durationTable5', '-append');

datetimeTable6_start = datetime;
save('Times.mat', 'datetimeTable6_start', '-append');
tic();
Table6;
timeTable6 = toc();
datetimeTable6_end = datetime;
load('Times.mat', 'datetimeTable6_start');
durationTable6 = datetimeTable6_end - datetimeTable6_start;
save('Times.mat', 'timeTable6', 'datetimeTable6_end', 'durationTable6', '-append');

datetimeTable7_start = datetime;
save('Times.mat', 'datetimeTable7_start', '-append');
tic();
Table7;
timeTable7 = toc();
datetimeTable7_end = datetime;
load('Times.mat', 'datetimeTable7_start');
durationTable7 = datetimeTable7_end - datetimeTable7_start;
save('Times.mat', 'timeTable7', 'datetimeTable7_end', 'durationTable7', '-append');

datetimeTable8_start = datetime;
save('Times.mat', 'datetimeTable8_start', '-append');
tic();
Table8;
timeTable8 = toc();
datetimeTable8_end = datetime;
load('Times.mat', 'datetimeTable8_start');
durationTable8 = datetimeTable8_end - datetimeTable8_start;
save('Times.mat', 'timeTable8', 'datetimeTable8_end', 'durationTable8', '-append');

datetimeTable9_start = datetime;
save('Times.mat', 'datetimeTable9_start', '-append');
tic();
Table9;
timeTable9 = toc();
datetimeTable9_end = datetime;
load('Times.mat', 'datetimeTable9_start');
durationTable9 = datetimeTable9_end - datetimeTable9_start;
save('Times.mat', 'timeTable9', 'datetimeTable9_end', 'durationTable9', '-append');

