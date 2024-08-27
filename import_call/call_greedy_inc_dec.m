function call_greedy_inc_dec(algorithm, instances, TP_targets, check_buffers, check_spares, dir_greedy_inc, allow_Q_min_adjustments, ENV)
rng(42);


%% parameters
if nargin == 0
    fprintf('No parameters given... using standard values...\n\n');
    algorithm = 'greedy-inc-dec';
    instances = {'validation_unbalanced_2'};
    TP_targets = [0.80 0.85 0.90];
    check_buffers = true;
    check_spares = true;
    dir_greedy_inc = 'greedy-inc';
    ENV = get_environment();
end


%% call
for i = 1:length(instances)
    instance = instances{i};
    load([ENV.INPUT_DIR, instance, '.mat'], 'machines', 'base_C_min', 'base_Q_min', 'base_C_max', 'base_Q_max', 'max_TP', 'number_instances_total', 'instances_mu', 'instances_p', 'instances_gamma', 'instances_costs_buffers', 'instances_costs_spares');
    
    % variables
    C_min = base_C_min * ones(1, machines-1);
    Q_min = base_Q_min * ones(1, machines);
    C_max = base_C_max * ones(1, machines-1);
    Q_max = base_Q_max * ones(1, machines);
    
    % prepare dir
    [~,~] = mkdir([ENV.DATA_DIR, instance]);
    
    for TP_target = TP_targets
        dir_local_inc = [ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '/', dir_greedy_inc, '/'];
        dir_local = [ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '/', algorithm, '/'];
        filename_global = [ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_', algorithm];
        [~,~] = mkdir(dir_local);
        delete_if_exists([filename_global, '.xlsx'])
        delete_if_exists([filename_global, '.mat'])
        
        table_values_rows = NaN(number_instances_total, 2*machines - 1 + 6);
        parfor k = 1:number_instances_total
            TP_inc = NaN;
            costs_inc = NaN;
            time_inc = NaN;
            c_inc = NaN(1, machines - 1);
            q_inc = NaN(1, machines);
            counter_iterations_inc = NaN;
            counter_evaluations_inc = NaN;
            counter_problems_inc = NaN;
            
            filename_local = num2str(k);
            if exist([dir_local, filename_local, '.mat'], 'file') == 2
                data = load_data_opt([dir_local, filename_local, '.mat']);
                table_values_rows(k, :) = data;
                continue;
            end
            if max_TP(k) < TP_target
                data = -1 * ones(1, 2*machines - 1 + 6);
                table_values_rows(k, :) = data;
                save_data_opt([dir_local, filename_local, '.mat'], data);
                continue;
            end
            
            c_dec = NaN(1, machines-1);
            q_dec = NaN(1, machines);
            
            % get general instance data
            mu = instances_mu(k, :);
            p = instances_p(k, :);
            gamma = instances_gamma(k, :);
            costs_buffers = instances_costs_buffers(k, :);
            costs_spares = instances_costs_spares(k, :);
            
            % get values from inc
            if exist([dir_local_inc, filename_local, '.mat'], 'file') == 2
                [status,message,messageId] = copyfile([dir_local_inc, filename_local, '.xlsx'], [dir_local, filename_local, '_inc', '.xlsx']);
                data = load_data_opt([dir_local_inc, filename_local, '.mat']);
                TP_inc = data(1);
                costs_inc = data(2);
                time_inc = data(3);
                c_inc = data(4:(4 + machines - 2));
                q_inc = data((4 + machines - 1):(4 + 2*machines - 2));
                counter_iterations_inc = data(end - 2);
                counter_evaluations_inc = data(end - 1);
                counter_problems_inc = data(end);
            else
                error(['Instance ', instance, ' with TP_target = ', num2str(TP_target, '%.2f'), ' and k = ', num2str(k), ' has no greedy-inc file.'])
            end
            
            % calculate
            tic();
            [c_dec(1, 1:machines-1), q_dec(1, 1:machines), TP_dec, costs_dec, counter_iterations_dec, counter_evaluations_dec, counter_problems_dec] = ...
                optimize_greedy_dec(machines, C_min, Q_min, C_max, Q_max, mu, p, gamma, TP_target, costs_buffers, costs_spares, c_inc, q_inc, dir_local, [filename_local, '_dec'], check_buffers, check_spares, allow_Q_min_adjustments, '');
            time_dec = toc();

            c = c_dec;
            q = q_dec;
            TP = TP_dec;
            costs = costs_dec;
            time = time_inc + time_dec;
            counter_iterations = counter_iterations_inc + counter_iterations_dec;
            counter_evaluations = counter_evaluations_inc + counter_evaluations_dec;
            counter_problems = counter_problems_inc + counter_problems_dec;
            
            % save values
            data = [TP costs time c(1, :) q(1, :) counter_iterations counter_evaluations counter_problems];
            table_values_rows(k, :) = data;
            save_data_opt([dir_local, filename_local, '.mat'], data);
        end
    
    % output
    table = array2table(table_values_rows, 'VariableNames', [{'TP', 'costs', 'time'}, [cellfun(@(c)['C_', num2str(c)], num2cell(1:(machines-1)),'UniformOutput',false) cellfun(@(c)['S_', num2str(c)], num2cell(1:(machines)),'UniformOutput',false)], {'iterations', 'evaluations', 'problems'}]);
    writetable(table,[filename_global, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
    disp(table);
    
    % save data
    save([filename_global, '.mat']);
    
    end
end
