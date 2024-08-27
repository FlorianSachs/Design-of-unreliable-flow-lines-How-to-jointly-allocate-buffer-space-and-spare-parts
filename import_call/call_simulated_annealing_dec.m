function call_simulated_annealing_dec(algorithm, instances, TP_targets, ENV)
rng(42);


%% parameters
if nargin == 0
    fprintf('No parameters given... using standard values...\n\n');
    algorithm = 'simulated-annealing-dec';
    instances = {'validation_unbalanced_2'};
    TP_targets = [0.80 0.85 0.90];
    ENV = get_environment();
end


%% algorithm parameters
alg.type = 'logarithmic';
alg.start_temperature = 20;
alg.param = NaN;
alg.T_min = 2;


%% call
for i = 1:length(instances)
    instance = instances{i};
    load([ENV.INPUT_DIR, instance, '.mat'], 'machines', 'base_C_min', 'base_Q_min', 'base_C_max', 'base_Q_max', 'max_TP', 'number_instances_total', 'instances_mu', 'instances_p', 'instances_gamma', 'instances_costs_buffers', 'instances_costs_spares');
    
    % variables
    C_min = base_C_min * ones(1, machines-1);
    Q_min = base_Q_min * ones(1, machines);
    C_max = base_C_max * ones(1, machines-1);
    Q_max = base_Q_max * ones(1, machines);
    alg.max_iterations = (2 * machines - 1) * 200;
    
    % prepare dir
    [~,~] = mkdir([ENV.DATA_DIR, instance]);
    
    for TP_target = TP_targets
        dir_local = [ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '/', algorithm, '/'];
        filename_global = [ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_', algorithm];
        [~,~] = mkdir(dir_local);
        delete_if_exists([filename_global, '.xlsx'])
        delete_if_exists([filename_global, '.mat'])
        
        table_values_rows = NaN(number_instances_total, 2*machines - 1 + 6);
        parfor k = 1:number_instances_total
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
            
            c = NaN(1, machines-1);
            q = NaN(1, machines);
            
            % get general instance data
            mu = instances_mu(k, :);
            p = instances_p(k, :);
            gamma = instances_gamma(k, :);
            costs_buffers = instances_costs_buffers(k, :);
            costs_spares = instances_costs_spares(k, :);
            
            % calculate
            tic();
            [c(1, 1:machines-1), q(1, 1:machines), TP, costs, counter_iterations, counter_evaluations, counter_problems] = ...
                optimize_simulated_annealing(machines, C_min, Q_min, C_max, Q_max, mu, p, gamma, TP_target, costs_buffers, costs_spares, C_max, Q_max, dir_local, filename_local, alg);
            time = toc();
            
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
