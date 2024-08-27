function call_genetic_algorithm(algorithm, instances, TP_targets, ENV)
rng(42);


%% parameters
if nargin == 0
    fprintf('No parameters given... using standard values...\n\n');
    algorithm = 'genetic-algorithm';
    instances = {'validation_unbalanced_2'};
    TP_targets = [0.80 0.85 0.90];
    ENV = get_environment();
end


%% python
%pyenv('Version','3.9')
%pyrun(["import platform", "print(platform.python_version())"])
code = [
    "import numpy as np"
    "import pygad"
    "import flowline"
    "import GA_flowline"
];
pyrun(code)


%% call
for i = 1:length(instances)
    instance = instances{i};
    load([ENV.INPUT_DIR, instance, '.mat'], 'machines', 'base_C_min', 'base_Q_min', 'base_C_max', 'base_Q_max', 'max_TP', 'number_instances_total', 'instances_mu', 'instances_p', 'instances_gamma', 'instances_costs_buffers', 'instances_costs_spares');
    
    % prepare dir
    [~,~] = mkdir([ENV.DATA_DIR, instance]);
    
    for TP_target = TP_targets
        dir_local = [ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '/', algorithm, '/'];
        filename_global = [ENV.DATA_DIR, instance, '/', num2str(TP_target, '%.2f'), '_', algorithm];
        [~,~] = mkdir(dir_local);
        delete_if_exists([filename_global, '.xlsx'])
        delete_if_exists([filename_global, '.mat'])
        
        line_general = py.flowline.flowline();
        line_general.fill(int64(machines), int64(base_C_max), int64(base_Q_max), 1, 0.01, 0.1, 1, 1);
        line_general.set_limits(int64(base_C_min), int64(base_C_max), int64(base_Q_min), int64(base_Q_max))
        
        table_values_rows = NaN(number_instances_total, 2*machines - 1 + 6);
        for k = 1:number_instances_total
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
            
            % get general instance data
            mu = instances_mu(k, :);
            p = instances_p(k, :);
            gamma = instances_gamma(k, :);
            costs_buffers = instances_costs_buffers(k, :);
            costs_spares = instances_costs_spares(k, :);
            
            % prepare python
            line = line_general;
            line.mu = pyrun("a = np.array(val, dtype=float)", "a", val = mu);
            line.p = pyrun("a = np.array(val, dtype=float)", "a", val = p);
            line.gamma = pyrun("a = np.array(val, dtype=float)", "a", val = gamma);
            if length(costs_buffers) == 1
                line.costs_buffers = pyrun("a = np.array([val], dtype=float)", "a", val = costs_buffers);
            else
                line.costs_buffers = pyrun("a = np.array(val, dtype=float)", "a", val = costs_buffers);
            end
            line.costs_spares = pyrun("a = np.array(val, dtype=float)", "a", val = costs_spares);
            line.check();
            
            % calculate
            tic();
            result = py.GA_flowline.genetic_algorithm(line, TP_target, false);
            time = toc();
            
            c = double(result.C);
            q = double(result.Q);
            costs = c * costs_buffers' + q * costs_spares';
            counter_iterations = double(result.counter_iterations);
            counter_evaluations = double(result.counter_evaluations);
            counter_problems = double(result.counter_problems);
            if (machines >= 3)
                TP = Spare_Decomposition(machines, c, q, mu, p, gamma, true);
            elseif (machines == 3)
                TP = Spare_3M(c(1), c(2), q(1), q(2), q(3), mu(1), mu(2), mu(3), p(1), p(2), p(3), gamma(1), gamma(2), gamma(3));
            else
                TP = Spare_2M(c, q(1), q(2), mu(1), mu(2), p(1), p(2), gamma(1), gamma(2));
            end
            
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
