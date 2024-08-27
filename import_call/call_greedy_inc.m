function call_greedy_inc(algorithm, instances, TP_targets, check_buffers, check_spares, allow_Q_min_adjustments, ENV)
rng(42);


%% parameters
if nargin == 0
    fprintf('No parameters given... using standard values...\n\n');
    algorithm = 'greedy-inc';
    instances = {'validation_unbalanced_2'};
    TP_targets = [0.80 0.85 0.90];
    check_buffers = true;
    check_spares = true;
    allow_Q_min_adjustments = true;
    ENV = get_environment();
end
TP_targets = sort(TP_targets);


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

    for TP_target_number = 1:length(TP_targets)
        TP_target = TP_targets(TP_target_number);
        if TP_target_number > 1
            dir_to_load = [ENV.DATA_DIR, instance, '/', num2str(TP_targets(TP_target_number - 1), '%.2f'), '/', algorithm, '/'];
        else
            dir_to_load = '';
        end

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

            TP = NaN;
            costs = NaN;
            c = NaN(1, machines-1);
            q = NaN(1, machines);

            % get general instance data
            mu = instances_mu(k, :);
            p = instances_p(k, :);
            gamma = instances_gamma(k, :);
            costs_buffers = instances_costs_buffers(k, :);
            costs_spares = instances_costs_spares(k, :);

            % start values
            Q_start = Q_min;
            if allow_Q_min_adjustments
                [~, ~, iso_start_spares] = get_availability_characteristics(machines, Q_max, mu, p, gamma, TP_target);
                if (any(iso_start_spares > Q_start))
                    Q_start = iso_start_spares;
                end
            end

            % calculate
            counter_iterations = 0;
            counter_evaluations = 0;
            counter_problems = 0;
            time = 0;
            for base_C_min_local = base_C_min:10
                C_start = base_C_min_local * ones(1, machines-1);

                time_load = 0;
                filename_to_load = '';
                if TP_target_number > 1 && exist([dir_to_load, filename_local, '_C', num2str(base_C_min_local), '.mat'], 'file') == 2
                    filename_to_load = [dir_to_load, filename_local, '_C', num2str(base_C_min_local), '_int', '.mat'];
                    if exist(filename_to_load, 'file') == 2
                        load_from_other = load_data_opt([dir_to_load, filename_local, '_C', num2str(base_C_min_local), '.mat']);
                        time_load = load_from_other(3);
                    else
                        filename_to_load = '';
                    end
                end

                tic();
                [c(1, 1:machines-1), q(1, 1:machines), TP, costs, counter_iterations_local, counter_evaluations_local, counter_problems_local] = ...
                    optimize_greedy_inc(machines, C_min, Q_min, C_max, Q_max, mu, p, gamma, TP_target, costs_buffers, costs_spares, C_start, Q_start, dir_local, [filename_local, '_C', num2str(base_C_min_local)], check_buffers, check_spares, allow_Q_min_adjustments, filename_to_load);
                time_local = toc();
                counter_iterations = counter_iterations + counter_iterations_local;
                counter_evaluations = counter_evaluations + counter_evaluations_local;
                counter_problems = counter_problems + counter_problems_local;
                time = time + time_local + time_load;
                data = [TP costs time c(1, :) q(1, :) counter_iterations counter_evaluations counter_problems];
                save_data_opt([dir_local, filename_local, '_C', num2str(base_C_min_local), '.mat'], data);
                if TP >= TP_target, break; end
                if ~check_buffers, break; end
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

end