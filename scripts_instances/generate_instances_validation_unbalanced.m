function generate_instances_validation_unbalanced(ENV)
rng(42);


%% parameters
if nargin == 0
    ENV = get_environment();
end


%% settings
filename = 'validation_unbalanced';
number_instances = 500;

deviation_mu_min = 0;
deviation_mu_max = 0.05;
deviation_mu_type = 'both';
deviation_p_min = 0;
deviation_p_max = 0.1;
deviation_p_type = 'both';
deviation_gamma_min = 0;
deviation_gamma_max = 0.1;
deviation_gamma_type = 'both';

base_mu = 1;
base_p = 0.005;
base_gamma = 0.05;

base_costs_buffers = 1;
base_costs_spares = 1;


%% create instances
number_instances_total = number_instances;
for machines = [2, 5]
    filename_spec = [filename, '_', num2str(machines)];
    delete_if_exists([ENV.INPUT_DIR, filename_spec, '.xlsx'])
    delete_if_exists([ENV.INPUT_DIR, filename_spec, '.mat'])

    switch machines
        case 2
            base_C_min = 1;
            base_Q_min = 1;
            base_C_max = 100;
            base_Q_max = 5;
        case 3
            base_C_min = 1;
            base_Q_min = 1;
            base_C_max = 50;
            base_Q_max = 5;
        case 5
            base_C_min = 1;
            base_Q_min = 1;
            base_C_max = 25;
            base_Q_max = 5;
        otherwise
            error('This case is not described.');
    end

    C = base_C_max * ones(1, machines-1);
    Q = base_Q_max * ones(1, machines);


    %% fill instances
    instances_costs_buffers = base_costs_buffers * ones(number_instances_total, machines-1);
    instances_costs_spares = base_costs_spares * ones(number_instances_total, machines);

    max_TP = NaN(number_instances_total, 1);
    time = NaN(number_instances_total, 1);
    iterations = NaN(number_instances_total, 1);
    epsilon = NaN(number_instances_total, 1);
    terminated_normally = NaN(number_instances_total, 1);


    %% random part

    % mu
    random_numbers = rand(number_instances, machines);
    switch deviation_mu_type
        case 'both'
            dev_rand = (random_numbers - .5) * 2;
        case 'worse'
            dev_rand = (random_numbers - 1) * 1;
        case 'better'
            dev_rand = (random_numbers - 0) * 1;
        otherwise
            error('Unvalid value of deviation_mu_type detected.');
    end
    temp_mu = base_mu * ones(1, machines);
    dev_rand = dev_rand * (deviation_mu_max - deviation_mu_min);
    dev_rand(dev_rand < 0) = dev_rand(dev_rand < 0) - deviation_mu_min;
    dev_rand(dev_rand > 0) = dev_rand(dev_rand > 0) + deviation_mu_min;
    dev_rand_mu = (1 + dev_rand) .* base_mu;

    % p
    random_numbers = rand(number_instances, machines);
    switch deviation_p_type
        case 'both'
            dev_rand = (random_numbers - .5) * 2;
        case 'worse'
            dev_rand = (random_numbers - 0) * 1;
        case 'better'
            dev_rand = (random_numbers - 1) * 1;
        otherwise
            error('Unvalid value of deviation_p_type detected.');
    end
    temp_p = base_p * ones(1, machines);
    dev_rand = dev_rand * (deviation_p_max - deviation_p_min);
    dev_rand(dev_rand < 0) = dev_rand(dev_rand < 0) - deviation_p_min;
    dev_rand(dev_rand > 0) = dev_rand(dev_rand > 0) + deviation_p_min;
    dev_rand_p = (1 + dev_rand) .* base_p;

    % gamma
    random_numbers = rand(number_instances, machines);
    switch deviation_gamma_type
        case 'both'
            dev_rand = (random_numbers - .5) * 2;
        case 'worse'
            dev_rand = (random_numbers - 1) * 1;
        case 'better'
            dev_rand = (random_numbers - 0) * 1;
        otherwise
            error('Unvalid value of deviation_gamma_type detected.');
    end
    temp_gamma = base_gamma * ones(1, machines);
    dev_rand = dev_rand * (deviation_gamma_max - deviation_gamma_min);
    dev_rand(dev_rand < 0) = dev_rand(dev_rand < 0) - deviation_gamma_min;
    dev_rand(dev_rand > 0) = dev_rand(dev_rand > 0) + deviation_gamma_min;
    dev_rand_gamma = (1 + dev_rand) .* base_gamma;


    %% combine everything
    instances_mu = dev_rand_mu .* ones(number_instances_total, machines);
    instances_p = dev_rand_p .* ones(number_instances_total, machines);
    instances_gamma = dev_rand_gamma .* ones(number_instances_total, machines);


    %% solve random instances
    parfor k = 1:number_instances_total
        switch machines
            case 2
                tic();
                max_TP(k) = Spare_2M(C, Q(1), Q(2), instances_mu(k, 1), instances_mu(k, 2), instances_p(k, 1), instances_p(k, 2), instances_gamma(k, 1), instances_gamma(k, 2));
                time(k) = toc();
            case 3
                tic();
                [max_TP(k), iterations(k), epsilon(k), terminated_normally(k)] = Spare_Decomposition(machines, C, Q, instances_mu(k, :), instances_p(k, :), instances_gamma(k, :), true);
                time(k) = toc();
            case 5
                tic();
                [max_TP(k), iterations(k), epsilon(k), terminated_normally(k)] = Spare_Decomposition(machines, C, Q, instances_mu(k, :), instances_p(k, :), instances_gamma(k, :), true);
                time(k) = toc();
            otherwise
                error('This case is not described.');
        end
    end


    %% save table
    table_parameters_rows = [instances_costs_buffers, instances_costs_spares, instances_mu, instances_p, instances_gamma, ...
        max_TP, time, iterations, epsilon, terminated_normally];

    table_parameters_columns = [cellfun(@(c)['c_b_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false) ...
        cellfun(@(c)['c_s_', num2str(c)],num2cell(1:machines),'UniformOutput',false) ...
        cellfun(@(c)['mu_', num2str(c)],num2cell(1:machines),'UniformOutput',false) ...
        cellfun(@(c)['p_', num2str(c)],num2cell(1:machines),'UniformOutput',false) ...
        cellfun(@(c)['gamma_', num2str(c)],num2cell(1:machines),'UniformOutput',false) ...
        'TP', 'time', 'iterations', 'epsilon', 'terminated_normally'];
    table_parameters = array2table(table_parameters_rows, 'VariableNames', table_parameters_columns);
    writetable(table_parameters, [ENV.INPUT_DIR, filename_spec, '.xlsx'], 'Sheet', 'Parameters', 'Range', 'A1', 'WriteMode','overwritesheet');
    disp(table_parameters);


    %% save data
    save([ENV.INPUT_DIR, filename_spec, '.mat']);


end

end
