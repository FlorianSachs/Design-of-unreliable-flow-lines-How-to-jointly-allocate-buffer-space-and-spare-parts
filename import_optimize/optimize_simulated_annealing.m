function [C_best, Q_best, TP_best, costs_best, counter_iterations, counter_evaluations, counter_problems] = optimize_simulated_annealing(machines, C_min, Q_min, C_max, Q_max, mu, p, gamma, TP_target, cost_buffers, cost_spares, C_start, Q_start, outputdir, filename, alg)
%optimize_simulated_annealing uses simulated annealing to reach target
%   Buffers or spare part inventories are increased to reach the target

if nargin == 0
    fprintf('No parameters given... using standard values...\n\n');
    addpath('import');
    machines = 2;
    TP_target = 0.85;
    C_min = 1 * ones(1, machines-1);
    Q_min = 1 * ones(1, machines);
    C_max = 100 * ones(1, machines-1);
    Q_max = 5 * ones(1, machines);
    mu = 1 * ones(1, machines);
    p = 0.005 * ones(1, machines);
    gamma = 0.05 * ones(1, machines);
    cost_buffers = 1 * ones(1, machines-1);
    cost_spares = 1 * ones(1, machines);
    C_start = C_min;
    Q_start = Q_min;
    outputdir = '';
    filename = '';
end

if ~exist('alg', 'var')
    %alg.type = 'linear';
    %alg.start_temperature = 20;
    %alg.param = 0.1;
    
    %alg.type = 'geometric';
    %alg.start_temperature = 20;
    %alg.param = 0.95;
    
    alg.type = 'logarithmic';
    alg.start_temperature = 20;
    alg.param = NaN;
    
    alg.max_iterations = (2 * machines - 1) * 200;
    alg.T_min = 2;
end


%% settings
big_step = 5;
number_specials = 4;
penalty = 1000;


%% init
rng(42);
table_values_rows = NaN(1, 2*machines + 10);

% get minimum for n and s
[~, ~, Q_min_iso] = get_availability_characteristics(machines, Q_max, mu, p, gamma, TP_target);
Q_min(Q_min < Q_min_iso) = Q_min_iso(Q_min < Q_min_iso);
%min_vector = [C_min Q_min];
%max_vector = [C_max Q_max];

% start values
C = C_start;
C(C_start < C_min) = C_min(C_start < C_min);
C(C_start > C_max) = C_max(C_start > C_max);

Q = Q_start;
Q(Q_start < Q_min) = Q_min(Q_start < Q_min);
Q(Q_start > Q_max) = Q_max(Q_start > Q_max);

% create output table columns
table_columns = [{'TP_old' 'TP' 'costs_old' 'costs' 'posN' 'posS' 'stepN' 'stepS' 'special' 'prob' 'random'} cellfun(@(c)['C_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false) cellfun(@(c)['Q_', num2str(c)],num2cell(1:machines),'UniformOutput',false)];


%% simulated annealing
if (machines >= 3)
    [TP_current, ~, ~, dec_terminated_normally] = Spare_Decomposition(machines, C, Q, mu, p, gamma, true);
    counter_problems = ~dec_terminated_normally ;
elseif (machines == 3)
    TP_current = Spare_3M(C(1), C(2), Q(1), Q(2), Q(3), mu(1), mu(2), mu(3), p(1), p(2), p(3), gamma(1), gamma(2), gamma(3));
    counter_problems = 0;
else
    TP_current = Spare_2M(C, Q(1), Q(2), mu(1), mu(2), p(1), p(2), gamma(1), gamma(2));
    counter_problems = 0;
end
counter_evaluations = 1;
counter_iterations = 0;
costs_current = sum(C .* cost_buffers) + sum(Q .* cost_spares);
if TP_current < TP_target
    costs_current = costs_current + penalty * ((TP_target - TP_current) / TP_target);
end
table_values_rows(counter_iterations + 1, 1:(2*machines + 10)) = [NaN TP_current NaN costs_current NaN NaN NaN NaN NaN NaN NaN C Q];

C_best = C;
Q_best = Q;
TP_best = TP_current;
costs_best = costs_current;

while true
    counter_iterations = counter_iterations + 1;
    
    switch alg.type
        case 'linear'
            T = alg.start_temperature - counter_iterations * alg.param;
        case 'geometric'
            T = alg.start_temperature - alg.param^counter_iterations;
        case 'logarithmic'
            T = (alg.start_temperature * log(2)) / log(1 + counter_iterations);
    end
    
    if (T <= alg.T_min || counter_iterations > alg.max_iterations), break; end
    
    % get next candidate
    success = false;
    while ~success
        random = randi(4*(machines - 1) + 2*machines + number_specials,1,1);
        n_tmp = C;
        s_tmp = Q;
        posN = NaN;
        posS = NaN;
        stepN = NaN;
        stepS = NaN;
        special = NaN;
        if (random <= machines - 1)
            % buffer_increase
            index_rand = random;
            posN = index_rand;
            stepN = 1;
            if (C(index_rand) + 1 <= C_max(index_rand))
                n_tmp(index_rand) = n_tmp(index_rand) + 1;
                success = true;
            end
        elseif (random <= 2*(machines - 1))
            % buffer_decrease
            index_rand = random - (machines - 1);
            posN = index_rand;
            stepN = -1;
            if (C(index_rand) - 1 >= C_min(index_rand))
                n_tmp(index_rand) = n_tmp(index_rand) - 1;
                success = true;
            end
        elseif (random <= 3*(machines - 1))
            % buffer_increase_big
            index_rand = random - 2*(machines - 1);
            posN = index_rand;
            stepN = big_step;
            if (C(index_rand) + big_step <= C_max(index_rand))
                n_tmp(index_rand) = n_tmp(index_rand) + big_step;
                success = true;
            end
        elseif (random <= 4*(machines - 1))
            % buffer_decrease_big
            index_rand = random - 3*(machines - 1);
            posN = index_rand;
            stepN = -big_step;
            if (C(index_rand) - big_step >= C_min(index_rand))
                n_tmp(index_rand) = n_tmp(index_rand) - big_step;
                success = true;
            end
        elseif (random <= 5*(machines - 1) + 1)
            % spare_increase
            index_rand = random - 4*(machines - 1);
            posS = index_rand;
            stepS = 1;
            if (Q(index_rand) + 1 <= Q_max(index_rand))
                s_tmp(index_rand) = s_tmp(index_rand) + 1;
                success = true;
            end
        elseif (random <= 6*(machines - 1) + 2)
            % spare_decrease
            index_rand = random - (5*(machines - 1) + 1);
            posS = index_rand;
            stepS = -1;
            if (Q(index_rand) - 1 >= Q_min(index_rand))
                s_tmp(index_rand) = s_tmp(index_rand) - 1;
                success = true;
            end
        else
            %different things
            index_rand = random - (6*(machines - 1) + 2);
            special = index_rand;
            switch index_rand
                case 1
                    if any(C + big_step <= C_max)
                        n_tmp(C + big_step <= C_max) = n_tmp(C + big_step <= C_max) + big_step;
                        success = true;
                    end
                case 2
                    if any(C - big_step >= C_min)
                        n_tmp(C - big_step >= C_min) = n_tmp(C - big_step >= C_min) - big_step;
                        success = true;
                    end
                case 3
                    if any(Q + 1 <= Q_max)
                        s_tmp(Q + 1 <= Q_max) = s_tmp(Q + 1 <= Q_max) + 1;
                        success = true;
                    end
                case 4
                    if any(Q - 1 >= Q_min)
                        s_tmp(Q - 1 >= Q_min) = s_tmp(Q - 1 >= Q_min) - 1;
                        success = true;
                    end
                otherwise
                    error('This should not happen. This special neighbor does not exist.');
            end
        end
        
        if (success)
            if (machines >= 3)
                [TP_new , ~, ~, dec_terminated_normally] = Spare_Decomposition(machines, n_tmp, s_tmp, mu, p, gamma, true);
                counter_problems = counter_problems + ~dec_terminated_normally ;
            elseif (machines == 3)
                TP_new = Spare_3M(n_tmp(1), n_tmp(2), s_tmp(1), s_tmp(2), s_tmp(3), mu(1), mu(2), mu(3), p(1), p(2), p(3), gamma(1), gamma(2), gamma(3));
            else
                TP_new = Spare_2M(n_tmp, s_tmp(1), s_tmp(2), mu(1), mu(2), p(1), p(2), gamma(1), gamma(2));
            end
            counter_evaluations = counter_evaluations + 1;
            costs_new = sum(n_tmp .* cost_buffers) + sum(s_tmp .* cost_spares);
            if TP_new < TP_target
                costs_new  = costs_new + penalty * ((TP_target - TP_new) / TP_target);
            end
        end
    end
    
    % check outcome
    TP_old = TP_current;
    costs_old = costs_current;
    Delta_costs = costs_new - costs_current;
    %Delta_TP = TP_new - TP_current;
    
    if Delta_costs <= 0
        C = n_tmp;
        Q = s_tmp;
        costs_current = costs_new;
        TP_current = TP_new;
        prob = NaN;
        random_event = NaN;
    else
        prob = min(1, exp(-Delta_costs / T));
        random_event = binornd(1, prob);
        if (random_event == 1)
            C = n_tmp;
            Q = s_tmp;
            costs_current = costs_new;
            TP_current = TP_new;
        end
    end
    
    if TP_best < TP_target && TP_current > TP_best || ...
            costs_current < costs_best && TP_current >= TP_target || ...
            (costs_current <= costs_best && TP_current > TP_best)
        C_best = C;
        Q_best = Q;
        TP_best = TP_current;
        costs_best = costs_current;
    end

    % output
    table_values_rows(counter_iterations+1, 1:(2*machines + 10)) = [TP_old TP_current costs_old costs_current posN posS stepN stepS special prob random_event C Q];
    %table_values = array2table(table_values_rows, 'VariableNames', table_columns);
    %writetable(table_values, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
    %disp(tail(table_values, 1));
end
counter_iterations = counter_iterations - 1;
table_values_rows_opt(1, :) = table_values_rows(end, :);
table_values_rows_opt(1, 1:(2*machines + 10)) = [NaN TP_best NaN costs_best NaN NaN NaN NaN NaN NaN NaN C_best Q_best];


%% output
table_values = array2table(table_values_rows, 'VariableNames', table_columns);
if ~isempty(filename)
    writetable(table_values, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
else
    disp(table_values);
end

table_opt = array2table(table_values_rows_opt, 'VariableNames', table_columns);
if ~isempty(filename)
    writetable(table_opt, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Best', 'Range', 'A1', 'WriteMode','overwritesheet');
else
    disp(table_opt);
end


end