function [C, Q, TP_current, costs_current, counter_iterations, counter_evaluations, counter_problems] = optimize_greedy_inc(machines, C_min, Q_min, C_max, Q_max, mu, p, gamma, TP_target, cost_buffers, cost_spares, C_start, Q_start, outputdir, filename, check_buffers, check_spares, allow_Q_min_adjustements, filename_to_load)
%optimize_greedy_inc uses greedy increase steps to reach target
%   Buffers or spare part inventories are increased to reach the target

if nargin == 0
    fprintf('No parameters given... using standard values...\n\n');
    addpath('import');
    machines = 5;
    TP_target = 0.85;
    C_min = 2 * ones(1, machines-1);
    Q_min = 1 * ones(1, machines);
    C_max = 100 * ones(1, machines-1);
    Q_max = 5 * ones(1, machines);
    mu = 1 * ones(1, machines);
    p = 0.005 * ones(1, machines);
    gamma = 0.05 * ones(1, machines);
    cost_buffers = 1 * ones(1, machines-1);
    cost_spares = 1 * ones(1, machines);
    C_start = 5 * ones(1, machines-1);
    Q_start = 1 * ones(1, machines);
    outputdir = '';
    filename = '';
    check_buffers = true;
    check_spares = true;
    allow_Q_min_adjustements = true;
    filename_to_load = '';
end


%% settings
max_iterations = 1000;


%% init
table_values_rows = NaN(1, 8*machines);

% get minimum for n and s
if allow_Q_min_adjustements
    [~, ~, Q_min_iso] = get_availability_characteristics(machines, Q_max, mu, p, gamma, TP_target);
    Q_min(Q_min < Q_min_iso) = Q_min_iso(Q_min < Q_min_iso);
end
min_vector = [C_min Q_min];
max_vector = [C_max Q_max];

% start values
C = C_start;
C(C_start < C_min) = C_min(C_start < C_min);
C(C_start > C_max) = C_max(C_start > C_max);

Q = Q_start;
Q(Q_start < Q_min) = Q_min(Q_start < Q_min);
Q(Q_start > Q_max) = Q_max(Q_start > Q_max);

% create output table columns
table_columns = [{'TP_old' 'TP' 'posN' 'posS'} cellfun(@(c)['n_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false) cellfun(@(c)['s_', num2str(c)],num2cell(1:machines),'UniformOutput',false) cellfun(@(c)['Delta_n_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false) cellfun(@(c)['Delta_s_', num2str(c)],num2cell(1:machines),'UniformOutput',false) cellfun(@(c)['TP_n_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false) cellfun(@(c)['TP_s_', num2str(c)],num2cell(1:machines),'UniformOutput',false) cellfun(@(c)['Inc_n_', num2str(c)],num2cell(1:(machines-1)),'UniformOutput',false) cellfun(@(c)['Inc_s_', num2str(c)],num2cell(1:machines),'UniformOutput',false)];


% check for feasibility of maximum values
if check_buffers
    C_max_start = C_max;
else
    C_max_start = C_start;
end
if check_spares
    Q_max_start = Q_max;
else
    Q_max_start = Q_start;
end

if (machines >= 3)
    [TP_current, ~, ~, dec_terminated_normally] = Spare_Decomposition(machines, C_max_start, Q_max_start, mu, p, gamma, true);
    counter_problems = ~dec_terminated_normally ;
elseif (machines == 3)
    TP_current = Spare_3M(C_max_start(1), C_max_start(2), Q_max_start(1), Q_max_start(2), Q_max_start(3), mu(1), mu(2), mu(3), p(1), p(2), p(3), gamma(1), gamma(2), gamma(3));
    counter_problems = 0;
else
    TP_current = Spare_2M(C_max_start, Q_max_start(1), Q_max_start(2), mu(1), mu(2), p(1), p(2), gamma(1), gamma(2));
    counter_problems = 0;
end
counter_evaluations = 1;
counter_iterations = 0;


%% increase greedy
if TP_current < TP_target
    costs_current = sum(C_max_start .* cost_buffers) + sum(Q_max_start .* cost_spares);
    C = C_max_start;
    Q = Q_max_start;
else
    % load if wanted
    if ~isempty(filename_to_load)
        load(filename_to_load, 'table_values_rows', 'C', 'Q', 'counter_iterations', 'counter_evaluations', 'counter_problems')
    end

    if (machines >= 3)
        [TP_current, ~, ~, dec_terminated_normally] = Spare_Decomposition(machines, C, Q, mu, p, gamma, true);
        counter_problems = counter_problems + ~dec_terminated_normally ;
    elseif (machines == 3)
        TP_current = Spare_3M(C(1), C(2), Q(1), Q(2), Q(3), mu(1), mu(2), mu(3), p(1), p(2), p(3), gamma(1), gamma(2), gamma(3));
    else
        TP_current = Spare_2M(C, Q(1), Q(2), mu(1), mu(2), p(1), p(2), gamma(1), gamma(2));
    end
    counter_evaluations = counter_evaluations + 1;
    costs_current = sum(C .* cost_buffers) + sum(Q .* cost_spares);
    while (TP_current < TP_target)
        TP = zeros(1, 2*machines-1);
        costs = zeros(1, 2*machines-1);
        counter_iterations = counter_iterations + 1;
        if (counter_iterations > max_iterations), break; end
        
        increase_steps = ones(1, 2*machines-1);
        
        % one buffer place or one spare part more
        %parfor i = 1:2*machines-1
        for i = 1:2*machines-1
            if (~check_buffers && i < machines)
                continue;
            end
            
            if (~check_spares && i >= machines)
                continue;
            end
            
            n_tmp = NaN;
            s_tmp = NaN;
            if (i < machines)
                increase_end = 5;
            else
                increase_end = 2;
            end
            
            dec_terminated_normally = true;
            for increase_step = 1:increase_end
                n_tmp = C;
                s_tmp = Q;
                
                % increase n or s
                increase_steps(i) = increase_step;
                if (i < machines)
                    n_tmp(i) = n_tmp(i) + increase_step;
                    if (n_tmp(i) > max_vector(i))
                        TP(i) = NaN;
                        costs(i) = NaN;
                        continue
                    end
                else 
                    s_tmp(i - machines + 1) = s_tmp(i - machines + 1) + increase_step;
                    if (s_tmp(i - machines + 1) > max_vector(i))
                        TP(i) = NaN;
                        costs(i) = NaN;
                        continue
                    end
                end
                
                if (machines >= 3)
                    [TP(i), ~, ~, dec_terminated_normally] = Spare_Decomposition(machines, n_tmp, s_tmp, mu, p, gamma, true);
                elseif (machines == 3)
                    TP(i) = Spare_3M(n_tmp(1), n_tmp(2), s_tmp(1), s_tmp(2), s_tmp(3), mu(1), mu(2), mu(3), p(1), p(2), p(3), gamma(1), gamma(2), gamma(3));
                else
                    TP(i)= Spare_2M(n_tmp, s_tmp(1), s_tmp(2), mu(1), mu(2), p(1), p(2), gamma(1), gamma(2));
                end
                counter_evaluations = counter_evaluations + 1;
                if (TP(i) > TP_current), break;
                else, TP(i) = NaN; end
            end
            counter_problems = counter_problems + ~dec_terminated_normally;
            
            costs(i) = sum(n_tmp .* cost_buffers) + sum(s_tmp .* cost_spares);
        end
        
        % increased buffer/spares but decreased TP
        TP(TP < TP_current) = NaN;
        
        % remove either buffers or spare parts from decision if wanted
        if (~check_buffers)
            TP(1:machines-1) = NaN;
        end
        if (~check_spares)
            TP(machines:2*machines-1) = NaN;
        end
        
        % compute differences and decision
        Delta_TP = TP - TP_current;
        Delta_costs = costs - costs_current;
        Delta = Delta_TP ./ Delta_costs;
        [val, pos] = max(Delta(:));
        
        % no feasible decisions left
        if isnan(val)
            table_values_rows(counter_iterations, 1:(8*machines)) = [TP_current NaN NaN NaN C Q Delta TP increase_steps];
            break;
        end
        
        % prepare for next iteration
        TP_old = TP_current;
        TP_current = TP(pos);
        costs_current = costs(pos);
        if (pos < machines)
            C(pos) = C(pos) + increase_steps(pos);
            posN = pos;
            posS = NaN;
        else
            Q(pos - machines + 1) = Q(pos - machines + 1) + increase_steps(pos);
            posN = NaN;
            posS = pos - machines + 1;
        end
        
        % output
        table_values_rows(counter_iterations, 1:(8*machines)) = [TP_old TP_current posN posS C Q Delta TP increase_steps];
        %table_values = array2table(table_values_rows, 'VariableNames', table_columns);
        %writetable(table_values, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
        %disp(tail(table_values, 1));
    end
end
table_values_rows_opt = table_values_rows(end, :);


%% output
table_values = array2table(table_values_rows, 'VariableNames', table_columns);
table_opt = array2table(table_values_rows_opt, 'VariableNames', table_columns);
if ~isempty(filename)
    writetable(table_values, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Results', 'Range', 'A1', 'WriteMode','overwritesheet');
    writetable(table_opt, [outputdir, '/', filename, '.xlsx'], 'Sheet', 'Best', 'Range', 'A1', 'WriteMode','overwritesheet');
    save([outputdir, '/', filename, '_int', '.mat'], 'table_values_rows', 'C', 'Q', 'counter_iterations', 'counter_evaluations', 'counter_problems')
else
    disp(table_values);
    disp(table_opt);
end

end