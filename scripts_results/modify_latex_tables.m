function modify_latex_tables(ENV)
rng(20);


%% parameters
if nargin == 0
    ENV = get_environment();
end


%% settings
TP_targets = {'0.60', '0.65', '0.70', '0.75', '0.80', '0.85', '0.90', '0.95'};
files = {'validation_unbalanced_2', 'insights_a_separate_3', ...
    'insights_b_costs_all_spares_5', 'insights_b_costs_all_spares_5_allocation', ...
    'insights_c_costs_one_spare_5', 'insights_c_costs_one_spare_5_allocation', ...
    'insights_d_bottlenecks_5', 'insights_d_bottlenecks_5_allocation', ...
    'validation_unbalanced_2', 'validation_unbalanced_2_best', 'validation_unbalanced_2_better', ...
    'validation_unbalanced_5', 'validation_unbalanced_5_best', 'validation_unbalanced_5_better'};


%% modify
for i = 1:length(files)
    file = files{i};
    line3 = '';
    switch file
        case {'validation_unbalanced_2'}
            line1 = '\begin{tabular}{LlRRRRRRRRRR}';
            line3 = ' &  & \multicolumn{1}{l}{\# infeasible} & \multicolumn{1}{l}{\# optimal} & \multicolumn{4}{l}{Computational time (s)} & \multicolumn{4}{l}{$\Delta_\text{costs}$ compared to ENUM} \\ ';
            line4 = 'TP^T & Algorithm & \multicolumn{1}{l}{solutions} & \multicolumn{1}{l}{solutions} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} \\ ';
        case {'validation_unbalanced_5'}
            line1 = '\begin{tabular}{LlRRRRRRRRRR}';
            line3 = ' &  & \multicolumn{1}{l}{\# infeasible} & \multicolumn{1}{l}{\# best} & \multicolumn{4}{l}{Computational time (s)} & \multicolumn{4}{l}{$\Delta_\text{costs}$ compared to GD} \\ ';
            line4 = 'TP^T & Algorithm & \multicolumn{1}{l}{solutions} & \multicolumn{1}{l}{solutions} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} \\ ';
        case 'validation_unbalanced_2_best'
            line1 = '\begin{tabular}{LlRRRRRRRR}';
            line3 = ' &  &  & \multicolumn{7}{l}{Proportion of also found optimal solutions} \\ ';
            line4 = 'TP^T & Algorithm & \multicolumn{1}{l}{\# optimal solutions} & \multicolumn{1}{l}{ENUM} & \multicolumn{1}{l}{GD} & \multicolumn{1}{l}{GI} & \multicolumn{1}{l}{GID} & \multicolumn{1}{l}{SAD} & \multicolumn{1}{l}{SAI} & \multicolumn{1}{l}{GA} \\ ';
        case {'validation_unbalanced_5_best'}
            line1 = '\begin{tabular}{LlRRRRRR}';
            line3 = ' &  &  & \multicolumn{5}{l}{Proportion of also found best solutions} \\ ';
            line4 = 'TP^T & Algorithm & \multicolumn{1}{l}{\# best solutions} & \multicolumn{1}{l}{GD} & \multicolumn{1}{l}{GI} & \multicolumn{1}{l}{GID} & \multicolumn{1}{l}{SAD} & \multicolumn{1}{l}{SAI} \\ ';
        case 'validation_unbalanced_2_better'
            line1 = '\begin{tabular}{LlRRRRRRRR}';
            line4 = 'TP^T & Algorithm & \multicolumn{1}{l}{\# solvable instances} & \multicolumn{1}{l}{ENUM} & \multicolumn{1}{l}{GD} & \multicolumn{1}{l}{GI} & \multicolumn{1}{l}{GID} & \multicolumn{1}{l}{SAD} & \multicolumn{1}{l}{SAI} & \multicolumn{1}{l}{GA} \\ ';
        case {'validation_unbalanced_5_better'}
            line1 = '\begin{tabular}{LlRRRRRR}';
            line4 = 'TP^T & Algorithm & \multicolumn{1}{l}{\# solvable instances} & \multicolumn{1}{l}{GD} & \multicolumn{1}{l}{GI} & \multicolumn{1}{l}{GID} & \multicolumn{1}{l}{SAD} & \multicolumn{1}{l}{SAI} \\ ';
        case {'insights_a_separate_3'}
            line1 = '\begin{tabular}{LlRRRRRRRRR}';
            line4 = 'TP^T & Case & \multicolumn{1}{l}{\# solvable} & \multicolumn{1}{l}{Mean TBC} & \multicolumn{1}{l}{Mean TCU} & \multicolumn{1}{l}{Mean TC} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} \\ ';
        case {'insights_b_costs_all_spares_5', 'insights_c_costs_one_spare_5'}
            line1 = '\begin{tabular}{LlLRRRRRRRRR}';
            line4 = 'TP^T & Case & \multicolumn{1}{l}{number solvable} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} \\ ';
        case 'insights_d_bottlenecks_5'
            line1 = '\begin{tabular}{LLLRRRRRRRRR}';
            line4 = 'TP^T & \multicolumn{1}{l}{Bottleneck} & \multicolumn{1}{l}{number solvable} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} & \multicolumn{1}{l}{Mean} & \multicolumn{1}{l}{Std. error} & \multicolumn{1}{l}{Min} & \multicolumn{1}{l}{Max} \\ ';
        case {'insights_b_costs_all_spares_5_allocation', 'insights_c_costs_one_spare_5_allocation'}
            line1 = '\begin{tabular}{LlRRRRRRRRR}';
            line4 = 'TP^T & \multicolumn{1}{l}{Case} & \multicolumn{1}{l}{$C_1$} & \multicolumn{1}{l}{$C_2$} & \multicolumn{1}{l}{$C_3$} & \multicolumn{1}{l}{$C_4$} & \multicolumn{1}{l}{$Q_1$} & \multicolumn{1}{l}{$Q_2$} & \multicolumn{1}{l}{$Q_3$} & \multicolumn{1}{l}{$Q_4$} & \multicolumn{1}{l}{$Q_5$} \\ ';
        case 'insights_d_bottlenecks_5_allocation'
            line1 = '\begin{tabular}{LLRRRRRRRRR}';
            line4 = 'TP^T & \multicolumn{1}{l}{Bottleneck} & \multicolumn{1}{l}{$C_1$} & \multicolumn{1}{l}{$C_2$} & \multicolumn{1}{l}{$C_3$} & \multicolumn{1}{l}{$C_4$} & \multicolumn{1}{l}{$Q_1$} & \multicolumn{1}{l}{$Q_2$} & \multicolumn{1}{l}{$Q_3$} & \multicolumn{1}{l}{$Q_4$} & \multicolumn{1}{l}{$Q_5$} \\ ';
        otherwise
            continue
    end

    if ~exist([ENV.OUTPUT_DIR_LATEX, file, '.tex'], 'file')
        continue
    end

    latex_file = readlines([ENV.OUTPUT_DIR_LATEX, file, '.tex']);
    latex_file{1} = line1;
    latex_file{4} = line4;
    tmp = latex_file{2};
    latex_file{2} = latex_file{3};
    latex_file{3} = tmp;
    if (line3 ~= ""), latex_file{3} = line3; end
    
    first_occurence_found = false;
    for k = 1:length(latex_file)
        for j = 1:length(TP_targets)
            TP_target = TP_targets{j};
            search_expression = ['^', TP_target, ' &'];
            if ~first_occurence_found
                index = regexp(latex_file{k}, search_expression, 'once');
                if ~isempty(index)
                    first_occurence_found = true;
                end
            else
                latex_file{k} = regexprep(latex_file{k}, search_expression, ['\\hline \n ', TP_target, ' &']);
            end
        end
    end

    [fid, msg] = fopen([ENV.OUTPUT_DIR_LATEX_MOD, file, '.tex'], 'w');
    if fid < 1
        error('Could not write output file: %s', msg);
    end
    fwrite(fid, strjoin(latex_file, '\n'));
    fclose(fid);
end

end