%% clear everything
clear variables


%% create important directories
[~,~] = mkdir('rawdata');
[~,~] = mkdir('import');
[~,~] = mkdir('import_general');
[~,~] = mkdir('import_optimize');
[~,~] = mkdir('input');
[~,~] = mkdir('output');
[~,~] = mkdir('scripts_instances');
[~,~] = mkdir('scripts_results');
[~,~] = mkdir('output/data');
[~,~] = mkdir('output/excel');
[~,~] = mkdir('output/latex');
[~,~] = mkdir('output/latex_mod');
[~,~] = mkdir('output/figures');
[~,~] = mkdir('paper');


%% add amendments to path
addpath('import');
addpath('import_general');
addpath('import_call');
addpath('import_optimize');
addpath('scripts_instances');
addpath('scripts_results');
addpath('.');


%% python
if count(py.sys.path, [pwd, '\\python_GA']) == 0
    insert(py.sys.path, int32(0), [pwd, '\\python_GA']);
end


%% format
format shortG
