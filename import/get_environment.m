function [ENV] = get_environment()

% output directory
tmp = what('output');
if isempty(tmp), tmp = what('../output'); end
ENV.OUTPUT_DIR = [tmp.path, '/'];
if isempty(tmp), error(['The output directory "output/" ', ENV.OUTPUT_DIR, ' does not exist.']); end


% output subdirectories
tmp = what([ENV.OUTPUT_DIR, 'data']);
ENV.OUTPUT_DIR_DATA = [tmp.path, '/'];
if isempty(tmp), error(['The output directory "data" ', ENV.OUTPUT_DIR_DATA, ' does not exist.']); end

tmp = what([ENV.OUTPUT_DIR, 'excel']);
ENV.OUTPUT_DIR_EXCEL = [tmp.path, '/'];
if isempty(tmp), error(['The output directory "output/excel" ', ENV.OUTPUT_DIR_EXCEL, ' does not exist.']); end

tmp = what([ENV.OUTPUT_DIR, 'latex']);
ENV.OUTPUT_DIR_LATEX = [tmp.path, '/'];
if isempty(tmp),     error(['The output directory "output/latex" ', ENV.OUTPUT_DIR_LATEX, ' does not exist.']); end

tmp = what([ENV.OUTPUT_DIR, 'latex_mod']);
ENV.OUTPUT_DIR_LATEX_MOD = [tmp.path, '/'];
if isempty(tmp),  error(['The output directory "output/latex_mod" ', ENV.OUTPUT_DIR_LATEX_MOD, ' does not exist.']); end

tmp = what([ENV.OUTPUT_DIR, 'figures']);
ENV.OUTPUT_DIR_FIGURES = [tmp.path, '/'];
if isempty(tmp), error(['The output directory "output/figures" ', ENV.OUTPUT_DIR_FIGURES, ' does not exist.']); end


% final output directory
tmp = what('../paper');
if isempty(tmp), tmp = what('paper'); end
ENV.PAPER_DIR = [tmp.path, '/'];
if isempty(tmp), error(['The paper directory ', ENV.PAPER_DIR, ' does not exist.']); end


% input directory
tmp = what('../input');
if isempty(tmp), tmp = what('input'); end
ENV.INPUT_DIR = [tmp.path, '/'];
if isempty(tmp),  error(['The input directory ', ENV.INPUT_DIR, ' does not exist.']); end


% input directory
tmp = what('../rawdata');
if isempty(tmp), tmp = what('rawdata'); end
ENV.DATA_DIR = [tmp.path, '/'];
if isempty(tmp),  error(['The data directory ', ENV.DATA_DIR, ' does not exist.']); end


end