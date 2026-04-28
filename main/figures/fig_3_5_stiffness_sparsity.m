clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'common'));

cfg = setup_explanatory_figures();
data = build_explanatory_case(cfg, false);

fig = figure('visible', 'off');

spy(data.Kff);
%title('Sparsity pattern of the reduced stiffness matrix K_{FF}');
xlabel('column index');
ylabel('row index');

axis square;

filename = fullfile(cfg.outdir, 'fig_3_5_stiffness_sparsity');
save_thesis_figure(fig, filename);
close(fig);