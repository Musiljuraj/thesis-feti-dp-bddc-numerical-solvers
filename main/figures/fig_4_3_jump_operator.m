%Title: Jump operator on the product interface space
%Caption: Schematic meaning and sparsity pattern of the jump operator \(B\). 
% Each row of \(B\) represents one equality constraint between two duplicated product-space copies of the same 
% assembled interface degree of freedom.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'common'));

cfg = setup_explanatory_figures();

% Keep the same decomposition size as the previous Chapter 4 figures.
cfg.n = 12;
cfg.nSubX = 4;
cfg.nSubY = 4;

data = build_explanatory_case(cfg, true);

if isfield(data, 'problem')
    prob = data.problem;
else
    prob = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, cfg.f);
end

% Build the actual jump operator from the current implementation.
[B, meta] = build_jump_operator_B(prob.prod);

fig = figure('visible', 'off');
set(fig, 'position', [100, 100, 1200, 520]);

% Colors consistent with previous figures
colInterior  = [0.82, 0.82, 0.82];
colInterface = [0.50, 0.70, 0.92];
colEdge      = [0.15, 0.15, 0.15];
colText      = [0.10, 0.10, 0.10];

% ============================================================
% Left panel: schematic meaning of one jump constraint
% ============================================================
subplot(1, 2, 1);
hold on;

% Small separated subdomains
xA = [0.05, 0.40, 0.40, 0.05, 0.05];
yA = [0.20, 0.20, 0.80, 0.80, 0.20];

xB = [0.60, 0.95, 0.95, 0.60, 0.60];
yB = yA;

plot(xA, yA, '-', 'Color', colEdge, 'LineWidth', 1.1);
plot(xB, yB, '-', 'Color', colEdge, 'LineWidth', 1.1);

% Interface-copy nodes
p1 = [0.40, 0.50];
p2 = [0.60, 0.50];

plot(p1(1), p1(2), 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colInterface, ...
    'MarkerEdgeColor', colEdge, ...
    'LineWidth', 0.8);

plot(p2(1), p2(2), 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colInterface, ...
    'MarkerEdgeColor', colEdge, ...
    'LineWidth', 0.8);

% Dashed line indicating that the two copies represent the same physical DOF
plot([p1(1), p2(1)], [p1(2), p2(2)], '--', ...
    'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 1.0);

% Labels
text(0.225, 0.86, '\Omega_i', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, ...
    'Color', colText);

text(0.775, 0.86, '\Omega_j', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, ...
    'Color', colText);

text(p1(1)-0.02, p1(2)-0.10, 'w_i', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, ...
    'Color', colText);

text(p2(1)+0.02, p2(2)-0.10, 'w_j', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, ...
    'Color', colText);

text(0.50, 0.08, 'constraint:  w_i - w_j = 0', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, ...
    'Color', colText);

axis equal;
axis off;
xlim([0, 1]);
ylim([0, 1]);

% ============================================================
% Right panel: sparsity pattern of actual B
% ============================================================
subplot(1, 2, 2);

spy(B);

xlabel('product interface index');
ylabel('constraint index');

% Avoid title inside final thesis figure if preferred.
% The LaTeX caption should carry the explanation.

axis square;

% Center both subplots reasonably inside the canvas
set(gcf, 'Color', 'w');

filename = fullfile(cfg.outdir, 'fig_4_3_jump_operator');
save_thesis_figure(fig, filename);
close(fig);