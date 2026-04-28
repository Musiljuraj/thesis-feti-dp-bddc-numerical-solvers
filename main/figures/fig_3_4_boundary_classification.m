%Title: Boundary classification on the finite element mesh
%Caption: Classification of the mesh nodes and boundary parts for the model Poisson problem on the unit square. The left and right sides form the Dirichlet boundary \(\Gamma_D\), while the top and bottom sides form the Neumann boundary \(\Gamma_N\). Interior nodes and Neumann boundary nodes remain active in the reduced discrete system, whereas the Dirichlet nodes are eliminated.

% \begin{figure}[htbp]
%   \centering
%   \includegraphics[width=0.78\textwidth]{output/figures/explanatory/fig_3_4_boundary_classification.pdf}
%   \caption{Classification of the mesh nodes and boundary parts for the model Poisson problem on the unit square. The left and right sides form the Dirichlet boundary \(\Gamma_D\), while the top and bottom sides form the Neumann boundary \(\Gamma_N\). Interior nodes and Neumann boundary nodes remain active in the reduced discrete system, whereas the Dirichlet nodes are eliminated.}
%   \label{fig:boundary-classification}
% \end{figure}

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'common'));

cfg = setup_explanatory_figures();

% Moderate mesh density for a clear classification figure.
cfg.n = 12;

data = build_explanatory_case(cfg, false);

p = data.p;
t = data.t;

nNodes = size(p, 1);

% Dirichlet nodes from boundary data.
dirichlet_nodes = unique(data.dirichlet(:));

% Identify all boundary nodes geometrically.
tol = 1e-12;
isBoundary = abs(p(:,1)) < tol | abs(p(:,1) - 1) < tol | ...
             abs(p(:,2)) < tol | abs(p(:,2) - 1) < tol;

% Dirichlet / Neumann / Interior classification.
isDirichlet = false(nNodes, 1);
isDirichlet(dirichlet_nodes) = true;

isNeumann  = isBoundary & ~isDirichlet;
isInterior = ~isBoundary;

% --- thesis-friendly colors ---
colMesh = [0.82, 0.82, 0.82];
colInt  = [0.55, 0.55, 0.55];
colNeu  = [0.50, 0.70, 0.92];   % brighter soft blue
colDir  = [0.92, 0.56, 0.56];   % brighter soft red

% Optional slightly darker colors for boundary labels
colNeuText = [0.35, 0.55, 0.82];
colDirText = [0.82, 0.42, 0.42];

fig = figure('visible', 'off');
set(fig, 'position', [100, 100, 900, 700]);

hold on;

% --- background mesh in gray only ---
triplot(t, p(:,1), p(:,2), '-', 'color', colMesh);

% --- nodes ---
hInterior = plot(p(isInterior,1), p(isInterior,2), 'o', ...
    'markersize', 4, ...
    'markerfacecolor', colInt, ...
    'markeredgecolor', colInt);

hNeumann = plot(p(isNeumann,1), p(isNeumann,2), 'o', ...
    'markersize', 6, ...
    'markerfacecolor', colNeu, ...
    'markeredgecolor', colNeu);

hDirichlet = plot(p(isDirichlet,1), p(isDirichlet,2), 's', ...
    'markersize', 6, ...
    'markerfacecolor', colDir, ...
    'markeredgecolor', colDir);

% --- boundary labels ---
text(-0.06, 0.50, '\Gamma_D', ...
    'fontsize', 12, ...
    'color', colDirText, ...
    'horizontalalignment', 'center');

text(1.06, 0.50, '\Gamma_D', ...
    'fontsize', 12, ...
    'color', colDirText, ...
    'horizontalalignment', 'center');

text(0.50, -0.07, '\Gamma_N', ...
    'fontsize', 12, ...
    'color', colNeuText, ...
    'horizontalalignment', 'center');

text(0.50, 1.07, '\Gamma_N', ...
    'fontsize', 12, ...
    'color', colNeuText, ...
    'horizontalalignment', 'center');

axis equal;
xlim([-0.10, 1.10]);
ylim([-0.10, 1.10]);
box on;

xlabel('x');
ylabel('y');

legend([hInterior, hNeumann, hDirichlet], ...
    {'interior nodes', 'Neumann boundary nodes', 'Dirichlet nodes'}, ...
    'location', 'northeastoutside');

filename = fullfile(cfg.outdir, 'fig_3_4_boundary_classification');
save_thesis_figure(fig, filename);
close(fig);