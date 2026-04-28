% Finite element approximation of the membrane deflection
% Finite element approximation \(u_h\) of the membrane displacement for the model Poisson problem on the unit square. The left and right boundary parts are subject to homogeneous Dirichlet conditions, while the remaining boundary parts satisfy homogeneous Neumann conditions.

% \begin{figure}[htbp]
%   \centering
%   \includegraphics[width=0.75\textwidth]{output/figures/explanatory/fig_3_3_membrane_surface.pdf}
%   \caption{Finite element approximation \(u_h\) of the membrane displacement for the model Poisson problem on the unit square. The left and right boundary parts are subject to homogeneous Dirichlet conditions, while the remaining boundary parts satisfy homogeneous Neumann conditions.}
%   \label{fig:membrane-surface}
% \end{figure}


clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'common'));

cfg = setup_explanatory_figures();

% Use a slightly finer mesh than the default conceptual figures,
% because the membrane surface should look smooth in the thesis.
cfg.n = 32;

% Constant unit load for the model Poisson problem.
% For -Delta u = f, this produces a positive displacement according
% to the chosen sign convention.
cfg.f = @(x,y) 1 + 0*x;

data = build_explanatory_case(cfg, false);

% Solve the reduced FEM system.
u_free = data.Kff \ data.ff;

% Reconstruct the full nodal vector, with homogeneous Dirichlet nodes equal to zero.
u = zeros(size(data.p, 1), 1);
u(data.free) = full(u_free);

fig = figure('visible', 'off');
set(fig, 'position', [100, 100, 900, 650]);

hSurf = trisurf(data.t, data.p(:,1), data.p(:,2), u);

% Clean surface appearance.
% If you want visible finite-element triangles, replace 'none' by 'k'.
set(hSurf, 'edgecolor', 'none');

shading interp;
grid on;
axis tight;
axis vis3d;

xlabel('x');
ylabel('y');
zlabel('u_h');

%title(sprintf('Finite element membrane deflection, n = %d', cfg.n));

colorbar;
view(45, 30);

filename = fullfile(cfg.outdir, 'fig_3_3_membrane_surface');
save_thesis_figure(fig, filename);
close(fig);