clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'common'));

cfg = setup_explanatory_figures();

% Finer mesh for a smooth membrane surface.
cfg.n = 24;

data = build_explanatory_case(cfg, false);

p = data.p;
t = data.t;

% Solve the reduced FEM system and reconstruct the full nodal vector.
u_free = data.Kff \ data.ff;

u = zeros(size(p, 1), 1);
u(data.free) = full(u_free);

% Plot membrane as downward deflection.
% Geometry uses z = -u, color uses positive u.
z = -u(:);

% --- thesis-friendly colors ---
colArrow  = [0.18, 0.18, 0.18];   % dark gray / almost black
colText   = [0.15, 0.15, 0.15];

fig = figure('visible', 'off');
set(fig, 'position', [100, 100, 1000, 760]);

hold on;

% --- main membrane surface ---
trisurf(t, p(:,1), p(:,2), z, u, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp');

shading interp;

% Slightly brighter pastel colormap
cmap = jet(256);
cmap = 0.70 * cmap + 0.30 * ones(size(cmap));
colormap(cmap);

cb = colorbar;
ylabel(cb, 'displacement magnitude');

% --- downward force arrows ---
% Manual arrows are used instead of quiver3 so the arrowheads are clearly visible.
xq = [0.22, 0.50, 0.78, 0.35, 0.65];
yq = [0.22, 0.25, 0.22, 0.72, 0.72];

z_top = 0.030;
z_bot = -0.035;

head_len = 0.018;
head_rad = 0.025;

for k = 1:length(xq)
    x = xq(k);
    y = yq(k);

    % Main arrow shaft
    plot3([x x], [y y], [z_top z_bot], '-', ...
        'Color', colArrow, ...
        'LineWidth', 2.2);

    % Arrowhead at the bottom
    plot3([x x + head_rad], [y y], [z_bot z_bot + head_len], '-', ...
        'Color', colArrow, ...
        'LineWidth', 2.2);

    plot3([x x - head_rad], [y y], [z_bot z_bot + head_len], '-', ...
        'Color', colArrow, ...
        'LineWidth', 2.2);

    plot3([x x], [y y + head_rad], [z_bot z_bot + head_len], '-', ...
        'Color', colArrow, ...
        'LineWidth', 2.2);

    plot3([x x], [y y - head_rad], [z_bot z_bot + head_len], '-', ...
        'Color', colArrow, ...
        'LineWidth', 2.2);
end

% Force label
text(0.52, 0.52, 0.040, 'f', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'Color', colArrow, ...
    'HorizontalAlignment', 'center');

xlabel('x');
ylabel('y');
zlabel('u_h');

box on;
grid on;

% Camera angle
view(35, 28);

% Less visually exaggerated vertical scale.
pbaspect([1 1 0.60]);

% Keep a little space above the arrows and below the membrane.
zlim([min(z) - 0.008, 0.055]);

% Keep the surface centered.
axis tight;

filename = fullfile(cfg.outdir, 'fig_2_1_model_problem');
save_thesis_figure(fig, filename);
close(fig);