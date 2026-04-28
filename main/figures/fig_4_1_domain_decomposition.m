% Title: Structured non-overlapping domain decomposition
% Caption: Structured non-overlapping decomposition of the computational domain into subdomains \(\Omega_i\). The underlying finite element mesh is shown together with the partition used for the domain decomposition framework.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'common'));

cfg = setup_explanatory_figures();

cfg.n = 12;
cfg.nSubX = 4;
cfg.nSubY = 4;

showSubLabels = false;

data = build_explanatory_case(cfg, true);

if isfield(data, 'problem')
    prob = data.problem;
else
    prob = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, cfg.f);
end

p   = prob.p;
t   = prob.t;
sub = prob.sub;
ddm = prob.ddm;

fig = figure('visible', 'off');
set(fig, 'position', [100, 100, 950, 950]);

hold on;

% --- smaller exploded layout ---
tileSize = 0.82;
gap      = 0.24;

% --- colors ---
colEdge     = [0.18, 0.18, 0.18];
colNodeFace = [0.82, 0.82, 0.82];
colNodeEdge = [0.15, 0.15, 0.15];
colLabel    = [0.20, 0.20, 0.20];

for s = 1:ddm.Nsub
    local_nodes = sub(s).nodes(:);

    tri_global = t(sub(s).elems, :);
    [tf, tri_local] = ismember(tri_global, local_nodes);
    if ~all(tf(:))
        error('fig_4_1_domain_decomposition: local triangle mapping failed for subdomain %d.', s);
    end

    bbox = sub(s).bbox;
    x0 = bbox(1); x1 = bbox(2);
    y0 = bbox(3); y1 = bbox(4);

    px = (p(local_nodes,1) - x0) / (x1 - x0);
    py = (p(local_nodes,2) - y0) / (y1 - y0);

    xShift = (sub(s).ix - 1) * (tileSize + gap);
    yShift = (cfg.nSubY - sub(s).iy) * (tileSize + gap);

    pDisp = [tileSize * px + xShift, ...
             tileSize * py + yShift];

    patch('Faces', tri_local, ...
          'Vertices', pDisp, ...
          'FaceColor', 'none', ...
          'EdgeColor', colEdge, ...
          'LineWidth', 0.75);

    plot(pDisp(:,1), pDisp(:,2), 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colNodeFace, ...
        'MarkerEdgeColor', colNodeEdge, ...
        'LineWidth', 0.75);

    if showSubLabels
        text(xShift + 0.5*tileSize, yShift + tileSize + 0.06, ...
            sprintf('\\Omega_{%d}', s), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 10, ...
            'Color', colLabel, ...
            'Interpreter', 'tex');
    end
end

axis equal;
axis off;

xmax = cfg.nSubX * tileSize + (cfg.nSubX - 1) * gap;
ymax = cfg.nSubY * tileSize + (cfg.nSubY - 1) * gap;

mx = 0.18;
my = 0.18;

xlim([-mx, xmax + mx]);
ylim([-my, ymax + my]);

ax = gca;
set(ax, 'units', 'normalized');
set(ax, 'position', [0.08, 0.08, 0.84, 0.84]);

filename = fullfile(cfg.outdir, 'fig_4_1_domain_decomposition');
save_thesis_figure(fig, filename);
close(fig);