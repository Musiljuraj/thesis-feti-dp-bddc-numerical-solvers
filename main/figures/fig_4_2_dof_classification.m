%Title: Classification of local degrees of freedom in the subdomain decomposition
%Caption: Exploded view of the structured domain decomposition with classification of local degrees of freedom. Shared interface degrees of freedom are shown in blue, while the selected primal (corner) degrees of freedom are shown in red. Ordinary non-interface nodes are shown in gray.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'common'));

cfg = setup_explanatory_figures();

cfg.n = 12;
cfg.nSubX = 4;
cfg.nSubY = 4;

data = build_explanatory_case(cfg, true);

if isfield(data, 'problem')
    prob = data.problem;
else
    prob = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, cfg.f);
end

p      = prob.p;
t      = prob.t;
sub    = prob.sub;
ddm    = prob.ddm;
primal = prob.primal;

if isfield(primal, 'glob_c')
    primal_glob = unique(primal.glob_c(:));
else
    error('fig_4_2_dof_classification: primal.glob_c not found.');
end

fig = figure('visible', 'off');
set(fig, 'position', [100, 100, 1200, 900]);

hold on;

% --- exploded layout ---
tileSize = 0.82;
gap      = 0.24;

% --- colors ---
colEdge       = [0.18, 0.18, 0.18];
colInterior   = [0.82, 0.82, 0.82];
colInterface  = [0.50, 0.70, 0.92];   % pastel blue
colPrimal     = [0.92, 0.56, 0.56];   % pastel red
colDirichlet  = [0.93, 0.82, 0.60];   % pastel gold / sand
colNodeEdge   = [0.15, 0.15, 0.15];

% Legend handles
hInterior  = [];
hInterface = [];
hPrimal    = [];
hDirichlet = [];

tol = 1e-12;

for s = 1:ddm.Nsub
    local_nodes = sub(s).nodes(:);

    % Local triangles in local numbering
    tri_global = t(sub(s).elems, :);
    [tf, tri_local] = ismember(tri_global, local_nodes);
    if ~all(tf(:))
        error('fig_4_2_dof_classification: local triangle mapping failed for subdomain %d.', s);
    end

    % Local bounding box
    bbox = sub(s).bbox;
    x0 = bbox(1); x1 = bbox(2);
    y0 = bbox(3); y1 = bbox(4);

    % Normalize local coordinates to [0,1] x [0,1]
    px = (p(local_nodes,1) - x0) / (x1 - x0);
    py = (p(local_nodes,2) - y0) / (y1 - y0);

    % Exploded placement
    xShift = (sub(s).ix - 1) * (tileSize + gap);
    yShift = (cfg.nSubY - sub(s).iy) * (tileSize + gap);

    pDisp = [tileSize * px + xShift, ...
             tileSize * py + yShift];

    % Draw local mesh
    patch('Faces', tri_local, ...
          'Vertices', pDisp, ...
          'FaceColor', 'none', ...
          'EdgeColor', colEdge, ...
          'LineWidth', 0.75);

    % --- geometric location on local boundary ---
    onLeft   = abs(px) < tol;
    onRight  = abs(px - 1) < tol;
    onBottom = abs(py) < tol;
    onTop    = abs(py - 1) < tol;

    isCornerGeom = (onLeft | onRight) & (onTop | onBottom);

    % --- Dirichlet nodes: global left/right boundaries ---
    isDirichlet = (sub(s).ix == 1         & onLeft)  | ...
                  (sub(s).ix == cfg.nSubX & onRight);

    % --- primal / corner nodes from actual primal set ---
    local_dofs = ddm.node2dof(local_nodes);
    isFree = (local_dofs > 0);

    isPrimal = false(size(local_nodes));
    isPrimal(isFree) = ismember(local_dofs(isFree), primal_glob);

    % Dirichlet overrides primal coloring
    isPrimal = isPrimal & ~isDirichlet;

    % --- vertical internal interfaces (non-corner only) ---
    isVerticalInterface = ...
        ((sub(s).ix > 1)         & onLeft  & ~isCornerGeom) | ...
        ((sub(s).ix < cfg.nSubX) & onRight & ~isCornerGeom);

    % --- horizontal internal interfaces (non-corner only) ---
    % We explicitly control top/bottom according to subdomain row:
    %   top row    -> bottom edge is interface, top edge is global Neumann
    %   bottom row -> top edge is interface, bottom edge is global Neumann
    %   middle rows -> both top and bottom are interface
    isHorizontalInterface = false(size(local_nodes));

    if sub(s).iy == 1
        % top row of subdomains
        isHorizontalInterface = onBottom & ~isCornerGeom;
    elseif sub(s).iy == cfg.nSubY
        % bottom row of subdomains
        isHorizontalInterface = onTop & ~isCornerGeom;
    else
        % middle rows
        isHorizontalInterface = (onTop | onBottom) & ~isCornerGeom;
    end

    % --- total interface nodes ---
    isInterface = (isVerticalInterface | isHorizontalInterface);

    % Exclusions: primal and Dirichlet should not be recolored blue
    isInterface = isInterface & ~isPrimal & ~isDirichlet;

    % --- everything else is interior (including outer Neumann non-corner nodes) ---
    isInterior = ~(isDirichlet | isPrimal | isInterface);

    % Plot nodes in layers
    if any(isInterior)
        hInterior = plot(pDisp(isInterior,1), pDisp(isInterior,2), 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', colInterior, ...
            'MarkerEdgeColor', colNodeEdge, ...
            'LineWidth', 0.75);
    end

    if any(isInterface)
        hInterface = plot(pDisp(isInterface,1), pDisp(isInterface,2), 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', colInterface, ...
            'MarkerEdgeColor', colNodeEdge, ...
            'LineWidth', 0.75);
    end

    if any(isPrimal)
        hPrimal = plot(pDisp(isPrimal,1), pDisp(isPrimal,2), 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', colPrimal, ...
            'MarkerEdgeColor', colNodeEdge, ...
            'LineWidth', 0.75);
    end

    if any(isDirichlet)
        hDirichlet = plot(pDisp(isDirichlet,1), pDisp(isDirichlet,2), 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', colDirichlet, ...
            'MarkerEdgeColor', colNodeEdge, ...
            'LineWidth', 0.75);
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

% Leave room for legend on the right
ax = gca;
set(ax, 'units', 'normalized');
set(ax, 'position', [0.06, 0.08, 0.72, 0.84]);

if ~isempty(hInterior) && ~isempty(hInterface) && ~isempty(hPrimal) && ~isempty(hDirichlet)
    legend([hInterior, hInterface, hPrimal, hDirichlet], ...
        {'interior nodes', 'interface nodes', 'corner / primal nodes', 'Dirichlet nodes'}, ...
        'location', 'northeastoutside');
end

filename = fullfile(cfg.outdir, 'fig_4_2_dof_classification');
save_thesis_figure(fig, filename);
close(fig);