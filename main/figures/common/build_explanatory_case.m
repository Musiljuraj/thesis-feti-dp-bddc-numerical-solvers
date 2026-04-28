function data = build_explanatory_case(cfg, need_ddm)
%BUILD_EXPLANATORY_CASE Build small FEM/DDM case for explanatory figures.

  if nargin < 2
    need_ddm = false;
  end

  data = struct();

  [p, t, bnd] = mesh_unit_square_P1(cfg.n);
  K = assemble_stiffness_P1(p, t);
  fvec = assemble_load_P1(p, t, cfg.f);

  [Kff, ff, free] = apply_dirichlet_elimination(K, fvec, bnd.dirichlet_nodes);

  data.p = p;
  data.t = t;
  data.bnd = bnd;
  data.K = K;
  data.fvec = fvec;
  data.Kff = Kff;
  data.ff = ff;
  data.free = free;
  data.dirichlet = bnd.dirichlet_nodes;

  if need_ddm
    data.problem = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, cfg.f);
  end
end