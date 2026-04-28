function cfg = setup_explanatory_figures()
%SETUP_EXPLANATORY_FIGURES Common setup for thesis explanatory figures.

  here        = fileparts(mfilename('fullpath'));  % .../main/figures/common
  figures_dir = fileparts(here);                   % .../main/figures
  main_dir    = fileparts(figures_dir);            % .../main
  root_dir    = fileparts(main_dir);               % project root

  addpath(here);
  addpath(figures_dir);
  addpath(main_dir);

  setup_paths();

  cfg = struct();

  cfg.root_dir    = root_dir;
  cfg.main_dir    = main_dir;
  cfg.figures_dir = figures_dir;

  cfg.n = 16;
  cfg.nSubX = 4;
  cfg.nSubY = 4;
  cfg.f = @(x,y) 1 + 0*x;

  cfg.outdir = fullfile(root_dir, 'output', 'figures', 'explanatory');

  if ~exist(cfg.outdir, 'dir')
    mkdir(cfg.outdir);
  end

  set(0, 'defaultaxesfontsize', 11);
  set(0, 'defaulttextfontsize', 11);
  set(0, 'defaultlinelinewidth', 1.1);
end