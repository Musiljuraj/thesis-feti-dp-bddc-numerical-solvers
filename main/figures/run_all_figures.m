clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));

fprintf('\nGenerating explanatory thesis figures...\n\n');

scripts = {
  'fig_3_3_membrane_surface.m'
  'fig_3_4_boundary_classification.m'
  'fig_3_5_stiffness_sparsity.m'
  'fig_4_1_domain_decomposition.m'
  'fig_4_2_dof_classification.m'
  
  
};

for k = 1:numel(scripts)
  script_path = fullfile(this_dir, scripts{k});

  if exist(script_path, 'file') ~= 2
    error('Figure script not found: %s', script_path);
  end

  fprintf('Running %s ...\n', scripts{k});
  run(script_path);
end

fprintf('\nDone.\n');