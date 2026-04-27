% ============================================================
% File: src/common/spectra/full_spectrum_operator.m
% ============================================================
function spec = full_spectrum_operator(applyA, n, opts)
%FULL_SPECTRUM_OPERATOR Compute the full spectrum of a matrix-free SPD operator.
% Thesis link: Chapter 6.3 (explicit spectral analysis from operator actions).
% The routine assembles an explicit matrix representation of an operator from
% its action on canonical basis vectors and computes its eigenvalues.
%
% This helper is intended for the additional diagnostic spectra of the
% non-preconditioned reduced operators A_FETI-DP and A_BDDC. In contrast to
% full_spectrum_precond.m, no preconditioner action is used here.
%
% Inputs:
%   applyA : function handle, y = applyA(x), x,y length n.
%   n      : positive integer, dimension of the operator domain/range.
%   opts   : optional struct with fields:
%            - symmetrize      (default true)  : replace A by its symmetric
%                                               part to suppress roundoff
%                                               asymmetry.
%            - store_matrices   (default false) : store A in output.
%            - verbose         (default false) : print progress messages.
%            - assemble_opts   (default struct): options passed to
%                                               assemble_from_apply.
%
% Output (struct spec):
%   spec.eigvals       : sorted eigenvalues of A.
%   spec.lmin          : minimum eigenvalue.
%   spec.lmax          : maximum eigenvalue.
%   spec.kappa         : lmax / lmin (Inf if lmin <= 0).
%   spec.symmA         : relative symmetry defect of assembled A before
%                        optional symmetrization.
%   spec.max_imag_eig  : maximum absolute imaginary part returned by eig(A).
%   spec.n_nonpositive : number of eigenvalues <= 0 after taking real parts.
%   (optional) spec.A if opts.store_matrices = true.
%
% Notes:
% - This routine is intended only for small n, because it explicitly assembles
%   a dense matrix and then computes the full eigenvalue decomposition.
% - In exact arithmetic, the operators used here should be symmetric in the
%   considered settings. Small numerical asymmetries may appear because the
%   matrix is assembled from operator actions, so symmetrization is enabled by
%   default.

  % ----------------------------
  % Input validation
  % ----------------------------
  if nargin < 2 || nargin > 3
    error('full_spectrum_operator:InvalidNargin: expected inputs (applyA, n [, opts]).');
  end

  if ~isa(applyA, 'function_handle')
    error('full_spectrum_operator:InvalidApplyA: applyA must be a function handle.');
  end

  if ~(isscalar(n) && isnumeric(n) && isreal(n) && isfinite(n) && n == floor(n) && n >= 1)
    error('full_spectrum_operator:InvalidN: n must be a positive integer scalar.');
  end

  if nargin < 3 || isempty(opts)
    opts = struct();
  end

  if ~isstruct(opts)
    error('full_spectrum_operator:InvalidOpts: opts must be a struct (or []).');
  end

  % ----------------------------
  % Options (defaults)
  % ----------------------------
  if ~isfield(opts, 'symmetrize') || isempty(opts.symmetrize)
    opts.symmetrize = true;
  end
  if ~islogical(opts.symmetrize) || ~isscalar(opts.symmetrize)
    error('full_spectrum_operator:InvalidSymmetrize: opts.symmetrize must be a logical scalar.');
  end

  if ~isfield(opts, 'store_matrices') || isempty(opts.store_matrices)
    opts.store_matrices = false;
  end
  if ~islogical(opts.store_matrices) || ~isscalar(opts.store_matrices)
    error('full_spectrum_operator:InvalidStoreMatrix: opts.store_matrices must be a logical scalar.');
  end

  if ~isfield(opts, 'verbose') || isempty(opts.verbose)
    opts.verbose = false;
  end
  if ~islogical(opts.verbose) || ~isscalar(opts.verbose)
    error('full_spectrum_operator:InvalidVerbose: opts.verbose must be a logical scalar.');
  end

  if ~isfield(opts, 'assemble_opts') || isempty(opts.assemble_opts)
    opts.assemble_opts = struct();
  end
  if ~isstruct(opts.assemble_opts)
    error('full_spectrum_operator:InvalidAssembleOpts: opts.assemble_opts must be a struct (or []).');
  end

  % If verbose is enabled here, enable coarse progress in assembly as well.
  if opts.verbose && ~isfield(opts.assemble_opts, 'verbose')
    opts.assemble_opts.verbose = true;
  end

  % ----------------------------
  % Initialize output structure
  % ----------------------------
  spec = struct();
  spec.eigvals       = [];
  spec.lmin          = NaN;
  spec.lmax          = NaN;
  spec.kappa         = NaN;

  spec.symmA         = NaN;
  spec.max_imag_eig  = NaN;
  spec.n_nonpositive = NaN;

  % ----------------------------
  % Assemble explicit matrix from matrix-free action
  % ----------------------------
  if opts.verbose
    fprintf('[full_spectrum_operator] assembling A (n=%d)\n', n);
  end

  A = assemble_from_apply(applyA, n, opts.assemble_opts);

  % ----------------------------
  % Symmetry diagnostic before optional symmetrization
  % ----------------------------
  spec.symmA = rel_symm_defect_(A);

  % ----------------------------
  % Optional symmetrization
  % ----------------------------
  % In exact arithmetic, the reduced operators considered here are symmetric.
  % Since the matrix is assembled numerically from operator applications, small
  % roundoff asymmetries can occur. We remove these artifacts before eig(A).
  if opts.symmetrize
    A = 0.5 * (A + A');
  end

  % ----------------------------
  % Eigenvalue computation
  % ----------------------------
  if opts.verbose
    fprintf('[full_spectrum_operator] eig(A)\n');
  end

  eigvals_raw = eig(A);

  if isempty(eigvals_raw)
    return;
  end

  % Eigenvalues should be real after symmetrization; record any imaginary part
  % for diagnostics and then keep the real parts, consistently with the
  % preconditioned spectral routine.
  spec.max_imag_eig = max(abs(imag(eigvals_raw(:))));

  eigvals = real(eigvals_raw(:));
  eigvals = sort(eigvals, 'ascend');

  spec.eigvals = eigvals;

  % ----------------------------
  % Derived scalars: lmin, lmax, kappa
  % ----------------------------
  spec.lmin = eigvals(1);
  spec.lmax = eigvals(end);

  spec.n_nonpositive = sum(eigvals <= 0);

  if spec.lmin > 0
    spec.kappa = spec.lmax / spec.lmin;
  else
    spec.kappa = Inf;
  end

  % ----------------------------
  % Optional storage of matrix
  % ----------------------------
  if opts.store_matrices
    spec.A = A;
  end
end


% ============================================================
% Local helper: relative symmetry defect
% ============================================================
function s = rel_symm_defect_(M)
%REL_SYMM_DEFECT_  Relative symmetry defect ||M-M^T||_F / ||M||_F.
  nrm = norm(M, 'fro');
  if nrm == 0
    s = 0;
  else
    s = norm(M - M', 'fro') / nrm;
  end
end