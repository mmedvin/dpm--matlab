%INCLUDE_FLAGS Flag variables defining current settings in MGLab.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

% Available Problems

   global problem_flag

   POISSON = 101;
   HELMHOLTZ = 102;
   CONVECT_DIFFUSE = 103;
   POISSON_BOLTZMAN = 104;
   CUT_SQUARE = 105;

% Available Interpolation operators

   global interp_flag

   CUBIC = 201;
   LINEAR = 202;
   EXPLICIT_BILINEAR = 203;
   OPERATOR_BASED = 204;

% Available Restriction operators

   global restrict_flag

   INJECTION = 301;
   HALF_WEIGHTING = 302;
   FULL_WEIGHTING = 303;
   BILINEAR_ADJOINT = 304;

% Available Smoothers

   global smooth_flag

   WEIGHTED_JACOBI = 401;
   GAUSS_SEIDEL = 402;
   RB_GAUSS_SEIDEL = 403;

% Available Preconditioners

   global precon_flag

   NONE = 500;
   JACOBI = 501;
   MG_CYCLE = 502;
   ILU = 505;
   SSOR = 506; 
   BLOCK_JACOBI = 507;

% Available Solvers

   global solver_flag

   DIRECT = 600;               % for Coarse-grid solve only
   VMG = 601;
   FMG = 602;
   PCG = 603;
   BICG_STAB = 604;
   GMRES = 605;
   SMOOTHER = 606;             % for Coarse-grid solve only
   CGS = 607;
   TFQMR = 608;
   SOR = 609;

% Available Coarsenings

   global coarsening_flag

   STANDARD = 701;
   GALERKIN = 702;
   AVERAGING = 703;

% Available Coarse-Grid Solves

   global coarse_solver_flag

% Right-hand Side

   global rhs_flag

% Available Plotting axes

  global x_axis_flag y_axis_flag

  ITERATIONS=801;
  TIME=802;
  MFLOPS=803;
  RESIDUAL=804;
  PRECON_RESIDUAL=805;
  ERROR=806;

% Available Cycles

   global cycle_flag

   V_CYCLE = 901;
   W_CYCLE = 902;
   HALF_V_CYCLE = 903;
