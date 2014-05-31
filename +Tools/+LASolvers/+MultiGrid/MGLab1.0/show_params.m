%SHOW_PARAMS Display the currently set parameters in a text window.
%
%       This M-File displays a summary of the currently defined global
%       parameters in the figure with handle "param_fig".
%
%       Accesses global variables in "include_globals"
%       Accesses global variables in "include_flags"
%       Accesses global variables in "include_figs"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

include_globals
include_flags
include_figs

% == Problem ============================================

p  =[sprintf(' Problem                                  ')];
if (problem_flag==POISSON)
p=[p;sprintf('    Poisson Equation                      ')];
elseif (problem_flag==HELMHOLTZ)
p=[p;sprintf('    Helmholtz Equation                    ')];
if imag(prob_args(1)) == 0
  p=[p;sprintf('        k = %6.2f                        ',prob_args(1))];
else
  p=[p;sprintf('        k = %6.2f + %6.2f * i           ',real(prob_args(1)), imag(prob_args(1)))];
end
elseif (problem_flag==CONVECT_DIFFUSE)
p=[p;sprintf('    Convection-Diffusion Equation         ')];
p=[p;sprintf('        lambda = %6d                   ',prob_args(1))];
p=[p;sprintf('        sigma  = %6d                   ',prob_args(2))];
elseif (problem_flag==CUT_SQUARE)
p=[p;sprintf('    Cut Square Equation                   ')];
p=[p;sprintf('        alpha  = %8.3f                 ',prob_args(1))];
end

p=[p;sprintf('    Problem size = (%3d,%3d)              ',nx1,ny1)];

% == Solver Menu Item ==============================================

s  =[sprintf(' Solver                                   ')];
if (solver_flag==VMG & cycle_flag==V_CYCLE)
  s=[s;sprintf('    V-Cycle                               ')];
elseif (solver_flag==VMG & cycle_flag==W_CYCLE)
  s=[s;sprintf('    W-Cycle                               ')];
elseif (solver_flag==PCG)
  s=[s;sprintf('    PCG                                   ')];
elseif (solver_flag==BICG_STAB)
  s=[s;sprintf('    Bi-CG(Stab)                           ')];
elseif (solver_flag==CGS)
  s=[s;sprintf('    CGS                                   ')];
elseif (solver_flag==GMRES)
  s=[s;sprintf('    GMRES (%3d)                           ',restart)];
elseif (solver_flag==FMG)
  s=[s;sprintf('    FMG                                   ')];
elseif (solver_flag==SOR)
  if SOR_omega == 1
      s=[s;sprintf('    Gauss-Seidel                          ')];
  else
      s=[s;sprintf('    SOR, omega = %8.2f                 ', SOR_omega)];
  end
end

% == Preconditoner ==============================================

r  =[sprintf(' Preconditioner                           ')];
if (precon_flag==MG_CYCLE & cycle_flag==V_CYCLE)
r=[r;sprintf('    V-Cycle                               ')];
elseif (precon_flag==MG_CYCLE & cycle_flag==W_CYCLE)
r=[r;sprintf('    W-Cycle                               ')];
elseif (precon_flag==JACOBI)
r=[r;sprintf('    Jacobi                                ')];
elseif (precon_flag==BLOCK_JACOBI)
r=[r;sprintf('    Block-Jacobi                          ')];
elseif (precon_flag==GAUSS_SEIDEL)
r=[r;sprintf('    Gauss-Seidel                          ')];
elseif (precon_flag==ILU)
r=[r;sprintf('    ILU                                   ')];
elseif (precon_flag==SSOR)
r=[r;sprintf('    SSOR                                  ')];
elseif (precon_flag==NONE)
r=[r;sprintf('    NONE                                  ')];
end

% == Stopping-Criteria ===========================================

t  =[sprintf(' Stopping Criteria                        ')];

if (prtol==0)
t=[t;sprintf('    (precon) Residual Tolerence = NONE    ')];
else
t=[t;sprintf('    (precon) Residual Tolerance = %4.0e   ',prtol)];
end
if (max_it==0)
t=[t;sprintf('    Iteration Limit = NONE                ')];
else
t=[t;sprintf('    Iteration Limit = %4d                ',max_it)];
end

% == MG Parameters ===========================================

m  =[sprintf(' MG Parameters                            ')];

m=[m;sprintf('    Number of Levels = %1d                  ',coarse_level)];
if (smooth_flag==WEIGHTED_JACOBI)
m=[m;sprintf('    Smoother = Weighted Jacobi (w=%4.2f)   ',wt)];
elseif (smooth_flag==GAUSS_SEIDEL)
m=[m;sprintf('    Smoother = Gauss-Seidel               ')];
elseif (precon_flag==RB_GAUSS_SEIDEL)
m=[m;sprintf('    Red/Black Gauss-Seidel                ')];
end
m=[m;sprintf('    Pre-smoothings =%3d                   ',nu1)];
m=[m;sprintf('    Post-smoothings =%3d                  ',nu2)];

if (restrict_flag==INJECTION)
m=[m;sprintf('    Restriction = Injection               ')];
elseif (restrict_flag==HALF_WEIGHTING)
m=[m;sprintf('    Restriction = Half Weighting          ')];
elseif (restrict_flag==FULL_WEIGHTING)
m=[m;sprintf('    Restriction = Full Weighting          ')];
elseif (restrict_flag==BILINEAR_ADJOINT)
m=[m;sprintf('    Restriction = Bilinear Adjoint        ')];
end

if (interp_flag==LINEAR)
m=[m;sprintf('    Prolongation = Linear                 ')];
elseif (interp_flag==CUBIC)
m=[m;sprintf('    Prolongation = Cubic                  ')];
elseif (interp_flag==OPERATOR_BASED)
m=[m;sprintf('    Prolongation = Operator-Based         ')];
elseif (interp_flag==EXPLICIT_BILINEAR)
m=[m;sprintf('    Prolongation = Explicit Bilinear      ')];
end

if (coarse_solver_flag==DIRECT)
m=[m;sprintf('    Coarse-grid Solver = Sparse GE        ')];
elseif (coarse_solver_flag==SMOOTHER)
m=[m;sprintf('    Prolongation = Smoother (30)          ')];
elseif (coarse_solver_flag==PCG)
m=[m;sprintf('    Prolongation = PCG                    ')];
elseif (coarse_solver_flag==BICG_STAB)
m=[m;sprintf('    Prolongation = Bi-CG(Stab)            ')];
elseif (coarse_solver_flag==GMRES)
m=[m;sprintf('    Prolongation = GMRES(%3d)             ',restart)];
end

if (coarsening_flag==STANDARD)
m=[m;sprintf('    Coarsener = Standard 5-point          ')];
elseif (coarsening_flag==GALERKIN)
m=[m;sprintf('    Coarsener = Galerkin                  ')];
elseif (coarsening_flag==AVERAGING)
m=[m;sprintf('    Coarsener = Coeff. Averaging          ')];
end

if (cycle_flag==V_CYCLE)
m=[m;sprintf('    MG Cycle Type = V-Cycle               ')];
elseif (cycle_flag==W_CYCLE)
m=[m;sprintf('    MG Cycle Type = W-Cycle               ')];
elseif (cycle_flag==HALF_V_CYCLE)
m=[m;sprintf('    MG Cycle Type = Half V-Cycle          ')];
end

page=[p;s;r;t;m];


if (param_fig==0)
    param_fig = figure('Position', param_position,...
               'Name', 'Parameters',...
               'NumberTitle', 'off', ...
               'Color','blue');
end

    figure(param_fig);
    [mm,nn] = size(page);
    subplot(1,1,1)
    hold off
    cla
    ht = 0.95;
    for j = 1:mm
       text(0.05, ht, page(j,:))
       axis('off')
       ht = ht - 0.05;
    end
   figure(main_fig);
