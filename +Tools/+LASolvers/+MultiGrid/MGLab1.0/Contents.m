% MGLab V1.0
% 
% The Graphical User Interface
% 
%     MGLab.m
%     show_params.m
%     set_defaults.m
%     run.m
%     multigrid_setup.m
%     menu_header.m
%     menu_item.m
%     version_info.m
%     Contents.m
% 
% Global Variables
% 
%     include_globals.m
%     extract_globals.m
%     include_flags.m
%     include_figs.m
% 
% Linear Systems
% 
%     get_matrix.m
%     sp_laplace.m
%     sp_cutsq2d.m
%     sp_convdiff.m
%     get_rhs.m
% 
% Solvers
% 
%     solve.m
%     converged.m
%     vmg.m
%     fmg.m
%     pcg.m
%     pcgs.m
%     pbicgstab.m
%     pgmres.m
%     sor.m
%     get_SOR_omega.m
% 
% Preconditioners
% 
%     precondition.m
%     precond_mg.m
% 
% Results
% 
%     update_results.m
% 
% Multigrid Routines: High Level
% 
%     mg_cycle.m
%     vmg_cycle.m
%     fmg_cycle.m
%     wmg_cycle.m
%     halfvmg_cycle.m
% 
% Multigrid Routines: Middle Level
% 
%     smooth.m
%     residual.m
%     restrict.m
%     coarse_grid_solve.m
%     interpolate.m
%     sp_prolong.m
% 
% Multigrid Routines: Low Level
% 
%     coarsest.m
%     max_level.m
% 
% Demos
% 
%     demo_globals.m
%     demo1.m
%     demo1_run.m
%     demo2.m
%     demo2_run.m
%     demo2_Vcycle.m 
%     demo3.m
%     demo3_run.m
%     demo3_vmg.m
%     sint.m
%     sint2.m
