%INCLUDE_GLOBALS Global matrices, vectors and scalars.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

% Global variables associated with the MG meshes

  global coarse_level

  global nx1 ny1 N1 A1 ARRAY1 X1 Y1 b1
  global nx2 ny2 N2 A2 ARRAY2 X2 Y2 
  global nx3 ny3 N3 A3 ARRAY3 X3 Y3
  global nx4 ny4 N4 A4 ARRAY4 X4 Y4
  global nx5 ny5 N5 A5 ARRAY5 X5 Y5

% Global parameters for the problems

  global prob_args

% Global parameters for the solver

  global rtol prtol max_it max_time max_mflop num_runs restart SOR_omega

% Global parameters for the smoother

  global nu1 nu2 wt 

% Global parameters for the linear system

  global generate_matrix generate_rhs matrix_type dimensions
