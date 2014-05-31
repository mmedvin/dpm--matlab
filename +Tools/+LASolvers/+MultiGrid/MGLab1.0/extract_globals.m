%EXTRACT_GLOBALS Copy global variables associated with a given grid level to 
%       local variables.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

% Get FINE and COARSE array grids

cmd_str = ['FINE = ARRAY',num2str(level),';']; eval(cmd_str)
cmd_str = ['COARSE = ARRAY',num2str(level+1),';']; eval(cmd_str)

% Get X_f,Y_f,X_c,Y_c mesh point locations

cmd_str = ['X_f = X',num2str(level),';'];  eval(cmd_str)
cmd_str = ['Y_f = Y',num2str(level),';'];  eval(cmd_str)
cmd_str = ['X_c = X',num2str(level+1),';'];  eval(cmd_str)
cmd_str = ['Y_c = Y',num2str(level+1),';'];  eval(cmd_str)

% Get nx_f,ny_f,nx_c,ny_c

cmd_str = ['nx_f = nx', num2str(level),';']; eval(cmd_str)
cmd_str = ['ny_f = ny', num2str(level),';']; eval(cmd_str)
cmd_str = ['nx_c = nx', num2str(level+1),';']; eval(cmd_str)
cmd_str = ['ny_c = ny', num2str(level+1),';']; eval(cmd_str)

% Get N_f,N_c

cmd_str = ['N_f = N', num2str(level),';'];   eval(cmd_str)
cmd_str = ['N_c = N', num2str(level+1),';']; eval(cmd_str)
