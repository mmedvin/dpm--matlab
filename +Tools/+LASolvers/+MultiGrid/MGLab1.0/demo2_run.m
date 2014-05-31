%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function demo2_run

include_globals
include_flags

set_defaults
coarse_level = 2;  % Two grid Algorithm

nx1 = 49; ny1 = nx1;

%smooth_flag = WEIGHTED_JACOBI;  % Select smoother
smooth_flag = GAUSS_SEIDEL;
nu1 = 0;   % No pre-smoothing
nu2 = 4;

prob_args=[10];

[A1,N1] = get_matrix(nx1,ny1);
generate_matrix = 0;

if (solver_flag==VMG | precon_flag==MG_CYCLE)
   multigrid_setup;
end

b1 = get_rhs(nx1,ny1);   
h = 1/(nx1+1);
[X,Y] = meshgrid([h:h:(1-h)]);
U_TRUE = sin(pi*X).*sin(pi*Y);
INIT_ERROR = sin(10*pi*X) .* sin(10*pi*Y) + ...
             sin(20*pi*X) .* sin(20*pi*Y) + ... 
             sin(30*pi*X) .* sin(30*pi*Y) + ...
             sin(40*pi*X) .* sin(40*pi*Y);
init_error = reshape(INIT_ERROR, N1, 1);
u_true = reshape(U_TRUE,N1,1);
b1 = A1 * u_true;
generate_rhs = 0;

level = 1;

u_rand = 2*rand(nx1*nx1,1)-1;   
u_in = u_true - init_error + 0.5 * u_rand;
u_out = u_in;

RELRES = [];

jj = 4;
for iter = 1:jj
   u_out = demo2_Vcycle(level, b1, u_out, u_true, nx1,iter);
   relres = norm(b1-A1*u_out)/norm(b1); RELRES = [RELRES; relres];
   fprintf('Relative residual = %g, at iteration %g\n', relres, iter);
end

logrho = (1/(jj-2)) * log10(RELRES(jj)/RELRES(2));
rho = 10^logrho;
fprintf('Average residual reduction factor = %g.\n', rho);

hold off

% Reset axes
% close(fig1); close; MGLab
