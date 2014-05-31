%PRECONDITION Precondition a vector.
%
%       Z = PRECONDITION(A, R) applies a preconditioner defined by the
%       global flag "precon_flag" to the vector R.
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function z = precondition(A, r)

include_flags 

if precon_flag == NONE

   z = r;

elseif precon_flag == JACOBI

   z = r ./ spdiags(A,[0]);

elseif precon_flag == GAUSS_SEIDEL

   z = tril(A) \ r;

elseif precon_flag == MG_CYCLE

   z = mg_cycle(1,r);

elseif precon_flag == BLOCK_JACOBI

        % ==========================
        % This assumes that nx1, N1 are available.
        % nlines is a parameter that specifies the block size in lines
        % nlines could be a global, set-able through the
        %     Solver->Precon->Block Jacobi->nlines submenu
        % ==========================
        nlines = 4;
        blk = nlines * nx1;
        nb = ceil(N1/blk);
 
        z = zeros(size(r));
 
        for j = 1:nb
           v = [(j-1)*blk+1 : min(N1, j*blk)]';
           z(v) = A(v,v) \ r(v);
        end
        % ==========================

end
