classdef ClassQinghaiMG < handle
    properties(Access = protected)
        % matrix;
        param;
        A; Ir; Ip;
        Matrix;
        
        Size;
        Levels;
    end
    
    methods (Access = public)
        function [Result,resRel, nIters] = Solve(obj,Rhs)
            Result = zeros(obj.Size^2,1);  % use zero as the initial guess.
            res0 = norm(Rhs,inf); % the initial residual
            
            jj=1;
            resRel = 1;
             while (resRel>= obj.param.TOL)  &&  (jj<=obj.param.maxIt) % for jj=1:obj.param.maxIt
                 
                if obj.param.cycleType=='v'
                    %        uS = vcycle_block(uS, rhs, obj.Levels);
                    Result = obj.vcycle(Result, Rhs, obj.Levels);
                elseif obj.param.cycleType=='f'
                    %        uS = fcycle_block(uS, rhs, obj.Levels);
                    Result = obj.fcycle(Result, Rhs, obj.Levels);
                else
                    error('unknown cycle type!');
                end
                
                tmp = Rhs-obj.Matrix*Result;
                res = norm(tmp,inf);
                resRel = res/res0;
                
                if obj.param.verbose == 1
                    info = ['  iter #', int2str(jj), ': max-norm of Residual=', num2str(resRel)];%, '; max-norm of error=', num2str(errRel)];
                    disp(info);
                end
                jj=jj+1;
                %    if (errRel<param.TOL)
%                 if (resRel<param.TOL)
%                     break
%                 end
            end
            nIters = jj-1;
        end
        
        %CONSTRUCTOR
        % mg_matrix_init.m
        %
        %   initialize the matrices used in the mulitgrid algorithm.
        %
        function obj = ClassQinghaiMG(Matrix,Size,tag, param,nVoffset)
            % tag can be 'u','v', or 'p'
            % nVoffset gives the flexibility of node-centering or cell-centering.
            %  e.g. nVoffset=-1 :  the nodes excluding the boundaries are variables
            %  e.g. nVoffset=+1 :  the nodes including the boundaries are variables
            
            obj.Size = Size;
            obj.Levels=log2(obj.Size-nVoffset);
            
            % initialize the data structure
            %
            obj.param = param;
            obj.A  = cell(obj.Levels,1);
            obj.Ir = cell(obj.Levels,1);
            obj.Ip = cell(obj.Levels,1);
            
            % initialize
            %
            obj.Matrix = Matrix;
            Ac=Matrix;
            
            for k=obj.Levels:-1:1
                N = 2^k+nVoffset;  % size of the current grid
                
                % set the operator
                %
                obj.A{k}  = Ac;
                
                % transfer matricies
                %
                if tag == 'u6'
                    interp  = obj.interp_matrix_u_6pt(N,N);
                elseif tag == 'u1'
                    interp  = obj.interp_matrix_u(N,N);
                elseif tag == 'v6'
                    interp  = obj.interp_matrix_v_6pt(N,N);
                elseif tag == 'v1'
                    interp  = obj.interp_matrix_v(N,N);
                elseif tag == 'p0'
                    interp  = obj.interp_matrix_p_const(N,N);
                elseif tag == 'p1'
                    interp  = obj.interp_matrix_p(N,N,nVoffset);
                else
                    error('unsupported variable tag!');
                end
                
                obj.Ip{k} = interp;
                % the denominator 4 follows from the dimensionality being 2
                %  and the coefficients of the bilinear interpolation.
                % In other words, we sum up the coefficients of fine values
                %  that are influenced from a single coarse value.
                obj.Ir{k} = interp'/4;
                
                % compute the new coarse grid matrix
                %
                Ac = obj.Ir{k} * Ac * obj.Ip{k};
                
            end
        end
        
        
        
        % fcycle.m
        %
        %   recursive fcycle code
        %
        function v = fcycle(obj,u,f,level)
            
            % data structure with the matricies, and parameters
            %
            %   global mg_matrix
            %   global mg_param
            
            % multigrid parameters
            %
            nu1      = obj.param.nu1;
            nu2      = obj.param.nu2;
            minlevel = obj.param.minlevel;
            
            % if on the last level solve the muthafucka
            %
            if( level == minlevel )
                B = obj.A{level};
                if (condest(B)>1.e12)
                    sz = size(B);
                    B = [B; ones(1,sz(2))];
                    f = [f;0];
                end
                v = B\f;
                %v = mg_matrix.A{level}\f;
            else
                
                % presmooth
                %
                v = obj.gs_lex(u,obj.A{level},f,nu1);
                
                % compute residual
                %
                r = f - obj.A{level}*v;
                
                % restrict
                %
                fc = obj.Ir{level}*r;
                
                % initialize
                %
                uc = 0.0*fc;
                
                % solve the error equation, recursively
                %
                ec = obj.fcycle(uc,fc,level-1);
                
                % correct
                %
                v = v + obj.Ip{level} * ec;
                
                % post smoothing
                %
                v = obj.gs_lex(v,obj.A{level},f,nu2);
                
                %%%%% now a vcycle
                
                % compute residual
                %
                r = f - obj.A{level}*v;
                
                % restrict
                %
                fc = obj.Ir{level}*r;
                
                % initialize
                %
                uc = 0.0*fc;
                
                % down the v
                %
                ec = obj.vcycle_sym(uc,fc,level-1);
                
                % correct
                %
                v = v + obj.Ip{level} * ec;
                
                % post smoothing
                %
                v = obj.gs_lex(v,obj.A{level},f,nu2);
                
            end
        end
        
        
        % fcycle_block
        %
        %   recursive fcycle code using block smoother
        %
        function v = fcycle_block(obj,u,f,level)
            
            % data structure with the matricies, and parameters
            %
            %global mg_matrix
            %global mg_param
            
            % multigrid parameters
            %
            nu1      = obj.param.nu1;
            nu2      = obj.param.nu2;
            minlevel = obj.param.minlevel;
            
            % block smoother parameters
            %
            bx = obj.param.bx;
            by = obj.param.by;
            sx = obj.param.sx;
            sy = obj.param.sy;
            
            
            % if on the last level solve the muthafucka
            %
            if( level == minlevel )
                B = obj.A{level};
                if (condest(B)>1.e12)
                    sz = size(B);
                    B = [B; ones(1,sz(2))];
                    f = [f;0];
                end
                v = B\f;
                %    v = mg_matrix.A{level}\f;
            else
                
                % presmooth
                %
                Nx = 2^level;  % size of the current grid
                Ny = Nx;
                v = obj.blocksmooth(u,obj.A{level},f,nu1,Nx,Ny,bx,by,sx,sy);
                
                % compute residual
                %
                r = f - obj.A{level}*v;
                
                % restrict
                %
                fc = obj.Ir{level}*r;
                
                % initialize
                %
                uc = 0.0*fc;
                
                % solve the error equation, recursively
                %
                ec = obj.fcycle_block(uc,fc,level-1);
                
                % correct
                %
                v = v + obj.Ip{level} * ec;
                
                % post smoothing
                %
                v = obj.blocksmooth(v,obj.A{level},f,nu2,Nx,Ny,bx,by,sx,sy);
                
                %%%%% now a vcycle
                
                % compute residual
                %
                r = f - obj.A{level}*v;
                
                % restrict
                %
                fc = obj.Ir{level}*r;
                
                % initialize
                %
                uc = 0.0*fc;
                
                % down the v
                %
                ec = obj.vcycle_block(uc,fc,level-1);
                
                % correct
                %
                v = v + obj.Ip{level} * ec;
                
                % post smoothing
                %
                v = obj.blocksmooth(v,obj.A{level},f,nu2,Nx,Ny,bx,by,sx,sy);
                
            end
        end
        
        
        % vcycle.m
        %
        %   recursive vcycle code
        %
        function v = vcycle(obj,u,f,level)
            
            % data structure with the matricies, and parameters
            %
            %   global mg_matrix
            %   global mg_param
            
            % multigrid parameters
            %
            nu1      = obj.param.nu1;
            nu2      = obj.param.nu2;
            minlevel = obj.param.minlevel;
            
            % if on the last level solve the muthafucka
            %
            if( level == minlevel )
                v = obj.A{level}\f;
            else
                
                % presmooth
                %
                v = obj.gs_lex(u,obj.A{level},f,nu1);
                
                % compute residual
                %
                r = f - obj.A{level}*v;
                
                % restrict
                %
                fc = obj.Ir{level}*r;
                
                % initialize
                %
                uc = 0.0*fc;
                
                % solve the error equation, recursively
                %
                ec = obj.vcycle(uc,fc,level-1);
                
                % correct
                %
                v = v + obj.Ip{level} * ec;
                
                % post smoothing
                %
                v = obj.gs_lex(v,obj.A{level},f,nu2);
                
            end
        end
        
        
        %
        % vcycle_block.m
        %
        %   recursive vcycle code using block smoother
        %
        function v = vcycle_block(obj,u,f,level)
            
            %disp(sprintf('level = %5i',level));
            
            % data structure with the matricies, and parameters
            %
            %   global mg_matrix
            %   global mg_param
            
            % multigrid parameters
            %
            nu1      = obj.param.nu1;
            nu2      = obj.param.nu2;
            minlevel = obj.param.minlevel;
            
            % block smoother parameters
            %
            bx = obj.param.bx;
            by = obj.param.by;
            sx = obj.param.sx;
            sy = obj.param.sy;
            
            
            % if on the last level solve the muthafucka
            %
            if( level == minlevel )
                B = obj.A{level};
                if (condest(B)>1.e12) % this is for pressure or phi
                    sz = size(B);
                    B = [B; ones(1,sz(2))];
                    f = [f;0];
                end
                v = B\f;
                %v = mg_matrix.A{level}\f;
            else
                
                % presmooth
                %
                Nx = 2^level;  % size of the current grid
                Ny = Nx;
                v = obj.blocksmooth(u,obj.A{level},f,nu1,Nx,Ny,bx,by,sx,sy);
                
                % compute residual
                %
                r = f - obj.A{level}*v;
                
                % restrict
                %
                fc = obj.Ir{level}*r;
                
                % initialize
                %
                uc = 0.0*fc;
                
                % solve the error equation, recursively
                %
                ec = obj.vcycle_block(uc,fc,level-1);
                
                % correct
                %
                v = v + obj.Ip{level} * ec;
                
                % post smoothing
                %
                v = obj.blocksmooth(v,obj.A{level},f,nu2,Nx,Ny,bx,by,sx,sy);
            end
        end
        
        %
        % vcycle_sym.m
        %
        %   recursive vcycle code with symmetric smoother
        %
        function v = vcycle_sym(obj,u,f,level)
            
            % data structure with the matricies, and parameters
            %
            %   global mg_matrix
            %   global mg_param
            
            % multigrid parameters
            %
            nu1      = obj.param.nu1;
            nu2      = obj.param.nu2;
            minlevel = obj.param.minlevel;
            
            % if on the last level solve the muthafucka
            %
            if( level == minlevel )
                B = obj.A{level};
                if (condest(B)>1.e12)
                    sz = size(B);
                    B = [B; ones(1,sz(2))];
                    f = [f;0];
                end
                v = B\f;
            else
                
                % presmooth
                %
                v = obj.gs_lex_sym(u,obj.A{level},f,nu1);
                
                
                % compute residual
                %
                r = f - obj.A{level}*v;
                
                % restrict
                %
                fc = obj.Ir{level}*r;
                
                % initialize
                %
                uc = 0.0*fc;
                
                % solve the error equation, recursively
                %
                ec = obj.vcycle_sym(uc,fc,level-1);
                
                % correct
                %
                v = v + obj.Ip{level} * ec;
                
                % post smoothing
                %
                v = obj.gs_lex_sym(v,obj.A{level},f,nu2);
                
            end
        end
        
    end
    
    methods (Access = protected) %aka core
        
        % blocksmooth.m -- perform block GS lex smoothing with overlapping blocks
        %
        %   input: u  -- current solution
        %          A  -- matrix (Nx*Ny) x (Nx*Ny)
        %          f  -- rhs
        %          Nx -- number of points in the x-direction
        %          Ny -- number of points in the y-direction
        %          bx -- block size in x direction
        %          by -- block size in y direciton
        %          sx -- skip size in x-direction, overlap size is bx-sx
        %          sy -- skip size in y-direction  overlap size is by-sy
        %          ns -- number of smoothing steps to apply
        %
        %   output: v -- smoothed solution
        %
        function v = blocksmooth(obj,u,A,f,ns,Nx,Ny,bx,by,sx,sy)
            
            % initialize
            %
            v = u;
            
            % apply some krylov iterations before smoothing
            %
            %v = pcg(A,f,1e-10,2,[],[],u);
            %v = gmres(A,f,2,1e-10,1,[],[],u);
            %v = bicgstab(A,f,1e-10,2,[],[],u);
            
            for n=1:ns
                
                for j1=1:sy:Ny
                    for i1=1:sx:Nx
                        
                        % form the restriction operator
                        %
                        i2 = min(i1 + (bx-1), Nx );
                        j2 = min(j1 + (by-1), Ny );
                        
                        bt = (i2-i1+1)*(j2-j1+1);  % total block size
                        I = zeros(bt,1);
                        J = I;
                        k = 1;
                        for i3=i1:i2
                            for j3=j1:j2
                                I(k) = i3;
                                J(k) = j3;
                                k=k+1;
                            end
                        end
                        
                        K = sub2ind([Nx,Ny],J,I); % double check whether I or J
                        % should come first
                        R = spalloc(bt,Nx*Ny,bt);
                        R(1:bt,K) = speye(bt);
                        
                        % compute the residual and restrict
                        %
                        RA = R*A;
                        rr = R*f - RA*v;
                        
                        
                        % form the restricted matrix, and compute the update
                        %
                        Ak = RA*R';
                        du = Ak\rr;
                        
                        % update the solution
                        %
                        v = v + R'*du;
                        
                    end
                end
            end
        end
        
        
        % cc_interp1.m
        %
        %  Form the 1D linear interpolation matrix for a cell centered grid
        %   with periodic boundary condtions.
        %
        % input -- Nf = number of cell centers on the fine grid
        % output --
        %
        function Ip = cc_interp1(obj,Nf)
            
            if Nf==2
                Ip = obj.cc_nm_interp1_const(Nf);
                return;
            end
            
            % number of centers on the coarse grid
            %
            Nc = Nf/2;
            
            % allocate space for the matrix
            %
            Ip = spalloc(Nf,Nc,4*Nc);
            
            % loop over columns
            %
            for j=1:Nc
                if j==1
                    Irange = [Nf,1:3];
                elseif j==Nc
                    Irange = [Nf-2:Nf,1];
                else
                    Irange = [2*j-2:2*j+1];
                end
                Ip(Irange,j) = [0.25; 0.75; 0.75; 0.25];
            end
        end
        
        
        % cc_nm_interp1.m
        %
        %  Form the 1D linear interpolation matrix for a cell centered grid.
        %
        % input -- Nf = number of cell centers on the fine grid
        % output --
        %
        function Ip = cc_nm_interp1(obj,Nf)
            
            if Nf<=2
                Ip = obj.cc_nm_interp1_const(Nf);
                return;
            end
            
            % number of centers on the coarse grid
            %
            Nc = Nf/2;
            
            % allocate space for the matrix
            %
            Ip = spalloc(Nf,Nc,3*Nc);
            
            % first column
            %
            Ip(1:3,1) = [0.75; 0.75; 0.25];
            % periodic BC, change needed for other BCs.
            Ip(Nf,1) = 0.25;
            
            % last column
            %
            Ip((Nf-2):Nf,Nc) = [0.25; 0.75; 0.75];
            % periodic BC, change needed for other BCs.
            Ip(1,Nc) = 0.25;
            
            % loop over columns
            %
            for j=2:Nc-1
                ifine = 2*j;
                Irange = [ifine-2:ifine+1];
                Ip(Irange,j) = [0.25; 0.75; 0.75; 0.25];
            end
        end
        
        %
        % cc_nm_interp1_const.m
        %
        %  Form the 1D constat interpolation matrix for a cell centered grid
        %
        % input -- Nf = number of cell centers on the fine grid
        % output --
        %
        function Ip = cc_nm_interp1_const(obj,Nf)
            
            % number of centers on the coarse grid
            %
            Nc = Nf/2;
            
            % allocate space for the matrix
            %
            Ip = spalloc(Nf,Nc,2*Nc);
            
            % copy value to cell above
            %
            for j=1:Nc
                ifine = 2*j;
                Irange = [ifine-1:ifine];  % the fine cells to be affected.
                Ip(Irange,j) = [1.0; 1.0];
            end
        end
        
        
        % ec_interp1.m
        %
        %  Form the 1D linear interpolation matrix for edge-centered variables.
        %
        % input  -- Nf = number of cell centers on the fine grid
        % output -- Ip
        %
        function Ip = ec_interp1(obj,Nf)
            
            % number of centers on the coarse grid
            %
            Nc = Nf/2;
            
            % numbers of fine and course unknowns
            %
            Nfe = Nf;
            Nce = Nc;
            
            % allocate space for the matrix
            %
            Ip = spalloc(Nfe,Nce,3*Nce);
            
            % loop over columns
            %
            for j=1:Nce
                jf = 2*j-1;  % the corresponding fine node.
                % the contribution of this coarse node to nearby fine nodes
                % this is for period MAC grids
                if j==1
                    rows = [Nfe, jf, jf+1];
                else
                    rows = [jf-1:jf+1];
                end
                Ip( rows, j ) = [1/2; 1; 1/2];
            end
        end
        
        
        % gs_lex.m
        %
        %   Perform n steps of gs lex.
        %
        %   input  -- u = initial value of solutions
        %             A = matrix
        %             f = right hand side
        %             n = number of sweeps to perfrom
        %   output -- v = result after gslex
        %
        function v = gs_lex(obj,u,A,f,n)
            
            D  = diag(diag(A));
            Ut = triu(A) - D;
            DL = A - Ut;
            
            v = u;
            for k=1:n
                v = DL\(f - Ut*v);
            end
        end
        
        
        % gs_lex_sym.m
        %
        %   Perform n steps of gs lex -- forwards and backwards, symmetrix
        %
        %   input  -- u = initial value of solutions
        %             A = matrix
        %             f = right hand side
        %             n = number of sweeps to perfrom
        %   output -- v = result after gslex
        %
        function v = gs_lex_sym(obj,u,A,f,n)
            
            D  = diag(diag(A));
            Ut = triu(A) - D;
            DL = A - Ut;
            
            Lt = tril(A) - D;
            DU = A - Lt;
            
            v = u;
            for k=1:n
                v = DL\(f - Ut*v);
                v = DU\(f - Lt*v);
            end
        end
        
        
        % interp_matrix_p.m
        %
        %   Form the interpolation matrix for the pressure
        %   which is stored at the cell centers
        %
        %   input -- Nx number of cells in x-direction, fine grid
        %            Ny number of cells in the y-direction, fine grid
        %   output
        %
        function Ip = interp_matrix_p(obj,Nx,Ny,nVoffset)
            
            % form two one-dimensional interpolation matricies
            %
            if abs(nVoffset)==1
                Ipx = obj.nc_nm_interp1(Nx, nVoffset);
                Ipy = obj.nc_nm_interp1(Ny, nVoffset);
            elseif nVoffset==0
                Ipx = obj.cc_nm_interp1(Nx);
                Ipy = obj.cc_nm_interp1(Ny);
            else
                error('unsupported nVoffset!');
            end
            
            % assemble the matrix for the two-dimensional problem
            %   note the x-y ordering in kron corresponds to the order in meshgrid
            %
            Ip = kron(Ipx,Ipy);            
        end
        
        
        % interp_matrix_u.m
        %
        %   Form the interpolation matrix for the horizontal velocity (u)
        %   which is stored at the vertical edges
        %
        %   input -- Nx number of cells in x-direction, fine grid
        %            Ny number of cells in the y-direction, fine grid
        %   output
        %
        %
        % interp_matrix_p.m
        %
        %   Form the interpolation matrix for the pressure
        %   which is stored at the cell centers; use constant
        %   interpolation
        %
        %   input -- Nx number of cells in x-direction, fine grid
        %            Ny number of cells in the y-direction, fine grid
        %   output
        %
        function Ip = interp_matrix_p_const(obj,Nx,Ny)
            
            % form two one-dimensional interpolation matricies
            %
            Ipx = obj.cc_nm_interp1_const(Nx);
            Ipy = obj.cc_nm_interp1_const(Ny);
            
            % assemble the matrix for the two-dimensional problem
            %   note the y-x ordering in kron, this correpsonds to
            %   the case when the first index is for the x-coordinate
            %
            Ip = kron(Ipy,Ipx);
        end
        
        function Iu = interp_matrix_u(obj,Nx,Ny)
            
            % form two one-dimensional interpolation matricies
            %
            Ipx = obj.ec_interp1(Nx);
            Ipy = obj.cc_interp1(Ny);
            
            % assemble the matrix for the two-dimensional problem
            %   note the y-x ordering in kron, this correpsonds to
            %   the case when the first index is for the x-coordinate
            %
            Iu = kron(Ipx,Ipy);
            
        end
        
        %
        % interp_matrix_u_6p.m
        %
        %   Form the interpolation matrix for the horizontal velocity (u)
        %   which is stored at the vertical edges -- uses a six point
        %   interpolation scheme which is the adjoint of restriction
        %
        %   input -- Nx number of cells in x-direction, fine grid
        %            Ny number of cells in the y-direction, fine grid
        %   output
        %
        function Iu = interp_matrix_u_6pt(obj,Nx,Ny)
            
            % form two one-dimensional interpolation matricies
            %
            Ipx = obj.ec_interp1(Nx);
            Ipy = obj.cc_nm_interp1_const(Ny);
            
            % assemble the matrix for the two-dimensional problem
            %   note the x-y ordering in kron correpsonds to the ordering in meshgrid
            %
            Iu = kron(Ipx,Ipy);                       
        end
        
        
        % interp_matrix_v.m
        %
        %   Form the interpolation matrix for the vertical velocity (v)
        %   which is stored at the horizontal edges
        %
        %   input -- Nx number of cells in x-direction, fine grid
        %            Ny number of cells in the y-direction, fine grid
        %   output
        %
        function Iv = interp_matrix_v(obj,Nx,Ny)
            
            % form two one-dimensional interpolation matricies
            %
            Ipx = obj.cc_interp1(Nx);
            Ipy = obj.ec_interp1(Ny);
            
            % assemble the matrix for the two-dimensional problem
            %   note the y-x ordering in kron, this correpsonds to
            %   the case when the first index is for the x-coordinate
            %
            Iv = kron(Ipx,Ipy);
            
        end
        
        
        % interp_matrix_v_6pt.m
        %
        %   Form the interpolation matrix for the vertical velocity (v)
        %   which is stored at the horizontal edges -- uses a six point
        %   interpolation scheme which is the adjoint of restriction
        %
        %   input -- Nx number of cells in x-direction, fine grid
        %            Ny number of cells in the y-direction, fine grid
        %   output
        %
        function Iv = interp_matrix_v_6pt(obj,Nx,Ny)
            
            % form two one-dimensional interpolation matricies
            %
            Ipx =  obj.cc_nm_interp1_const(Nx);
            Ipy = obj.ec_interp1(Ny);
            
            % assemble the matrix for the two-dimensional problem
            %   note the x-y ordering in kron correpsonds to the ordering in meshgrid
            %
            Iv = kron(Ipx,Ipy);
            
        end
        
        
        %
        % nc_nm_interp1.m
        %
        %  Form the 1D linear interpolation matrix for a node centered grid.
        %
        % input -- Nf = number of nodes on the fine grid
        % output --
        %
        function Ip = nc_nm_interp1(obj,Nf, nVoffset)
            
            if Nf==3
                if nVoffset==1 % this
                    Ip = [1 0; 1/2 1/2; 0 1];
                elseif nVoffset==-1
                    Ip = [1/2; 1; 1/2];
                end
            else
                
                Nc = (Nf-nVoffset)/2+nVoffset; % number of variables on the coarse grid
                
                Ip = spalloc(Nf,Nc,3*Nc); % allocate space for the matrix
                
                if nVoffset==1
                    
                    Ip(1:2,1) = [1; 0.5]; % first column,  assume homogeneous Dirichlet BC
                    
                    Ip(Nf-1:Nf,Nc) = [0.5; 1]; % last column, assume homogeneous Dirichlet BC
                    
                    for j=2:Nc-1 % loop over columns to set the interpolation operator.
                        Ip(2*j-2:2*j,j) = [0.5; 1; 0.5];
                    end
                elseif nVoffset==-1
                    for j=1:Nc % assume homogeneous Dirichlet BC
                        Ip(2*j-1:2*j+1,j) = [0.5; 1; 0.5];
                    end
                end
            end
        end
        
    end
end
