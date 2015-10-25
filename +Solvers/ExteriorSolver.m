classdef ExteriorSolver < Solvers.SuperHomoSolver

    properties(AbortSet = true)%Access = protected, AbortSet = true)
        A1;
        A2;
        k;
  %      w0;
  %      w1;
  %      f;
    end
    
    methods
        function obj = ExteriorSolver(Arguments)                      
            obj = obj@Solvers.SuperHomoSolver(Arguments);
            
            obj.k=Arguments.CoeffsParams.k;%supposed to be constant obj.Coeffs.k;
            
            obj.CreateDirectOperatorA1();
            obj.CreateSolverA2();
            
        end
        
        function u = P_Omega(obj,xi_gamma,Uinc)
            
            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Nm));
            
			w = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.GridGamma));
			
            if ~exist('Uinc','var')
                %Uinc=zeros(obj.Grid.Nx,obj.Grid.Ny);%(size(obj.w0));
				 Uinc=spalloc(obj.Grid.Nx,obj.Grid.Ny,0);
				%do nothing
			else
				w(obj.GridGamma) = Uinc(obj.GridGamma);
            end
            
            %obj.w0(obj.GridGamma) = Uinc(obj.GridGamma);
            
            %GLW = obj.Solve(xi_gamma - obj.w0);
			GLW = obj.Solve(xi_gamma - w);
            
            u(obj.Scatterer.Nm)=GLW(obj.Scatterer.Nm) + Uinc(obj.Scatterer.Nm);
        end
        
        function Qj = Qcol2(obj,w,dw)
            
            if isobject(w)%numel(w) == numel(obj.Scatterer.BasisArg)
                %w, dw expected to be uinc and it's normal derivative
                tmp = spalloc(obj.Grid.Nx,obj.Grid.Ny,length(obj.GridGamma));
                tmp(obj.GridGamma)=obj.Scatterer.Expansion(w,dw,obj.NoSource,obj.Coeffs);
                w=tmp;
            elseif numel(w) == numel(obj.GridGamma)%~= numel(obj.w0)
                tmp = spalloc(obj.Grid.Nx,obj.Grid.Ny,length(obj.GridGamma));
                tmp(obj.GridGamma)=w;
                w=tmp;
            end
            
            GLW = obj.Solve(w);            
            Qj = obj.Qcol(GLW,w);
        end
    end
    
    methods(Access = protected)
        
        function f = Lu(obj,u,msk)
            f = obj.A1*u;
            
            if exist('msk','var');
                f = f(msk,:);            
            end
%             if exist('msk','var');
%                 f = obj.A1(msk,:)*u;
%             else
%                 f = obj.A1*u;
%             end
        end
                        
        function u = Gf(obj,f)
            rhs = fft(full(f),[],2);
            
            y = obj.A2\rhs(:);
            y = reshape(y,obj.Grid.Nx,obj.Grid.Ny);
            
            u = ifft(y,[],2);  % Backward Fourier transform
        end
        
        function Qj = Qcol(obj,GLW,w)
            Qj = GLW(obj.GridGamma)-w(obj.GridGamma);
        end
        
        
%         function u = Solve(obj,x)            
%              obj.f(obj.Scatterer.Mp) = obj.Lu(x(:),obj.Scatterer.Mp);
%              %obj.Truncate(f);
%              u = obj.Gf(obj.f);                        
%         end
        
%         function calc_QnW(obj)
%             for j = 1:obj.Basis.NBss
%                
%                 [obj.w0(obj.GridGamma),obj.w1(obj.GridGamma)] = obj.ExpandedBasis(obj.Basis.Indices(j)) ;
%                 
%                 obj.myW0(obj.GridGamma,j) = obj.w0(obj.GridGamma);
%                 obj.myW1(obj.GridGamma,j) = obj.w1(obj.GridGamma);
%                 
% %                 GLW0 = obj.Solve(obj.w0);
% %                 GLW1 = obj.Solve(obj.w1);
% %                 
% %                 obj.myQ0(:,j) = obj.Qcol(GLW0,obj.w0);
% %                 obj.myQ1(:,j) = obj.Qcol(GLW1,obj.w1);
% 
%                 obj.myQ0(:,j) = obj.Qcol2(obj.w0);
%                 obj.myQ1(:,j) = obj.Qcol2(obj.w1);
% 
% 
%             end
%             obj.IsReadyQnW = true;
%         end
        
        
        
        
        function CreateDirectOperatorA1(obj)
                hr  = obj.Grid.dx;
                hth = obj.Grid.dy;
                r = obj.Grid.x;
                nr=obj.Grid.Nx;
                nth=obj.Grid.Ny;
                
                hr2=hr^2;
                hth2=hth^2;
                k2=obj.k^2;
                
                
                cp = 1+hr/2./r;
                cm = 1-hr/2./r;
                b  = (5/3+hr2/6./r.^2) ;
                                                
                rip1 = r + hr;
                rim1 = r - hr;
                %     d  = 1/hth2./rip1 - 1./r.^2 - k2/2;
                
                Amn     = sparse( - 5/3/hr2 - b/hth2./r.^2 +  (2/3+hr2/12./r.^2)*k2 ) ;
                Amp1n   = sparse( 5*cp/6/hr2 + k2/12 - 1/6/hth2./rip1.^2 - hr*(1/hth2./rip1.^2 - 1./r.^2 - k2/2)./12./r ) ;
                Amm1n   = sparse( 5*cm/6/hr2 + k2/12 - 1/6/hth2./rim1.^2 + hr*(1/hth2./rim1.^2 - 1./r.^2 - k2/2)./12./r ) ;
                Amnpm1  = sparse( b/2/hth2./r.^2 + k2/12 - 1/6/hr2 ) ;
                Amp1npm1= sparse( cp.*(1/hth2./rip1.^2 + 1./hr2)/12 ) ;
                Amm1npm1= sparse( cm.*(1/hth2./rim1.^2 + 1./hr2)/12 ) ;
                
                
                
                
                %     Rm1 = sparse(2:nr,1:nr-1,1,nr,nr);
                %     Rp1=Rm1';
                %     R = Rm1+Rp1;
                
%                 Tm1 = sparse(2:nth,1:nth-1,1,nth,nth);
%                 Tp1=Tm1';
%                 T = Tm1+Tp1;
%                 T(1,nth)=1;
%                 T(nth,1)=1;
                
                %         Ir = speye(nr);
                It = speye(nth);
                
                nn = nr*nth;
                
                
                Amm1npm1 =repmat([0 Amm1npm1(2:end)],1,nth-1);
                Amp1npm1 =repmat([Amp1npm1(1:end-1) 0],1,nth-1);
                
                obj.A1 = kron(It,diag(Amn)) + kron(It,diag(Amp1n(1:end-1),1) + diag(Amm1n(2:end),-1) ) ...
                    + diag( repmat(Amnpm1,1,nth-1),-nr ) + diag( repmat(Amnpm1,1,nth-1),nr ) ...
                    + diag( Amnpm1,nn-nr) + diag( Amnpm1,-(nn-nr)) ... % periodicity
                    + diag( Amm1npm1(2:end), -nr-1) + diag( [ Amm1npm1 0], nr-1) ...
                    + diag( Amm1npm1(2:nr), -nn+nr-1) + diag( [ Amm1npm1(1:nr) 0], nn-nr-1) ...% periodicity
                    + diag( [0,Amp1npm1], -nr+1) + diag( Amp1npm1(1:end-1), nr+1) ...
                    + diag( [0,Amp1npm1(1:nr)], -(nn-nr-1)) + diag( Amp1npm1(1:nr-1), nn-nr+1) ... % periodicity
                    ;
                
        end
            
        function CreateSolverA2(obj)
            r = obj.Grid.x;
            theta = obj.Grid.y;
            ak = obj.k;
            
            ak2 = ak^2;
            rho = r(2:end-1).';
            r_out=r(end);
            
            hr = abs(r(1)-r(2));
            htheta = abs(theta(1)-theta(2));
            nr = numel(r);
            ntheta = numel(theta);
            
            eig= - (2*sin(theta/2)/htheta).^2;
            
            anu2=(-eig./(1 + eig*htheta^2/12)).';
            anu=sqrt(anu2);
            han = zeros(numel(anu),3);
%             ierr = zeros(numel(anu),3);
            
            
            
            hr2 = hr^2;
            
            % The sub-diagonal
            a = repmat((rho - hr/2)./rho/hr2 - hr2*(ak2*(1./rho./hr./2 - 1/hr2) + 1./hr./rho.^3 )/12,1,ntheta) ...
                + ((1 - hr./2./rho).*(htheta^2/12/hr2 + 1./12./((rho - hr).^2)))*eig;
            
            % The super-diagonal
            c = repmat((rho + hr/2)./rho/hr2 + hr2*(ak2*(1./rho./hr./2 + 1/hr2) + 1./hr./rho.^3 )/12,1,ntheta) ...
                + ((1 + hr./2./rho).*(htheta^2/12/hr2 + 1./12./((rho + hr).^2)))*eig;
            
            % The diagonal
            b = repmat(ak2*(hr2./2./rho.^2 + 5)/6 - 2/hr2,1,ntheta) ...
                + (hr2./2./rho.^4 + 5./rho.^2 + (ak2/2 - 1/hr2)*htheta^2 )*(eig/6);
            
            asympt=sqrt(anu./(anu+1)).*(anu+1.)/((r_out-hr/2)*ak);
            han(1:2,:)          = besselh([[-1;0],[0;1],[1;2]],2,ak*(r_out-hr/2));
            han(3:end,:)    = besselh([anu(3:end)-1,anu(3:end),anu(3:end)+1],2,ak*(r_out-hr/2));
            ierr=isinf(han);
            tmperr = sum(ierr,2);
            
            alpha = zeros(ntheta-1,1);
            m=tmperr==0;
            alpha(m) = 0.5*ak*(han(m,1)-han(m,3)) ./han(m,2) ;
            
            m = find(tmperr>0 & ak*asympt>30);
            alpha(m)=ak*(-asympt(m));
            
            if any((tmperr>0 & ak*asympt<=30)),error('grrr'),end
            
            c1=1./hr*(1 - hr2/24.*( (anu2+2)/(r_out-hr/2)^2-ak2) - alpha*hr2/8/(r_out-hr/2));
            c2=1./2.*(alpha+hr^2/24.*(-3.*anu2/(r_out-hr/2)^3+ak2/(r_out-hr/2))-alpha.*hr2/8.*(anu2/(r_out-hr/2)^2-ak2));
            beta=(c1+c2)./(c1-c2);
            
            a=[zeros(1,ntheta); a; -beta.'];
            b=[ones(1,ntheta) ; b; ones(1,ntheta)];
            c=[zeros(1,ntheta) ; c ; zeros(1,ntheta)];
            
            n=nr*(ntheta);
            
            obj.A2 = sparse(2:n,1:n-1,a(2:end),n,n,3*n) + sparse(1:n,1:n,b,n,n,3*n) + sparse(1:n-1,2:n,c(1:end-1),n,n,3*n);
            
        end
        
        
        
    end
    
    
end
