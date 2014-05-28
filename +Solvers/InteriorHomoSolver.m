classdef InteriorHomoSolver < Solvers.SuperHomoSolver
   
    properties(Access = protected, AbortSet = true)
        A;
        k;
    end
    
    methods
        function obj = InteriorHomoSolver( ...
                BasisIndices,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams,CollectRhs)
            obj = obj@Solvers.SuperHomoSolver( ...
                BasisIndices,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams,CollectRhs);
            
            if  numel(obj.Coeffs.k)>1
                GridK = Tools.Grid.CartesianGrid( ...
                    obj.Grid.x1 - obj.Grid.dx , ...
                    obj.Grid.xn + obj.Grid.dx , ...
                    obj.Grid.Nx + 2           , ...
                    obj.Grid.y1 - obj.Grid.dy , ...
                    obj.Grid.yn + obj.Grid.dy , ...
                    obj.Grid.Ny + 2         ) ;
                
                [X,Y] = GridK.mesh();
                
                ScattK = struct('r',abs(X+1i.*Y));%ScattererClsHandle(GridK,obj.ScattererAddParams);
                WNPlr= Tools.Coeffs.WaveNumberPolarR(ScattK,WaveNumberAddParams);
                obj.k = sparse(WNPlr.k);
            else
                obj.k = obj.Coeffs.k.*ones(obj.Grid.Size+2);
            end
            
            obj.HlmSemA();
            
        end
        
        function u = P_Omega(obj,xi_gamma)
            
            rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
            
            rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
            GLW = obj.Gf(rhs(:));
            
            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));
            u(obj.Scatterer.Np)=xi_gamma(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).';                        
        end
    end
    
    methods(Access = protected)
      
        function f = Lu(obj,u,msk)
            if exist('msk','var');
                f = obj.A(msk,:)*u;
            else
                f = obj.A*u;
            end
        end
        
        function u = Gf(obj,f)
           u = obj.A\f;           
        end
                
        function Qj = Qcol(obj,GLW,~)
            Qj = -GLW(obj.GridGamma,:);
        end
        
        function HlmSemA(obj)%(x,y,k)
            
            %  global   rc Ic Ir Ym1 Yp1 %Xm1 Xp1
            %
            %             dx=x(2)-x(1);
            %             dy=y(2)-y(1);
            
            dx = obj.Grid.dx;
            dy = obj.Grid.dy;
            
            dx2=dx^2;
            dy2=dy^2;
            %             rows  = length(x);
            %             cols = length(y);
            rows = obj.Grid.Nx;
            cols = obj.Grid.Ny;
            
            rc=rows*cols;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  verify  input
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %k  = VerifyK(k, cols,rows);
            
            
            [Ic,Ir,Ym1,Yp1] = InitialDerivationMatrixA();%(k,dx,dy,rows,cols);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %clean up
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % A(1:cols,:)=0;
            % A((rc-cols+1):rc,:) = 0;
           % obj.A(1:cols:rc,:) = 0;
           % obj.A(cols:cols:rc,:) = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %external boundary
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            c1 = (1/dx2+1/dy2);
            cx = 1/dx2-c1/6;
            cy = 1/dy2-c1/6;
            
            % % % % % % % % % % % % % % % %
            %1,j
            % % % % % % % % % % % % % % % %
            k0j = sparse(diag(obj.k(2:end-1,1)));
            k1j = sparse(diag(obj.k(2:end-1,2)));
            k2j = sparse(diag(obj.k(2:end-1,3)));
            
            alpha =12*Ic+dx2.*(k0j.^2);
            alpha = alpha/( dy2.* (6.*Ic+dx2.*(k1j.^2) ) );
            
            a=1i.*k1j;%sparse(ones(size()));
            b=sparse(0*(1./(2*diag(a))));
            
            Ahat1j = 2.*diag(b)./dx-(a.*dy2)/dx + dy2.*(k1j.^2).*(k2j-k0j)./6;
            A1j = (2/3).*(k1j.^2) -(5/3).*c1.*Ic ;
            A1jm1 = Ym1.*cy + diag((obj.k(2:end-2,1).^2)./12  ,-1)  - diag( (alpha(2:end,2:end)*b(2:end)/dx) , -1);
            A1jp1 = Yp1.*cy + diag((obj.k(3:end-1,1).^2)./12,1) - diag((alpha(1:end-1,1:end-1)*b(1:end-1)/dx),1);
            Aipm1jpm1 = c1/12;
            
            obj.A(1:cols,1:cols)= A1j + alpha.*Ahat1j  + A1jm1  + A1jp1 ;	 %(1,j-1), (1,j),(1,j+1)
            obj.A(1:cols,(cols+1):2*cols)= 2.*cx.*Ic+ (k0j.^2+k2j.^2)./12 + 2.*Aipm1jpm1.*(Ym1+Yp1);%(i+1=2,j)
            
            % A(1:cols,1:(2*cols)) = [A1j + alpha.*Ahat1j  + A1jm1  + A1jp1,2.*cx.*Ic+ (k0j.^2+k2j.^2)./12 + 2.*Aipm1jpm1.*(Ym1+Yp1)];
            
            
            obj.A([1,cols],1:2*cols)=0;%clean up of corners, we put there Dirichlet condition atm
            
            clear k0j k1j k2j tmp
            
            
            % % % % % % % % % % % % % % % %
            %   i,1
            % % % % % % % % % % % % % % % %
           % sparse(1:cols:rc,1:cols:rc,Ir,cols,cols)
            
            %obj.A(1:cols:rc,1:cols:rc) = Ir;
            obj.A(1:cols:rc,:) = sparse(1:cols,1:cols:rc,1,rows,rc);
            
            % % % % % % % % % % % % % % % %
            % n,j
            % % % % % % % % % % % % % % % %
            knp1j = sparse(diag(obj.k(2:end-1,end)));
            knj		= sparse(diag(obj.k(2:end-1,end-1)));
            knm1j = sparse(diag(obj.k(2:end-1,end-2)));
            
            alpha =12*Ic+dx2.*(knp1j.^2);
            alpha = alpha/( dy2.* (6.*Ic+dx2.*(knj.^2) ) );
            
            a=-1i.*knj;%sparse(ones(size(knj)));
            b=sparse(0*(1./(2*diag(a))));
            
            Ahatnj = 2.*diag(b)./dx-(a.*dy2)/dx + dy2.*(knj.^2).*(knp1j-knm1j)./6;
            Anj = (2/3).*(knj.^2) -(5/3).*c1.*Ic ;
            Anjm1 = Ym1.*cy + diag((obj.k(2:end-2,end).^2)./12  ,-1)  + diag( (alpha(2:end,2:end)*b(2:end)/dx) , -1);
            Anjp1 = Yp1.*cy + diag((obj.k(3:end-1,end).^2)./12,1) + diag((alpha(1:end-1,1:end-1)*b(1:end-1)/dx),1);
            Aipm1jpm1 = c1/12;
            
            obj.A((rc-cols+1):rc,(rc-cols+1):rc)= Anj - alpha.*Ahatnj  + Anjm1  + Anjp1 ;	 %(1,j-1), (1,j),(1,j+1)
            obj.A((rc-cols+1):rc,(rc-2*cols+1):(rc-cols))= 2.*cx.*Ic+ (knp1j.^2+knm1j.^2)./12 + 2.*Aipm1jpm1.*(Ym1+Yp1);%(i+1=2,j)
            obj.A([(rc-cols+1),rc],(rc-2*cols+1):rc)=0;%clean up of corners, we put there Dirichlet condition atm
            
            clear knp1j knj knm1j tmp
            
            
            % % % % % % % % % % % % % % % %
            %  i,n
            % % % % % % % % % % % % % % % %
%             obj.A(cols:cols:rc,cols:cols:rc) = Ir;
            obj.A(cols:cols:rc,:) = sparse(1:cols,cols:cols:rc,1,rows,rc);
            
            % corners
            obj.A([1,rc-cols+1],[1,rc-cols+1]) =Ir([1,end],[1,end]);
            obj.A([cols,rc],[cols,rc]) = Ir([1,end],[1,end]);
            
            
            function [Ic,Ir,Ym1,Yp1] = InitialDerivationMatrixA()%(k,dx,dy,rows,cols)
                % This function return  matrix A - inner scheme for Helmholtz equation
                
                %  global rc  Ic Ir Ym1 Yp1 Xm1 Xp1
                
                dx2 = dx.^2;
                dy2 = dy.^2;
                
                Ir = speye(rows);
                Ic = speye(cols);
                
                Ym1 = sparse(2:cols,1:cols-1,1,cols,cols);
                Yp1=Ym1';
                Y = Ym1+Yp1;
                
                Xm1 = sparse(2:rows,1:rows-1,1,rows,rows);
                Xp1=Xm1';
                X = Xm1+Xp1;
                
                II = speye(rc);
                
                k2 = sparse(obj.k).^2;
                
                c1 = (1/dx2+1/dy2);
                cx = 1/dx2-c1/6;
                cy = 1/dy2-c1/6;
                
                k2ij=k2(2:end-1,2:end-1);
                Aij = 2/3.*diag(k2ij(:))-(5/3)*c1*II;
                
                k2im1j=k2(2:end-1,1:end-2);
                Aim1j = cx.*kron(Xm1,Ic) + kron(Xm1,Ic) .*diag(k2im1j(cols+1:rc),-cols)/12;
                
                k2ip1j=k2(2:end-1,3:end);
                Aip1j = cx.*kron(Xp1,Ic) + kron(Xp1,Ic).*diag(k2ip1j(1:rc-cols),cols)/12;
                
                k2ijm1=k2(1:end-2,2:end-1);
                Aijm1 = cy.* kron(Ir,Ym1) + kron(Ir,Ym1).*diag(k2ijm1(2:end),-1)/12;
                
                k2ijp1=k2(3:end,2:end-1);
                Aijp1 = cy.* kron(Ir,Yp1) + kron(Ir,Yp1) .*diag(k2ijp1(1:end-1),1)/12;
                
                Aipm1jpm1 = c1/12;
                
                obj.A = Aim1j+Aijm1+Aij+Aijp1+Aip1j+kron(X,Y)*Aipm1jpm1;
            end
            
        end
        
    end
    
    
end
