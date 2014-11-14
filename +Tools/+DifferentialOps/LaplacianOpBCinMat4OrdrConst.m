classdef LaplacianOpBCinMat4OrdrConst <Tools.DifferentialOps.SuperLaplacianOp
    
    properties
        BC_x1=0;
        BC_xn=0;
        BC_y1=0;
        BC_yn=0;
    end
    
        methods (Access=public)
        
        function obj = LaplacianOpBCinMat4OrdrConst(ParamsStruct)
			%%%% Create Coeffs
            obj.Grid=ParamsStruct.Grid;
            obj.Coeffs = ParamsStruct.CoeffsHandle([],ParamsStruct.CoeffsParams);
            
            obj.BC_x1 = ParamsStruct.BC_x1;
            obj.BC_xn = ParamsStruct.BC_xn;
            obj.BC_y1 = ParamsStruct.BC_y1;
            obj.BC_yn = ParamsStruct.BC_yn;
            
            obj.LaplacianA();
        end
		
        function Rhs = AdjustRhs(obj,Rhs,Exact)
			if exist('Exact','var')			
                Nx = obj.Grid.Nx;
                Ny = obj.Grid.Ny;
                
                Rhs(1:Nx,  1      ) = Exact(1:Nx,  1      );
                Rhs(1:Nx,  Ny     ) = Exact(1:Nx,  Ny     );
                Rhs(1   ,  1: Ny  ) = Exact(1   ,  1:Ny   );
                Rhs(Nx  ,  1: Ny  ) = Exact(Nx  ,  1:Ny   );
            else 
                %warning('not verified part')
                Rhs(:	, 1		) = obj.BC_y1;%m
                Rhs(:   , end	) = obj.BC_yn;%p
                Rhs(1	, :     ) = obj.BC_x1; %m
                Rhs(end	, :     ) = obj.BC_xn; %p
			end
        end
        
        
		function b = ApplyOp(obj,x,mask)
			% returns matrix-vector multiplication result, i.e b = A*x
			% if parameter 'mask' is esixts, a mask 'mask' is applied on the returned value b,
			% more precisely b(mask) is returned, i.e. only the indices in 'mask'
			
			if exist('mask','var');
				b = obj.A(mask,:)*x;
			else
				b = obj.A*x;
			end
		end
		
		function u = Solve(obj,f)
			u = obj.A\f;
        end
        
        
        function rhs = Bf(obj,F)
            
            dx = obj.Grid.dx;
            dy = obj.Grid.dy;
            dx2=dx^2;
            dy2=dy^2;
            rows = obj.Grid.Nx;
            cols = obj.Grid.Ny;
            rc=rows*cols;
            
            [Ic,B]...%,Dx1,Dxn]
            = InitialDerivationMatrix();
            
            %%%%%
            Fxeq1 = F(1:end ,   1       );%F(2:2:(Rows-1)  ,   2           );
            Fxeqn = F(1:end ,   end-1   );%F(2:2:(Rows-1)  ,   end-1       );
            Fyeq1 = F(1       , end-1   );%F(2             ,   2:2:(Cols-1));
            Fyeqn = F(end-1   , end-1   );%F(end-1         ,   2:2:(Cols-1));
            %%%%%
            
            %F=sparse(F(:));
            rhs = sparse(B*F(:)); return
            %Fx1 = Dx1*F;
            %Fxn = Dxn*F;
            % Fy1=Dy1*F;
            % Fyn=Dyn*F;
            %clear B Dx1 Dxn Dy1 Dyn
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %external boundary
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % % % % % % % % % % % % % % % %
            %1,j
            % % % % % % % % % % % % % % % %
            
            rhs(1:cols) = Fyeq1;
            
            % % % % % % % % % % % % % % % %
            %   i,1
            % % % % % % % % % % % % % % % %
            rhs(1:cols:rc)=Fxeq1;
            
            % % % % % % % % % % % % % % % %
            % n,j
            % % % % % % % % % % % % % % % %
            
            rhs((rc-cols+1):rc) = Fyeqn;
            
            % % % % % % % % % % % % % % % %
            %  i,n
            % % % % % % % % % % % % % % % %

            rhs(cols:cols:rc)=Fxeqn;
            
            
            function [Ic,B]...%,Dx1,Dxn,Dy1,Dyn]
                = InitialDerivationMatrix()%(dx,dy,rows,cols)
                % This function return 3 matrices,
                % matrix A is inner scheme for Helmholtz equation
                % matrix B is inner scheme for inhomogeneous part F of eq  Hu=F
                % matrix C is additional derivative on F for BC
                % we assuming here that F is given on finer grid
                
                
                %%%%%%%%%%%%%%%%%%%%%
                % Create A
                %%%%%%%%%%%%%%%%%%%%%
                
                % Ir = speye(rows);
                Ic = speye(cols);
                
                %%%%%%%%%%%%%%%%%%%%%
                % Create B
                %%%%%%%%%%%%%%%%%%%%%
                
                %Cols=2*cols+1;
                %Rows=2*rows+1;
                %RC=Rows*Cols;
                
                Cols=cols;
                Rows=rows;
                RC=rc;
                
                BIr = speye(Rows);
                BIc = speye(Cols);
                
                BYm1 = sparse(2:Cols,1:Cols-1,1,Cols,Cols);
                BY = BYm1+BYm1';
                
                BXm1 = sparse(2:Rows,1:Rows-1,1,Rows,Rows);
                BX = BXm1+BXm1';
                
                BII = speye(RC);
                B=(kron(BX,BIc) + kron(BIr,BY) + 8*BII)./12;
                
                % find indices related to equations for i,j
                % i.e. eliminate equations related to i+-1/2
%                 MM=[];
%                 for indx=Cols:2*Cols:(RC-2*Cols)
%                     MM=[MM, (indx+2):2:(indx+Cols)];
%                 end
%                 B=B(MM,:);% eliminate
                
%                 %%%%%%%%%%%%%%%%%%%%%
%                 % Create Dx1,Dxn,Dy1,Dyn
%                 %%%%%%%%%%%%%%%%%%%%%
%                 T=kron(BXm1'-BXm1,BIc);
%                 Dx1=T(MM(1:cols),:)./dx;% eliminate
%                 Dxn=T(MM(rc-cols+1:rc),:)./dx;
%                 
%                 T=kron(BIr,BYm1'-BYm1);
%                 Dy1=(T(MM(1:cols:rc),:)./dy);% eliminate
%                 Dyn=(T(MM(cols:cols:rc),:)./dy);
                
            end
            
            
        end

        
        function rhs = OldBf(obj,F)
            
            dx = obj.Grid.dx;
            dy = obj.Grid.dy;
            dx2=dx^2;
            dy2=dy^2;
            rows = obj.Grid.Nx;
            cols = obj.Grid.Ny;
            rc=rows*cols;
            
            [Ic,B,Dx1,Dxn] = InitialDerivationMatrix();
            
            %%%%%
            Fxeq1 = F(2:2:(Rows-1)  ,   2           );
            Fxeqn = F(2:2:(Rows-1)  ,   end-1       );
            Fyeq1 = F(2             ,   2:2:(Cols-1));
            Fyeqn = F(end-1         ,   2:2:(Cols-1));
            %%%%%
            
            %F=sparse(F(:));
            rhs = sparse(B*F(:));
            %Fx1 = Dx1*F;
            %Fxn = Dxn*F;
            % Fy1=Dy1*F;
            % Fyn=Dyn*F;
            %clear B Dx1 Dxn Dy1 Dyn
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %external boundary
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % % % % % % % % % % % % % % % %
            %1,j
            % % % % % % % % % % % % % % % %
            
            rhs(1:cols) = Fxeq1;
            
            % % % % % % % % % % % % % % % %
            %   i,1
            % % % % % % % % % % % % % % % %
            rhs(1:cols:rc)=Fyeq1;
            
            % % % % % % % % % % % % % % % %
            % n,j
            % % % % % % % % % % % % % % % %
            
            rhs((rc-cols+1):rc) = Fxeqn;
            
            % % % % % % % % % % % % % % % %
            %  i,n
            % % % % % % % % % % % % % % % %

            rhs(cols:cols:rc)=Fyeqn;
            
            
            function [Ic,B,Dx1,Dxn,Dy1,Dyn] = InitialDerivationMatrix()%(dx,dy,rows,cols)
                % This function return 3 matrices,
                % matrix A is inner scheme for Helmholtz equation
                % matrix B is inner scheme for inhomogeneous part F of eq  Hu=F
                % matrix C is additional derivative on F for BC
                % we assuming here that F is given on finer grid
                
                
                %%%%%%%%%%%%%%%%%%%%%
                % Create A
                %%%%%%%%%%%%%%%%%%%%%
                
                % Ir = speye(rows);
                Ic = speye(cols);
                
                %%%%%%%%%%%%%%%%%%%%%
                % Create B
                %%%%%%%%%%%%%%%%%%%%%
                
                Cols=2*cols+1;
                Rows=2*rows+1;
                RC=Rows*Cols;
                
                BIr = speye(Rows);
                BIc = speye(Cols);
                
                BYm1 = sparse(2:Cols,1:Cols-1,1,Cols,Cols);
                BY = BYm1+BYm1';
                
                BXm1 = sparse(2:Rows,1:Rows-1,1,Rows,Rows);
                BX = BXm1+BXm1';
                
                BII = speye(RC);
                B=(kron(BX,BIc) + kron(BIr,BY) - BII)./3;
                
                % find indices related to equations for i,j
                % i.e. eliminate equations related to i+-1/2
                MM=[];
                for indx=Cols:2*Cols:(RC-2*Cols)
                    MM=[MM, (indx+2):2:(indx+Cols)];
                end
                B=B(MM,:);% eliminate
                
                %%%%%%%%%%%%%%%%%%%%%
                % Create Dx1,Dxn,Dy1,Dyn
                %%%%%%%%%%%%%%%%%%%%%
                T=kron(BXm1'-BXm1,BIc);
                Dx1=T(MM(1:cols),:)./dx;% eliminate
                Dxn=T(MM(rc-cols+1:rc),:)./dx;
                
                T=kron(BIr,BYm1'-BYm1);
                Dy1=(T(MM(1:cols:rc),:)./dy);% eliminate
                Dyn=(T(MM(cols:cols:rc),:)./dy);
                
            end
            
            
        end
       
    end
    methods (Access=protected)
        function LaplacianA(obj)
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
            
            
            [Ic,Ir,Ym1,Yp1] = InitialDerivationMatrixA();
            obj.A = obj.A.*obj.Coeffs.a;
            
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

           % A1j =  -(5/3).*c1.*Ic ;
           % A1jm1 = Ym1.*cy ;
           % A1jp1 = Yp1.*cy ;
           % Aipm1jpm1 = c1/12;
           
            obj.A(1:cols,1:2*cols)=sparse(1:cols,1:cols,1,cols,2*cols);
            %obj.A(1:cols,1:cols)= A1j + A1jm1  + A1jp1 ;	 %(1,j-1), (1,j),(1,j+1)
            %obj.A(1:cols,(cols+1):2*cols)= 2.*cx.*Ic + 2.*Aipm1jpm1.*(Ym1+Yp1);%(i+1=2,j)
            
            % A(1:cols,1:(2*cols)) = [A1j + alpha.*Ahat1j  + A1jm1  + A1jp1,2.*cx.*Ic+ (k0j.^2+k2j.^2)./12 + 2.*Aipm1jpm1.*(Ym1+Yp1)];
            
            
           % obj.A([1,cols],1:2*cols)=0;%clean up of corners, we put there Dirichlet condition atm
            
            
            % % % % % % % % % % % % % % % %
            %   i,1
            % % % % % % % % % % % % % % % %
           % sparse(1:cols:rc,1:cols:rc,Ir,cols,cols)
            
            %obj.A(1:cols:rc,1:cols:rc) = Ir;
            obj.A(1:cols:rc,:) = sparse(1:cols,1:cols:rc,1,rows,rc);
            
            % % % % % % % % % % % % % % % %
            % n,j
            % % % % % % % % % % % % % % % %
            
           % Anj =  -(5/3).*c1.*Ic ;
           % Anjm1 = Ym1.*cy ;
           % Anjp1 = Yp1.*cy ;
           % Aipm1jpm1 = c1/12;
            
           obj.A((rc-cols+1):rc,(rc-2*cols+1):rc)= sparse(1:cols,cols+1:2*cols,1,cols,2*cols);%sparse((rc-cols+1):rc,(rc-cols+1):rc,1,cols,2*cols);
            %obj.A((rc-cols+1):rc,(rc-cols+1):rc)=  Anj + Anjm1  + Anjp1 ;	 %(1,j-1), (1,j),(1,j+1)
            %obj.A((rc-cols+1):rc,(rc-2*cols+1):(rc-cols))= 2.*cx.*Ic + 2.*Aipm1jpm1.*(Ym1+Yp1);%(i+1=2,j)
            %obj.A([(rc-cols+1),rc],(rc-2*cols+1):rc)=0;%clean up of corners, we put there Dirichlet condition atm
            
            % % % % % % % % % % % % % % % %
            %  i,n
            % % % % % % % % % % % % % % % %
%             obj.A(cols:cols:rc,cols:cols:rc) = Ir;
            obj.A(cols:cols:rc,:) = sparse(1:cols,cols:cols:rc,1,rows,rc);
            
            % corners
           % obj.A([1,rc-cols+1],[1,rc-cols+1]) =Ir([1,end],[1,end]);
           % obj.A([cols,rc],[cols,rc]) = Ir([1,end],[1,end]);
            
            
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

                c1 = (1/dx2+1/dy2);
                cx = 1/dx2-c1/6;
                cy = 1/dy2-c1/6;
                
                Aij = -(5/3)*c1*II;
                
                Aim1j = cx.*kron(Xm1,Ic) ;
                
                Aip1j = cx.*kron(Xp1,Ic) ;
                
                Aijm1 = cy.* kron(Ir,Ym1) ;
                
                Aijp1 = cy.* kron(Ir,Yp1) ;
                
                Aipm1jpm1 = c1/12;
                
                obj.A = Aim1j+Aijm1+Aij+Aijp1+Aip1j+kron(X,Y)*Aipm1jpm1;
            end
            
        end
 
       
        
    end
end