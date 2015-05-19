classdef SuperLaplacianOp<Tools.DifferentialOps.SuperDiffOp
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
    properties(Access = public)
        A;
    end
    properties(Access = protected)
        Coeffs;
        Grid;
    end
    
    methods (Access=public)
        
        function obj = SuperLaplacianOp(ParamsStruct)

            if nargin>0
                obj.Grid=ParamsStruct.Grid;
			%%%% Create Coeffs
            
            if ~isfield(ParamsStruct,'Order')
                ParamsStruct.Order=2;
            end
            
            if ParamsStruct.Order==4
                GridK = Tools.Grid.CartesianGrid( ...
                    obj.Grid.x1 - 2*obj.Grid.dx , ...
                    obj.Grid.xn + 2*obj.Grid.dx , ...
                    obj.Grid.Nx + 4           , ...
                    obj.Grid.y1 - 2*obj.Grid.dy , ...
                    obj.Grid.yn + 2*obj.Grid.dy , ...
                    obj.Grid.Ny + 4         ) ;
            else
			GridK = Tools.Grid.CartesianGrid( ...
				obj.Grid.x1 - obj.Grid.dx/2 , ...
				obj.Grid.xn + obj.Grid.dx/2 , ...
				2*obj.Grid.Nx + 1           , ...
				obj.Grid.y1 - obj.Grid.dy/2 , ...
				obj.Grid.yn + obj.Grid.dy/2 , ...
				2*obj.Grid.Ny + 1         ) ;
            end
            
			
            if isfield(ParamsStruct.CoeffsParams,'FocalDistance')

                [Eta,Phi] = GridK.ToElliptical(ParamsStruct.CoeffsParams.FocalDistance);
                ScattK = struct('FocalDistance',ParamsStruct.CoeffsParams.FocalDistance, 'Eta',Eta,'Phi',Phi);
            else
                ScattK = struct('r',GridK.R);%,'th',GridK.Theta);
            end
			 coeffs = ParamsStruct.CoeffsHandle(ScattK,ParamsStruct.CoeffsParams);
			%%%%
			
            obj.VerifyCoeffs(coeffs);
            obj.LaplacianA();
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
            rhs = F(:);
            % do nothing
        end
    end
    methods (Access=protected)
        function LaplacianA(obj)%(x,y,k)
            dx = obj.Grid.dx;
            dy = obj.Grid.dy;
			
			dx2 = dx.^2;
			dy2 = dy.^2;
			
            rows = obj.Grid.Nx;
            cols = obj.Grid.Ny;
           
            rc=rows*cols;
			
			a_iphalf_j = obj.Coeffs.a(2:2:end-1,3:2:end)./dx2;
			a_imhalf_j = obj.Coeffs.a(2:2:end-1,1:2:end-1)./dx2;
			
			b_i_jphalf = obj.Coeffs.b(3:2:end,2:2:end-1)./dy2;
			b_i_jmhalf = obj.Coeffs.b(1:2:end-1,2:2:end-1)./dy2;
			
			
			Aij = -( (a_iphalf_j(:) + a_imhalf_j(:))+ (b_i_jphalf(:) + b_i_jmhalf(:)) + obj.Coeffs.Sigma(:) );
			
			Aip1j = [zeros(1,cols)          ,a_iphalf_j(1:end-cols) ].';
			Aim1j = [a_imhalf_j(cols+1:end) ,zeros(1,cols)          ].';
			Aijp1 = [0                      ,b_i_jphalf(1:end-1)    ].';
			Aijm1 = [b_i_jmhalf(2:end)      ,0                      ].';
			
			Aijm1(cols:cols:rc)=0;
			Aijp1(cols+1:cols:rc)=0;
			
			obj.A = spdiags([Aim1j,Aijm1,Aij,Aijp1,Aip1j],[-cols,-1,0,1,cols], rc,rc);
		end
          
        function VerifyCoeffs(obj,coeffs)
            if numel(coeffs.a)==1
                obj.Coeffs.a=  coeffs.a*ones( 2*obj.Grid.Nx + 1, 2*obj.Grid.Ny + 1 ) ;
            else
                obj.Coeffs.a = coeffs.a;
            end
            if numel(coeffs.b)==1
                obj.Coeffs.b=  coeffs.b*ones( 2*obj.Grid.Nx + 1, 2*obj.Grid.Ny + 1 ) ;
            else
                obj.Coeffs.b = coeffs.b;
            end
            
            if numel(coeffs.sigma)==1
                obj.Coeffs.Sigma=  coeffs.sigma*ones( obj.Grid.Nx, obj.Grid.Ny) ;
            else 
                obj.Coeffs.Sigma  = coeffs.sigma;
            end
            
            
            if any(size(obj.Coeffs.Sigma) ~= [obj.Grid.Nx, obj.Grid.Ny])
                error('non constant sigma should be evaluated at x1:xn, y1:yn')
            end
                                    
            if     any(size(obj.Coeffs.a) ~= [2*obj.Grid.Nx+1 , 2*obj.Grid.Ny+1]) || ...
                   any(size(obj.Coeffs.b) ~= [2*obj.Grid.Nx+1 , 2*obj.Grid.Ny+1])                   
                error('non constant coefficients a and b should be evaluated at x1-dx/2:dx:xn+dx/2, y1-dy/2:dy:yn+dy/2')
            end
            
            
%             obj.Coeffs.a = sparse(obj.Coeffs.a);
%             obj.Coeffs.b = sparse(obj.Coeffs.b);
%             obj.Coeffs.Sigma = sparse(obj.Coeffs.Sigma);
        end
        
    end
    
end

