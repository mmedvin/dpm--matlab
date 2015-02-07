classdef LapOp4OrdrVarCoeffBCinRhs<Tools.DifferentialOps.LaplacianOpBCinRhs
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
    
    methods (Access=public)
        function obj = LapOp4OrdrVarCoeffBCinRhs(ParamsStruct) 
            obj = obj@Tools.DifferentialOps.LaplacianOpBCinRhs(ParamsStruct);	
        end
        				
        function Rhs = AdjustRhs(obj,Rhs,Exact,BC_yn,BC_x1,BC_xn)
            
            %assuming zero BC
            assert(obj.BC_y1); assert(obj.BC_yn); assert(obj.BC_x1); assert(obj.BC_xn); 
            return
            error('TBD')           
            Nx = obj.Grid.Nx;
            Ny = obj.Grid.Ny;
			
			dx2 = obj.Grid.dx.^2;
			dy2 = obj.Grid.dy.^2;
			
            
            % -(((axm2h+axmh) uxm3h)/(36 hx^2))
            % -((axp2h+axph) uxp3h)/(36 hx^2) ...            
            % -(uxm2h (axm2h-8 axmh+a[x,y]))/(24 hx^2) ...
            % -(uxp2h (axp2h-8 axph+a[x,y]))/(24 hx^2) ...
            % +(uxmh (axm2h+2 axmh+axph+2 a[x,y]))/(12 hx^2) ...
            % +(uxph (axmh+axp2h+2 axph+2 a[x,y]))/(12 hx^2) ...
            % -((axm2h+40 axmh+axp2h+40 axph+18 a[x,y]) u[x,y])/(72 hx^2)

									
			Rhs(2:end-1	, 2			) = Rhs(2:end-1	, 2			) - a_iphalf_j.* obj.BC_y1./dx2;%m
			Rhs(2:end-1	, end-1		) = Rhs(2:end-1	, end-1		) - a_imhalf_j.* obj.BC_yn./dx2;%p
			Rhs(2		,  2:end-1	) = Rhs(2		, 2:end-1	) - b_i_jphalf.* obj.BC_x1./dy2; %m
			Rhs(end-1	,  2:end-1	) = Rhs(end-1	, 2:end-1	) - b_i_jmhalf.* obj.BC_xn./dy2; %p
            
                        
			%a_iphalf_j = obj.Coeffs.a(3:2:end,2);
			%a_imhalf_j = obj.Coeffs.a(3:2:end,end-1);
			
			%b_i_jphalf = obj.Coeffs.b(2    ,3:2:end);
			%b_i_jmhalf = obj.Coeffs.b(end-1,3:2:end);
									
			%Rhs(2:end-1	, 2			) = Rhs(2:end-1	, 2			) - a_iphalf_j.* obj.BC_y1./dx2;%m
			%Rhs(2:end-1	, end-1		) = Rhs(2:end-1	, end-1		) - a_imhalf_j.* obj.BC_yn./dx2;%p
			%Rhs(2		,  2:end-1      ) = Rhs(2		, 2:end-1	) - b_i_jphalf.* obj.BC_x1./dy2; %m
			%Rhs(end-1	,  2:end-1      ) = Rhs(end-1	, 2:end-1	) - b_i_jmhalf.* obj.BC_xn./dy2; %p

			


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
            
            % -(((axm2h+axmh) uxm3h)/(36 hx^2))
            % -((axp2h+axph) uxp3h)/(36 hx^2) ...            
            % -(uxm2h (axm2h-8 axmh+a[x,y]))/(24 hx^2) ...
            % -(uxp2h (axp2h-8 axph+a[x,y]))/(24 hx^2) ...
            % +(uxmh (axm2h+2 axmh+axph+2 a[x,y]))/(12 hx^2) ...
            % +(uxph (axmh+axp2h+2 axph+2 a[x,y]))/(12 hx^2) ...
            % -((axm2h+40 axmh+axp2h+40 axph+18 a[x,y]) u[x,y])/(72 hx^2)
             
            axm2h = obj.Coeffs.a(3:end-2, 1:end-4);
            axmh  = obj.Coeffs.a(3:end-2, 2:end-3);
            aij   = obj.Coeffs.a(3:end-2, 3:end-2);  
            axph  = obj.Coeffs.a(3:end-2, 4:end-1); 
            axp2h = obj.Coeffs.a(3:end-2, 5:end  );

            bym2h = obj.Coeffs.b(1:end-4, 3:end-2);
            bymh  = obj.Coeffs.b(2:end-3, 3:end-2);
            bij   = obj.Coeffs.b(3:end-2, 3:end-2);
            byph  = obj.Coeffs.b(4:end-1, 3:end-2);
            byp2h = obj.Coeffs.b(5:end  , 3:end-2);

            
            Aij = -( (axm2h(:) + 40*axmh(:) + axp2h(:) + 40*axph(:) + 18*aij(:))./(72*dx2) + (bym2h(:) + 40*bymh(:) + byp2h(:) + 40*byph(:) + 18*bij(:))./(72*dy2)  + obj.Coeffs.Sigma(:) );
            
            Aip3j = [ zeros(1, 3*cols)                                                                          ,-(axp2h(1:end-3*cols) + axph(1:end-3*cols))/(36*dx2)                                       ].'; %(1:end-3*cols)
            Aip2j = [ zeros(1, 2*cols)                                                                          ,-(axp2h(1:end-2*cols) - 8*axph(1:end-2*cols) + aij(1:end-2*cols))/(24*dx2)                 ].'; %(1:end-2*cols)
			Aip1j = [ zeros(1,   cols)                                                                          , (axmh(1:end-cols) + axp2h(1:end-cols) + 2*axph(1:end-cols) + 2*aij(1:end-cols))./(12*dx2) ].'; %(1:end-cols)
            Aim1j = [ (axm2h(cols+1:end) + 2*axmh(cols+1:end) + axph(cols+1:end) + 2*aij(cols+1:end))./(12*dx2) , zeros(1,   cols)                                                                          ].'; %(cols+1:end)
            Aim2j = [-(axm2h(2*cols+1:end) - 8*axmh(2*cols+1:end) + aij(2*cols+1:end))/(24*dx2)                 , zeros(1, 2*cols)                                                                          ].'; %(2*cols+1:end)
            Aim3j = [-(axm2h(3*cols+1:end) + axmh(3*cols+1:end))/(36*dx2)                                       , zeros(1, 3*cols)                                                                          ].'; %(3*cols+1:end)
			
            Aijp3 = [ 0 0 0                                                                 ,-(byp2h(1:end-3) + byph(1:end-3))/(36*dy2)                                     ].'; %(1:end-3)
            Aijp2 = [ 0 0                                                                   ,-(byp2h(1:end-2) - 8*byph(1:end-2) + bij(1:end-2))/(24*dy2)                    ].'; %(1:end-2)
            Aijp1 = [ 0                                                                     , (bymh(1:end-1) + byp2h(1:end-1) + 2*byph(1:end-1) + 2*bij(1:end-1))./(12*dy2) ].'; %(1:end-1)
            Aijm1 = [ (bym2h(2:end) + 2*bymh(2:end) + byph(2:end) + 2*bij(2:end))./(12*dy2) , 0                                                                             ].'; %(2:end)
            Aijm2 = [-(bym2h(3:end) - 8*bymh(3:end) + bij(3:end))/(24*dy2)                  , 0 0                                                                           ].'; %(3:end)
            Aijm3 = [-(bym2h(4:end) + bymh(4:end))/(36*dy2)                                 , 0 0 0                                                                         ].'; %(4:end) 

%             Aijp3 = [ 0 0 0                                                                 ,-(byp2h(4:end) + byph(4:end))/(36*dy2)                                     ].'; %(1:end-3)
%             Aijp2 = [ 0 0                                                                   ,-(byp2h(3:end) - 8*byph(3:end) + bij(3:end))/(24*dy2)                    ].'; %(1:end-2)
%             Aijp1 = [ 0                                                                     , (bymh(2:end) + byp2h(2:end) + 2*byph(2:end) + 2*bij(2:end))./(12*dy2) ].'; %(1:end-1)
%             Aijm1 = [ (bym2h(1:end-1) + 2*bymh(1:end-1) + byph(1:end-1) + 2*bij(1:end-1))./(12*dy2) , 0                                                                             ].'; %(2:end)
%             Aijm2 = [-(bym2h(1:end-2) - 8*bymh(1:end-2) + bij(1:end-2))/(24*dy2)                  , 0 0                                                                           ].'; %(3:end)
%             Aijm3 = [-(bym2h(1:end-3) + bymh(1:end-3))/(36*dy2)                                 , 0 0 0                                                                         ].'; %(4:end) 


            
            Aijm3(cols  :cols:rc)=0;
            Aijm3(cols-1:cols:rc)=0;
            Aijm3(cols-2:cols:rc)=0;
            %
            Aijm2(cols-1:cols  :rc)=0;
            Aijm2(cols  :cols  :rc)=0;
            %
            Aijm1(cols  :cols:rc)=0;
            Aijp1(cols+1:cols:rc)=0;
            %
            Aijp2(cols+1:cols:rc)=0;
            Aijp2(cols+2:cols:rc)=0;
            %
            Aijp3(cols+1:cols:rc)=0;
            Aijp3(cols+2:cols:rc)=0;
            Aijp3(cols+3:cols:rc)=0;
            
            
            obj.A = spdiags([Aim3j,Aim2j,Aim1j,Aijm3,Aijm2,Aijm1,Aij,Aijp1,Aijp2,Aijp3,Aip1j,Aip2j,Aip3j],[-3*cols,-2*cols,-cols,-3,-2,-1,0,1,2,3,cols,2*cols,3*cols], rc,rc);
            
            

            
            
			%a_iphalf_j = obj.Coeffs.a(2:2:end-1,3:2:end)./dx2;
			%a_imhalf_j = obj.Coeffs.a(2:2:end-1,1:2:end-1)./dx2;
			
			%b_i_jphalf = obj.Coeffs.b(3:2:end,2:2:end-1)./dy2;
			%b_i_jmhalf = obj.Coeffs.b(1:2:end-1,2:2:end-1)./dy2;           
            
			%Aip1j = [zeros(1,cols)          ,a_iphalf_j(1:end-cols) ].';
			%Aim1j = [a_imhalf_j(cols+1:end) ,zeros(1,cols)          ].';
			%Aijp1 = [0                      ,b_i_jphalf(1:end-1)    ].';
			%Aijm1 = [b_i_jmhalf(2:end)      ,0                      ].';
			
			%Aijm1(cols:cols:rc)=0;
			%Aijp1(cols+1:cols:rc)=0;
			
			%obj.A = spdiags([Aim1j,Aijm1,Aij,Aijp1,Aip1j],[-cols,-1,0,1,cols], rc,rc);
            
            
            %B=full(obj.A); %dbg
		end
          
        function VerifyCoeffs(obj,coeffs)
            if numel(coeffs.a)==1
                obj.Coeffs.a = coeffs.a*ones( obj.Grid.Nx + 4, obj.Grid.Ny + 4 );
            else
                obj.Coeffs.a = coeffs.a;
            end
            if numel(coeffs.b)==1
                obj.Coeffs.b =  coeffs.b*ones( obj.Grid.Nx + 4, obj.Grid.Ny + 4 ) ;
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
                                    
            if     any(size(obj.Coeffs.a) ~= [obj.Grid.Nx+4 , obj.Grid.Ny+4]) || ...
                   any(size(obj.Coeffs.b) ~= [obj.Grid.Nx+4 , obj.Grid.Ny+4])                   
                error('non constant coefficients a and b should be evaluated at x1-2dx:dx:xn+2dx, y1-2dy:dy:yn+2dy')
            end
            
            
            obj.Coeffs.a = sparse(obj.Coeffs.a);
            obj.Coeffs.b = sparse(obj.Coeffs.b);
            obj.Coeffs.Sigma = sparse(obj.Coeffs.Sigma);
        end
        
    end
 
    
end

