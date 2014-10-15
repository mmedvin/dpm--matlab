classdef LaplacianOp4thOrder<Tools.DifferentialOps.SuperDiffOp
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A;
    end
    
    methods
        
        function LaplacianOp4thOrder(ParamsStruct)
            obj.Grid=ParamsStruct.Grid;

			%%%% Create Coeffs
			GridK = Tools.Grid.CartesianGrid( ...
				obj.Grid.x1 - obj.Grid.dx/2 , ...
				obj.Grid.xn + obj.Grid.dx/2 , ...
				2*obj.Grid.Nx + 1           , ...
				obj.Grid.y1 - obj.Grid.dy/2 , ...
				obj.Grid.yn + obj.Grid.dy/2 , ...
				2*obj.Grid.Ny + 1         ) ;
			
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
        
        function ApplyBC(obj)
            
            rows = obj.Grid.Nx;
            cols = obj.Grid.Ny;
            
            rc=rows*cols;
            
            obj.A(1:cols,:)         = sparse(1:cols,1:cols,1,rows,rc);			% Ic; %1,j
            obj.A(1:cols:rc,:)      = sparse(1:cols,1:cols:rc,1,rows,rc);		% i,1
            obj.A((rc-cols+1):rc,:) = sparse(1:cols,(rc-cols+1):rc,1,rows,rc);	% n,j
            obj.A(cols:cols:rc,:)   = sparse(1:cols,cols:cols:rc,1,rows,rc);	% i,n
            
        end
        
        function rhs = Bf(obj,F)
            
            row = obj.Grid.Nx;
            column = obj.Grid.Ny;
            
            rhs = zeros(row*column);
            n=1;
            for j=2:row+1
                for i=2:column+1
                    
                    rhs(n)=F(j,i)+h^2/12*((1/h^2-obj.a_x(x(i),y(j))/(2*obj.a(x(i),y(j))*h))*F(j,i+1) ...
                          +(1/h^2+obj.a_x(x(i),y(j))/(2*obj.a(x(i),y(j))*h))*F(j,i-1)+(1/h^2-obj.b_y(x(i),y(j))/(2*obj.b(x(i),y(j))*h))*F(j+1,i) ...
                          +(1/h^2+obj.b_y(x(i),y(j))/(2*obj.b(x(i),y(j))*h))*F(j-1,i)+(2*obj.b_y(x(i),y(j))^2/obj.b(x(i),y(j))^2-3*obj.b_yy(x(i),y(j))/(2*obj.b(x(i),y(j))) ...
                          +2*obj.a_x(x(i),y(j))^2/obj.a(x(i),y(j))^2-3*obj.a_xx(x(i),y(j))/(2*obj.a(x(i),y(j)))-4/h^2)*F(j,i));
                    n=n+1;
                end
            end       
        end
        
        function LaplacianA(obj)
            
            %%%%%%%%%%%%%
            h = obj.Grid.dx;
            assert(h == obj.Grid.dy);
						
            row = obj.Grid.Nx;
            column = obj.Grid.Ny;
           
            rc=row*column;
            %%%%%%%%%%%%%            
            
            L=zeros(rc,rc);
            k=0;
            
            n=1;
            for j=2:row+1
                for i=2:column+1
                    
                    L(n,n)=-(obj.a(x(i)+h/2,y(j))+obj.a(x(i)-h/2,y(j))+obj.b(x(i),y(j)+h/2)+obj.b(x(i),y(j)-h/2))/(h^2) ...
                        +k^2-h^2/12*(-4*obj.b(x(i),y(j))/h^4+2*k^2/h^2-4*obj.a(x(i),y(j))/h^4+2*k^2/h^2 ...
                        +(2*obj.b_xx(x(i),y(j))+2*obj.a_yy(x(i),y(j))-2*obj.a_x(x(i),y(j))*obj.b_x(x(i),y(j))/obj.a(x(i),y(j)) ...
                        -2*obj.b_y(x(i),y(j))*obj.a_y(x(i),y(j))/obj.b(x(i),y(j)))/h^2+(2*obj.a_x(x(i),y(j))^2/obj.a(x(i),y(j))^2 ...
                        -3*obj.a_xx(x(i),y(j))/(2*obj.a(x(i),y(j))))*(-k^2+2*obj.b(x(i),y(j))/h^2)+(2*obj.b_y(x(i),y(j))^2/obj.b(x(i),y(j))^2 ...
                        -3*obj.b_yy(x(i),y(j))/(2*obj.b(x(i),y(j))))*(-k^2+2*obj.a(x(i),y(j))/h^2));
                    
                    if mod(n,row)~=1 % excludes left column of interior points
                        
                        L(n,n-1)=obj.a(x(i)-h/2,y(j))/h^2-h^2/12*(-obj.a_x(x(i),y(j))/h^3+2*obj.a(x(i),y(j))/h^4 ...
                            +2*obj.b(x(i),y(j))/h^4-k^2/h^2+obj.b(x(i),y(j))*obj.a_x(x(i),y(j))/(obj.a(x(i),y(j))*h^3) ...
                            -obj.a_x(x(i),y(j))*k^2/(2*obj.a(x(i),y(j))*h)+obj.a_xxx(x(i),y(j))/(4*h) ...
                            -5*obj.a_xx(x(i),y(j))*obj.a_x(x(i),y(j))/(4*obj.a(x(i),y(j))*h) ...
                            +obj.a_x(x(i),y(j))^3/(obj.a(x(i),y(j))^2*h)+obj.a_xyy(x(i),y(j))/(2*h)-obj.a_yy(x(i),y(j))/h^2-2*obj.b_x(x(i),y(j))/h^3 ...
                            -obj.b_y(x(i),y(j))*obj.a_xy(x(i),y(j))/(2*obj.b(x(i),y(j))*h)+obj.b_y(x(i),y(j))*obj.a_y(x(i),y(j))/(obj.b(x(i),y(j))*h^2) ...
                            +(2*obj.b_y(x(i),y(j))^2/obj.b(x(i),y(j))^2-3*obj.b_yy(x(i),y(j))/(2*obj.b(x(i),y(j))))*(obj.a_x(x(i),y(j))/(2*h) ...
                            -obj.a(x(i),y(j))/h^2));
                        
                        if n>row
                            L(n,n-column-1)=-h^2/12*(obj.b_y(x(i),y(j))/(2*h^3)-obj.b(x(i),y(j))/h^4 ...
                                +obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.a(x(i),y(j))*h^2) ...
                                -obj.a_x(x(i),y(j))*obj.b(x(i),y(j))/(2*obj.a(x(i),y(j))*h^3)+obj.a_x(x(i),y(j))/(2*h^3)-obj.a(x(i),y(j))/h^4 ...
                                +obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.b(x(i),y(j))*h^2) ...
                                -obj.a(x(i),y(j))*obj.b_y(x(i),y(j))/(2*obj.b(x(i),y(j))*h^3) ...
                                -(obj.b_xy(x(i),y(j))+obj.a_xy(x(i),y(j)))/(2*h^2)+(obj.b_x(x(i),y(j))+obj.a_y(x(i),y(j)))/h^3);
                        end
                        
                        if n<=row*(column-1)
                            L(n,n+column-1)=-h^2/12*(-obj.b_y(x(i),y(j))/(2*h^3)-obj.b(x(i),y(j))/h^4 ...
                                -obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.a(x(i),y(j))*h^2) ...
                                -obj.a_x(x(i),y(j))*obj.b(x(i),y(j))/(2*obj.a(x(i),y(j))*h^3)+obj.a_x(x(i),y(j))/(2*h^3)-obj.a(x(i),y(j))/h^4 ...
                                -obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.b(x(i),y(j))*h^2) ...
                                +obj.a(x(i),y(j))*obj.b_y(x(i),y(j))/(2*obj.b(x(i),y(j))*h^3) ...
                                +(obj.b_xy(x(i),y(j))+obj.a_xy(x(i),y(j)))/(2*h^2)-(-obj.b_x(x(i),y(j))+obj.a_y(x(i),y(j)))/h^3);
                        end
                        
                    end
                    
                    if mod(n,row)~=0 % exludes right column of interior points
                        
                        L(n,n+1)=obj.a(x(i)+h/2,y(j))/h^2-h^2/12*(obj.a_x(x(i),y(j))/h^3+2*obj.a(x(i),y(j))/h^4+2*obj.b(x(i),y(j))/h^4 ...
                            -k^2/h^2-obj.b(x(i),y(j))*obj.a_x(x(i),y(j))/(obj.a(x(i),y(j))*h^3)+obj.a_x(x(i),y(j))*k^2/(2*obj.a(x(i),y(j))*h) ...
                            -obj.a_xxx(x(i),y(j))/(4*h)+5*obj.a_xx(x(i),y(j))*obj.a_x(x(i),y(j))/(4*obj.a(x(i),y(j))*h)-obj.a_x(x(i),y(j))^3/(obj.a(x(i),y(j))^2*h) ...
                            -obj.a_xyy(x(i),y(j))/(2*h)-obj.a_yy(x(i),y(j))/h^2+2*obj.b_x(x(i),y(j))/h^3+obj.b_y(x(i),y(j))*obj.a_xy(x(i),y(j))/(2*obj.b(x(i),y(j))*h) ...
                            +obj.b_y(x(i),y(j))*obj.a_y(x(i),y(j))/(obj.b(x(i),y(j))*h^2)+(2*obj.b_y(x(i),y(j))^2/obj.b(x(i),y(j))^2 ...
                            -3*obj.b_yy(x(i),y(j))/(2*obj.b(x(i),y(j))))*(-obj.a_x(x(i),y(j))/(2*h)-obj.a(x(i),y(j))/h^2));
                        
                        if n>row
                            L(n,n-column+1)=-h^2/12*(obj.b_y(x(i),y(j))/(2*h^3)-obj.b(x(i),y(j))/h^4-obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.a(x(i),y(j))*h^2) ...
                                +obj.a_x(x(i),y(j))*obj.b(x(i),y(j))/(2*obj.a(x(i),y(j))*h^3)-obj.a_x(x(i),y(j))/(2*h^3)-obj.a(x(i),y(j))/h^4 ...
                                -obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.b(x(i),y(j))*h^2)-obj.a(x(i),y(j))*obj.b_y(x(i),y(j))/(2*obj.b(x(i),y(j))*h^3) ...
                                +(obj.b_xy(x(i),y(j))+obj.a_xy(x(i),y(j)))/(2*h^2)-(obj.b_x(x(i),y(j))-obj.a_y(x(i),y(j)))/h^3);
                        end
                        
                        if n<=row*(column-1)
                            L(n,n+column+1)=-h^2/12*(-obj.b_y(x(i),y(j))/(2*h^3)-obj.b(x(i),y(j))/h^4+obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.a(x(i),y(j))*h^2) ...
                                +obj.a_x(x(i),y(j))*obj.b(x(i),y(j))/(2*obj.a(x(i),y(j))*h^3)-obj.a_x(x(i),y(j))/(2*h^3)-obj.a(x(i),y(j))/h^4 ...
                                +obj.a_x(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.b(x(i),y(j))*h^2)+obj.a(x(i),y(j))*obj.b_y(x(i),y(j))/(2*obj.b(x(i),y(j))*h^3) ...
                                -(obj.b_xy(x(i),y(j))+obj.a_xy(x(i),y(j)))/(2*h^2)-(obj.b_x(x(i),y(j))+obj.a_y(x(i),y(j)))/h^3);
                        end
                        
                    end
                    
                    if n>row % excludes bottom row of interior points
                        L(n,n-column)=obj.b(x(i),y(j)-h/2)/h^2-h^2/12*(-obj.b_y(x(i),y(j))/h^3+2*obj.b(x(i),y(j))/h^4+2*obj.a(x(i),y(j))/h^4-k^2/h^2 ...
                            +obj.a(x(i),y(j))*obj.b_y(x(i),y(j))/(obj.b(x(i),y(j))*h^3)-obj.b_y(x(i),y(j))*k^2/(2*obj.b(x(i),y(j))*h)+obj.b_yyy(x(i),y(j))/(4*h) ...
                            -5*obj.b_yy(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.b(x(i),y(j))*h)+obj.b_y(x(i),y(j))^3/(obj.b(x(i),y(j))^2*h)+obj.b_yxx(x(i),y(j))/(2*h) ...
                            -obj.b_xx(x(i),y(j))/h^2-2*obj.a_y(x(i),y(j))/h^3-obj.a_x(x(i),y(j))*obj.b_xy(x(i),y(j))/(2*obj.a(x(i),y(j))*h) ...
                            +obj.a_x(x(i),y(j))*obj.b_x(x(i),y(j))/(obj.a(x(i),y(j))*h^2)+(2*obj.a_x(x(i),y(j))^2/obj.a(x(i),y(j))^2 ...
                            -3*obj.a_xx(x(i),y(j))/(2*obj.a(x(i),y(j))))*(obj.b_y(x(i),y(j))/(2*h)-obj.b(x(i),y(j))/h^2));
                    end
                    
                    if n<=row*(column-1) % excludes top row of interior points
                        L(n,n+column)=obj.b(x(i),y(j)+h/2)/h^2-h^2/12*(obj.b_y(x(i),y(j))/h^3+2*obj.b(x(i),y(j))/h^4+2*obj.a(x(i),y(j))/h^4-k^2/h^2 ...
                            -obj.a(x(i),y(j))*obj.b_y(x(i),y(j))/(obj.b(x(i),y(j))*h^3)+obj.b_y(x(i),y(j))*k^2/(2*obj.b(x(i),y(j))*h)-obj.b_yyy(x(i),y(j))/(4*h) ...
                            +5*obj.b_yy(x(i),y(j))*obj.b_y(x(i),y(j))/(4*obj.b(x(i),y(j))*h)-obj.b_y(x(i),y(j))^3/(obj.b(x(i),y(j))^2*h)-obj.b_yxx(x(i),y(j))/(2*h) ...
                            -obj.b_xx(x(i),y(j))/h^2+2*obj.a_y(x(i),y(j))/h^3+obj.a_x(x(i),y(j))*obj.b_xy(x(i),y(j))/(2*obj.a(x(i),y(j))*h) ...
                            +obj.a_x(x(i),y(j))*obj.b_x(x(i),y(j))/(obj.a(x(i),y(j))*h^2)+(2*obj.a_x(x(i),y(j))^2/obj.a(x(i),y(j))^2 ...
                            -3*obj.a_xx(x(i),y(j))/(2*obj.a(x(i),y(j))))*(-obj.b_y(x(i),y(j))/(2*h)-obj.b(x(i),y(j))/h^2));
                    end
                    
                    n=n+1;
                end
            end  
            
            [i,j,s] = find(L);
            obj.A = sparse(i,j,s,rc,rc);            
        end
        
    end
    
end

