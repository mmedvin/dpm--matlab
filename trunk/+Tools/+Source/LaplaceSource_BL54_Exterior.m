classdef LaplaceSource_BL54_Exterior < Tools.Source.SuperSource
    properties (Dependent = true)
        Source;%Fn;Ff;Fnn;Fff;    % WN;
    end
    
    properties(Access = protected)       
        Scatterer; 
		CoeffsClsrHndl;
		CoeffsParams;
		ExParams;
    end
    
    methods
        function [F,Fn,Ff,Fnn,Fff] = Derivatives(obj)
            
            if obj.IsDummy
                F=0;
                Fn=0;
                Ff=0;
                Fnn=0;
                Fff=0;
			else
				try
					eta = obj.Scatterer.Eta;
					phi = obj.Scatterer.Phi;
					fd  = obj.Scatterer.FocalDistance;
				catch exception
					if strcmp(exception.identifier,'MATLAB:nonExistentField')
						eta  = obj.Scatterer.eta;
						phi = obj.Scatterer.phi;
						fd  = obj.Scatterer.FocalDistance;						
					else
						rethrow(exception);
					end
				end
				
				x  = fd*cosh(eta).*cos(phi);
                y  = fd*sinh(eta).*sin(phi);
                xn = fd*sinh(eta).*cos(phi);
				xnn = x;
                yn = fd*cosh(eta).*sin(phi);	
				ynn=y;
                xf =-fd*cosh(eta).*sin(phi);
				xff=-x;
                yf = fd*sinh(eta).*cos(phi);
				yff=-y;

                c = obj.ExParams.c;
                d = obj.ExParams.d;

                F   = 56*(x.^c).*(y.^(d-2)) + 72*(x.^(c-2)).*(y.^d);
				
                if nargout>1, 
                    Fx = c*(d-1)*d*(x.^(c-1)).*(y.^(d-2)) + (c-2)*(c-1)*c*(x.^(c-3)).*(y.^d);
                    Fy = (d-2)*(d-1)*d*(x.^c).*(y.^(d-3)) + (c-1)*c*d*(x.^(c-2)).*(y.^(d-1));
                    Fn = Fx.*xn + Fy.*yn;
                end
                if nargout>2, Ff = Fx.*xf + Fy.*yf;  end
                %if nargout>3, Fnn =  end
                %if nargout>4, Fff =  end
            end
            
        end
            
        function obj = LaplaceSource_BL54_Exterior(Scatterer, CoeffsClsrHndl,CoeffsParams,ExParams)
			obj.Scatterer = Scatterer;
			obj.CoeffsClsrHndl = CoeffsClsrHndl;
			obj.CoeffsParams = CoeffsParams;
			obj.ExParams = ExParams;
			obj.IsDummy = false;
		end
        
		function S = get.Source(obj)
 			S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
			
			%[F,Fn,~,Fnn] = obj.Derivatives();
            [F,Fn] = obj.Derivatives();
            
            try
                eta = obj.Scatterer.Eta;
                phi = obj.Scatterer.Phi;
                fd  = obj.Scatterer.FocalDistance;
            catch exception
                if strcmp(exception.identifier,'MATLAB:nonExistentField')
                    eta  = obj.Scatterer.eta;
                    phi = obj.Scatterer.phi;
                    fd  = obj.Scatterer.FocalDistance;
                else
                    rethrow(exception);
                end
            end
            
            x  = fd*cosh(eta).*cos(phi);
            y  = fd*sinh(eta).*sin(phi);
						
			u = (x.^obj.ExParams.c) .* (y.*obj.ExParams.d);
			
			S(obj.Scatterer.Outside) = F(obj.Scatterer.Outside);
            S(1:end,1)=	u(1:end,1);
            S(1,1:end)= u(1,1:end);
            S(1:end,end)= u(1:end,end);
            S(end,1:end)= u(end,1:end);
			
			S(obj.Scatterer.GridGamma)	= F(obj.Scatterer.GridGamma) ...
										+ obj.Scatterer.deta.*Fn(obj.Scatterer.GridGamma) ; %...  
										%+ (obj.Scatterer.deta.^2).*Fnn(obj.Scatterer.GridGamma)/2;%taylor
			
			%%tmp = obj.Derivatives();
			%S(obj.Scatterer.Inside) = F(obj.Scatterer.Inside);   %was obj.Source(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
		end
	end
end

