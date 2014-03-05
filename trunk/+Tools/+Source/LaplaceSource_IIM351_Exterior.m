classdef LaplaceSource_IIM351_Exterior < Tools.Source.SuperSource
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

				coeffs = obj.CoeffsClsrHndl(obj.Scatterer,obj.CoeffsParams);
				B = coeffs.Derivatives(obj,'a');
                F   = -2*B*sin(x).*cos(y);
				
                if nargout>1, Fn  =  -2*B*cos(x).*cos(y).*xn +  2*B*sin(x).*sin(y).*yn;  end
                if nargout>2, Ff  =  -2*B*cos(x).*cos(y).*xf +  2*B*sin(x).*sin(y).*yf;  end
                if nargout>3, Fnn = 2*B*(cos(x).*(2*sin(y).*xn.*yn - cos(y).*xnn) + sin(x).*(cos(y).*(xn.^2+yn.^2) + sin(y).*ynn));  end
                if nargout>4, Fff = 2*B*(cos(x).*(2*sin(y).*xf.*yf - cos(y).*xff) + sin(x).*(cos(y).*(xf.^2+yf.^2) + sin(y).*yff));  end
            end
            
        end
            
        function obj = LaplaceSource_IIM351_Exterior(Scatterer, CoeffsClsrHndl,CoeffsParams,ExParams)
			obj.Scatterer = Scatterer;
			obj.CoeffsClsrHndl = CoeffsClsrHndl;
			obj.CoeffsParams = CoeffsParams;
			obj.ExParams = ExParams;
			obj.IsDummy = false;
		end
        
		function S = get.Source(obj)
 			S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
			
			[F,Fn,~,Fnn] = obj.Derivatives();
			
% 			Coeffs	= obj.CoeffsClsrHndl(obj.Scatterer.TheScatterer,obj.CoeffsParams);
			Exact	= Tools.Exact.ExLapElps351(obj.Scatterer, obj.ExParams);
			
			
			
			S(obj.Scatterer.Outside) = F(obj.Scatterer.Outside);
            S(1:end,1)=	Exact.u(1:end,1);
            S(1,1:end)= Exact.u(1,1:end);
            S(1:end,end)= Exact.u(1:end,end);
            S(end,1:end)= Exact.u(end,1:end);
			
			S(obj.Scatterer.GridGamma)	= F(obj.Scatterer.GridGamma) ...
										+ obj.Scatterer.dr.*Fn(obj.Scatterer.GridGamma) ...  
										+ (obj.Scatterer.dr.^2).*Fnn(obj.Scatterer.GridGamma)/2;%taylor
			
			%%tmp = obj.Derivatives();
			%S(obj.Scatterer.Inside) = F(obj.Scatterer.Inside);   %was obj.Source(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
		end
	end
end

