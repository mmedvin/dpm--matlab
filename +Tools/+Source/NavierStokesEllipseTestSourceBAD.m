classdef NavierStokesTestSource < Tools.Source.SuperSource
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
				
                Coeff = obj.CoeffsClsrHndl(obj.Scatterer,obj.CoeffsParams);
                
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

                cd = c^2 - d^2;
            
                F   = cd*(s + cd)*exp(c*x).*sin(d*y);
				
%                 if nargout>1, 
%                     Fx = ;
%                     Fy = ;
%                     Fn = Fx.*xn + Fy.*yn;
%                 end
%                 if nargout>2, Ff = Fx.*xf + Fy.*yf;  end
%                 if nargout>3, 
%                     Fxx = ;
%                     Fyy = ;
%                     Fxy = ;
%                     Fnn = Fxx.*(xn.^2) + Fyy.*(yn.^2) + 2*Fxy.*xn.*yn + Fx.*xnn + Fy.*ynn;
%                 end
%                 if nargout>4,                     
%                     Fff =  Fxx.*(xf.^2) + Fyy.*(yf.^2) + 2*Fxy.*xf.*yf + Fx.*xff + Fy.*yff;
%                 end
            end
            
        end
            
        function obj = NavierStokesTestSource(Scatterer, CoeffsClsrHndl,CoeffsParams,ExParams)
			obj.Scatterer = Scatterer;
			obj.CoeffsClsrHndl = CoeffsClsrHndl;
			obj.CoeffsParams = CoeffsParams;
			obj.ExParams = ExParams;
			obj.IsDummy = false;
		end
        
		function S = get.Source(obj)
 			S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
			            
            eta  = obj.Scatterer.Eta0;
            phi = obj.Scatterer.phi;
            fd  = obj.Scatterer.FocalDistance;
            
            c = obj.ExParams.c;
            d = obj.ExParams.d;
            
            s = Coeff.sigma;
            
            cd = c^2 - d^2;
            
            F   = cd*(s + cd)*exp(c*x).*sin(d*y);
            S(obj.Scatterer.Np) = F(obj.Scatterer.Np);
            
%  			
% 			S(obj.Scatterer.GridGamma)	= F + obj.Scatterer.deta.*Fn + (obj.Scatterer.deta.^2).*Fnn/2;%taylor
% 			                                   
% 			Fin = obj.Derivatives();
% 			S(obj.Scatterer.Inside) = Fin(obj.Scatterer.Inside);   %was obj.Source(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
		end
	end
end

