classdef LaplaceSource < Tools.Source.SuperSource
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
					r = obj.Scatterer.R;
				catch exception
					if strcmp(exception.identifier,'MATLAB:nonExistentField')
						r = obj.Scatterer.r;
					else
						rethrow(exception);
					end
				end
                F   = 8*r.^2 + 4;
                
                if nargout>1, Fn  = 16*r;  end
                if nargout>2, Ff  = 0   ;  end
                if nargout>3, Fnn = 16.*ones(size(r))  ;  end
                if nargout>4, Fff = 0   ;  end
            end
            
        end
            
        function obj = LaplaceSource(Scatterer, CoeffsClsrHndl,CoeffsParams,ExParams)
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
			Exact	= Tools.Exact.ExLapCrclVarCoeffs(obj.Scatterer, obj.ExParams);
			
			
			
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

