classdef ConstLapCoeffs < Tools.Coeffs.AbstractCoeffs
	properties
		a;
		b;
		sigma;
	end
	
	%     methods(Static = true)
	%         function [k,kr] = kkr(r,r0,k0)
	%             k=obj.k;
	%             kr=0;
	%         end
	%     end
	
	
	methods
		
		function varargout = Derivatives(obj,WhichOne)
			switch WhichOne
				case 'a'
					varargout(1)={obj.a};
				case 'b'
					varargout(1)={obj.b};
				case {'sigma','s'}
					varargout(1)={obj.sigma};
				otherwise
					error('unknown coeff')
			end
			
			nout = nargout;
			varargout(2:nout)={0};
		end
		
		function obj=ConstLapCoeffs(~,Params) %(k0,r,r0)
			obj.a       = Params.a;
			obj.b       = Params.b;
			obj.sigma   = Params.sigma;
		end
	end
end
