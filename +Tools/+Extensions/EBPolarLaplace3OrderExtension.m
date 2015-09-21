classdef EBPolarLaplace3OrderExtension < Tools.Extensions.SuperPolarTwoTupleExtension
methods

  function obj = EBPolarLaplace3OrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
        end
        
        function res = Expansion(obj,Xi0,Xi1,Src)
            [xi0,xi0t,xi0tt] = Xi0.Derivatives();
            [xi1,xi1t,xi1tt] = Xi1.Derivatives();
             
			 f      = Src.Derivatives();
			 [a,ar] = obj.Coeffs.Derivatives('ar');
             [~,at] = obj.Coeffs.Derivatives('at');
             [b,br] = obj.Coeffs.Derivatives('br');
             [~,bt] = obj.Coeffs.Derivatives('bt');
			 sigma	= obj.Coeffs.Derivatives('sigma');
			 
			% urr2 = (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2))./a - xi1./obj.r0 ;
			 
             sint = sin(obj.th);
             cost = cos(obj.th);
             sincost = sint.*cost;
             sint2 = sint.^2; 
             cost2 = cost.^2;
             
             urr = (f + sigma.*xi0 ... 
             -  (xi0t.*( (obj.r0.*(br-ar) + 2*(a-b) ).*sincost   + bt.*cost2  +at.*sint2)./(obj.r0.^2)  ... 
            + xi1.* ( (bt-at).*sincost  +(obj.r0.*ar+b).*cost2 + (a +obj.r0.* br).*sint2 )./obj.r0 ...
            + 2*sincost./obj.r0 .* (b-a).* xi1t  ...
            + (b.*cost2 + a.*sint2).*xi0tt./(obj.r0.^2)) )./(a.*cost2 + b.*sint2);
	 
		 
             res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*urr/2;% + (obj.dr.^3).*u3r/6;
        end
end
end