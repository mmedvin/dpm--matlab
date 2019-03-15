classdef EBPolarLaplace5OrderExtension < Tools.Extensions.SuperPolarTwoTupleExtension
methods

  function obj = EBPolarLaplace5OrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
        end
        
        function res = Expansion(obj,Xi0,Xi1,Src)
                 [xi0,xi0t,xi0tt,xi0ttt,xi0tttt,xi0tttttt] = Xi0.Derivatives();
             [xi1,xi1t,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
            % assert('tbd')
             
			 
			 % 			 %this is temporary part
			 % 			 %assuming sigma constant or at least doesn't depends on r
             
             [f,fr,frr,ft,ftt]      = Src.Derivatives();
             
             [a,ar,arr,a3r] = obj.Coeffs.Derivatives('ar');
             [~,at,att, art, artt] = obj.Coeffs.Derivatives('at');
             [b,br,brr] = obj.Coeffs.Derivatives('br');
             [~,bt,btt,b3t,brt,brrt] = obj.Coeffs.Derivatives('bt');
             sigma	= obj.Coeffs.Derivatives('sigma');
             
             r=obj.r0;
             r2 = r^2;
             r3=r^3;
             r4=r^4;
             
             urr = (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2))./a - xi1./obj.r0 ;

             %              u3r = (fr + sigma.*xi1 - arr.*xi1 - ar.*urr ...
             %                  -  (btr.*xi0t + bt.*xi1t + br.*xi0tt + b.*xi1tt)./(obj.r0.^2)...
             %                  + 2*(bt.*xi0t + b.*xi0tt)./(obj.r0.^3) )./a ...
             %                  + (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2)).*ar./a./a ...
             %                  - urr./obj.r0 + xi1./obj.r0.^2 ;
			 
			u3r = -(urr./r) + xi1./r2 - (ar.*(f + sigma.*xi0 - (bt.*xi0t + b.*xi0tt)/r2 - ar.*xi1))./(a^.2) ...
                + (fr-(ar.*urr) + (2*(bt.*xi0t + b.*xi0tt))./r3 - arr.*xi1 + sigma*xi1 ...
                - (brt.*xi0t + br.*xi0tt + bt.*xi1t + b.*xi1tt)/r2)./a;
             
            urrt = -((at.*(f + sigma.*xi0 - (bt.*xi0t + b.*xi0tt)/r2 - ar.*xi1))./(a^.2)) ...
                 - xi1t/r + (ft+sigma*xi0t - (b.*xi0ttt + btt.*xi0t + 2*bt.*xi0tt)/r2 - art.*xi1 - ar.*xi1t)./a;
      
            urrtt = ((2*(at.^2))./(a.^3) - att./(a.^2)).*(f + sigma*xi0 - (bt.*xi0t + b.*xi0tt)/r2 - ar.*xi1) ...
                  - (2*at.*(ft+sigma*xi0t - (b.*xi0ttt + btt.*xi0t + 2*bt.*xi0tt)/r2 - art.*xi1 ...
                  - ar.*xi1t))./(a.^2) - xi1tt/r + (ftt+sigma*xi0tt ...
                  - (3*bt.*xi0ttt + b3t.*xi0t + 3*btt.*xi0tt + b.*xi0tttt)/r2 - artt.*xi1 - 2*art.*xi1t - ar.*xi1tt)./a;
             
			u4r = -(u3r/r) + (2*urr)/r2 - (2*xi1)/r3 ...
                +((2*(ar.^2))./(a.^3) - arr./(a.^2)).*(f + sigma*xi0 - (bt.*xi0t + b.*xi0tt)/r2 - ar.*xi1) ...
                - (2*ar.*(fr - (ar.*urr) + (2*(bt.*xi0t + b.*xi0tt))/r3 - arr.*xi1 + sigma*xi1 ...
                -(brt.*xi0t + br.*xi0tt + bt.*xi1t + b.*xi1tt)/r2))./(a.^2) ...
                + (frr-(ar.*u3r) - 2*arr.*urr + sigma*urr - (6*(bt.*xi0t + b.*xi0tt))/r4...
                - a3r.*xi1 + (4*(brt.*xi0t + br.*xi0tt + bt.*xi1t + b.*xi1tt))/r3 ...
                - (bt.*urrt + b.*urrtt + brrt.*xi0t + brr.*xi0tt + 2*brt.*xi1t + 2*br.*xi1tt)/r2)./a; 
            
            
            
             res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*urr/2 + (obj.dr.^3).*u3r/6 + (obj.dr.^4).*u4r/24 ;
   
        end
end
end