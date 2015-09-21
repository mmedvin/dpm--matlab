classdef EBPolarLaplace5OrderExtension < Tools.Extensions.SuperPolarTwoTupleExtension
methods

  function obj = EBPolarLaplace5OrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
        end
        
        function res = Expansion(obj,Xi0,Xi1,Src)
                 [xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
             [xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
             assert('tbd')
             
			 
			 % 			 %this is temporary part
			 % 			 %assuming sigma constant or at least doesn't depends on r
			 % 			 br = LapCoeffs.br;
			 % 			 btr= LapCoeffs.btr;
			 % 			 btrr= LapCoeffs.btrr;
			 % 			 u3r = (fr + sigma.*xi1 - arr.*xi1 - ar.*urr -  (btr.*xi0t + bt.*xi1t + br.*xi0tt + b.*xi1tt)./(obj.r0.^2) + 2*(bt.*xi0t + b.*xi0tt)./(obj.r0.^3) )./a ...
			 % 				 + (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2)).*ar./a./a ...
			 % 				 - urr./obj.r0 + xi1./obj.r0.^2 ;
			 %
			 %
			 % 			 u4r =   (frr + sigma.*urr - a3r.*xi1  - 2*arr.*urr - ar.*u3r ...
			 % 					-(btrr.*xi0t + btr.*xi1t + btr.*xi1t +bt.*urrt
			 % 				+ br.*xi0tt + b.*xi1tt)./(obj.r0.^2) +  2*(btr.*xi0t + bt.*xi1t + br.*xi0tt + b.*xi1tt)./(obj.r0.^3)
			 %
			 % 			 + 2*(bt.*xi0t + b.*xi0tt)./(obj.r0.^3) )./a ...
			 % 				 + (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2)).*ar./a./a ...
			 % 				 - urr./obj.r0 + xi1./obj.r0.^2 ;
			 
             res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 ;
   
        end
end
end