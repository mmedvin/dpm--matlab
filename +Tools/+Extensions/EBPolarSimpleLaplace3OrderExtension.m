classdef EBPolarSimpleLaplace3OrderExtension < Tools.Extensions.SuperPolarTwoTupleExtension
methods

  function obj = EBPolarSimpleLaplace3OrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
  end
        
  function res = Expansion(obj,Xi0,Xi1,Src)
      [xi0,~,xi0tt] = Xi0.Derivatives();
      xi1 = Xi1.Derivatives();
      
      f      = Src.Derivatives();

      sigma	= obj.Coeffs.Derivatives('sigma');
      
      urr = f + sigma.*xi0  - xi1./obj.r0 - xi0tt./(obj.r0.^2);
      
      res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*urr/2;% + (obj.dr.^3).*u3r/6;
  end
end
end