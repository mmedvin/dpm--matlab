classdef HelmholtzSourceStarShaped < Tools.Source.HelmholtzSourceR
    properties (Dependent = true)
        %Source;%Fn;Ff;Fnn;Fff;    % WN;
    end
    
    properties(Access = protected)       
        %         Scatterer;
        %         Exact;
        %         IsExactReady=false;
        %
        %         WaveNumberHndl;
        %         WaveNumberAddParams;
    end
    
    methods
        function [F,Fn,Fnn,Fs,Fss] = Derivatives(obj,Scatterer)
            [F,Fr,Frr,Ft,Ftt,Frt] = Derivatives@Tools.Source.HelmholtzSourceR(obj);
            if nargout>1 && exist('Scatterer','var')
                
                [x,dx] = Scatterer.XHandle.Derivatives(Scatterer.BasisArg);
                [y,dy] = Scatterer.YHandle.Derivatives(Scatterer.BasisArg);
                r_ = sqrt(x.^2 + y.^2);                                
                
                Fx = Fr.*x./r_ - Ft.*y./(r_.^2);
                Fy = Fr.*y./r_ + Ft.*x./(r_.^2);
                
                Fxx = Ftt.*y.^2./(r_.^4) + Frr.*(x.^2)./(r_.^2) - 2*Frt.*x.*y./(r_.^3) ...
                    + 2*Ft.*x.*y./(r_.^4) + Fr.*(y.^2)./(r_.^3) ;
                
                Fyy = Ftt.*x.^2./(r_.^4) + Frr.*(y.^2)./(r_.^2) + 2*Frt.*x.*y./(r_.^3) ...
                    - 2*Ft.*x.*y./(r_.^4) + Fr.*(x.^2)./(r_.^3) ;
                
                Fxy =-Ftt.*x.*y./(r_.^4) + Frr.*x.*y./(r_.^2) + Frt.*(x.^2 - y.^2)./(r_.^3) ...
                    + Ft.*(y.^2-x.^2)./(r_.^4) - (Fr.*x.*y)./(r_.^3) ;
                
                
                curv = Scatterer.Curvature(Scatterer.BasisArg);
                h = Scatterer.Metrics(Scatterer.BasisArg);
                
                Fn  = (Fx.*dy  - Fy.*dx)./h;
                Fnn = (Fxx.*dy.^2 - 2*Fxy.*dx.*dy + Fyy.*dx.^2 )./(h.^2);
                
                Fs  = (Fx.*dx + Fy.*dy)./h;
                Fss = (Fxx.*dx.^2 + 2*Fxy.*dx.*dy + Fyy.*dy.^2)./(h.^2) + curv.*Fn;
            end
        end
            
        function obj = HelmholtzSourceStarShaped(Scatterer, WaveNumberHndl,WaveNumberAddParams,~)
            obj = obj@Tools.Source.HelmholtzSourceR(Scatterer, WaveNumberHndl,WaveNumberAddParams);
        end        
        
        function S = getSource(obj)
            S=obj.Derivatives();
        end    
    end
    
end

