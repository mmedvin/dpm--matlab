classdef HelmholtzSource < Tools.Source.SuperHelmholtzSource
    properties (Dependent = true)
        Source;%Fn;Ff;Fnn;Fff;    % WN;
    end
    
    properties(Access = protected)       
        Scatterer; 
        Exact;
        IsExactReady=false;
        
        WaveNumberHndl;
        WaveNumberAddParams;        
    end
    
    methods
        function [F,Fn,Ff,Fnn,Fff] = Derivatives(obj,Ex)
            
            if obj.IsDummy
                F=0;
                Fn=0;
                Ff=0;
                Fnn=0;
                Fff=0;
            else
                if ~exist('Ex','var') 
                   if ~obj.IsExactReady
                       
                       WN = obj.WaveNumberHndl(obj.Scatterer,obj.WaveNumberAddParams);
                       
                       obj.Exact = Tools.Exact.ExactElpsVarKr(obj.Scatterer, WN);
                       obj.IsExactReady = true;
                   end
                    Ex = obj.Exact;
                end
                               
                [u,un,unn,uf,uff] = Ex.Derivatives();
                
                if nargout==1,  s                = Ex.calc_s();  end
                if nargout==2, [s,sn]            = Ex.calc_s();  end
                if nargout==3, [s,sn,snn]        = Ex.calc_s();  end
                if nargout==4, [s,sn,snn,sf]     = Ex.calc_s();  end
                if nargout==5, [s,sn,snn,sf,sff] = Ex.calc_s();  end
                
                F   = s.*u;%urr + ur./r + utt./(r.^2) + k.^2.*u;%
                
                if nargout>1, Fn  = un.*s + u.*sn;               end
                if nargout>2, Ff  = uf.*s + u.*sf;               end
                if nargout>3, Fnn = unn.*s + 2*un.*sn + u.*snn;  end
                if nargout>4, Fff = uff.*s + 2*uf.*sf + u.*sff;  end
            end
            
        end
            
        function obj = HelmholtzSource(Scatterer, WaveNumberHndl,WaveNumberAddParams)%(Exact)%(FocalDist,eta,phi,k0,r0)
                obj.Scatterer = Scatterer;
                obj.WaveNumberHndl = WaveNumberHndl;
                obj.WaveNumberAddParams = WaveNumberAddParams;
                obj.IsDummy = false;                
        end
        
        
        function S = get.Source(obj)
            S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
            
%             ScatGamma = struct('FocalDistance',obj.Scatterer.FocalDistance,'Eta', obj.Scatterer.Eta0,'Phi', obj.Scatterer.Phi(obj.Scatterer.GridGamma));
            
            WN            = obj.WaveNumberHndl(obj.Scatterer.TheScatterer,obj.WaveNumberAddParams);
            ExactOnGamma  = Tools.Exact.ExactElpsVarKr(obj.Scatterer.TheScatterer, WN);
            
            [F,Fn,~,Fnn] = obj.Derivatives(ExactOnGamma);
            
            S(obj.Scatterer.GridGamma) = F + obj.Scatterer.deta.*Fn  + (obj.Scatterer.deta.^2).*Fnn/2;%taylor
            
            
           
            tmp = obj.Derivatives();
            S(obj.Scatterer.Inside) = tmp(obj.Scatterer.Inside);   %was obj.Source(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
        end
        
        
%         function [UV,UVx,UVy,UVxx,UVyy] = ChainRule(U,V)
%             [u,ux,uy,uxx,uyy] = U.Derivatives();
%             [v,vx,vy,vxx,vyy] = V.Derivatives();
%             
%             UV   = u.*v;
%             UVx  = ux.*v + u.*vx;
%             UVy  = uy.*v + u.*vy;
%             UVxx = uxx.*v + 2*ux.*vx + u.*vxx;
%             UVyy = uyy.*v + 2*uy.*vy + u.*vyy;
%         end
%         
        
        

    
    end
    
end

