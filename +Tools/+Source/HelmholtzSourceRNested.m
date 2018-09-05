classdef HelmholtzSourceRNested < Tools.Source.HelmholtzSourceR
    
    methods
        function obj = HelmholtzSourceRNested(Scatterer, WaveNumberHndl,WaveNumberAddParams,~)
            obj.Scatterer = Scatterer;
            obj.WaveNumberHndl = WaveNumberHndl;
            obj.WaveNumberAddParams = WaveNumberAddParams;
            obj.IsDummy = false;
        end
        
        function S = getSource(obj)
            S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
            
            %             ScatGamma = struct('FocalDistance',obj.Scatterer.FocalDistance,'Eta', obj.Scatterer.Eta0,'Phi', obj.Scatterer.Phi(obj.Scatterer.GridGamma));
            
            WN            = obj.WaveNumberHndl(obj.Scatterer.TheScatterer,obj.WaveNumberAddParams);
            ExactOnGamma  = Tools.Exact.ExactVarKr(obj.Scatterer.TheScatterer, WN);
            
            [F,Fr,Frr] = obj.Derivatives(ExactOnGamma);
            
            %S(obj.Scatterer.GridGamma) = F + obj.Scatterer.dn.*Fr  + (obj.Scatterer.dn.^2).*Frr/2;%taylor
            S(obj.Scatterer.GridGamma) = F + obj.Scatterer.dr.*Fr  + (obj.Scatterer.dr.^2).*Frr/2;%taylor
            
            tmp = obj.Derivatives();
            S(obj.Scatterer.Inside) = tmp(obj.Scatterer.Inside);   %was obj.Source(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
            S=tmp;
        end
    end
    
end
