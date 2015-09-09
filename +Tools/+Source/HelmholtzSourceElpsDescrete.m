classdef HelmholtzSourceElpsDescrete < Tools.Source.DescreteSource
    
    
    
    
    
    methods
        function obj = HelmholtzSourceElpsDescrete(Scatterer, WaveNumberHndl,WaveNumberParams,Params)
            
            [X,Y] = Params.Grid.mesh();
                
            Scatt.FocalDistance = Scatterer.FocalDistance;
            [Scatt.Eta,Scatt.Phi] = Params.Grid.ToElliptical(Scatt.FocalDistance);
            
            
            Scatt.r=abs(X+1i.*Y);
            WNPlr=WaveNumberHndl(Scatt,WaveNumberParams);
            
            Exact  = Tools.Exact.ExactElpsVarKr(Scatt, WNPlr);
            
            %TmpSrc =  Tools.Source.HelmholtzSource(Scatterer, WaveNumberHndl,WaveNumberParams);
            u = Exact.Derivatives();
            s = Exact.calc_s();
            
            Params.Src = s.*u;
            obj =  obj@Tools.Source.DescreteSource(Scatterer, [],[],Params);
            
        end
        
    end
    
    
end
