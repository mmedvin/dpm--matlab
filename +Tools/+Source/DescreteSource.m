classdef DescreteSource < Tools.Source.SuperSource
    properties (Dependent = true)
        Source;%Fn;Ff;Fnn;Fff;    % WN;
    end
    
    properties(Access = protected)       
        Scatterer; 
        %Exact;
        %IsExactReady=false;
        
        %WaveNumberHndl;
        %WaveNumberAddParams;
        Src;
        Der;
        
        mFx;mFy;mFxx;mFyy;mFxy;
    end
    
    methods
        function [F,Fx,Fy,Fxx,Fyy,Fxy] = Derivatives(obj)
            
            if obj.IsDummy
                F=0;
                Fx=0;
                Fy=0;
                Fxx=0;
                Fyy=0;
                Fxy=0;
            else
                if isfield(obj.Scatterer,'GridGamma')
                    Fx = obj.mFx;
                    Fy = obj.mFy;
                    Fxx= obj.mFxx;
                    Fyy= obj.mFyy;
                    Fxy= obj.mFxy;
                    F=obj.Src(obj.Scatterer.GridGamma);
                else
                    [Fx,Fy,Fxx,Fyy,Fxy]       = obj.Der.CartesianDerivatives(obj.Src);
                    F=obj.Src;
                end
                
            end
            
        end
            
        function obj = DescreteSource(Scatterer, ~,~,Params)
                obj.Scatterer = Scatterer;
                %obj.WaveNumberHndl = WaveNumberHndl;
                %obj.WaveNumberAddParams = WaveNumberAddParams;
                obj.Der = Tools.Common.SecondDerivative(Params.Grid.Nx,Params.Grid.Ny,Params.Grid.dx,Params.Grid.dy);
                obj.Src = Params.Src;
                
            if isfield(obj.Scatterer,'GridGamma')
                [obj.mFx,obj.mFy,obj.mFxx,obj.mFyy,obj.mFxy]       = obj.Der.CartesianDerivatives(obj.Src);
                obj.mFx = obj.mFx(obj.Scatterer.GridGamma); 
                obj.mFy = obj.mFy(obj.Scatterer.GridGamma);
                obj.mFxx= obj.mFxx(obj.Scatterer.GridGamma);
                obj.mFyy= obj.mFyy(obj.Scatterer.GridGamma);
                obj.mFxy= obj.mFxy(obj.Scatterer.GridGamma);
            end
                
                
                obj.IsDummy = false;                
        end
        
        
        function S = get.Source(obj)
            S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
                        
            S(obj.Scatterer.Inside) = obj.Src(obj.Scatterer.Inside);
        end    
    end
    
end

