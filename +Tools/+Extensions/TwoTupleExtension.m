classdef TwoTupleExtension < Tools.Extensions.SuperExtension
    properties
        Basis;
        Scatterer;
        Coeffs;
        Grid;
    end

    properties(Access=protected)
        NoXi;
        
        GridGamma;
        BasisArg;
        
        CoeffsHandle;
        CoeffsParams
    end
    methods
        function obj = TwoTupleExtension(Arguments)
            obj.Basis           = Arguments.Basis;
            %obj.Coeffs         = Arguments.Coeffs;
            obj.Scatterer       = Arguments.Scatterer;
            obj.CoeffsHandle 	= Arguments.CoeffsHandle;
            obj.CoeffsParams    = Arguments.CoeffsParams;
            obj.Coeffs          = Arguments.CoeffsHandle(obj.Scatterer.TheScatterer,Arguments.CoeffsParams);%TODO: Consider to remove TheScatterer
            obj.Grid            = Arguments.Grid;
            
            obj.W         = cell(1,2);
            obj.W{1}      = spalloc( Arguments.Grid.Nx*Arguments.Grid.Ny, obj.Basis.NBss0, numel(obj.GridGamma)*obj.Basis.NBss0);
            obj.W{2}      = spalloc( Arguments.Grid.Nx*Arguments.Grid.Ny, obj.Basis.NBss1, numel(obj.GridGamma)*obj.Basis.NBss1);
            
            obj.NoXi = obj.Basis.Handle();
        end
        function Expand(obj)
            if 0
                tmp1 = zeros(numel(obj.GridGamma),numel(obj.Basis.Indices0));
                parfor n=1:numel(obj.Basis.Indices0)
                    tmp1(:,n) = ExpandedBasis0(obj,obj.Basis.Indices0(n));
                end
                obj.W{1}(obj.GridGamma,:) = tmp1;
                
                tmp2 = zeros(numel(obj.GridGamma),numel(obj.Basis.Indices1));
                parfor n=1:numel(obj.Basis.Indices1)
                    tmp2(:,n) = ExpandedBasis1(obj,obj.Basis.Indices1(n));
                end
                obj.W{2}(obj.GridGamma,:) = tmp2;
            else
                tmp1                      = arrayfun(@(n) obj.ExpandedBasis0(n), obj.Basis.Indices0, 'UniformOutput', false);
                tmp2                      = arrayfun(@(n) obj.ExpandedBasis1(n), obj.Basis.Indices1, 'UniformOutput', false);
                obj.W{1}(obj.GridGamma,:) = cell2mat(tmp1);
                obj.W{2}(obj.GridGamma,:) = cell2mat(tmp2);
            end
        end
        
        function xi0j = ExpandedBasis0(obj,n)
            
            Xi   = obj.Basis.Handle(obj.BasisArg, n, obj.Basis.MoreParams);            
            
            xi0j = obj.Expansion(Xi      ,    obj.NoXi,   obj.NoSource);            
        end
        
        function xi1j = ExpandedBasis1(obj,n)
            
            Xi   = obj.Basis.Handle(obj.BasisArg, n, obj.Basis.MoreParams);
            
            xi1j = obj.Expansion(obj.NoXi,    Xi      ,   obj.NoSource);       
        end
        
        function ExpandSource(obj,SourceHandle,SourceParams)
            if isequal(SourceHandle , @Tools.Source.NavierStokesDescreteSrc)%TODO remove it and make a subclass
                SourceParams.Np = obj.Scatterer.Np;
                SourceParams.Mp = obj.Scatterer.Mp;
                SourceParams.GridGamma = obj.Scatterer.GridGamma;
                SourceParams.Grid = obj.Scatterer.Grid;
                SourceParams.GridLinesIntersection = obj.Scatterer.GridLinesIntersection;
            end
            
            Source = SourceHandle(obj.Scatterer.TheScatterer,obj.CoeffsHandle,obj.CoeffsParams,SourceParams);

            obj.Wf                             = {spalloc(obj.Grid.Nx*obj.Grid.Ny,1,numel(obj.Scatterer.GridGamma))};
            obj.Wf{1}(obj.Scatterer.GridGamma) = obj.Expansion(obj.NoXi,obj.NoXi,Source);
        end

        function val = Expansion(obj,xi0,xi1,Src)
            val = obj.Scatterer.Expansion(xi0, xi1, Src, obj.Coeffs);
        end
        function val = get.BasisArg(obj)
            val = obj.Scatterer.BasisArg; 
        end
        function GG = get.GridGamma(obj)
            GG = obj.Scatterer.GridGamma;
        end
    end
    
end