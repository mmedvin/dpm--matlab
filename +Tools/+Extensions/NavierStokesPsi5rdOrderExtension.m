classdef NavierStokesPsi5rdOrderExtension < Tools.Extensions.NavierStokesPsi4rdOrderExtension

    methods
        function obj = NavierStokesPsi5rdOrderExtension(Arguments)
                    obj = obj@Tools.Extensions.NavierStokesPsi4rdOrderExtension(Arguments);
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
            obj.Wf{1}(obj.Scatterer.GridGamma) = obj.Expansion(obj.NoXi,obj.NoXi,obj.noXiPsi,Source);
        end

    end


end