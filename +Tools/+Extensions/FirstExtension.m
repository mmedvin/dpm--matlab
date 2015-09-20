classdef FirstExtension < Tools.Extensions.SuperExtension
    properties
        Basis;
        Scatterer;
        Coeffs;
        Grid;
    end

    properties(Access=protected)
        NoXi;
        
        GridGamma;
    end
    methods
        function obj = FirstExtension(Arguments)
            obj.Basis     = Arguments.Basis;
            obj.Coeffs    = Arguments.Coeffs;
            obj.Scatterer = Arguments.Scatterer;
            obj.Grid      = Arguments.Grid;
            
            obj.W         = cell(1,2);
            obj.W{1}      = spalloc( Arguments.Grid.Nx*Arguments.Grid.Ny, obj.Basis.NBss, numel(obj.Scatterer.GridGamma)*obj.Basis.NBss);
            obj.W{2}      = spalloc( Arguments.Grid.Nx*Arguments.Grid.Ny, obj.Basis.NBss, numel(obj.Scatterer.GridGamma)*obj.Basis.NBss);
            
            obj.NoXi = obj.Basis.Handle();
        end
        function Expand(obj)
            [tmp1,tmp2]                         = arrayfun(@(n) obj.ExpandedBasis(n), obj.Basis.Indices, 'UniformOutput', false);
            obj.W{1}(obj.Scatterer.GridGamma,:) = cell2mat(tmp1);
            obj.W{2}(obj.Scatterer.GridGamma,:) = cell2mat(tmp2);
        end
        
        function [xi0j,xi1j] = ExpandedBasis(obj,n)
            
            Xi   = obj.Basis.Handle(obj.Scatterer.BasisArg, n, obj.Basis.MoreParams);            
            
            xi0j = obj.Scatterer.Expansion(Xi      ,    obj.NoXi,   obj.NoSource,   obj.Coeffs);
            xi1j = obj.Scatterer.Expansion(obj.NoXi,    Xi      ,   obj.NoSource,   obj.Coeffs);
            
        end
        
        function ExpandSource(obj,Source)
            obj.Wf                          = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.GridGamma));
            obj.Wf(obj.Scatterer.GridGamma) = obj.Scatterer.Expansion(obj.NoXi,obj.NoXi,Source,obj.Coeffs);
        end

        function GG = get.GridGamma(obj)
            GG = obj.Scatterer.GridGamma;
        end
    end
    
end