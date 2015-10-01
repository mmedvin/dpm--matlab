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
    end
    methods
        function obj = TwoTupleExtension(Arguments)
            obj.Basis     = Arguments.Basis;
            obj.Coeffs    = Arguments.Coeffs;
            obj.Scatterer = Arguments.Scatterer;
            obj.Grid      = Arguments.Grid;
            
            obj.W         = cell(1,2);
            obj.W{1}      = struct('GridGamma',obj.GridGamma, 'msk',obj.Scatterer.Mp,'W',...
                                    spalloc( Arguments.Grid.Nx*Arguments.Grid.Ny, obj.Basis.NBss, numel(obj.GridGamma)*obj.Basis.NBss));
            obj.W{2}      = struct('GridGamma',obj.GridGamma, 'msk',obj.Scatterer.Mp, 'W',...
                                    spalloc( Arguments.Grid.Nx*Arguments.Grid.Ny, obj.Basis.NBss, numel(obj.GridGamma)*obj.Basis.NBss));
            
            obj.NoXi = obj.Basis.Handle();
        end
        function Expand(obj)
            [tmp1,tmp2]                         = arrayfun(@(n) obj.ExpandedBasis(n), obj.Basis.Indices, 'UniformOutput', false);
            obj.W{1}.W(obj.GridGamma,:) = cell2mat(tmp1);
            obj.W{2}.W(obj.GridGamma,:) = cell2mat(tmp2);
        end
        
        function [xi0j,xi1j] = ExpandedBasis(obj,n)
            
            Xi   = obj.Basis.Handle(obj.BasisArg, n, obj.Basis.MoreParams);            
            
            xi0j = obj.Expansion(Xi      ,    obj.NoXi,   obj.NoSource);
            xi1j = obj.Expansion(obj.NoXi,    Xi      ,   obj.NoSource);
            
        end
        
        function ExpandSource(obj,Source)
            obj.Wf                          = struct('GridGamma',obj.GridGamma, 'msk',obj.Scatterer.Mp,'W',...
                                                        spalloc(obj.Grid.Nx*obj.Grid.Ny,1,numel(obj.Scatterer.GridGamma)));
            obj.Wf.W(obj.Scatterer.GridGamma) = obj.Expansion(obj.NoXi,obj.NoXi,Source);
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