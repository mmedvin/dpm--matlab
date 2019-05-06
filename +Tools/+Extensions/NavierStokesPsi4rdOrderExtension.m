classdef NavierStokesPsi4rdOrderExtension < Tools.Extensions.SuperPolarTwoTupleExtension
    %NAVIERSTOKESPSI3RDORDEREXTENSION 
    
    properties(Access =public)
        XiPsi;
        noXiPsi;
        Wpsi;
    end
    
    properties(Access= protected)
        ArgsXiPsi;
    end
    
    methods
        function obj = NavierStokesPsi4rdOrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
            
            % if isfield(Arguments,'PsiBC')
            
            obj.Wpsi = {spalloc(obj.Grid.Nx*obj.Grid.Ny,1,numel(obj.Scatterer.GridGamma))};

            
            obj.ArgsXiPsi.PsiBC = Arguments.PsiBC;
            obj.ArgsXiPsi.Scatterer = obj.Scatterer.TheScatterer();
            obj.XiPsi = Tools.Misc.XiPsi(obj.ArgsXiPsi);
            obj.noXiPsi = Tools.Misc.XiPsi();
            
            % end

        end
        
        function Update(obj,PsiBC)
            obj.ArgsXiPsi.PsiBC = PsiBC;
            obj.XiPsi = Tools.Misc.XiPsi(obj.ArgsXiPsi);
            obj.ExpandPsi();
        end
        
        function Expand(obj)
            Expand@Tools.Extensions.SuperPolarTwoTupleExtension(obj);
            obj.ExpandPsi();
        end
        
        function xi0j = ExpandedBasis0(obj,n)
            
            Xi   = obj.Basis.Handle(obj.BasisArg, n, obj.Basis.MoreParams);
            
            xi0j = obj.Expansion(Xi      ,    obj.NoXi, obj.noXiPsi,  obj.NoSource);
        end
        
        function xi1j = ExpandedBasis1(obj,n)
            
            Xi   = obj.Basis.Handle(obj.BasisArg, n, obj.Basis.MoreParams);
            
            xi1j = obj.Expansion(obj.NoXi,    Xi      , obj.noXiPsi,   obj.NoSource);
        end
        
        function ExpandPsi(obj)
            obj.Wpsi{1}(obj.Scatterer.GridGamma) = obj.Expansion(obj.NoXi,obj.NoXi,obj.XiPsi,obj.NoSource);
        end
        
        function ExpandSource(obj,SourceHandle,SourceParams)
% do I need it here?
% if isequal(SourceHandle , @Tools.Source.NavierStokesDescreteSrc)%TODO remove it and make a subclass
%     SourceParams.Np = obj.Scatterer.Np;
%     SourceParams.Mp = obj.Scatterer.Mp;
%     SourceParams.GridGamma = obj.Scatterer.GridGamma;
%     SourceParams.Grid = obj.Scatterer.Grid;
% end
            
            Source = SourceHandle(obj.Scatterer.TheScatterer,obj.CoeffsHandle,obj.CoeffsParams,SourceParams);
            
            obj.Wf                             = {spalloc(obj.Grid.Nx*obj.Grid.Ny,1,numel(obj.Scatterer.GridGamma))};
            obj.Wf{1}(obj.Scatterer.GridGamma) = obj.Expansion(obj.NoXi,obj.NoXi,obj.noXiPsi,Source);
        end
        
         function XiP = Expansion(obj,Omega,Omega_r,XiPsi,Source)
             r = obj.r0;
             [xi0P,xi1P,xi0Ptt,xi1Ptt,xi0Ptttt] = XiPsi.Derivatives();
             [O,~,Ott] = Omega.Derivatives();
             Or = Omega_r.Derivatives();
             Src = Source.Derivatives();
             
             dr = obj.Scatterer.dr;
             dr2=dr.^2;
             dr3=dr2.*dr;
             dr4=dr3.*dr;
             
             r2=r.*r;
             r3=r2.*r;
             r4=r3.*r;
             
             %Orr = Src - Or./r - Ott./r2 + obj.CoeffsParams.sigma*O;
             
             Prr = O - xi1P./r - xi0Ptt./r2;
             %Prrtt = Ott - xi1Ptt./r - xi0Ptttt./r2;
             
             P3r = Or + xi1P./r2 - Prr./r + 2*xi0Ptt./r3 - xi1Ptt./r2;
             %P4r = Orr - 2*xi1P./r3 + 2*Prr./r2 - P3r./r - 6*xi0Ptt./r4 + 4*xi1Ptt./r3 - Prrtt./r2;
             
             XiP = xi0P + xi1P.*dr + Prr.*(dr2/2)+ P3r.*(dr3/6);% + P4r.*(dr4/24);
         end

    end
end

