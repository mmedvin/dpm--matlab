classdef EllipticScatterer < Tools.Scatterer.SingleScatterer
    properties (Access = public)
        
        Eta;
        Phi;
        
        FocalDistance;
        Eta0;
        
        eta;
        phi;
        deta;
        
        MetricsAtScatterer;
        
        BasisArg;
        TheScatterer;
%         ScattererForSource;
                
        %consider to remove next 2
        R;
        TH;
        
        Inside;
        Outside;
    end
    
    properties (Access = private)
        
        AddParams;
    end
    
    methods
        function obj = EllipticScatterer(Grid,AddParams)
            
            if ~exist('Grid','var') || ~exist('AddParams','var')
                error('Costructor called without args');
            end
            
            obj = obj@Tools.Scatterer.SingleScatterer(Grid);
            
            obj.AddParams = AddParams;
%             obj.Eta0 = AddParams.Eta0;
%             obj.FocalDistance = AddParams.FocalDistance;
            
            obj.MyGrid();
            obj.SplitGrid();%(obj.Eta,obj.Eta0);
            
            obj.GridGamma = intersect(obj.Np,obj.Nm)';
            
            obj.GridOnScatterer();
            
            obj.MetricsAtScatterer = Tools.Metrics.EllipticalMetrics(obj.FocalDistance,obj.Eta0,obj.phi);
        end
  
        function eta0 = get.Eta0(obj)
            eta0 = obj.AddParams.Eta0;
        end
        
        function FD = get.FocalDistance(obj)
            FD = obj.AddParams.FocalDistance;
        end

        
        function BA = get.BasisArg(obj)
            BA = obj.phi;
        end
        
        function TS = get.TheScatterer(obj)
            TS = struct('FocalDistance',obj.FocalDistance,'Eta',obj.Eta0,'Phi',obj.phi);
        end
        
        function I = get.Inside(obj)            
            I = obj.Eta <= obj.Eta0;            
        end
        
        function O = get.Outside(obj)            
            O= obj.Eta > obj.Eta0;            
        end
        
%         function S4S = get.ScattererForSource(obj)
%             GridF = Grids( ...
%                         obj.Grid.x1 - obj.Grid.dx/2 , ...
%                         obj.Grid.xn + obj.Grid.dx/2 , ...
%                         obj.Grid.Nx * 2 + 1         , ...
%                         obj.Grid.y1 - obj.Grid.dy/2 , ...
%                         obj.Grid.yn + obj.Grid.dy/2 , ...
%                         obj.Grid.Ny * 2 + 1         ) ;
%             
%            % AddP = struct('Eta0',obj.Eta0,'FocalDistance',obj.FocalDistance);            
%            % S4S = EllipticScatterer(GridF,AddP);
%            S4S = EllipticScatterer(GridF,obj.AddParams);
%         end
        
        function res = Expansion(obj,Xi0,Xi1,Source,WaveNumber)
            
           % [F,Fn,Ff,Fnn,Fff] = Source.Derivatives();
            [F,Fn,Ff,Fnn,Fff] = Source.Derivatives();
            
            [xi0,xi0f,xi0ff,xi0ffff] = Xi0.Derivatives();
            [xi1,~,xi1ff] = Xi1.Derivatives();
            
            [k,kn,kf,knn,kff] = WaveNumber.Derivatives();
            
            [h,hn,hnn,~,~,hf,hff] = obj.MetricsAtScatterer.metrics();
            
            
            xinn = (h.^2).*F - xi0ff - xi0.*(h.*k).^2;
            xi3n = 2*h.*hn.*F + Fn.*h.^2 - 2*(h.*hn.*k.^2 + k.*kn.*h.^2).*xi0 - xi1.*(h.*k).^2 - xi1ff;
            xinnff = 2*(hf.^2 + h.*hff).*F + 4*h.*hf.*Ff +  (h.^2).*Fff - xi0ffff ...
                - 2*((hf.*k).^2 + h.*hff.*k.^2 + 4*h.*hf.*k.*kf + (h.*kf).^2 + k.*kff.*h.^2).*xi0 ...
                - 4*(h.*hf.*k.^2 + k.*kf.*h.^2).*xi0f - xi0ff.*(h.*k).^2;
            xi4n = 2*(hn.^2 + h.*hnn).*F + 4*h.*hn.*Fn +  (h.^2).*Fnn - xinnff ...
                - 2*((hn.*k).^2 + h.*hnn.*k.^2 + 4*h.*hn.*k.*kn + (h.*kn).^2 + k.*knn.*h.^2).*xi0 ...
                - 4*(h.*hn.*k.^2 + k.*kn.*h.^2).*xi1 - xinn.*(h.*k).^2;
            
            res = xi0 + obj.deta.*xi1 + (obj.deta.^2).*xinn/2 + (obj.deta.^3).*xi3n/6 + (obj.deta.^4).*xi4n/24 ;
            
        end
    end
    
    methods(Access = protected)        
        function MyGrid(obj)
              Z=obj.Grid.Z();
            
            obj.Eta = real(acosh(Z/obj.FocalDistance));
            obj.Phi = imag(acosh(Z/obj.FocalDistance));                       
        end 
        
        function GridOnScatterer(obj)
            obj.eta = obj.Eta(obj.GridGamma);
            obj.deta = obj.eta - obj.Eta0;
            obj.phi = obj.Phi(obj.GridGamma);
        end
    end    
end
