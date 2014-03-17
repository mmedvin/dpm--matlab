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
        ExpansionType;
        AddParams;
    end
    
    methods
        function obj = EllipticScatterer(Grid,AddParams)
            
            if ~exist('Grid','var') || ~exist('AddParams','var')
                error('Costructor called without args');
            end
            
            obj = obj@Tools.Scatterer.SingleScatterer(Grid);
            
            obj.AddParams = AddParams;
			if isfield(AddParams,'ExpansionType')
				obj.ExpansionType = AddParams.ExpansionType;
			else
				obj.ExpansionType = 25;
			end
%             obj.Eta0 = AddParams.Eta0;
%             obj.FocalDistance = AddParams.FocalDistance;
            
            obj.MyGrid();
			
			if obj.ExpansionType == 33
				obj.SecondOrderSplitGrid();
			else
				obj.SplitGrid();
			end
            
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
        
		function res = Expansion5thOrdrHelm(obj,Xi0,Xi1,Source,WaveNumber)
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
        function res = Expansion3thOrdrLap(obj,Xi0,Xi1,F,LapCoeffs)
             [xi0,xi0f,xi0ff] = Xi0.Derivatives();
             [xi1] = Xi1.Derivatives();
             
			 f = F.Derivatives();
			 [a,an] = LapCoeffs.Derivatives('a');
			 [b,bf] = LapCoeffs.Derivatives('b');
			 sigma	= LapCoeffs.Derivatives('sigma');
			 h2 = obj.MetricsAtScatterer.metrics().^2;

			 unn = ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 )./a ;
             
             res = xi0 + obj.deta.*xi1 + (obj.deta.^2).*unn/2 ;
             
         end
		function res = Expansion5thOrdrHomoLap(obj,Xi0,Xi1,F,LapCoeffs)
             [xi0,xi0t,xi0tt,xi0tttt,xi0tttttt] = Xi0.Derivatives();
             [xi1,~,xi1tt,xi1tttt,~] = Xi1.Derivatives();
             assert('tbd')
             res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 ;
         end
		
        function res = Expansion(obj,Xi0,Xi1,Source,Coeffs)
			
			switch obj.ExpansionType
				case 25
					res = Expansion5thOrdrHelm(obj,Xi0,Xi1,Source,Coeffs);
				case 33
					res = Expansion3thOrdrLap(obj,Xi0,Xi1,Source,Coeffs);
				case 35
					res = Expansion5thOrdrHomoLap(obj,Xi0,Xi1,Source,Coeffs);
			end
			
			
            
                      
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
