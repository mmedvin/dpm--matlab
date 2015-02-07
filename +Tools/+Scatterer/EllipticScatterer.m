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
        ExpansionType = 25;
        Stencil=5
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
			end
%             obj.Eta0 = AddParams.Eta0;
%             obj.FocalDistance = AddParams.FocalDistance;
            
            obj.MyGrid();
			
            if isfield(AddParams,'Stencil')
                obj.Stencil = AddParams.Stencil;
			end
            
            obj.SplitGrid(obj.Stencil);
            
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
            
            [xi0,xi0f,xi0ff,~,xi0ffff] = Xi0.Derivatives();
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
             [xi1,xi1f,xi1ff] = Xi1.Derivatives();
             
			 [f,fn] = F.Derivatives();
			 [a,an] = LapCoeffs.Derivatives('an');
             [~,af] = LapCoeffs.Derivatives('af');
             [b,bn] = LapCoeffs.Derivatives('bn');
			 [~,bf] = LapCoeffs.Derivatives('bf');
			 sigma	= LapCoeffs.Derivatives('sigma');
             [h,hn] = obj.MetricsAtScatterer.metrics();
             
			 h2 = h.^2;

			 unn1 = ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 )./a ;
             
             fd = obj.FocalDistance;
             fd4 = fd^4;

             x  = fd*cosh(obj.Eta0).*cos(obj.phi);
             y  = fd*sinh(obj.Eta0).*sin(obj.phi);
             xn = fd*sinh(obj.Eta0).*cos(obj.phi);            
             xf =-fd*cosh(obj.Eta0).*sin(obj.phi);
             
             cn = fd4*(b-a).*(sinh(2*obj.Eta0).*(cos(4*obj.phi) - 3) + sinh(4*obj.Eta0).*cos(2.*obj.phi))./(8*h2) + x.*y.*(bf-af) + an.*xn.^2 + bn.*xf.^2; 
            
             cf = fd4*(a-b).*((cosh(4*obj.Eta0) - 3).*sin(2*obj.phi) + cosh(2*obj.Eta0).*sin(4.*obj.phi))./(8*h2) + x.*y.*(bn-an) + bf.*xn.^2 + an.*xf.^2; 
            
             unn = ((f + sigma.*xi0).*(h.^4) - (b.*(xn.^2) + a.*(xf.^2)).*xi0ff - 2*x.*y.*(b-a).*xi1f - cn.*xi1 - cf.*xi0f)./(a.*(xn.^2) + b.*(xf.^2));
            
             %(a.*(xn.^2) + b.*(xf.^2)).*unn + (b.*(xn.^2) + a.*(xf.^2)).*xi0ff + 2*x.*y.*xi1f + c1.*xi1 + c2.*xi0f;
             
             
             
             res = xi0 + obj.deta.*xi1 + (obj.deta.^2).*unn/2 ;%+ (obj.deta.^3).*unnn/6;
             
         end
		function res = Expansion5thOrdrHomoLap(obj,Xi0,Xi1,F,LapCoeffs)
             %[xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
             %[xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
             
             [xi0,xi0f,xi0ff,xi0fff,xi0ffff] = Xi0.Derivatives();
             [xi1,xi1f,xi1ff] = Xi1.Derivatives();
             
			 [f,fn,fnn,f_f,f_ff] = F.Derivatives();
			 [a,an,ann,a3n,af,aff,anf,anff] = LapCoeffs.Derivatives('a');
			 [b,bf,bff,b3f,bn,bnn,bfn,bfnn] = LapCoeffs.Derivatives('b');
			 sigma	= LapCoeffs.Derivatives('sigma');
             [h,hn,hnn,h3n,h4n,hf,hff,h3f,h4f] = obj.MetricsAtScatterer.metrics();
             
			 h2 = h.^2;
             %assume a = b = const, sigma =0;
             %unn   = f.*h2./a - xi0ff;
             %unnn  = (h2.*fn + 2*f.*h.*hn)./a - xi1ff;
             %unnff = (2*f.*(hf.^2) + h2.*f_ff + 2*h.*(2.*f_f.*hf + f.*hff))./a - xi0ffff;  
             %unnnn = (2*f.*(hn.^2) + h2.*fnn  + 2*h.*(2*fn.*hn    + f.*hnn))./a - unnff;
             

             			 unn  = ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 )./a ;
             
                          %atm assumimg sigma is constant
                         
                         unnn =  (fn.*h2 + 2*f.*h.*hn - an.*unn - bfn.*xi0f - bn.*xi0ff - ann.*xi1 - bf.*xi1f - b.*xi1ff)./a ...
                                ... ( (fn + sigma.*xi1).*h2 + (f + sigma.*xi0).*(2*h.*hn) -  bfn.*xi0f  -  bf.*xi1f - bn.*xi0ff - b.*xi1ff - ann.*xi1 - an.*unn )./a ...
                               - ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(an./(a.^2));
                           
                          
                          unnf = ( (f_f + sigma.*xi0f).*h2 + (f + sigma.*xi0).*(2*h.*hf) - bff.*xi0f - 2*bf.*xi0ff - b.*xi0fff - anf.*xi1  - an.*xi1f)./a ...
                               - ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*af./(a.^2);
             
                          unnff= ( (f_ff + sigma.*xi0f).*h2 + (f_f + sigma.*xi0f).*(3*h.*hf) + (f + sigma.*xi0).*(2*hf.^2 + 2*h.*hff) ...
                          - b3f.*xi0f - 3*bff.*xi0ff - 3*bf.*xi0fff - b.*xi0ffff - anff.*xi1 - 2*anf.*xi1f - an.*xi1ff)./a ...
                               - 2*( (f_f + sigma.*xi0f).*h2 + (f + sigma.*xi0).*(2*h.*hf) - bff.*xi0f - bf.*xi0ff - bf.*xi0ff - b.*xi0fff - anf.*xi1  - an.*xi1f).*af./(a.^2) ...
                               - ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*aff./(a.^2) ...
                               + 2*( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(af.^2)./(a.^3);
             
                          unnnn= ( (fnn + sigma.*unn).*h2 + (fn + sigma.*xi1).*(3*h.*hn) + (f + sigma.*xi0).*(2*hn.^2 + 2*h.*hnn) ...
                                    -  bfnn.*xi0f - 2*bfn.*xi1f - bf.*unnf - bnn.*xi0ff - 2*bn.*xi1ff - b.*unnff - a3n.*xi1 - 2*ann.*unn - an.*unnn )./a ...
                               - 2*( (fn + sigma.*xi1).*h2 + (f + sigma.*xi0).*(2*h.*hn) -  bfn.*xi0f  -  bf.*xi1f - bn.*xi0ff - b.*xi1ff - ann.*xi1 - an.*unn ).*an./(a.^2)...
                               - ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*ann./(a.^2) ...
                               + 2*( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(an.^2)./(a.^3);
              
              
             res = xi0 + obj.deta.*xi1 + (obj.deta.^2).*unn/2 + (obj.deta.^3).*unnn/6 + (obj.deta.^4).*unnnn/24;

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
