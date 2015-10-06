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
			
                obj.Stencil = AddParams.Stencil;
            
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
            TS = struct('FocalDistance',obj.FocalDistance,'Eta',obj.Eta0,'Phi',obj.phi,'GridGamma',obj.GridGamma);
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
        function res = Expansion3thOrdrLap(obj,Xi0,Xi1,Src,LapCoeffs)
%              [xi0,xi0f,xi0ff] = Xi0.Derivatives();
%              [xi1,xi1f,xi1ff] = Xi1.Derivatives();
%              
% 			 [f,fn] = Src.Derivatives();
% 			 [a,an] = LapCoeffs.Derivatives('an');
%              [~,af] = LapCoeffs.Derivatives('af');
%              [b,bn] = LapCoeffs.Derivatives('bn');
% 			 [~,bf] = LapCoeffs.Derivatives('bf');
% 			 sigma	= LapCoeffs.Derivatives('sigma');
%              [h,hn] = obj.MetricsAtScatterer.metrics();
%              
% 			 h2 = h.^2;
% 
% 			 %unn1 = ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 )./a ;
%              
%              fd = obj.FocalDistance;
%              fd4 = fd^4;
% 
%              x  = fd*cosh(obj.Eta0).*cos(obj.phi);
%              y  = fd*sinh(obj.Eta0).*sin(obj.phi);
%              xn = fd*sinh(obj.Eta0).*cos(obj.phi);            
%              xf =-fd*cosh(obj.Eta0).*sin(obj.phi);
%              
%              cn = fd4*(b-a).*(sinh(2*obj.Eta0).*(cos(4*obj.phi) - 3) + sinh(4*obj.Eta0).*cos(2.*obj.phi))./(8*h2) + x.*y.*(bf-af) + an.*xn.^2 + bn.*xf.^2; 
%             
%              cf = fd4*(a-b).*((cosh(4*obj.Eta0) - 3).*sin(2*obj.phi) + cosh(2*obj.Eta0).*sin(4.*obj.phi))./(8*h2) + x.*y.*(bn-an) + bf.*xn.^2 + an.*xf.^2; 
%             
%              unn = ((f + sigma.*xi0).*(h.^4) - (b.*(xn.^2) + a.*(xf.^2)).*xi0ff - 2*x.*y.*(b-a).*xi1f - cn.*xi1 - cf.*xi0f)./(a.*(xn.^2) + b.*(xf.^2));
             
             %(a.*(xn.^2) + b.*(xf.^2)).*unn + (b.*(xn.^2) + a.*(xf.^2)).*xi0ff + 2*x.*y.*xi1f + c1.*xi1 + c2.*xi0f;
             
        
            [xi0,xi0f,xi0ff,xi0fff,xi0ffff] = Xi0.Derivatives();
             [xi1,xi1f,xi1ff,xi1fff] = Xi1.Derivatives();
             
			 %[F,Fn,Fnn,Ff,Fff] = Src.Derivatives(); Fnf=0;
             F = Src.Derivatives(); Fnf=0;
			 %[a,an,ann,a3n,af,aff,a3f,anf,anff,annf] = LapCoeffs.Derivatives('a');
			 %[b,bn,bnn,b3n,bf,bff,b3f,bnf,bnff,bnnf] = LapCoeffs.Derivatives('b');
             [a,an] = LapCoeffs.Derivatives('an');
             [~,af] = LapCoeffs.Derivatives('af');
             [b,bn] = LapCoeffs.Derivatives('bn');
             [~,bf] = LapCoeffs.Derivatives('bf');

             sigma	= LapCoeffs.Derivatives('sigma');
             [h,hn,hnn,h3n,h4n,hf,hff,h3f,h4f,hnf] = obj.MetricsAtScatterer.metrics();
             
             sigma_f=0;  sigma_ff=0; sigma_n=0; sigma_nf=0; sigma_nn=0;
             

			 h2 = h.^2;
             h3 = h.^3;
             h4 = h.^4;
             
             fd = obj.FocalDistance;
             %fd4 = fd^4;

             x  = fd*cosh(obj.Eta0).*cos(obj.phi);
             y  = fd*sinh(obj.Eta0).*sin(obj.phi);
             xn = fd*sinh(obj.Eta0).*cos(obj.phi);            
             yn = fd*cosh(obj.Eta0).*sin(obj.phi);
             %xf =-fd*cosh(obj.Eta0).*sin(obj.phi);
             %yf = fd*sinh(obj.Eta0).*cos(obj.phi);
             
             xf =-yn;%fd*cosh(eta).*sin(phi);
             yf = xn;%fd*sinh(eta).*cos(phi);
            
%            ynf = x;
%            xnf=-y; 
%            xnn = x;
%             x3n=xn;
%             ynn=y;           
%             y3n=yn;           
%             xnnf = xf;
%             xff=-x;
%             x3f=-xf;
%             xffn=-xn;             
%             ynnf=yf;
%             yff=-y;
%             y3f=-yf;
%             yffn=-ynf;
            
            
             C   = a.*xn.^2 + b.*xf.^2;
             iC  = 1./C;

             P   = a.*xf.^2 + b.*xn.^2;
             
             curv = (fd.^2).*cosh(obj.Eta0).*sinh(obj.Eta0)./(h3);
                          
             Xixx    = ((fd.^2).*sin(2.*obj.phi).^2 - h2.*cos(2.*obj.phi)).*curv./(h.^3);             
             
             Phxx    = ((fd.*xn./h).^2 - fd.^2/2).*sin(2*obj.phi)./h4 + sin(2*obj.phi).*curv.^2;
            
             
             unn = h.^4.*iC.*(F+sigma.*xi0)+(-1).*iC.*P.*xi0ff+2.*(a+(-1).*b).*iC.*x.*xi1f.*y+iC.*xi1.*((-1).*bn.* ...
                    xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+iC.*xi0f.*((-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y);
             
             
             res = xi0 + obj.deta.*xi1 + (obj.deta.^2).*unn/2 ;%+ (obj.deta.^3).*unnn/6;
             
         end
		function res = Expansion5thOrdrLap(obj,Xi0,Xi1,Src,LapCoeffs)
             %[xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
             %[xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
             
             [xi0,xi0f,xi0ff,xi0fff,xi0ffff] = Xi0.Derivatives();
             [xi1,xi1f,xi1ff,xi1fff] = Xi1.Derivatives();
             
			 [F,Fn,Fnn,Ff,Fff] = Src.Derivatives(); Fnf=0;
			 [a,an,ann,a3n,af,aff,a3f,anf,anff,annf] = LapCoeffs.Derivatives('a');
			 [b,bn,bnn,b3n,bf,bff,b3f,bnf,bnff,bnnf] = LapCoeffs.Derivatives('b');
			 sigma	= LapCoeffs.Derivatives('sigma');
             [h,hn,hnn,h3n,h4n,hf,hff,h3f,h4f,hnf] = obj.MetricsAtScatterer.metrics();
             
             sigma_f=0;  sigma_ff=0; sigma_n=0; sigma_nf=0; sigma_nn=0;
                   
             
			 h2 = h.^2;
             h3 = h.^3;
             h4 = h.^4;
             
             fd = obj.FocalDistance;
             %fd4 = fd^4;
             
             x  = fd*cosh(obj.Eta0).*cos(obj.phi);
             y  = fd*sinh(obj.Eta0).*sin(obj.phi);
             xn = fd*sinh(obj.Eta0).*cos(obj.phi);            
             yn = fd*cosh(obj.Eta0).*sin(obj.phi);
             %xf =-fd*cosh(obj.Eta0).*sin(obj.phi);
             %yf = fd*sinh(obj.Eta0).*cos(obj.phi);
             
             xf =-yn;%fd*cosh(eta).*sin(phi);
             yf = xn;%fd*sinh(eta).*cos(phi);
            
%            ynf = x;
%            xnf=-y; 
%            xnn = x;
%             x3n=xn;
%             ynn=y;           
%             y3n=yn;           
%             xnnf = xf;
%             xff=-x;
%             x3f=-xf;
%             xffn=-xn;             
%             ynnf=yf;
%             yff=-y;
%             y3f=-yf;
%             yffn=-ynf;

             
             C   = a.*xn.^2 + b.*xf.^2;
             Cn  = an.*xn.^2  + 2*a.*xn.*x + bn.*xf.^2 - 2.*b.*xf.*y;
             Cnn = ann.*xn.^2 + 4*an.*xn.*x + 2*a.*x.^2 + 2*a.*xn.^2 - 2*bn.*xf.*y + 2*b.*y.^2 + 2.*b.*xf.^2; 
             Cf  = af.*xn.^2  - 2*a.*xn.*y + bf.*xf.^2 - 2.*b.*xf.*x; 
             Cff = aff.*xn.^2 - 4*af.*xn.*y + 2*a.*y.^2 - 2*a.*xn.*yf + bff.*xf.^2 - 2*bf.*x.*xf.^2 - 2.*bf.*xf.*x + 2.*b.*x.^2 - 2.*b.*xf.^2;
             Cnf = anf.*xn.^2 - 2*an.*xn.*y + 2*af.*xn.*x - 2*a.*y.*x + 2*a.*xn.*xf + bnf.*xf.^2 - 2*bn.*xf.*x - 2.*bf.*xf.*y + 2.*b.*xf.*x - 2.*b.*xf.*yf;
             
             iC  = 1./C;
             iCf =-Cf./(C.^2);
             iCff=-(Cff.*C - 2*Cf.^2)./(C.^3);
             iCn =-Cn./(C.^2);
             iCnn=-(Cnn.*C - 2*Cn.^2)./(C.^3);
             iCnf=-(Cnf.*C - 2*Cn.*Cf)./(C.^3);

             P   = a.*xf.^2 + b.*xn.^2;
             Pn  = an.*xf.^2  - 2*a.*xf.*y + bn.*xn.^2 + 2.*b.*xn.*x;
             Pnn = ann.*xf.^2 - 2*an.*xf.*y - 2*an.*xf.*y + 2*a.*y.^2 - 2*a.*xf.*yn + bnn.*xn.^2 + 2*bn.*xn.*x + 2.*bn.*xn.*x + 2.*b.*x.^2 + 2.*b.*xn.^2;              
             Pf  = af.*xf.^2  - 2*a.*xf.*x + bf.*xn.^2 - 2.*b.*xn.*y;
             Pff = aff.*xf.^2 - 2*af.*xf.*x - 2*af.*xf.*x + 2*a.*x.^2 - 2*a.*xf.^2 + bff.*xn.^2 - 2*bf.*xn.*y - 2.*bf.*xn.*y + 2.*b.*y.^2 - 2.*b.*xn.*yf ;             
             Pnf = anf.*xf.^2 - 2*an.*xf.*x - 2*af.*xf.*y + 2*a.*x.*y - 2*a.*xf.*yf + bnf.*xn.^2 - 2*bn.*xn.*y + 2.*bf.*xn.*x - 2.*b.*y.*x + 2.*b.*xn.*xf;
             
             curv = (fd.^2).*cosh(obj.Eta0).*sinh(obj.Eta0)./(h3);
             curv_n = fd.^2.*cosh(obj.Eta0).^2.*h.^(-3) + fd.^2.*h.^(-3).*sinh(obj.Eta0).^2 - 3.*fd.^2.*cosh(obj.Eta0).*h.^(-4).*sinh(obj.Eta0).*hn;
             curv_f = (-3).*fd.^2.*cosh(obj.Eta0).*h.^(-4).*sinh(obj.Eta0).*hf;

             curv_nn = fd.^2.*cosh(obj.Eta0).*h.^(-3).*sinh(obj.Eta0) + 2.*fd.^2.*sinh(obj.Eta0).*(cosh(obj.Eta0).*h.^(-3) ...
                     - 3.*h.^(-4).*sinh(obj.Eta0).*hn)+fd.^2.*cosh(obj.Eta0).*(h.^(-3).*sinh(obj.Eta0)+(-6).*cosh(obj.Eta0).*h.^(-4).*hn ...
                     + sinh(obj.Eta0).*(12.*h.^(-5).*hn.^2+(-3).*h.^(-4).*hnn));

             curv_ff = fd.^2.*cosh(obj.Eta0).*sinh(obj.Eta0).*(12.*h.^(-5).*hf.^2+(-3).*h.^(-4).*hff);

             curv_nf =-3.*fd.^2.*cosh(obj.Eta0).^2.*h.^(-4).*hf - 3.*fd.^2.*h.^(-4).*sinh(obj.Eta0).^2.*hf+12.*fd.^2.*cosh(obj.Eta0).*h.^(-5).*sinh(obj.Eta0).*hf.*hn ...
                     - 3.*fd.^2.*cosh(obj.Eta0).*h.^(-4).*sinh(obj.Eta0).*hnf;
                          
             Xixx    = ((fd.^2).*sin(2.*obj.phi).^2 - h2.*cos(2.*obj.phi)).*curv./(h.^3);
             Xixx_n  = h.^(-4).*(cos(2.*obj.phi).*h.^2.*curv.*hn - 3.*fd.^2.*sin(2.*obj.phi).^2.*curv.*hn+(-1).*cos(2.*obj.phi).*h.^3.*curv_n+fd.^2.*h.*sin(2.*obj.phi).^2.*curv_n);

             Xixx_f  = h.^(-4).*(cos(2.*obj.phi).*h.^2.*curv.*hf+(-3).*fd.^2.*sin(2.*obj.phi).^2.*curv.*hf+h.^3.*(2.*sin(2.*obj.phi).*curv ...
                     - cos(2.*obj.phi).*curv_f) + fd.^2.*h.*sin(2.*obj.phi).*(4.*cos(2.*obj.phi).*curv+sin(2.*obj.phi).*curv_f));

             Xixx_nn = h.^(-5).*(12.*fd.^2.*sin(2.*obj.phi).^2.*curv.*hn.^2+cos(2.*obj.phi).*h.^3.*(2.*hn.*curv_n+curv.*hnn) ...
                     - 3.*fd.^2.*h.*sin(2.*obj.phi).^2.*(2.*hn.*curv_n+curv.*hnn) - cos(2.*obj.phi).*h.^4.*curv_nn ...
                     + h.^2.*((-2).*cos(2.*obj.phi).*curv.*hn.^2+fd.^2.*sin(2.*obj.phi).^2.*curv_nn));
             
             Xixx_ff = h.^(-5).*(12.*fd.^2.*sin(2.*obj.phi).^2.*curv.*hf.^2+h.^3.*(2.*cos(2.*obj.phi).* hf.*curv_f+curv.*((-4).*sin(2.*obj.phi).*hf ...
                     + cos(2.*obj.phi).*hff))+(-3).*fd.^2.*h.*sin(2.*obj.phi).*(2.*sin(2.*obj.phi).*hf.*curv_f+curv.*(8.*cos(2.*obj.phi).*hf+sin(2.*obj.phi).*hff)) ...
                     + h.^4.*(4.*cos(2.*obj.phi).*curv+4.*sin(2.*obj.phi).*curv_f - cos(2.*obj.phi).*curv_ff)+h.^2.*(curv.*(8.*fd.^2.*cos(4.*obj.phi) ...
                     - 2.*cos(2.*obj.phi).*hf.^2)+fd.^2.*(4.*sin(4.*obj.phi).*curv_f+sin(2.*obj.phi).^2.*curv_ff)));
             
             Xixx_nf = h.^(-5).*(12.*fd.^2.*sin(2.*obj.phi).^2.*curv.*hf.*hn+h.^3.*(cos(2.*obj.phi).*(curv_f.*hn+hf.*curv_n)+curv.*((-2).*sin(2.*obj.phi).*hn ...
                     + cos(2.*obj.phi).*hnf)) - 3.*fd.^2.*h.*sin(2.*obj.phi).*(sin(2.*obj.phi).*(curv_f.*hn+hf.*curv_n)+curv.*(4.*cos(2.*obj.phi).*hn ...
                     + sin(2.*obj.phi).*hnf))+h.^4.*(2.*sin(2.*obj.phi).*curv_n+(-1).*cos(2.*obj.phi).*curv_nf)+h.^2.*((-2).*cos(2.*obj.phi).*curv.*hf.*hn ...
                     + fd.^2.*(2.*sin(4.*obj.phi).*curv_n+sin(2.*obj.phi).^2.*curv_nf)));
             
             
             Phxx    = ((fd.*xn./h).^2 - fd.^2/2).*sin(2*obj.phi)./h4 + sin(2*obj.phi).*curv.^2;
             Phxx_n  = 2.*h.^(-7).*sin(2.*obj.phi).*(fd.^2.*h.^2.*hn+(-3).*fd.^2.*hn.*xn.^2+h.^7.*curv.*curv_n+fd.^2.*h.*xn.*x);
             
             Phxx_f  = h.^(-7).*((-1).*fd.^2.*cos(2.*obj.phi).*h.^3+2.*fd.^2.*h.^2.*sin(2.*obj.phi).*hf + 2.*h.^7.*curv.*(cos(2.*obj.phi).*curv+sin(2.*obj.phi).*curv_f) ...
                     - 6.*fd.^2.*sin(2.*obj.phi).*hf.*xn.^2+2.*fd.^2.*h.*xn.*(cos(2.*obj.phi).*xn - sin(2.*obj.phi).*y));
             
             Phxx_nn = 2.*h.^(-8).*sin(2.*obj.phi).*(21.*fd.^2.*hn.^2.*xn.^2+fd.^2.*h.^3.*hnn+(-3).*fd.^2.*h.*xn.*(xn.*hnn+4.*hn.*x) ...
                     + h.^8.*(curv_n.^2+curv.*curv_nn)+fd.^2.*h.^2.*((-5).*hn.^2+x.^2+xn.^2));

             
             Phxx_ff = 2.*h.^(-8).*(fd.^2.*h.^4.*sin(2.*obj.phi)+fd.^2.*h.^3.*(4.*cos(2.*obj.phi).*hf +sin(2.*obj.phi).*hff)+h.^8.*((-2).*sin(2.*obj.phi).*curv.^2 ...
                     + sin(2.*obj.phi).*curv_f.^2+curv.*(4.*cos(2.*obj.phi).*curv_f+sin(2.*obj.phi).*curv_ff))+ 21.*fd.^2.*sin(2.*obj.phi).*hf.^2.*xn.^2 ...
                     - 3.*fd.^2.*h.*xn.*(sin(2.*obj.phi).*hff.*xn+4.*hf.*(cos(2.*obj.phi).*xn- sin(2.*obj.phi).*y)) ...
                     - fd.^2.*h.^2.*(5.*sin(2.*obj.phi).*hf.^2+2.*sin(2.*obj.phi).*xn.^2 + sin(2.*obj.phi).*y.^2 ...
                     - xn.*(-4.*cos(2.*obj.phi).*y - sin(2.*obj.phi).*xn)));
             
             Phxx_nf = 2.*h.^(-8).*(21.*fd.^2.*sin(2.*obj.phi).*hf.*hn.*xn.^2+fd.^2.*h.^3.*(2.*cos(2.*obj.phi).*hn+sin(2.*obj.phi).*hnf) ...
                     + h.^8.*(sin(2.*obj.phi).*curv_f.*curv_n+curv.*(2.*cos(2.*obj.phi).*curv_n+sin(2.*obj.phi).*curv_nf)) - 3.*fd.^2.*h.*xn.*(2.*hn.*(cos(2.*obj.phi).*xn ...
                     - sin(2.*obj.phi).*y)+sin(2.*obj.phi).*(xn.*hnf+2.*hf.*x)) - fd.^2.*h.^2.*(5.*sin(2.*obj.phi).*hf.*hn + sin(2.*obj.phi).*y.*x ...
                     - xn.*(2.*cos(2.*obj.phi).*x+sin(2.*obj.phi).*xf))   );
             

             
             unn = h.^4.*iC.*(F+sigma.*xi0)+(-1).*iC.*P.*xi0ff+2.*(a+(-1).*b).*iC.*x.*xi1f.*y+iC.*xi1.*((-1).*bn.* ...
                    xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+iC.*xi0f.*((-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y);
              
              
             unnf =  h.^4.*iCf.*(F+sigma.*xi0)+iC.*(h.^4.*(Ff+sigma_f.*xi0)+4.*h.^3.*hf.*(F+sigma.*xi0)+(-1).*P.*xi0fff+ ...
                     2.*(a+(-1).*b).*x.*xi1ff.*y)+xi0ff.*((-1).*iCf.*P+iC.*((-1).*Pf+(-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).* ...
                     af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y))+xi1f.*(2.*(a+(-1).*b).*iCf.*x.*y+iC.*((-1).*bn.* ...
                     xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+3.*(af+(-1).*bf).*x.*y+2.*(a+(-1).*b).*xf.*y+2.* ...
                     (a+(-1).*b).*x.*yf))+xi1.*(iCf.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+( ...
                     -1).*bf).*x.*y)+iC.*((-1).*Xixx_f.*(a+(-1).*b).*h.^4+(-1).*bnf.*xf.^2+(-2).*bn.*xf.*(-x)+(-1).*(af+( ...
                     -1).*bf).*h.^4.*Xixx+(-4).*(a+(-1).*b).*h.^3.*hf.*Xixx+(-1).*anf.*xn.^2+(-2).*an.*xn.*(-y)+(aff+(-1) ...
                     .*bff).*x.*y+(af+(-1).*bf).*xf.*y+(af+(-1).*bf).*x.*yf))+xi0f.*(iCf.*((-1).*(a+(-1).*b).*h.^4.*Phxx+ ...
                     (-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y)+iC.*((-1).*Phxx_f.*(a+(-1).*b).*h.^4+(-1).*(af+ ...
                     (-1).*bf).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.*hf.*Phxx+h.^4.*sigma+(-1).*aff.*xf.^2+(-2).*af.*xf.* ...
                     (-x)+(-1).*bff.*xn.^2+(-2).*bf.*xn.*(-y)+(anf+(-1).*bnf).*x.*y+(an+(-1).*bn).*xf.*y+(an+(-1).*bn).*x.*yf));
               
               
            unnn = h.^4.*iCn.*(F+sigma.*xi0)+((-1).*iCn.*P+(-1).*iC.*Pn).*xi0ff+iC.*unn.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx ...
                   +(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+iC.*(h.^4.*(Fn+sigma_n.*xi0)+4.*h.^3.*hn.*(F+sigma.*xi0)+(-1).*P.*xi1ff+2.*(a+(-1).*b).*unnf.*x.*y) ...
                   +xi1f.*(2.*(a+(-1).*b).*iCn.*x.*y+iC.*((-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+3.*(an+(-1).*bn).*x.*y+2.*(a+(-1).*b).*xn.*y+ ...
                   2.*(a+(-1).*b).*x.*yn))+xi1.*(iCn.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+( ...
                   af+(-1).*bf).*x.*y)+iC.*((-1).*Xixx_n.*(a+(-1).*b).*h.^4+h.^4.*sigma+(-1).*bnn.*xf.^2+(-1).*(an+(-1) ...
                   .*bn).*h.^4.*Xixx+(-4).*(a+(-1).*b).*h.^3.*hn.*Xixx+(-1).*ann.*xn.^2+(-2).*bn.*xf.*(-y)+(-2).*an.* ...
                   xn.*x+(anf+(-1).*bnf).*x.*y+(af+(-1).*bf).*xn.*y+(af+(-1).*bf).*x.*yn))+xi0f.*(iCn.*((-1).*(a+(-1) ...
                   .*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y)+iC.*((-1).*Phxx_n.*(a+(-1).*b) ...
                   .*h.^4+(-1).*(an+(-1).*bn).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.*hn.*Phxx+(-1).*anf.*xf.^2+(-1).* ...
                   bnf.*xn.^2+(-2).*af.*xf.*(-y)+(-2).*bf.*xn.*x+(ann+(-1).*bnn).*x.*y+(an+(-1).*bn).*xn.*y+(an+(-1).*bn).*x.*yn));
             

            unnff = 8.*h.^3.*hf.*iCf.*(F+sigma.*xi0)+h.^4.*(2.*Ff.*iCf+F.*iCff+(2.*sigma_f.*iCf+iCff.*sigma).*xi0)+(-2) ...
                    .*iCf.*(P.*xi0fff+2.*((-1).*a+b).*x.*xi1ff.*y)+iC.*(Fff.*h.^4+8.*Ff.*h.^3.*hf+sigma_ff.*h.^4.*xi0+ ...
                    8.*sigma_f.*h.^3.*hf.*xi0+4.*h.^2.*(3.*hf.^2+h.*hff).*(F+sigma.*xi0)+(-2).*Pf.*xi0fff+(-1).*P.* ...
                    xi0ffff+4.*(af+(-1).*bf).*x.*xi1ff.*y+2.*(a+(-1).*b).*x.*xi1fff.*y+xi1ff.*((-1).*bn.*xf.^2+(-1).*(a+ ...
                    (-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+xi0fff.*((-1).*(a+(-1).*b).*h.^4.*Phxx+( ...
                    -1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y)+4.*(a+(-1).*b).*xi1ff.*(xf.*y+x.*yf))+xi0ff.*(( ...
                    -1).*iCff.*P+(-2).*iCf.*Pf+2.*iCf.*((-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+( ...
                    an+(-1).*bn).*x.*y)+iC.*((-1).*Pff+h.^4.*sigma+2.*((-1).*Phxx_f.*(a+(-1).*b).*h.^4+(-1).*(af+(-1).* ...
                    bf).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.*hf.*Phxx+2.*af.*x.*xf+(-1).*aff.*xf.^2+(-1).*bff.*xn.^2+( ...
                    -2).*bf.*xn.*(-y)+(anf+(-1).*bnf).*x.*y+(an+(-1).*bn).*xf.*y+(an+(-1).*bn).*x.*yf)))+xi1f.*(4.*(af+( ...
                    -1).*bf).*iCf.*x.*y+2.*(a+(-1).*b).*iCff.*x.*y+4.*(a+(-1).*b).*iCf.*xf.*y+2.*iCf.*((-1).*bn.*xf.^2+( ...
                    -1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+4.*(a+(-1).*b).*iCf.*x.*yf+2.*iC.* ...
                    ((-1).*Xixx_f.*(a+(-1).*b).*h.^4+2.*bn.*x.*xf+(-1).*bnf.*xf.^2+(-1).*(af+(-1).*bf).*h.^4.*Xixx+(-4) ...
                    .*(a+(-1).*b).*h.^3.*hf.*Xixx+(-1).*anf.*xn.^2+(-2).*an.*xn.*(-y)+2.*(aff+(-1).*bff).*x.*y+3.*(af+( ...
                    -1).*bf).*xf.*y+3.*(af+(-1).*bf).*x.*yf+(a+(-1).*b).*((-1).*x.*y+2.*xf.*yf+x.*(-y))))+xi1.*(iCff.*(( ...
                    -1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+2.*iCf.*((-1).* ...
                    Xixx_f.*(a+(-1).*b).*h.^4+2.*bn.*x.*xf+(-1).*bnf.*xf.^2+(-1).*(af+(-1).*bf).*h.^4.*Xixx+(-4).*(a+( ...
                    -1).*b).*h.^3.*hf.*Xixx+(-1).*anf.*xn.^2+(-2).*an.*xn.*(-y)+(aff+(-1).*bff).*x.*y+(af+(-1).*bf).*xf.* ...
                    y+(af+(-1).*bf).*x.*yf)+iC.*(4.*bnf.*x.*xf+(-1).*bnff.*xf.^2+(-2).*bn.*(x.^2+(-1).*xf.^2)+((-1).* ...
                    aff+bff).*h.^4.*Xixx+2.*((-1).*af+bf).*h.^3.*(Xixx_f.*h+4.*hf.*Xixx)+(-1).*(a+(-1).*b).*h.^2.*(h.*( ...
                    Xixx_ff.*h+8.*Xixx_f.*hf)+4.*(3.*hf.^2+h.*hff).*Xixx)+(-1).*anff.*xn.^2+(-4).*anf.*xn.*(-y)+(-2).* ...
                    an.*((-1).*xn.^2+(-y).^2)+(a3f+(-1).*b3f).*x.*y+2.*(aff+(-1).*bff).*(xf.*y+x.*yf)+(af+(-1).*bf).*(( ...
                    -1).*x.*y+2.*xf.*yf+x.*(-y))))+xi0f.*(2.*h.^4.*iCf.*sigma+iCff.*((-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).* ...
                    af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y)+2.*iCf.*((-1).*Phxx_f.*(a+(-1).*b).*h.^4+(-1).*(af+( ...
                    -1).*bf).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.*hf.*Phxx+2.*af.*x.*xf+(-1).*aff.*xf.^2+(-1).*bff.* ...
                    xn.^2+(-2).*bf.*xn.*(-y)+(anf+(-1).*bnf).*x.*y+(an+(-1).*bn).*xf.*y+(an+(-1).*bn).*x.*yf)+iC.*(2.* ...
                    sigma_f.*h.^4+((-1).*aff+bff).*h.^4.*Phxx+2.*((-1).*af+bf).*h.^3.*(Phxx_f.*h+4.*hf.*Phxx)+(-1).*(a+( ...
                    -1).*b).*h.^2.*(h.*(Phxx_ff.*h+8.*Phxx_f.*hf)+4.*(3.*hf.^2+h.*hff).*Phxx)+8.*h.^3.*hf.*sigma+4.* ...
                    aff.*x.*xf+(-1).*a3f.*xf.^2+(-2).*af.*(x.^2+(-1).*xf.^2)+(-1).*b3f.*xn.^2+(-4).*bff.*xn.*(-y)+(-2).* ...
                    bf.*((-1).*xn.^2+(-y).^2)+(anff+(-1).*bnff).*x.*y+2.*(anf+(-1).*bnf).*(xf.*y+x.*yf)+(an+(-1).*bn).*(( ...
                    -1).*x.*y+2.*xf.*yf+x.*(-y))));
             
             
             
             
           unnnf = 4.*h.^3.*(hn.*iCf+hf.*iCn).*(F+sigma.*xi0)+h.^4.*(Ff.*iCn+F.*iCnf+sigma_f.*iCn.*xi0+iCnf.*sigma.* ...
                   xi0+iCf.*(Fn+sigma_n.*xi0))+(-1).*P.*(iCn.*xi0fff+iCf.*xi1ff)+2.*(a+(-1).*b).*x.*(iCf.*unnf+iCn.* ...
                   xi1ff).*y+unn.*(iCf.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).* ...
                   x.*y)+iC.*((-1).*Xixx_f.*(a+(-1).*b).*h.^4+2.*bn.*x.*xf+(-1).*bnf.*xf.^2+(-1).*(af+(-1).*bf).*h.^4.* ...
                   Xixx+(-4).*(a+(-1).*b).*h.^3.*hf.*Xixx+(-1).*anf.*xn.^2+(-2).*an.*xn.*(-y)+(aff+(-1).*bff).*x.*y+(af+ ...
                   (-1).*bf).*xf.*y+(af+(-1).*bf).*x.*yf))+iC.*(Ff.*h.^4+4.*Fn.*h.^3.*hf+4.*Ff.*h.^3.*hn+sigma_nf.* ...
                   h.^4.*xi0+4.*sigma_n.*h.^3.*hf.*xi0+4.*sigma_f.*h.^3.*hn.*xi0+12.*h.^2.*hf.*hn.*(F+sigma.*xi0)+4.* ...
                   h.^3.*hnf.*(F+sigma.*xi0)+(-1).*Pn.*xi0fff+(-1).*Pf.*xi1ff+(-1).*P.*xi1fff+2.*(af+(-1).*bf).*unnf.* ...
                   x.*y+2.*(a+(-1).*b).*unnff.*x.*y+2.*(a+(-1).*b).*unnf.*xf.*y+2.*(an+(-1).*bn).*x.*xi1ff.*y+2.*(a+( ...
                   -1).*b).*xi1ff.*xn.*y+unnf.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1) ...
                   .*bf).*x.*y)+xi1ff.*((-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).* ...
                   x.*y)+2.*(a+(-1).*b).*unnf.*x.*yf+2.*(a+(-1).*b).*x.*xi1ff.*yn)+xi0ff.*((-1).*iCnf.*P+(-1).*iCn.*Pf+ ...
                   (-1).*iCf.*Pn+iCn.*((-1).*(a+(-1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.* ...
                   y)+iC.*((-1).*Phxx_n.*(a+(-1).*b).*h.^4+(-1).*(an+(-1).*bn).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.* ...
                   hn.*Phxx+(-1).*Pnf+(-1).*anf.*xf.^2+(-2).*bf.*x.*xn+(-1).*bnf.*xn.^2+(-2).*af.*xf.*(-y)+(ann+(-1).* ...
                   bnn).*x.*y+(an+(-1).*bn).*xn.*y+(an+(-1).*bn).*x.*yn))+xi1f.*(2.*(an+(-1).*bn).*iCf.*x.*y+2.*(af+( ...
                   -1).*bf).*iCn.*x.*y+2.*(a+(-1).*b).*iCnf.*x.*y+2.*(a+(-1).*b).*iCn.*xf.*y+2.*(a+(-1).*b).*iCf.*xn.* ...
                   y+iCn.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+iCf.*(( ...
                   -1).*(a+(-1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y)+2.*(a+(-1).*b).* ...
                   iCn.*x.*yf+2.*(a+(-1).*b).*iCf.*x.*yn+iC.*((-1).*Phxx_f.*(a+(-1).*b).*h.^4+(-1).*Xixx_n.*(a+(-1).*b) ...
                   .*h.^4+(-1).*(af+(-1).*bf).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.*hf.*Phxx+h.^4.*sigma+2.*(a+(-1).*b) ...
                   .*x.^2+2.*af.*x.*xf+(-1).*aff.*xf.^2+(-1).*bnn.*xf.^2+(-1).*(an+(-1).*bn).*h.^4.*Xixx+(-4).*(a+(-1) ...
                   .*b).*h.^3.*hn.*Xixx+(-2).*an.*x.*xn+(-1).*ann.*xn.^2+(-1).*bff.*xn.^2+(-2).*bn.*xf.*(-y)+(-2).*bf.* ...
                   xn.*(-y)+4.*(anf+(-1).*bnf).*x.*y+3.*(an+(-1).*bn).*xf.*y+3.*(af+(-1).*bf).*xn.*y+2.*(a+(-1).*b).* ...
                   (-y).*y+3.*(an+(-1).*bn).*x.*yf+2.*(a+(-1).*b).*xn.*yf+3.*(af+(-1).*bf).*x.*yn+2.*(a+(-1).*b).*xf.* ...
                   yn))+xi1.*(h.^4.*iCf.*sigma+iCnf.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+( ...
                   af+(-1).*bf).*x.*y)+iCn.*((-1).*Xixx_f.*(a+(-1).*b).*h.^4+2.*bn.*x.*xf+(-1).*bnf.*xf.^2+(-1).*(af+( ...
                   -1).*bf).*h.^4.*Xixx+(-4).*(a+(-1).*b).*h.^3.*hf.*Xixx+(-1).*anf.*xn.^2+(-2).*an.*xn.*(-y)+(aff+(-1) ...
                   .*bff).*x.*y+(af+(-1).*bf).*xf.*y+(af+(-1).*bf).*x.*yf)+iCf.*((-1).*Xixx_n.*(a+(-1).*b).*h.^4+(-1).* ...
                   bnn.*xf.^2+(-1).*(an+(-1).*bn).*h.^4.*Xixx+(-4).*(a+(-1).*b).*h.^3.*hn.*Xixx+(-2).*an.*x.*xn+(-1).* ...
                   ann.*xn.^2+(-2).*bn.*xf.*(-y)+(anf+(-1).*bnf).*x.*y+(af+(-1).*bf).*xn.*y+(af+(-1).*bf).*x.*yn)+iC.*( ...
                   sigma_f.*h.^4+(-1).*Xixx_nf.*(a+(-1).*b).*h.^4+(-1).*Xixx_n.*(af+(-1).*bf).*h.^4+(-1).*Xixx_f.*(an+( ...
                   -1).*bn).*h.^4+(-4).*Xixx_n.*(a+(-1).*b).*h.^3.*hf+(-4).*Xixx_f.*(a+(-1).*b).*h.^3.*hn+4.*h.^3.*hf.* ...
                   sigma+(af+(-1).*bf).*x.^2+2.*bnn.*x.*xf+(-1).*bnnf.*xf.^2+(-1).*(anf+(-1).*bnf).*h.^4.*Xixx+(-4).*( ...
                   an+(-1).*bn).*h.^3.*hf.*Xixx+(-4).*(af+(-1).*bf).*h.^3.*hn.*Xixx+(-12).*(a+(-1).*b).*h.^2.*hf.*hn.* ...
                   Xixx+(-4).*(a+(-1).*b).*h.^3.*hnf.*Xixx+(-2).*anf.*x.*xn+(-2).*an.*xf.*xn+2.*bn.*xf.*xn+(-1).*annf.* ...
                   xn.^2+(-2).*an.*x.*(-y)+2.*bn.*x.*(-y)+(-2).*bnf.*xf.*(-y)+(-2).*ann.*xn.*(-y)+(anff+(-1).*bnff).*x.*y+( ...
                   anf+(-1).*bnf).*xf.*y+(aff+(-1).*bff).*xn.*y+(af+(-1).*bf).*(-y).*y+(anf+(-1).*bnf).*x.*yf+(af+(-1).* ...
                   bf).*xn.*yf+(aff+(-1).*bff).*x.*yn+(af+(-1).*bf).*xf.*yn))+xi0f.*(h.^4.*iCn.*sigma+iCnf.*((-1).*(a+( ...
                   -1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y)+iCn.*((-1).*Phxx_f.*(a+(-1) ...
                   .*b).*h.^4+(-1).*(af+(-1).*bf).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.*hf.*Phxx+2.*af.*x.*xf+(-1).* ...
                   aff.*xf.^2+(-1).*bff.*xn.^2+(-2).*bf.*xn.*(-y)+(anf+(-1).*bnf).*x.*y+(an+(-1).*bn).*xf.*y+(an+(-1).* ...
                   bn).*x.*yf)+iCf.*((-1).*Phxx_n.*(a+(-1).*b).*h.^4+(-1).*(an+(-1).*bn).*h.^4.*Phxx+(-4).*(a+(-1).*b) ...
                   .*h.^3.*hn.*Phxx+(-1).*anf.*xf.^2+(-2).*bf.*x.*xn+(-1).*bnf.*xn.^2+(-2).*af.*xf.*(-y)+(ann+(-1).*bnn) ...
                   .*x.*y+(an+(-1).*bn).*xn.*y+(an+(-1).*bn).*x.*yn)+iC.*(sigma_n.*h.^4+(-1).*Phxx_nf.*(a+(-1).*b).* ...
                   h.^4+(-1).*Phxx_n.*(af+(-1).*bf).*h.^4+(-1).*Phxx_f.*(an+(-1).*bn).*h.^4+(-4).*Phxx_n.*(a+(-1).*b).* ...
                   h.^3.*hf+(-4).*Phxx_f.*(a+(-1).*b).*h.^3.*hn+(-1).*(anf+(-1).*bnf).*h.^4.*Phxx+(-4).*(an+(-1).*bn).* ...
                   h.^3.*hf.*Phxx+(-4).*(af+(-1).*bf).*h.^3.*hn.*Phxx+(-12).*(a+(-1).*b).*h.^2.*hf.*hn.*Phxx+(-4).*(a+( ...
                   -1).*b).*h.^3.*hnf.*Phxx+4.*h.^3.*hn.*sigma+(an+(-1).*bn).*x.^2+2.*anf.*x.*xf+(-1).*anff.*xf.^2+(-2) ...
                   .*bff.*x.*xn+2.*af.*xf.*xn+(-2).*bf.*xf.*xn+(-1).*bnff.*xn.^2+2.*af.*x.*(-y)+(-2).*bf.*x.*(-y)+(-2).* ...
                   aff.*xf.*(-y)+(-2).*bnf.*xn.*(-y)+(annf+(-1).*bnnf).*x.*y+(ann+(-1).*bnn).*xf.*y+(anf+(-1).*bnf).*xn.* ...
                   y+(an+(-1).*bn).*(-y).*y+(ann+(-1).*bnn).*x.*yf+(an+(-1).*bn).*xn.*yf+(anf+(-1).*bnf).*x.*yn+(an+(-1) ...
                   .*bn).*xf.*yn));



             
           unnnn = h.^4.*iCnn.*(F+sigma.*xi0)+2.*h.^3.*iCn.*(h.*(Fn+sigma_n.*xi0)+4.*hn.*(F+sigma.*xi0))+((-1).*iCnn.* ...
                   P+(-2).*iCn.*Pn+(-1).*iC.*Pnn).*xi0ff+((-2).*iCn.*P+(-2).*iC.*Pn).*xi1ff+iC.*unnn.*((-1).*bn.*xf.^2+ ...
                   (-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+iC.*(Fnn.*h.^4+8.*Fn.*h.^3.*hn+( ...
                   -1).*P.*unnff+sigma_nn.*h.^4.*xi0+8.*sigma_n.*h.^3.*hn.*xi0+4.*h.^2.*(3.*hn.^2+h.*hnn).*(F+sigma.* ...
                   xi0)+2.*(a+(-1).*b).*unnnf.*x.*y)+unnf.*(4.*(a+(-1).*b).*iCn.*x.*y+iC.*((-1).*(a+(-1).*b).*h.^4.* ...
                   Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+5.*(an+(-1).*bn).*x.*y+4.*(a+(-1).*b).*(xn.*y+x.*yn)))+unn.*( ...
                   2.*iCn.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+iC.*( ...
                   h.^4.*sigma+2.*((-1).*Xixx_n.*(a+(-1).*b).*h.^4+(-1).*bnn.*xf.^2+(-1).*(an+(-1).*bn).*h.^4.*Xixx+( ...
                   -4).*(a+(-1).*b).*h.^3.*hn.*Xixx+(-2).*an.*x.*xn+(-1).*ann.*xn.^2+(anf+(-1).*bnf).*x.*y+2.*bn.*xf.* ...
                   y+(af+(-1).*bf).*xn.*y+(af+(-1).*bf).*x.*yn)))+xi1f.*(2.*(a+(-1).*b).*iCnn.*x.*y+2.*iCn.*((-1).*(a+( ...
                   -1).*b).*h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+3.*(an+(-1).*bn).*x.*y+2.*(a+(-1).*b).*xn.*y+2.* ...
                   (a+(-1).*b).*x.*yn)+2.*iC.*((-1).*Phxx_n.*(a+(-1).*b).*h.^4+(-1).*(an+(-1).*bn).*h.^4.*Phxx+(-4).*( ...
                   a+(-1).*b).*h.^3.*hn.*Phxx+(-1).*anf.*xf.^2+(-2).*bf.*x.*xn+(-1).*bnf.*xn.^2+2.*(ann+(-1).*bnn).*x.* ...
                   y+2.*af.*xf.*y+3.*(an+(-1).*bn).*xn.*y+3.*(an+(-1).*bn).*x.*yn+(a+(-1).*b).*(2.*x.*y+2.*xn.*yn)))+ ...
                   xi1.*(iCnn.*((-1).*bn.*xf.^2+(-1).*(a+(-1).*b).*h.^4.*Xixx+(-1).*an.*xn.^2+(af+(-1).*bf).*x.*y)+2.* ...
                   iCn.*((-1).*Xixx_n.*(a+(-1).*b).*h.^4+h.^4.*sigma+(-1).*bnn.*xf.^2+(-1).*(an+(-1).*bn).*h.^4.*Xixx+( ...
                   -4).*(a+(-1).*b).*h.^3.*hn.*Xixx+(-2).*an.*x.*xn+(-1).*ann.*xn.^2+(anf+(-1).*bnf).*x.*y+2.*bn.*xf.* ...
                   y+(af+(-1).*bf).*xn.*y+(af+(-1).*bf).*x.*yn)+iC.*(2.*sigma_n.*h.^4+8.*h.^3.*hn.*sigma+(-1).*b3n.* ...
                   xf.^2+((-1).*ann+bnn).*h.^4.*Xixx+2.*((-1).*an+bn).*h.^3.*(Xixx_n.*h+4.*hn.*Xixx)+(-1).*(a+(-1).*b) ...
                   .*h.^2.*(h.*(Xixx_nn.*h+8.*Xixx_n.*hn)+4.*(3.*hn.^2+h.*hnn).*Xixx)+(-4).*ann.*x.*xn+(-1).*a3n.* ...
                   xn.^2+(-2).*an.*(x.^2+xn.^2)+(annf+(-1).*bnnf).*x.*y+4.*bnn.*xf.*y+(-2).*bn.*(xf.^2+y.^2)+2.*(anf+( ...
                   -1).*bnf).*(xn.*y+x.*yn)+(af+(-1).*bf).*(2.*x.*y+2.*xn.*yn)))+xi0f.*(iCnn.*((-1).*(a+(-1).*b).* ...
                   h.^4.*Phxx+(-1).*af.*xf.^2+(-1).*bf.*xn.^2+(an+(-1).*bn).*x.*y)+2.*iCn.*((-1).*Phxx_n.*(a+(-1).*b).* ...
                   h.^4+(-1).*(an+(-1).*bn).*h.^4.*Phxx+(-4).*(a+(-1).*b).*h.^3.*hn.*Phxx+(-1).*anf.*xf.^2+(-2).*bf.* ...
                   x.*xn+(-1).*bnf.*xn.^2+(ann+(-1).*bnn).*x.*y+2.*af.*xf.*y+(an+(-1).*bn).*xn.*y+(an+(-1).*bn).*x.*yn) ...
                   +iC.*(((-1).*ann+bnn).*h.^4.*Phxx+2.*((-1).*an+bn).*h.^3.*(Phxx_n.*h+4.*hn.*Phxx)+(-1).*(a+(-1).*b) ...
                   .*h.^2.*(h.*(Phxx_nn.*h+8.*Phxx_n.*hn)+4.*(3.*hn.^2+h.*hnn).*Phxx)+(-1).*annf.*xf.^2+(-4).*bnf.*x.* ...
                   xn+(-1).*bnnf.*xn.^2+(-2).*bf.*(x.^2+xn.^2)+(a3n+(-1).*b3n).*x.*y+4.*anf.*xf.*y+(-2).*af.*(xf.^2+ ...
                   y.^2)+2.*(ann+(-1).*bnn).*(xn.*y+x.*yn)+(an+(-1).*bn).*(2.*x.*y+2.*xn.*yn)));
             
                         
                           
                          
             
             
%              %assume a = b = const, sigma =0;
%              %unn   = f.*h2./a - xi0ff;
%              %unnn  = (h2.*fn + 2*f.*h.*hn)./a - xi1ff;
%              %unnff = (2*f.*(hf.^2) + h2.*f_ff + 2*h.*(2.*f_f.*hf + f.*hff))./a - xi0ffff;  
%              %unnnn = (2*f.*(hn.^2) + h2.*fnn  + 2*h.*(2*fn.*hn    + f.*hnn))./a - unnff;
             
% 
%                           			 unn  = ( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 )./a ;
%              
%                                        %atm assumimg sigma is constant
%              
%                                       unnn =  (Fn.*h2 + 2*f.*h.*hn - an.*unn - bfn.*xi0f - bn.*xi0ff - ann.*xi1 - bf.*xi1f - b.*xi1ff)./a ...
%                                              ... ( (fn + sigma.*xi1).*h2 + (f + sigma.*xi0).*(2*h.*hn) -  bfn.*xi0f  -  bf.*xi1f - bn.*xi0ff - b.*xi1ff - ann.*xi1 - an.*unn )./a ...
%                                             - ( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(an./(a.^2));
%              
%              
%                                        unnf = ( (Ff + sigma.*xi0f).*h2 + (f + sigma.*xi0).*(2*h.*hf) - bff.*xi0f - 2*bf.*xi0ff - b.*xi0fff - anf.*xi1  - an.*xi1f)./a ...
%                                             - ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*af./(a.^2);
%              
%                                        unnff= ( (Fff + sigma.*xi0f).*h2 + (f_f + sigma.*xi0f).*(3*h.*hf) + (f + sigma.*xi0).*(2*hf.^2 + 2*h.*hff) ...
%                                        - b3f.*xi0f - 3*bff.*xi0ff - 3*bf.*xi0fff - b.*xi0ffff - anff.*xi1 - 2*anf.*xi1f - an.*xi1ff)./a ...
%                                             - 2*( (f_f + sigma.*xi0f).*h2 + (f + sigma.*xi0).*(2*h.*hf) - bff.*xi0f - bf.*xi0ff - bf.*xi0ff - b.*xi0fff - anf.*xi1  - an.*xi1f).*af./(a.^2) ...
%                                             - ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*aff./(a.^2) ...
%                                             + 2*( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(af.^2)./(a.^3);
%              
%                                        unnnn= ( (Fnn + sigma.*unn).*h2 + (fn + sigma.*xi1).*(3*h.*hn) + (f + sigma.*xi0).*(2*hn.^2 + 2*h.*hnn) ...
%                                                  -  bfnn.*xi0f - 2*bfn.*xi1f - bf.*unnf - bnn.*xi0ff - 2*bn.*xi1ff - b.*unnff - a3n.*xi1 - 2*ann.*unn - an.*unnn )./a ...
%                                             - 2*( (fn + sigma.*xi1).*h2 + (f + sigma.*xi0).*(2*h.*hn) -  bfn.*xi0f  -  bf.*xi1f - bn.*xi0ff - b.*xi1ff - ann.*xi1 - an.*unn ).*an./(a.^2)...
%                                             - ( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*ann./(a.^2) ...
%                                             + 2*( (f + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(an.^2)./(a.^3);
              
              
             res = xi0 + obj.deta.*xi1 + (obj.deta.^2).*unn/2 + (obj.deta.^3).*unnn/6 + (obj.deta.^4).*unnnn/24;

        end
        
        function res = Expansion5thOrdrLap2(obj,Xi0,Xi1,Src,LapCoeffs)
            %[xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
            %[xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
            
            [xi0,xi0f,xi0ff,xi0fff,xi0ffff] = Xi0.Derivatives();
            [xi1,xi1f,xi1ff,xi1fff] = Xi1.Derivatives();
            
            [F,Fn,Fnn,Ff,Fff] = Src.Derivatives(); Fnf=0;
            [a,an,ann,a3n,af,aff,a3f,anf,anff,annf] = LapCoeffs.Derivatives('a');
            [b,bn,bnn,b3n,bf,bff,b3f,bnf,bnff,bnnf] = LapCoeffs.Derivatives('b');
            sigma	= LapCoeffs.Derivatives('sigma');
            [h,hn,hnn,h3n,h4n,hf,hff,h3f,h4f,hnf] = obj.MetricsAtScatterer.metrics();
            
            sigma_f=0;  sigma_ff=0; sigma_n=0; sigma_nf=0; sigma_nn=0;
            
            
            h2 = h.^2;
            h3 = h.^3;
            h4 = h.^4;
            
            fd = obj.FocalDistance;
            %fd4 = fd^4;
            
            x  = fd*cosh(obj.Eta0).*cos(obj.phi);
            y  = fd*sinh(obj.Eta0).*sin(obj.phi);
            xn = fd*sinh(obj.Eta0).*cos(obj.phi);
            yn = fd*cosh(obj.Eta0).*sin(obj.phi);

  
            
            
            %assume a = b = const, sigma =0;
            unn   = F.*h2./a - xi0ff;
            unnn  = (h2.*Fn + 2*F.*h.*hn)./a - xi1ff;
            unnff = (2*F.*(hf.^2) + h2.*Fff + 2*h.*(2.*Ff.*hf + F.*hff))./a - xi0ffff;
            unnnn = (2*F.*(hn.^2) + h2.*Fnn  + 2*h.*(2*Fn.*hn    + F.*hnn))./a - unnff;
            
            
%             unn  = ( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 )./a ;
%             
%             %atm assumimg sigma is constant
%             
%             unnn =  (Fn.*h2 + 2*F.*h.*hn - an.*unn - bnf.*xi0f - bn.*xi0ff - ann.*xi1 - bf.*xi1f - b.*xi1ff)./a ...
%                 ... ( (fn + sigma.*xi1).*h2 + (f + sigma.*xi0).*(2*h.*hn) -  bnf.*xi0f  -  bf.*xi1f - bn.*xi0ff - b.*xi1ff - ann.*xi1 - an.*unn )./a ...
%                 - ( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(an./(a.^2));
%             
%             
%             unnf = ( (Ff + sigma.*xi0f).*h2 + (F + sigma.*xi0).*(2*h.*hf) - bff.*xi0f - 2*bf.*xi0ff - b.*xi0fff - anf.*xi1  - an.*xi1f)./a ...
%                 - ( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*af./(a.^2);
% 
%             unnff= ( (Fff + sigma.*xi0f).*h2 + (Ff + sigma.*xi0f).*(3*h.*hf) + (F + sigma.*xi0).*(2*hf.^2 + 2*h.*hff) ...
%                 - b3f.*xi0f - 3*bff.*xi0ff - 3*bf.*xi0fff - b.*xi0ffff - anff.*xi1 - 2*anf.*xi1f - an.*xi1ff)./a ...
%                 - 2*( (Ff + sigma.*xi0f).*h2 + (F + sigma.*xi0).*(2*h.*hf) - bff.*xi0f - bf.*xi0ff - bf.*xi0ff - b.*xi0fff - anf.*xi1  - an.*xi1f).*af./(a.^2) ...
%                 - ( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*aff./(a.^2) ...
%                 + 2*( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(af.^2)./(a.^3);
% 
%             unnnn= ( (Fnn + sigma.*unn).*h2 + (Fn + sigma.*xi1).*(3*h.*hn) + (F + sigma.*xi0).*(2*hn.^2 + 2*h.*hnn) ...
%                 -  bnnf.*xi0f - 2*bnf.*xi1f - bf.*unnf - bnn.*xi0ff - 2*bn.*xi1ff - b.*unnff - a3n.*xi1 - 2*ann.*unn - an.*unnn )./a ...
%                 - 2*( (Fn + sigma.*xi1).*h2 + (F + sigma.*xi0).*(2*h.*hn) -  bnf.*xi0f  -  bf.*xi1f - bn.*xi0ff - b.*xi1ff - ann.*xi1 - an.*unn ).*an./(a.^2)...
%                 - ( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*ann./(a.^2) ...
%                 + 2*( (F + sigma.*xi0).*h2 -  bf.*xi0f - b.*xi0ff - an.*xi1 ).*(an.^2)./(a.^3);
%                          
%             
            res = xi0 + obj.deta.*xi1 + (obj.deta.^2).*unn/2 + (obj.deta.^3).*unnn/6 + (obj.deta.^4).*unnnn/24;
            
        end

        
        function res = Expansion(obj,Xi0,Xi1,Source,Coeffs)
			
			switch obj.ExpansionType
				case 25
					res = Expansion5thOrdrHelm(obj,Xi0,Xi1,Source,Coeffs);
				case 33
					res = Expansion3thOrdrLap(obj,Xi0,Xi1,Source,Coeffs);
				case 35
					res = Expansion5thOrdrLap(obj,Xi0,Xi1,Source,Coeffs);
                case 36
                    res = Expansion5thOrdrLap2(obj,Xi0,Xi1,Source,Coeffs);
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
