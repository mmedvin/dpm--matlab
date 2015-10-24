classdef PolarScatterer < Tools.Scatterer.SingleScatterer
    properties (Access = public)
        R;
        Th;
        r;
        th;
        
        dr;
        
        r0;
        
        Inside;
        Outside;
%         ScattererForSource;
        BasisArg;
        TheScatterer;
            
        MetricsAtScatterer = Tools.Metrics.PolarMetrics;
    end
    
%     properties(Access = private)
%         ExpansionType;
%     end
    
    methods
        
        function L = get.Inside(obj)
            L=obj.R <= obj.r0;
        end
        function L = get.Outside(obj)
            L = obj.R>obj.r0;
        end
        
        function BA = get.BasisArg(obj)
            BA = obj.th;
        end
        
        function TS = get.TheScatterer(obj)
            TS = struct('r',obj.r0,'th',obj.th);%*ones(size(obj.th))
        end
                
%         function S4S = get.ScattererForSource(obj)
%             GridF = Grids( ...
%                 obj.Grid.x1 - obj.Grid.dx/2 , ...
%                 obj.Grid.xn + obj.Grid.dx/2 , ...
%                 obj.Grid.Nx * 2 + 1         , ...
%                 obj.Grid.y1 - obj.Grid.dy/2 , ...
%                 obj.Grid.yn + obj.Grid.dy/2 , ...
%                 obj.Grid.Ny * 2 + 1         ) ;
%             
%             AddP = struct('r0',obj.r0);
%             S4S = PolarScatterer(GridF,AddP);
%         end      
        
        function obj = PolarScatterer(Grid,Params)
            obj = obj@Tools.Scatterer.SingleScatterer(Grid);
            obj.r0 = Params.r0;
            
            %obj.ExpansionType = Params.ExpansionType;            
            
            obj.MyGrid();
			
            obj.SplitGrid(Params.Stencil);
            	
            obj.GridGamma = intersect(obj.Np,obj.Nm)';
            
            obj.GridOnScatterer();
            
        end     
        
%          function res = Expansion5thOrdrHomoHelm(obj,Xi0,Xi1,F,WaveNumber)
%              [xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
%              [xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
%              
%              dr2=obj.dr.^2;
%              dr3=obj.dr.^3;
%              dr4=obj.dr.^4;
%              k2=WaveNumber.k.^2;
%              
%              R2=obj.r0.^2;
%              % R3=R.^3;
%              R4=obj.r0.^4;
%              
%              res = (1 - dr2.*k2./2 + dr3.*k2./obj.r0./6  - dr4.*k2.*(3./R2 - k2)/24).*xi0 ...
%                  +  (obj.dr - dr2./2./obj.r0 + dr3.*(2./R2-k2)./6 - dr4.*(3./R2 - k2)./12./obj.r0).*xi1 ...
%                  + xi0tt.*dr2.*(-1 + obj.dr./obj.r0 - dr2.*(11./R2 - 2.*k2)./12)./R2./2 ...
%                  + xi1tt.*dr3.*( obj.dr/obj.r0/2 - 1/3)./R2./2 ...
%                  + xi0tttt.*dr4./R4./24;
%          end
%          
%          function res = Expansion7thOrdrHomoHelm(obj,Xi0,Xi1,F,WaveNumber)
%              [xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
%              [xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
%              
%              %[k,kr,krr,k3r,k4r,k5r] = calc_k(k0,r,r0);
%              [k,kr,krr,k3r] = WaveNumber.Derivatives();
%              %k=0; kr=0; krr=0; k3r=0; %k4r=0; k5r=0;
%              
%              %  [f,fr,frtt,frr,ftt] = F.Derivatives();%verify!!!
%              f=0; fr=0; frr=0; f3r=0; ftt=0; frtt=0; ftttt=0;fttrr=0;
%              
%              xirr = f - xi1./obj.r0 - xi0tt./(obj.r0.^2) - xi0.*(k.^2);
%              xi3r = fr + xi1./(obj.r0.^2) - xirr./obj.r0 + 2.*xi0tt./(obj.r0.^3) - xi1tt./(obj.r0.^2) - (k.^2).*xi1 - 2.*k.*kr.*xi0;
%              xittrr = ftt - xi1tt./obj.r0 - xi0tttt./(obj.r0.^2) - (k.^2).*xi0tt;
%              xi4r = frr - 2.*(kr.^2 + k.*krr).*xi0 - (2./(obj.r0.^3)+ 4.*k.*kr).*xi1 + (2./(obj.r0.^2)-k.^2).*xirr - xi3r./obj.r0 - 6.*xi0tt./(obj.r0.^4) + 4.*xi1tt./(obj.r0.^3) - xittrr./(obj.r0.^2);
%              
%              xitt3r = frtt + xi1tt./(obj.r0.^2) - xittrr./obj.r0 + 2.*xi0tttt./(obj.r0.^3) - xi1tttt./(obj.r0.^2)-...
%                  (k.^2).*xi1tt - 2.*k.*kr.*xi0tt;
%              
%              xi5r = f3r - 2.*(3.*kr.*krr+k.*k3r).*xi0 - 6.*(-1./(obj.r0./4)+kr.^2+k.*krr).*xi1-...
%                  6.*(1./(obj.r0.^3) + k.*kr).*xirr + (3./(obj.r0.^2) - k.^2).*xi3r - xi4r./obj.r0 +...
%                  24.*xi0tt./(obj.r0.^5) - 18.*xi1tt./(obj.r0.^4)+ 6.*xittrr./(obj.r0.^3)  - xitt3r./(obj.r0.^2);
%              
%              %                res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 + (obj.dr.^5).*xi5r/120;
%                           
%              xi4trr = ftttt - xi1tttt./obj.r0 - xi0tttttt./(obj.r0.^2) - xi0tttt.*(k.^2);
%              xitt4r = fttrr - 2.*xi1tt./(obj.r0.^3) + (2./(obj.r0.^2)-k.^2).*xittrr - xitt3r./obj.r0 - 6*xi0tttt./(obj.r0.^4) + 4*xi1tttt./(obj.r0.^3) - xi4trr./(obj.r0.^2);
%              
%              xi6r = 24*xirr./(obj.r0.^4) - 24*xi1./(obj.r0.^5) - 12*xi3r./(obj.r0.^3) + (4./(obj.r0.^2) - k^2).*xi4r - xi5r./obj.r0 ...
%                  - 120*xi0tt./(obj.r0.^6) + 96*xi1tt./(obj.r0.^5) - 36*xittrr./(obj.r0.^4) + 8*xitt3r./(obj.r0.^3) - xitt4r./(obj.r0.^2);
%                           
%              res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 + (obj.dr.^5).*xi5r/120 + (obj.dr.^6).*xi6r/720;
%          end
%          
%          function res = Expansion5thOrdrHelm(obj,Xi0,Xi1,F,WaveNumber)
%              [xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
%              [xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
%              
%              [f,fr,frr,ftt] = F.Derivatives();
%              [k,kr,krr] = WaveNumber.Derivatives('r');
%              [~,kt,ktt] = WaveNumber.Derivatives('t');
%              
%              
%              %     function xg = xigamma(th,xi0,xi1,xi0t,xi0tt,xi1tt,xi0tttt,R,dr,f,fr,frr,ftt)
%              
%              %   global k0 NHR
%              %    [k,kr,krr,~,~,kt,ktt] = calc_k(k0,R,th,NHR);
%              
%              xirr = f - xi1./obj.r0 - xi0tt./(obj.r0.^2) - xi0.*(k.^2);
%              xi3r = fr + xi1./(obj.r0.^2) - xirr./obj.r0 + 2*xi0tt./(obj.r0.^3) - xi1tt./(obj.r0.^2) ...
%                  - (k.^2).*xi1 - 2.*k.*kr.*xi0;
%              
%              k2xi0_tt = (k.^2).*xi0tt + 4*k.*kt.*xi0t + (2*(kt.^2)+2*k.*ktt).*xi0 ;
%              xittrr = ftt - xi1tt./obj.r0 - xi0tttt./(obj.r0.^2) - k2xi0_tt;%(k.^2).*xi0tt;
%              
%              xi4r = frr - 2*(kr.^2 + k.*krr).*xi0 - (2./(obj.r0.^3)+ 4*k.*kr).*xi1 ...
%                  + (2/(obj.r0.^2)-k.^2).*xirr - xi3r./obj.r0 - 6*xi0tt./(obj.r0.^4) ...
%                  + 4*xi1tt./(obj.r0.^3) - xittrr./(obj.r0.^2);
%              
%              res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 ;
% 
%          end
%          
%          function res = Expansion3thOrdrLap(obj,Xi0,Xi1,F,LapCoeffs)
%              [xi0,xi0t,xi0tt] = Xi0.Derivatives();
%              [xi1,xi1t,xi1tt] = Xi1.Derivatives();
%              
% 			 f      = F.Derivatives();
% 			 [a,ar] = LapCoeffs.Derivatives('ar');
%              [~,at] = LapCoeffs.Derivatives('at');
%              [b,br] = LapCoeffs.Derivatives('br');
%              [~,bt] = LapCoeffs.Derivatives('bt');
% 			 sigma	= LapCoeffs.Derivatives('sigma');
% 			 
% 			% urr2 = (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2))./a - xi1./obj.r0 ;
% 			 
%              sint = sin(obj.th);
%              cost = cos(obj.th);
%              sincost = sint.*cost;
%              sint2 = sint.^2; 
%              cost2 = cost.^2;
%              
%              urr = (f + sigma.*xi0 ... 
%              -  (xi0t.*( (obj.r0.*(br-ar) + 2*(a-b) ).*sincost   + bt.*cost2  +at.*sint2)./(obj.r0.^2)  ... 
%             + xi1.* ( (bt-at).*sincost  +(obj.r0.*ar+b).*cost2 + (a +obj.r0.* br).*sint2 )./obj.r0 ...
%             + 2*sincost./obj.r0 .* (b-a).* xi1t  ...
%             + (b.*cost2 + a.*sint2).*xi0tt./(obj.r0.^2)) )./(a.*cost2 + b.*sint2);
% 	 
% 		 
%              res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*urr/2;% + (obj.dr.^3).*u3r/6;
%              
%          end
%          
%          function res = Expansion5thOrdrHomoLap(obj,Xi0,Xi1,F,LapCoeffs)
%              [xi0,xi0t,xi0tt,~,xi0tttt,xi0tttttt] = Xi0.Derivatives();
%              [xi1,~,xi1tt,~,xi1tttt,~] = Xi1.Derivatives();
%              assert('tbd')
%              
% 			 
% 			 % 			 %this is temporary part
% 			 % 			 %assuming sigma constant or at least doesn't depends on r
% 			 % 			 br = LapCoeffs.br;
% 			 % 			 btr= LapCoeffs.btr;
% 			 % 			 btrr= LapCoeffs.btrr;
% 			 % 			 u3r = (fr + sigma.*xi1 - arr.*xi1 - ar.*urr -  (btr.*xi0t + bt.*xi1t + br.*xi0tt + b.*xi1tt)./(obj.r0.^2) + 2*(bt.*xi0t + b.*xi0tt)./(obj.r0.^3) )./a ...
% 			 % 				 + (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2)).*ar./a./a ...
% 			 % 				 - urr./obj.r0 + xi1./obj.r0.^2 ;
% 			 %
% 			 %
% 			 % 			 u4r =   (frr + sigma.*urr - a3r.*xi1  - 2*arr.*urr - ar.*u3r ...
% 			 % 					-(btrr.*xi0t + btr.*xi1t + btr.*xi1t +bt.*urrt
% 			 % 				+ br.*xi0tt + b.*xi1tt)./(obj.r0.^2) +  2*(btr.*xi0t + bt.*xi1t + br.*xi0tt + b.*xi1tt)./(obj.r0.^3)
% 			 %
% 			 % 			 + 2*(bt.*xi0t + b.*xi0tt)./(obj.r0.^3) )./a ...
% 			 % 				 + (f + sigma.*xi0 - ar.*xi1 -  (bt.*xi0t + b.*xi0tt)./(obj.r0.^2)).*ar./a./a ...
% 			 % 				 - urr./obj.r0 + xi1./obj.r0.^2 ;
% 			 
%              res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 ;
%          end
%         
          function res = Expansion(obj,Xi0,Xi1,F,Coeffs)
%              
% 			 switch obj.ExpansionType
% 				 case 15 
% 					 res = Expansion5thOrdrHomoHelm(obj,Xi0,Xi1,F,Coeffs);
% 				 case 17
% 					 res = Expansion7thOrdrHomoHelm(obj,Xi0,Xi1,F,Coeffs);
% 				 case 25
% 					 res = Expansion5thOrdrHelm(obj,Xi0,Xi1,F,Coeffs);
% 				 case 33
% 					 res = Expansion3thOrdrLap(obj,Xi0,Xi1,F,Coeffs);
% 				 case 35
% 					 res = Expansion5thOrdrHomoLap(obj,Xi0,Xi1,F,Coeffs);
% 			 end
          end
    end
    
    
    methods(Access = protected)
        
        function GridOnScatterer(obj)
            obj.r  = obj.R(obj.GridGamma);
            obj.dr = obj.r - obj.r0;
            obj.th = obj.Th(obj.GridGamma);
        end
        
        function MyGrid(obj)
            Z = obj.Grid.Z();
            obj.R=abs(Z);
            obj.Th=angle(Z);
        end	
    end
end
