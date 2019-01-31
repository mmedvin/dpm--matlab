classdef NavierStokesSourceR < Tools.Source.SuperSource
    properties (Dependent = true)
        Source;%Fn;Ff;Fnn;Fff;    % WN;
    end
    
    properties%(Access = protected)       
        Scatterer; 
        Exact;
        IsExactReady=false;
        
		CoeffsClsrHndl;
		CoeffsParams;
        ExParams;
    end
    
    methods
        function [F,Fr,Frr,Ft,Ftt] = Derivatives(obj,Ex)
            
            if obj.IsDummy
                F=0;
                Fr=0;
                Frr=0;
                %F3r=0;
                Ft=0;
                Ftt=0;
            else
                
                try
                    r = obj.Scatterer.R;
                    th = obj.Scatterer.Th;
                catch exception
                    if strcmp(exception.identifier,'MATLAB:nonExistentField')
                        r  = obj.Scatterer.r;
                        th = obj.Scatterer.th;
                    else
                        rethrow(exception);
                    end
                end
                
                if numel(r)==1
                    r=ones(size(th))*r;
                end
                
            end
            
            coeffs = obj.CoeffsClsrHndl(obj.Scatterer,obj.CoeffsParams);
            sigma = coeffs.Derivatives('sigma');
            r0 = obj.ExParams.r0;
            p = obj.ExParams.p;
            
            [F,Fr,Frr] = obj.CalcF(th,r,r0,sigma,p);
            Ft=0;Ftt=0;
        end
            
        
            
        function obj = NavierStokesSourceR(Scatterer, CoeffsClsrHndl,CoeffsParams,ExParams)
            obj.Scatterer = Scatterer;
            obj.CoeffsClsrHndl = CoeffsClsrHndl;
            obj.CoeffsParams = CoeffsParams;
            obj.ExParams = ExParams;
            obj.IsDummy = false;
        end
        
        
        function S = get.Source(obj)
            S=obj.getSource();
        end
        function S = getSource(obj)
            S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
            
            %ScatGamma = struct('FocalDistance',obj.Scatterer.FocalDistance,'Eta', obj.Scatterer.Eta0,'Phi', obj.Scatterer.Phi(obj.Scatterer.GridGamma));
            
            Fin = obj.Derivatives();
            
            r  = obj.Scatterer.r;
            th = obj.Scatterer.th;
            
            if numel(r)==1
                r=ones(size(th))*r;
            end


            coeffs = obj.CoeffsClsrHndl(obj.Scatterer,obj.CoeffsParams);
            sigma = coeffs.Derivatives('sigma');
            r0 = obj.ExParams.r0;
            p = obj.ExParams.p;
            
            [F,Fr,Frr] = obj.CalcF(th,r,r0,sigma,p);
            
            %S(obj.Scatterer.GridGamma) = F + obj.Scatterer.dn.*Fr  + (obj.Scatterer.dn.^2).*Frr/2;%taylor
            S(obj.Scatterer.GridGamma) = F + obj.Scatterer.dr.*Fr  + (obj.Scatterer.dr.^2).*Frr/2;%taylor
            
            tmp = obj.Derivatives();
            S(obj.Scatterer.Inside) = Fin(obj.Scatterer.Inside);   %was obj.Source(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
            S=tmp;
        end    
    end
    
    methods(Access = protected)
        function [F,Fr,Frr] = CalcF(obj,theta, r,r0,sigma,p)

            if p==2
                F = 64 - 8*(2*r.^2 - r0.^2)*sigma;
                Fr = -32*r*sigma;
                Frr = -32*sigma*ones(size(r));


            elseif p>6
                c1 = (r.^2 - r0^2);
                
                F = 4*p*(c1.^(p-4)) .* (4*(p-1)*((p-1)*p*r.^4 - 4*(p-1)*(r.^2)*(r0^2) + 2*r0^4) + (c1.^2).*(p*r.^2 - r0^2) *sigma);
                
                Fr = 8*(p-1)*p*r.*(c1.^(p-5)).*(4*(p^3)*r.^4 - 12*(p^2)*(r.^4 + 2*(r.^2)*(r0^2) ) ...
                    - 2*(r0^2)*(sigma*r.^4  + (r.^2)*(24 - 2*sigma*r0^2) + (r0^2)*(24 + sigma*r0^2)) ...
                    + p*(24*r0^4 + sigma*r.^6 + (r.^4)*(8 - 2*sigma*r0^2 ) + (r.^2)*(r0^2)*(72 + sigma*r0^2 )));
                
                Frr = 8*(p-1)*p*(c1.^(p-6)).* (4*(p-2)*((p-1)*p*(2*p-5)*r.^6 - (p-1)*(17*p-42)*(r.^4)*(r0^2)...
                    + 6*(5*p-12)*(r.^2)*(r0^4) - 6*r0^6) + (c1.^2).*(p*(2*p-3)*r.^4 + (10 - 7*p)*(r.^2)*(r0^2) + 2*r0^4)*sigma);
            else
                assert(true);
            end
           
        end
    end
    
end

