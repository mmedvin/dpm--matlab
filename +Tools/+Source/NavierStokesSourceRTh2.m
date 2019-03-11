classdef NavierStokesSourceRTh2 < Tools.Source.SuperSource
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
        %function [F,Fr,Frr,Ft,Ftt,Frt] = Derivatives(obj,Ex)
            function [F,Fr,Frr] = Derivatives(obj,Ex)
            if obj.IsDummy
                F=0;
                Fr=0;
                 Frr=0;
%                 Ft=0;
%                 Ftt=0;
% 
%                 Ftt=0;
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
            
            
            %[F,Fr,Frr,Ft,Ftt] = obj.CalcF(th,r,r0,sigma);
            [F,Fr,Frr] = obj.CalcF(th,r,r0,sigma);
        end
            
        
            
        function obj = NavierStokesSourceRTh2(Scatterer, CoeffsClsrHndl,CoeffsParams,ExParams)
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
           % S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
            S=zeros(obj.Scatterer.Size(1),obj.Scatterer.Size(2));
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
            
            
            [F,Fr,Frr] = obj.CalcF(th,r,r0,sigma);
            
            %S(obj.Scatterer.GridGamma) = F + obj.Scatterer.dn.*Fr  + (obj.Scatterer.dn.^2).*Frr/2;%taylor
            S(obj.Scatterer.GridGamma) = F + obj.Scatterer.dr.*Fr  + (obj.Scatterer.dr.^2).*Frr/2;%taylor
            
            tmp = obj.Derivatives();
            S(obj.Scatterer.Inside) = Fin(obj.Scatterer.Inside);   %was obj.Source(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
            %S=tmp;
        end    
    end
    
    methods(Access = protected)
        function [F,Fr,Frr,Ft,Ftt] = CalcF(obj,theta, r,r0,sigma)
            
            F = 2*(2 + sigma).*sin(r.*cos(theta)).*sin(r.*sin(theta));
            Fr = 2*(2 + sigma).*(cos(r.*sin(theta)).*sin(theta).*sin(r.*cos(theta)) + cos(theta).*cos(r.*cos(theta)).*sin(r.*sin(theta)));
            Frr = 2*(2 + sigma).*(cos(r.*cos(theta)).*cos(r.*sin(theta)).*sin(2*theta) - sin(r.*cos(theta)).*sin(r.*sin(theta)));
            
            Ft = 2*r.*(2 + sigma).*(cos(theta).*cos(r.*sin(theta)).*sin(r.*cos(theta)) - cos(r.*cos(theta)).*sin(theta).*sin(r.*sin(theta)));
            Ftt = -2*r.*(2 + sigma).*(sin(r.*cos(theta)).*(cos(r.*sin(theta)).*sin(theta) + r.*sin(r.*sin(theta))) + cos(r.*cos(theta)).*(r.*cos(r.*sin(theta)).*sin(2*theta) + cos(theta).*sin(r.*sin(theta))));
            
        end
    end
end

