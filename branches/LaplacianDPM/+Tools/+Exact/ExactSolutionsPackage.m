classdef ExactSolutionsPackage
    %ExactSolutionsPackage to be expanded, but the main idea is to keep all
    % exact solutions in one class.
    %   Detailed explanation goes here
    
    properties
        u;
        dudn;
    end
    
    methods(Static=false)
        function obj = ExactSolutionsPackage(Type,Params)
            HankOrPlane = Type;
            if strcmpi(HankOrPlane,'PlaneWave')                        
            obj.u          = @(phi) -ExactSolutionsPackage.Uinc(Params.ExParams,phi,Params.IncAng,Params.k);
            obj.dudn    = @(phi) -ExactSolutionsPackage.detaUinc(Params.ExParams,phi,Params.IncAng,Params.k);           
        elseif strcmpi(HankOrPlane,'Hankel')
            f1          = @(phi) ExactSolutionsPackage.ExactHank(Params.ExParams,phi,Params.k);
            dfdn    = @(phi) ExactSolutionsPackage.dnExactHank(Params.ExParams,phi,Params.k);
        end
        end
    end
    methods(Static=true)
        function uinc = Uinc(Params,phi,IncAng,k)
            
            if strcmpi(Params.ScattererType,'ellipse')
                x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
                y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
            elseif strcmpi(Params.ScattererType,'circle')
                x = Params.r .* cos(phi);
                y = Params.r .* sin(phi);
            elseif strcmpi(Params.ScattererType,'submarine')
                e = (cos(phi).^2/Params.ellipse.a^2 + sin(phi).^2/Params.ellipse.b^2);
                r = sqrt((1 + Params.tower.c*sin(phi).^Params.tower.p) ./ e);
                
                x = r .* cos(phi);
                y = r .* sin(phi);
                
            end
            
            uinc = exp( 1i.* k .* (x.*cos(IncAng) + y.*sin(IncAng)) );
        end
        
        function duinc = detaUinc(Params,phi,IncAng,k)
            
            if strcmpi(Params.ScattererType,'ellipse')
                dx = Params.FocalDistance * sinh(Params.eta) .* cos(phi);
                dy = Params.FocalDistance * cosh(Params.eta) .* sin(phi);
            elseif strcmpi(Params.ScattererType,'circle')
                dx = cos(phi);
                dy = sin(phi);
            elseif strcmpi(Params.ScattererType,'submarine')
                e = (cos(phi).^2/Params.ellipse.a^2 + sin(phi).^2/Params.ellipse.b^2);
                r = sqrt((1 + Params.tower.c*sin(phi).^Params.tower.p) ./ e);
                %         dr = Params.tower.c * Params.tower.p * cos(phi).*sin(phi).^(Params.tower.p-1)./e ...
                %             - 2*(-1/Params.ellipse.a^2 + 1/Params.ellipse.b^2).*cos(phi).* sin(phi).*r.^2 ./e;
                %         dr = dr/2/r;
                de = (-1/Params.ellipse.a^2 + 1/Params.ellipse.b^2).* sin(2*phi);
                dg = Params.tower.c*Params.tower.p*cos(phi).*sin(phi).^(-1 + Params.tower.p);
                dr = (dg - de.*r.^2)/2./r./e;
                
                dx = cos(phi).*dr;
                dy = sin(phi).*dr;
                
            end
            
            uinc = ExactSolutionsPackage.Uinc(Params,phi,IncAng,k)  ;
            
            duinc = 1i .* k .*  uinc .* (dx.*cos(IncAng) + dy.*sin(IncAng));
            %   h = FocalDist*sqrt(sinh(eta).^2 + sin(phi).^2);
            % duinc = duinc./h;
            
        end
        
        
        function u = ExactHank(Params,phi,k)
            if strcmpi(Params.ScattererType,'ellipse')
                x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
                y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
            elseif strcmpi(Params.ScattererType,'circle')
                x = Params.r .* cos(phi);
                y = Params.r .* sin(phi);
                %         r= Params.r;
                %         th = phi;
            end
            
            Z=x+1i*y;
            r=abs(Z);
            th=angle(Z);
            
            
            
            bn=Params.HankelIndex;
            u = besselh(bn,Params.HankelType,k*r).*exp(bn*1i*th);
        end
        
        function un = dnExactHank(Params,phi,k)
            
            if strcmpi(Params.ScattererType,'ellipse')
                x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
                y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
            elseif strcmpi(Params.ScattererType,'circle')
                x = Params.r .* cos(phi);
                y = Params.r .* sin(phi);
                %         r= Params.r;
                %         th = phi;
            end
            
            Z=x+1i*y;
            r=abs(Z);
            th=angle(Z);
            
            bn=Params.HankelIndex;
            
            dudr  = k*0.5*(besselh(bn-1,Params.HankelType,k*r) - besselh(bn+1,Params.HankelType,k*r)).*exp(bn*1i*th);
            
            if strcmpi(Params.ScattererType,'ellipse')
                
                dudth  = 1i*bn*besselh(bn,Params.HankelType,k*r).*exp(bn*1i*th);
                
                drdeta = Params.FocalDistance^2 * cosh(Params.eta).*sinh(Params.eta) ./r;
                dthdeta = Params.FocalDistance^2 * cos(phi).*sin(phi) ./(r.^2);
                un = dudr.*drdeta + dudth.*dthdeta;
                
                %         SM = EllipticalMetrics(Params.FocalDistance,Params.eta,phi);
                % un = un.*SM.h;
            elseif strcmpi(Params.ScattererType,'circle')
                un = dudr;
            end
        end
        
        
    end
    
end

