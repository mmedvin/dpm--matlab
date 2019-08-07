classdef NSExact2 < Tools.Exact.SuperNavierStokesExact
    % Exact Functions for Navier Stokes
    % Psi and Omega and their derivatives
    % Psi = sin(x).*sin(y)
    % Omega = Laplacian Psi
    
    methods(Static)
        
        function str = toString()
            str = 'sin(x).*sin(y)';
        end
        
        function P = Psi(theta,Params)
            r = Params.r;
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            P = sin(x).*sin(y);
        end
        
        function DP = DrDPsi(theta,Params)
            r = Params.r;
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            dxdr = cos(theta);
            dydr = sin(theta);
            
            DP = cos(x).*sin(y).*dxdr + cos(y).*sin(x).*dydr;
        end
        
        function DP = DPsiDThetaTheta(theta,Params)
            r = Params.r;
            r2 = r.*r;
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            dxdtheta = -y; %-r.*sin(theta);
            dydtheta =  x; %r.*cos(theta);
            
            %DP = -(r.*(sin(x).*(r.*sin(y) + cos(y).*sin(theta)) + cos(x).*(cos(theta).*sin(y) + r.*cos(y).*sin(2*theta))));
            DP = -x.^2.*sin(x).*sin(y)-y.^2.*sin(x).*sin(y)-x.*cos(x).*sin(y)-y.*cos(y).*sin(x)-x.*y.*cos(x).*cos(y).*2.0;
        end
        
        function DP = DPsiD4Theta(theta,Params)
            r = Params.r;
            r2 = r.*r;
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            DP = x.^2.*sin(x).*sin(y)+x.^4.*sin(x).*sin(y)+y.^2.*sin(x).*sin(y)+y.^4.*sin(x).*sin(y)+x.*cos(x).*sin(y)+y.*cos(y).*sin(x)+x.^3.*cos(x).*sin(y).*6.0+y.^3.*cos(y).*sin(x).*6.0+x.*y.*cos(x).*cos(y).*1.4e1+x.*y.^3.*cos(x).*cos(y).*4.0+x.^3.*y.*cos(x).*cos(y).*4.0-x.*y.^2.*cos(x).*sin(y).*6.0-x.^2.*y.*cos(y).*sin(x).*6.0+x.^2.*y.^2.*sin(x).*sin(y).*6.0;
            
        end
        
        function DP = DPsiDrThetaTheta(theta,Params)
            r = Params.r;
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            DP = -cos(y).*sin(theta).*sin(x)-cos(theta).*cos(x).*sin(y)-x.*cos(x).*cos(y).*sin(theta).*6.0-x.*cos(theta).*sin(x).*sin(y)-y.*sin(theta).*sin(x).*sin(y)-x.^2.*cos(theta).*cos(x).*sin(y)+x.^2.*cos(y).*sin(theta).*sin(x)-y.^2.*cos(y).*sin(theta).*sin(x)+x.*y.*cos(x).*sin(theta).*sin(y);
            
        end
        
        function O = Omega(theta,Params)
            
            O = -2*Tools.Exact.NSExact2.Psi(theta,Params);
            
        end
        
        function DO = DrDOmega(theta,Params)
            DO = -2*Tools.Exact.NSExact2.DrDPsi(theta,Params);
        end
        
        function LP = LaplacianOmega(theta,Params)
            r = Params.r;
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            LP = -2*Tools.Exact.NSExact2.Omega(theta,Params);
            
        end
    end
    
end