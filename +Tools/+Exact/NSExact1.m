classdef NSExact1 < Tools.Exact.SuperNavierStokesExact
    % Exact Functions for Navier Stokes
    % Psi and Omega and their derivatives
    % Psi = x*y*(r^2 - r0^2)^2
    % Omega = Laplacian Psi
    
    methods(Static)
        
        function str = toString()
            str = 'x*y*(r^2 - r0^2)^2';
        end
        function P = Psi(theta,Params)
            r0=Params.r0;
            if numel(Params.r)==1
                r=ones(size(theta))*Params.r;
            else
                r = Params.r;
            end
            
            
            P = (r.^2).*cos(theta).*sin(theta).*(r.^2-r0.^2).^2;
        end
        
        function DP = DrDPsi(theta,Params)
            r0=Params.r0;
            if numel(Params.r)==1
                r=ones(size(theta))*Params.r;
            else
                r = Params.r;
            end
            
            DP = r.*(r.^2-r0.^2).*(3*r.^2-r0.^2).*sin(2*theta);
        end
        
        function DP = DPsiDThetaTheta(theta,Params)
            %         x = r .* cos(theta);
            %         y = r .* sin(theta);
            
            DP = -4*Tools.Exact.NSExact1.Psi(theta,Params);
            
        end
        
        function DP = DPsiD4Theta(theta,Params)
            %         x = r .* cos(theta);
            %         y = r .* sin(theta);
            
            DP = 16*Tools.Exact.NSExact1.Psi(theta,Params);
            
        end
        
        function DP = DPsiDrThetaTheta(theta,Params)
            %         x = r .* cos(theta);
            %         y = r .* sin(theta);
            
            DP = -4* Tools.Exact.NSExact1.DrDPsi(theta,Params);
            
        end
        
        function O = Omega(theta,Params)
            
            r0=Params.r0;
            if numel(Params.r)==1
                r=ones(size(theta))*Params.r;
            else
                r = Params.r;
            end
            
            O =  4*(r.^2).*(4*r.^2-3*r0^2).*sin(2*theta);
            
        end
        
        function DO = DrDOmega(theta,Params)
            
            r0=Params.r0;
            if numel(Params.r)==1
                r=ones(size(theta))*Params.r;
            else
                r = Params.r;
            end
            
            DO=8*r.*(8*r.^2-3*r0^2).*sin(2*theta);
        end
        
        function LO = LaplacianOmega(theta,Params)
            LO = 384*(Params.r.^2).*cos(theta).*sin(theta);
        end
        
    end
    
end