classdef IncidentWave
    % this module is taken from Radar project, 
    % and probably is similar to Uinc functions inside various Helmholtz runners
    % however it needs to be verified before replaced there
    % TODO: use it in Helmholtz runners
    
    methods(Static)
        
        function [uinc,uinc_t,uinc_tt,uinc_3t,uinc_4t] = PlaneWave(Params,phi,IncAng,k)
            
            IsStarshaped = false;
            
            if strcmpi(Params.ScattererType,'ellipse')
                x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
                y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
            elseif strcmpi(Params.ScattererType,'circle')
                x = Params.r .* cos(phi);
                y = Params.r .* sin(phi);
            elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
                IsStarshaped = true;
                try
                    x = Params.Parameterization.XHandle.Derivatives(phi);
                    y = Params.Parameterization.YHandle.Derivatives(phi);
                catch
                    x = Params.r.*cos(phi);
                    y = Params.r.*sin(phi);
                end
            end
            
            uinc = exp( 1i.* k .* (x.*cos(IncAng) + y.*sin(IncAng)) );
            
            if nargout > 1 && IsStarshaped
                %try
                [x,xt,xtt,x3t,x4t] = Params.Parameterization.XHandle.Derivatives(phi);
                [y,yt,ytt,y3t,y4t] = Params.Parameterization.YHandle.Derivatives(phi);
                %catch
                %x = Params.r.*cos(phi);
                %y = Params.r.*sin(phi);
                %end
                
                uinc_t  = 1i .* k .*  uinc    .* (xt .*cos(IncAng) + yt .*sin(IncAng));
                uinc_tt = 1i .* k .* (uinc_t  .* (xt .*cos(IncAng) + yt .*sin(IncAng)) +     uinc    .* (xtt.*cos(IncAng) + ytt.*sin(IncAng)) );
                uinc_3t = 1i .* k .* (uinc_tt .* (xt .*cos(IncAng) + yt .*sin(IncAng)) + 2 * uinc_t  .* (xtt.*cos(IncAng) + ytt.*sin(IncAng)) + uinc .* (x3t.*cos(IncAng) + y3t.*sin(IncAng)));
                uinc_4t = 1i .* k .* (uinc_3t .* (xt .*cos(IncAng) + yt .*sin(IncAng)) + 3 * uinc_tt .* (xtt.*cos(IncAng) + ytt.*sin(IncAng)) ...
                    +       3  *  uinc_t  .* (x3t.*cos(IncAng) + y3t.*sin(IncAng)) +     uinc    .* (x4t.*cos(IncAng) + y4t.*sin(IncAng)));
            end
            
        end
        
        
        function [duinc,duinc_t,duinc_tt] = dPlaneWave(Params,phi,IncAng,k)
            % calculate normal derivative of a planewave
            % in case of ellipse it returns derivative with respect to elliptical radii which is not exactly normal derivative
            
            IsStarshaped = false;
            if strcmpi(Params.ScattererType,'ellipse')
                dx = Params.FocalDistance * sinh(Params.eta) .* cos(phi);
                dy = Params.FocalDistance * cosh(Params.eta) .* sin(phi);
            elseif strcmpi(Params.ScattererType,'circle')
                dx = cos(phi);
                dy = sin(phi);
            elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
                IsStarshaped = true;
                [x,dx] = Params.Parameterization.XHandle.Derivatives(phi);
                [y,dy] = Params.Parameterization.YHandle.Derivatives(phi);
            end
            
            uinc = Tools.Common.IncidentWave.PlaneWave(Params,phi,IncAng,k)  ;
            
            duinc = 1i .* k .*  uinc .* (dx.*cos(IncAng) + dy.*sin(IncAng));
            %   h = FocalDist*sqrt(sinh(eta).^2 + sin(phi).^2);
            % duinc = duinc./h;
            
            if IsStarshaped
                
                h = sqrt(dx.^2 + dy.^2);
                duinc = 1i .* k .*  uinc .* (dy.*cos(IncAng) - dx.*sin(IncAng))./h;
                
                if nargout > 1
                    [x,xt,xtt,x3t,x4t] = Params.Parameterization.XHandle.Derivatives(phi);
                    [y,yt,ytt,y3t,y4t] = Params.Parameterization.YHandle.Derivatives(phi);
                    
                    [uinc,uinc_t,uinc_tt] = Uinc(Params,phi,IncAng,k);
                    
                    ht  = (xt.*xtt + yt.*ytt)./h;
                    htt = (xtt.^2 + ytt.^2 + xt.*x3t + yt.*y3t - ht.^2)./h;
                    h3t = (3*xtt.*x3t + 3*ytt.*y3t + xt.*x4t + yt.*y4t +  - 3*ht.*htt)./h;
                    
                    duinc_t = 1i .* k .* ( uinc_t  .* (yt .*cos(IncAng) - xt .*sin(IncAng))./h  + uinc   .* (ytt.*cos(IncAng) - xtt.*sin(IncAng))./h + uinc   .* (yt .*cos(IncAng) - xt .*sin(IncAng)).*(-ht./(h.^2)) );
                    duinc_tt = 1i .* k .*( uinc_tt .* (yt .*cos(IncAng) - xt .*sin(IncAng))./h  + uinc_t .* (ytt.*cos(IncAng) - xtt.*sin(IncAng))./h + uinc_t .* (yt .*cos(IncAng) - xt .*sin(IncAng)).*(-ht./(h.^2)) ...
                        +uinc_t  .* (ytt.*cos(IncAng) - xtt.*sin(IncAng))./h  + uinc   .* (y3t.*cos(IncAng) - x3t.*sin(IncAng))./h + uinc   .* (ytt.*cos(IncAng) - xtt.*sin(IncAng)).*(-ht./(h.^2)) ...
                        +uinc_t  .* (yt .*cos(IncAng) - xt .*sin(IncAng)).*(-ht./(h.^2))  + uinc .* (ytt.*cos(IncAng) - xtt.*sin(IncAng)).*(-ht./(h.^2))  + uinc .* (yt.*cos(IncAng) - xt.*sin(IncAng)).*((2*ht.^2)./(h.^3) - htt./(h.^2)) ...
                        );
                    
                    
                    
                    
                end
                
            end
            
        end
        
    end
end