classdef SubmarineBody < Tools.Body.SuperStarShapedBody
    properties(Access = protected)
        ellipse;
        tower;
        LastTheta=NaN;
        LastNargOut = 0 ;
    end
    methods        
        function obj = SubmarineBody(ellipse,tower)
            obj.ellipse = ellipse;
            obj.tower = tower;                        
        end
        
        function [R, DR, DRR, D3R, D4R] = Derivatives(obj,theta)

            NotTheSameTheta = true; % Guilty until proven otherwise
            szTh     = size(theta); 
            szLTh   =  size(obj.LastTheta);
            if szTh(:) == szLTh(:)
                NotTheSameTheta = any(obj.LastTheta(:) ~= theta(:));
            end
                
            if (obj.LastNargOut < nargout) || NotTheSameTheta % used to reduce computation effort
                obj.LastNargOut = nargout;
                obj.LastTheta   = theta;
                
                e = (cos(theta).^2/obj.ellipse.a^2 + sin(theta).^2/obj.ellipse.b^2);
                g = 1 + obj.tower.c*sin(theta).^obj.tower.p;
                obj.r   = sqrt(g./e);
                
                if nargout > 1
                    de      = (-1/obj.ellipse.a^2 + 1/obj.ellipse.b^2).* sin(2*theta);
                    dg      = obj.tower.c*obj.tower.p*cos(theta).*sin(theta).^(-1 + obj.tower.p);
                    obj.dr  = (dg - de.*obj.r.^2)/2./obj.r./e;
                end
                
                if nargout > 2
                    dee = 2*(-1/obj.ellipse.a^2 + 1/obj.ellipse.b^2).*cos(2*theta);
                    %d3e = -4*de;
                    %d4e = -4*dee;
                    dgg = 0.5*obj.tower.c*obj.tower.p.*(obj.tower.p+obj.tower.p*cos(2*theta) -2).*sin(theta).^(-2 + obj.tower.p);
                    obj.drr = (dgg - dee.*obj.r.^2 - 4*obj.r.*obj.dr.*de - 2*e.*obj.dr.^2)/2./obj.r./e;
                end
                
                if nargout > 3
                    d3g     = 0.5*obj.tower.c*obj.tower.p.*cos(theta).*(4-6*obj.tower.p + obj.tower.p^2 ...
                            + obj.tower.p^2*cos(2*theta) -2).*sin(theta).^(-3 + obj.tower.p);
                    obj.d3r = (d3g + 4*de.*obj.r.^2 - 6*obj.r.*obj.dr.*dee - 6*de.*obj.dr.^2 ...
                            - 6*obj.r.*de.*obj.drr - 6*e.*obj.dr.*obj.drr)/2./obj.r./e;
                end
                
                if nargout > 4
                    d4g = (1/8)*obj.tower.c*obj.tower.p.*(...
                        -32+56*obj.tower.p -24*obj.tower.p^2 + 3*obj.tower.p^3 ...
                        + 4*(-4 + 8*obj.tower.p - 6*obj.tower.p^2 + obj.tower.p^3).*cos(2*theta) ...
                        + (obj.tower.p^3).*cos(4*theta)).*sin(theta).^(-4 + obj.tower.p);
                    obj.d4r = (d4g + 4*dee.*obj.r.^2 - 12*dee.*obj.dr.^2 - 12*obj.r.*dee.*obj.drr + 32*obj.r.*obj.dr.*de ...
                            - 24*de.*obj.dr.*obj.drr - 6*e.*obj.drr.^2 - 8*e.*obj.dr.*obj.d3r - 8*obj.r.*de.*obj.d3r)/2./obj.r./e;
                end
                
%                 if nargout > 5
%                     d5g = (1/8)*obj.tower.c*obj.tower.p.*cos(theta).*(...
%                        -160 - 320*obj.tower.p + 200*obj.tower.p^2 - 40*obj.tower.p^3 + 3*obj.tower.p^4 ...
%                        + 4*(8 - 20*obj.tower.p + 20*obj.tower.p^2 - 10*obj.tower.p^3 + obj.tower.p^4).*cos(2*theta) ...
%                        + (obj.tower.p^4).*cos(4*theta)).*sin(theta).^(-5 + obj.tower.p);
%                     
%                     obj.d5r    = (d5g + 40*obj.r.*obj.dr.*dee - 16*de.*obj.r.^2 + 80*de.*obj.dr.^2 + 80*obj.r.*de.*obj.drr ...
%                                - 60*obj.dr.*dee.*obj.drr - 30*de.*obj.drr.^2 - 20*e.*obj.drr.*obj.d3r - 40*de.*obj.dr.*obj.d3r ...
%                                - 20*obj.r.*dee.*obj.d3r - 10*obj.r.*de.*obj.d4r - 10*e.*obj.dr.*obj.d4r)/2./obj.r./e;
%                 end
            end
            
            %[R, DR, DRR, D3R, D4R, D5R] = Derivatives@Tools.Body.SuperStarShapedBody(obj);
			
			R   = obj.r;
			if nargout>1, DR  = obj.dr; end
			if nargout>2, DRR  = obj.drr; end
			if nargout>3, D3R  = obj.d3r; end
            if nargout>4, D4R  = obj.d4r; end
			%if nargout>5, D5R  = obj.d5r; end
            
        end

        
    end
end
