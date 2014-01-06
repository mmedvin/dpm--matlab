classdef SubmarineScatterer < Tools.Scatterer.SingleScatterer
    
    properties (Access = public)
         
        % grid
        R;
        Th;
         
        % grid on scatterer
        r;
        th;        
        dn;
        nrml_th;
        
        %scatterer params
        ellipse;
        tower;        
        
        Inside;
        Outside;
       % ScattererForSource;
        BasisArg;
        TheScatterer;
        
        Submarine;
        
        MetricsAtScatterer = Tools.Metrics.PolarMetrics;
    end
    
    properties(Access = private)
        ExpansionType;
    end
    
    methods
        
        function L = get.Inside(obj)
            L=obj.R <= obj.Submarine.Derivatives(obj.Th); %obj.MyShape(obj.Th);
        end
        function L = get.Outside(obj)
            L = obj.R > obj.Submarine.Derivatives(obj.Th); %obj.MyShape(obj.Th);
        end
        
        function BA = get.BasisArg(obj)
            BA = obj.nrml_th;
        end
        
        function TS = get.TheScatterer(obj)
            TS = 'TBD';%struct('R',obj.r0,'Th',obj.th);
           % error('TBD');
		end
		
		function obj = SubmarineScatterer(Grid,AddParams)
            obj = obj@Tools.Scatterer.SingleScatterer(Grid);
            obj.ellipse = AddParams.ellipse;
            obj.tower = AddParams.tower;
            
            obj.Submarine = Tools.Body.SubmarineBody(obj.ellipse,obj.tower);
            
            obj.ExpansionType = AddParams.ExpansionType;
            
            obj.MyGrid();
            obj.SplitGrid();%(obj.Eta,obj.Eta0);
            obj.GridGamma = intersect(obj.Np,obj.Nm)';
            obj.GridOnScatterer();
            
            assert(max(  [abs(obj.Grid.yn),abs(obj.Grid.y1)] ) >  max(abs(obj.Submarine.Derivatives(obj.nrml_th(:)).*sin(obj.nrml_th(:)) )) );
            assert(max([abs(obj.Grid.xn),abs(obj.Grid.x1)] ) > max(abs(obj.Submarine.Derivatives(obj.nrml_th(:)).*cos(obj.nrml_th(:)))));
        end
        
        function res = Expansion(obj,Xi0,Xi1,~,WaveNumber) % ~ was F
            
            % Xi1.Derivatives() are evaluated at obj.nrml_th rather than obj.th??
            [xi0,xi0t,xi0tt,xi0ttt,xi0tttt] = Xi0.Derivatives();
            [xi1,xi1t,xi1tt,xi1ttt,xi1tttt] = Xi1.Derivatives();
            [r0,dr,d2r,d3r,d4r]= obj.Submarine.Derivatives(obj.nrml_th); %obj.MyShape(obj.nrml_th); % r and its derivatives evaluated at the normal points                        
            %[k,kn,kf,knn,kff] = WaveNumber.Derivatives();
            k = WaveNumber.Derivatives();

            [h,ht,htt,h3t] = obj.Metrics(r0,dr,d2r,d3r,d4r);			
            [curv,ds_curv,dsds_curv] = obj.curvature(obj.nrml_th);
            
            xi0s    = xi0t./h;
            xi0ss   = (xi0tt./h - ht.*xi0t./h.^2)./h;
			%xi03s   = ((xi0ttt./h - 2*ht.*xi0tt./h.^2 - htt.*xi0t./h.^2 + 2*ht.^2.*xi0t./h.^3)./h ...
            %             - ht.*(xi0tt./h - ht.*xi0t./h.^2)./h.^2 )./h;
            %% bad/old xi03s   = (2*xi0t.*ht.^2./h.^3 - (2*xi0tt.*ht + xi0t.*htt)./h.^2 + xi0ttt./h)./h;
            xi0sss  = xi0ttt./(h.^3) - (htt.*xi0t + 3*ht.*xi0tt)./(h.^4) + 3*(ht.^2).*xi0t./(h.^5) ;
            
            xi0ssss = (xi0tttt./(h.^3) - 6*ht.*xi0ttt./(h.^4) ...
                    - (h3t.*xi0t + 4*htt.*xi0tt)./(h.^4) + (10*htt.*ht.*xi0t + 15*(ht.^2).*xi0tt)./(h.^5) ...
                    - 15*(ht.^3).*xi0t./(h.^6))./h;
            
            xi1s = xi1t./h;
            xi1ss= (xi1tt./h - ht.*xi1t./h.^2)./h;
            
            u_nn    = xi1.*curv - xi0ss - xi0.*k.^2;
            u_nnn   = u_nn.*curv +(curv.^2-k.^2).*xi1 - ds_curv.*xi0s - xi1ss - 2*curv.*xi0ss;
            
            %u_nns   = xi1s.*curv + xi1.*ds_curv -  xi0sss - xi0s.*k.^2;
            u_nnss  = xi1.*dsds_curv + 2*xi1s.*ds_curv + xi1ss.*curv - xi0ss.*k.^2 - xi0ssss;
            
            u_nnnn  = curv.*u_nnn + (2*curv.^2 - k.^2).*u_nn + 2*curv.^3.*xi1 ...
                    - 2*ds_curv.*xi1s - 4*curv.*xi1ss - u_nnss - 6*curv.*ds_curv.*xi0s - 6*curv.^2.*xi0ss  ;


            % Complete the extension by plugging normal derivatives into Taylor along  with the lengths of the normals:
            res = xi0 + obj.dn.*xi1 + ((obj.dn.^2)/2).*u_nn + obj.dn.^3./6.*u_nnn + obj.dn.^4./24.*u_nnnn;
        end
        
        
        
        
    end
    
    
    methods(Access = protected)
        
       
%         function zero = FindMyZeros(obj,theta,x1,y1)
%         % for given given (x1,y1) we want to find closest (x0,y0) on the
%         % scatterer i.e. the line trough these 2 point will be the normal
%         % such line should satisfy (y0-y1) = -1/(dr/dt) (x0-x1)
%         
%             [r_,dr_] = obj.MyShape(theta);
%                                                                         
%             x0 = r_.*cos(theta);
%             y0 = r_.*sin(theta);
%             
%             zero = (y0-y1).*dr_ - (x0-x1);
%             if abs(theta) < 10^10
%                 zero=0;
%             end
%             assert(zero<10)
%         
%         end 
        function dDist2dTheta = FindMyZeros(obj,theta,x1,y1)
            % for given given (x1,y1) we want to find closest (x0,y0) on the
            % scatterer i.e. the line trough these 2 point will be the normal
            % the distanse^2 is given by (x0-x1)^2+(y0-y1)^2 and we looking for
            % an extremum of this equation, i.e. we have to find the root\zero
            % of the first derivative of the distance equation....we also change
            % unknow (x0,y0) with it's parametrization on the curve....
            
            [r_,dr_] = obj.Submarine.Derivatives(theta); %obj.MyShape(theta);
                                                                        
            x0 = r_.*cos(theta);
            y0 = r_.*sin(theta);
            dx0dtheta = dr_.*cos(theta) - r_.*sin(theta);
            dy0dtheta = dr_.*sin(theta) + r_.*cos(theta);
            
            dDist2dTheta  = (x0 - x1).*dx0dtheta + (y0-y1).*dy0dtheta;                                   
        end
        
        function GridOnScatterer(obj)
            obj.r  = obj.R(obj.GridGamma);
            obj.th = obj.Th(obj.GridGamma);
            
            x1 = obj.r .* cos(obj.th);
            y1 = obj.r .* sin(obj.th);
            
            options = [];%optimset('TolX',10^-10);
            
            [obj.nrml_th,fval,exitflag,output] ...
             = arrayfun(@(indx) ...
                           fzero(@(arg) obj.FindMyZeros(arg,x1(indx),y1(indx)), obj.th(indx),options), ... % obj.th(indx) is an initial guess of root finding algorithm of 'fzero'
                           1:numel(obj.th));%,'UniformOutput',false);
            
            obj.nrml_th = obj.nrml_th.';
            [r0,dr0] = obj.Submarine.Derivatives(obj.nrml_th); %obj.MyShape(obj.nrml_th);% points of submarine normal to grid gamma
            
            if 0
                % a debug test
                x0=r0.*cos(obj.nrml_th);
                y0=r0.*sin(obj.nrml_th);
                %	zero = (y0-y1).*dr0 - (x0-x1);
                
                dx0   = dr0.*cos(obj.nrml_th) - r0.*sin(obj.nrml_th);
                dy0   = dr0.*sin(obj.nrml_th) + r0.*cos(obj.nrml_th);
            
                nrm = sqrt(dx0.^2 + dy0.^2);
            
                mdn = abs((x1+1i*y1)-(x0+1i*y0)).*sign(obj.r - r0);
                
                mx1 = x0 + mdn.*dy0./nrm;
                my1 = y0 - mdn.*dx0./nrm;
                norm(mx1-x1);
                norm(my1-y1);
                
                %another test
%                 mth0=pi/4;
%                 [mr0,mdr0] = obj.Submarine.Derivatives(mth0); %obj.MyShape(mth0);
%                 mx0=mr0.*cos(mth0);
%                 my0=mr0.*sin(mth0);
%                 
%                 delta=0.001;
%                 x1=mx0+delta;
%                 y1=my0 - delta/mdr0;
%                 
%                 nth0 = fzero(@(arg) obj.FindMyZeros(arg,x1,y1), angle(x1+1i*y1));
%                 [nr0,ndr0] = obj.Submarine.Derivatives(nth0); %obj.MyShape(nth0);
%                 
%                 nx0=nr0.*cos(nth0);
%                 ny0=nr0.*sin(nth0);
%                 
%                 nx0-mx0
%                 ny0-my0
            end
%             obj.dr = obj.r - r0;
                        
            %distance between (obj.r, obj.th) and (MyShape(obj.nrml_th),obj.nrml_th)
            obj.dn = sqrt((obj.r).^2+r0.^2-2*(obj.r).*r0.*cos(obj.th-obj.nrml_th)); % distance along normal to grid bdy nodes

            % sanity check 1:  check that no normal is longer than the diagonal of the grid
          %  assert(max(obj.dn) <= sqrt(obj.Grid.dx.^2 + obj.Grid.dy.^2));
            
            % sanity check 2: check that no normal is longer than curvature
%             c = obj.curvature(0:0.00001:2*pi);
%             if max(obj.dn) > min(abs(c))
%              warning('max(obj.dn)=%d <= min(c)=%d', max(obj.dn) ,min(abs(c)))
%             end
            
            
            %now we add direction to dn (delta n)
            sgn = sign(obj.r - r0);
            obj.dn = obj.dn.*sgn;
                        
        end
        
        function [curv,ds_curv,dsds_curv] = curvature(obj,theta)
            
            [r0,dr,d2r,d3r,d4r] = obj.Submarine.Derivatives(theta); %obj.MyShape(theta);
            h = sqrt(r0.^2+dr.^2);
            
            curv = - (r0.^2 + 2*dr.^2 -r0.*d2r)./(h.^3);
            ds_curv   =  (d3r.*r0.*(dr.^2 + r0.^2) + dr.*(3*d2r.*dr.^2 + (-3 *d2r.^2 + 4*dr.^2).*r0 - 3*d2r.*r0.^2 + r0.^3))./h.^6;

                                    dsds_curv = -((2*dr.*(d3r + dr) - d4r.*r0 + d2r.*(3*d2r + 2*r0))./(h.^5) ...
                                  - (4*dr.*(d2r + r0).*(-d3r.*r0 + dr.*(3*d2r + 2*r0)))./(h.^7) ...
                                  - (-3*curv.*d2r.*(d2r + r0))./(h.^4) ...
                                  - (-3*dr.*(d2r + r0).*ds_curv)./(h.^3) ...
                                  - 3*(-curv).*dr.*((d3r + dr)./(h.^4) - (3*dr.*(d2r + r0).^2)./(h.^6)));
            
            
% %             x0 = r_.*cos(theta);
% %             y0 = r_.*sin(theta);
%             dx0dtheta = dr_.*cos(theta) - r_.*sin(theta);
%             dy0dtheta = dr_.*sin(theta) + r_.*cos(theta);
%             
%             c=abs( ...
%                 dx0dtheta + 1i*dy0dtheta + dr_.*(cos(theta) + 1i.*sin(theta)) ...
%                 + drr_.*(sin(theta) - 1i.*cos(theta)) ...
%             );
                        
        end
        
        function [h,ht,htt,h3t] = Metrics(obj,r0,dr,d2r,d3r,d4r)
            h    = sqrt(r0.^2+dr.^2); %sqrt(dx0.^2 + dy0.^2);
            ht   = (dr.*d2r + r0.*dr)./h;
            htt   = (dr.^2 + d2r.*r0 + d2r.^2 + dr.*d3r - ht.^2)./h;
            h3t   = (3*dr.*d2r + d3r.*r0 + 3*d2r.*d3r + dr.*d4r - 3*ht.*htt)./h ;
        end
        
        
        function MyGrid(obj)
            Z = obj.Grid.Z();
            obj.R=abs(Z);
            obj.Th=angle(Z);
        end
        
    end
    
end
