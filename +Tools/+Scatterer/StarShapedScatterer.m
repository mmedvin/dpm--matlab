classdef StarShapedScatterer < Tools.Scatterer.SingleScatterer
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
   properties(SetAccess = protected)%(AbortSet = true) commented because this make problem with save/load
       
       MetricsAtScatterer;
       
       BasisArg;
       TheScatterer;    % need for non homogenious problem
       
       Inside;
       Outside;
       
 
        
        
     % following should be a protected-properties, but Matlab's debug won't allow to see it, so...grrrr   
        % give grid point (x0,y0), drop a normal to the scatterer, the point
        % on scatterer is defined by the parameter t = nrml_t, i.e. (x(t),y(t))

        a; 
        b; 
        
        nrml_t; 
        dn;      % the distance between (x0,y0) and (x(t),y(t))
        
        XHandle;
        YHandle;
        
        FZeroAlg=1;
        
        R;
        Th;
        r;
        th;
        
        ScatParamToGridAng;
        RAtScatParamToGridAng;
        
    end
    

    
    methods
        
        function obj = StarShapedScatterer(Grid,Params)
            
            obj = obj@Tools.Scatterer.SingleScatterer(Grid);

            obj.XHandle = Params.Parameterization.XHandle;
            obj.YHandle = Params.Parameterization.YHandle;

            obj.MyGrid();
            
            HashInfo = obj.TryToLoad();
            
            if HashInfo.Loaded
                
                obj.ScatParamToGridAng      = HashInfo.ScatParamToGridAng;
                obj.RAtScatParamToGridAng   = HashInfo.RAtScatParamToGridAng;
                
            else
                obj.CreateScatParamToGridAng();                
            end
            
            obj.SplitGrid(Params.Stencil);
            obj.GridGamma = intersect(obj.Np,obj.Nm)';
            
            obj.GridOnScatterer(HashInfo);            
            
            if ~HashInfo.Loaded
                obj.Save(HashInfo,obj.ScatParamToGridAng, obj.RAtScatParamToGridAng, obj.nrml_t);
            end
            
            obj.MetricsAtScatterer = Tools.Metrics.StarShapedParametricMetric(obj.XHandle, obj.YHandle);
            
        end
        
        function I = get.Inside(obj)
            I = obj.R <= obj.RAtScatParamToGridAng;
        end
        
        function O = get.Outside(obj)
            O  = obj.R > obj.RAtScatParamToGridAng;
        end

        function BA = get.BasisArg(obj)
            BA = obj.nrml_t;
        end
        
        function TS = get.TheScatterer(obj)
            
            x = obj.XHandle.Derivatives(obj.nrml_t);
            y = obj.YHandle.Derivatives(obj.nrml_t);
            
            %x = obj.XHandle.Derivatives(obj.ScatParamToGridAng);
            %y = obj.YHandle.Derivatives(obj.ScatParamToGridAng);
            
            z=x+1i*y;
            TS = struct('r',abs(z),'th',angle(z));
            %TS = struct('r',obj.RAtScatParamToGridAng(obj.GridGamma));%,'Th',obj.th);
            
        end
    
        function [curv,ds_curv,dsds_curv] = Curvature(obj,angl)
            [curv,ds_curv,dsds_curv] = obj.MetricsAtScatterer.curvature(angl);
        end
        
        function [h,ht,htt,h3t] = Metrics(obj,angl)
            [h,ht,htt,h3t] = obj.MetricsAtScatterer.metrics(angl);
        end
        function res = Expansion(obj,Xi0,Xi1,Source,WaveNumber) 
            
            % Xi1.Derivatives() are evaluated at obj.nrml_t rather than obj.th??
            [xi0,xi0t,xi0tt,xi0ttt,xi0tttt] = Xi0.Derivatives();
            [xi1,xi1t,xi1tt,xi1ttt,xi1tttt] = Xi1.Derivatives();            
            [h,ht,htt,h3t] = obj.Metrics(obj.nrml_t); %obj.MetricsAtScatterer.metrics(obj.nrml_t);			
            [curv,ds_curv,dsds_curv] = obj.Curvature(obj.nrml_t); %obj.MetricsAtScatterer.curvature(obj.nrml_t);
            
            
            [x,dx] = obj.XHandle.Derivatives(obj.nrml_t);
            [y,dy] = obj.YHandle.Derivatives(obj.nrml_t);
            r_ = sqrt(x.^2 + y.^2);
            
            [k,kn,knn,ks,kss] = WaveNumber.Derivatives(obj);
            %%%%%%%%%%%%%%%%%%%%%%
%             [k,kr,krr] = WaveNumber.Derivatives();
%                         
%             kx = kr.*x./r_;%kr.*cos(obj.nrml_t); % + kt*(-sin(obj.nrml_t)/r;
%             ky = kr.*y./r_;%kr.*sin(obj.nrml_t); % + kt*(cos(obj.nrml_t)/r;
%             
%             kxx = (kr.*y.^2)./(r_.^3) + (krr.*x.^2)./(r_.^2);
%             kyy = (kr.*x.^2)./(r_.^3) + (krr.*y.^2)./(r_.^2);
%             kxy = -(kr.*x.*y)./(r_.^3) + (krr.*x.*y)./(r_.^2);
% 
%             %%%%%%%%%%%%%%%%%%%%%%                        
%             %%%%%%%%%%%%%%%%%%%%%%
%             kn  = (kx.*dy  - ky.*dx)./h;
%             knn = (kxx.*dy.^2 - 2*kxy.*dx.*dy + kyy.*dx.^2 )./(h.^2); ...     + kx.*dyy.*dy + ky.*dxx.*dx          - kn.* ht./(h.^2);
% 
%             ks  = (kx.*dx + ky.*dy)./h;
%             kss = (kxx.*dx.^2 + 2*kxy.*dx.*dy + kyy.*dy.^2)./(h.^2) + curv.*kn;

            [F,Fn,Fnn,Fs,Fss] = Source.Derivatives(obj);
            %%%%%%%%%%%%%%%%%%%%%%
%             [F,Fr,Frr,Ft,Ftt,Frt] = Source.Derivatives();
%             
%             Fx = Fr.*x./r_ - Ft.*y./(r_.^2);
%             Fy = Fr.*y./r_ + Ft.*x./(r_.^2);
%             
%             Fxx = Ftt.*y.^2./(r_.^4) + Frr.*(x.^2)./(r_.^2) - 2*Frt.*x.*y./(r_.^3) ...
%                 + 2*Ft.*x.*y./(r_.^4) + Fr.*(y.^2)./(r_.^3) ;
%             
%             Fyy = Ftt.*x.^2./(r_.^4) + Frr.*(y.^2)./(r_.^2) + 2*Frt.*x.*y./(r_.^3) ...
%                 - 2*Ft.*x.*y./(r_.^4) + Fr.*(x.^2)./(r_.^3) ;
%             
%             Fxy =-Ftt.*x.*y./(r_.^4) + Frr.*x.*y./(r_.^2) + Frt.*(x.^2 - y.^2)./(r_.^3) ...
%                 + Ft.*(y.^2-x.^2)./(r_.^4) - (Fr.*x.*y)./(r_.^3) ;
% 
% 
%             Fn  = (Fx.*dy  - Fy.*dx)./h;
%             Fnn = (Fxx.*dy.^2 - 2*Fxy.*dx.*dy + Fyy.*dx.^2 )./(h.^2); ...     + kx.*dyy.*dy + ky.*dxx.*dx          - kn.* ht./(h.^2);
% 
%             Fs  = (Fx.*dx + Fy.*dy)./h;
%             Fss = (Fxx.*dx.^2 + 2*Fxy.*dx.*dy + Fyy.*dy.^2)./(h.^2) + curv.*Fn;
            
            %%%%%%%%%%%%%%%%%%%%%%
            
            
            xi0s    = xi0t./h;
            xi0ss   = (xi0tt./h - ht.*xi0t./h.^2)./h;
            xi0sss  = xi0ttt./(h.^3) - (htt.*xi0t + 3*ht.*xi0tt)./(h.^4) + 3*(ht.^2).*xi0t./(h.^5) ;
            
            xi0ssss = (xi0tttt./(h.^3) - 6*ht.*xi0ttt./(h.^4) ...
                    - (h3t.*xi0t + 4*htt.*xi0tt)./(h.^4) + (10*htt.*ht.*xi0t + 15*(ht.^2).*xi0tt)./(h.^5) ...
                    - 15*(ht.^3).*xi0t./(h.^6))./h;
            
            xi1s = xi1t./h;
            xi1ss= (xi1tt./h - ht.*xi1t./h.^2)./h;
            
            u_nn    = F + xi1.*curv - xi0ss - xi0.*k.^2;
            u_nnn   = Fn + u_nn.*curv +(curv.^2-k.^2).*xi1 - ds_curv.*xi0s - xi1ss - 2*curv.*xi0ss - 2*k.*kn.*xi0;%- 2*k.*kn.*xi0
            
            %u_nns   = xi1s.*curv + xi1.*ds_curv -  xi0sss - xi0s.*k.^2 - 2*k.*ks.*xi0; %- 2*k.*ks.*xi0
            u_nnss  = Fss + xi1.*dsds_curv + 2*xi1s.*ds_curv + xi1ss.*curv - xi0ss.*k.^2 - xi0ssss - 4*k.*ks.*xi0s -2*(ks.^2 + k.*kss).*xi0; % - 4*k.*ks.*xi0s -2*(ks.^2 + k.*kss).*xi0
            
            u_nnnn  = Fnn + curv.*u_nnn + (2*curv.^2 - k.^2).*u_nn + (2*curv.^3 - 4*k.*kn).*xi1 ...
                    - 2*ds_curv.*xi1s - 4*curv.*xi1ss - u_nnss - 6*curv.*ds_curv.*xi0s - 6*curv.^2.*xi0ss -2*(kn.^2 + k.*knn).*xi0 ; %- 4*k.*kn.*xi1 -2*(kn.^2 + k.*knn).*xi0


            % Complete the extension by plugging normal derivatives into Taylor along  with the lengths of the normals:
            res = xi0 + obj.dn.*xi1 + ((obj.dn.^2)/2).*u_nn + obj.dn.^3./6.*u_nnn + obj.dn.^4./24.*u_nnnn;
        end
        
    end
    
    methods(Access = protected)
        
        function g = FindMyZeros(obj,t,theta)

            x0 = obj.XHandle.Derivatives(t);
            y0 = obj.YHandle.Derivatives(t);
            
            arg = angle(x0+1i*y0);
            if sign(arg) < sign(t) && abs(theta) > pi/2
                arg = arg + 2*pi;
            elseif sign(arg) > sign(t) && abs(theta) > pi/2
                arg = arg - 2*pi;
            end
            
            g = arg - theta;                
        end
        
        function val = ExtractAllPropertiesVals(obj,A)
            Ameta = metaclass(A);
            val=[];
            
            if numel(Ameta.Properties)==0
                val = 1;
            end
            
            for i = 1:numel(Ameta.Properties)
                nm1 = Ameta.Properties{i}.Name;
                val1 = eval(['A.',nm1]);
                if isobject(val1)
                    val = [val obj.ExtractAllPropertiesVals(val1)]; %#ok<AGROW>
                else
                    val = [val val1]; %#ok<AGROW>
                end
            end
        end
        
        function HashInfo = hash_func(obj)
			P=11;
			Ppow = P;
            HashInfo.Loaded = 0;
            HashInfo.X = obj.XHandle;
            HashInfo.Y = obj.YHandle;
            HashInfo.Xmeta = metaclass(HashInfo.X);
            HashInfo.Ymeta = metaclass(HashInfo.Y);
            
			HashInfo.HashVal = abs(obj.Grid.x1*obj.Grid.xn*obj.Grid.Nx) + abs(obj.Grid.y1*obj.Grid.yn*obj.Grid.Ny);
			
			Str=strsplit(HashInfo.Xmeta.Name,'.');
			Str=Str{end};
			for i=1:numel(Str)
				 HashInfo.HashVal = int64(HashInfo.HashVal + double(Str(i)) *Ppow);
				 Ppow=Ppow*P;
			end
			
			Str=strsplit(HashInfo.Ymeta.Name,'.');
			Str=Str{end};
			Ppow = P;
			for i=1:numel(Str)
				 HashInfo.HashVal = int64(HashInfo.HashVal + double(Str(i))*Ppow);
				Ppow=Ppow*P;
            end
            
           Xvals= obj.ExtractAllPropertiesVals(HashInfo.X);
           Yvals= obj.ExtractAllPropertiesVals(HashInfo.Y);
           vals = conv(Xvals,Yvals);
           Ppow = P;
           
           for i = 1:numel(vals)
   				 HashInfo.HashVal = int64(HashInfo.HashVal + vals(i)*Ppow);
				Ppow=Ppow*P;
           end
                        
            
        end
        
        function Gstr = TypeOfGrid(obj)
            m=metaclass(obj.Grid);
            if isempty(strfind(m.Name, 'Polar')), 
                Gstr = 'Cartesian'; 
            else
                Gstr = 'Polar'; 
            end
            %isempty(strfind(m.Name, 'Cartesian'))
        end
        
        function HashInfo = TryToLoad(obj)
            
            HashInfo = obj.hash_func();
            
            path = [ pwd filesep 'SavedData' filesep 'StarShapedBody' filesep obj.TypeOfGrid filesep ,...
                'h' num2str(HashInfo.HashVal),...
                'x' num2str(obj.Grid.x1) '_' num2str(obj.Grid.xn) 'y' num2str(obj.Grid.y1) '_' num2str(obj.Grid.yn), ...
                'g' num2str(obj.Grid.Nx) 'x' num2str(obj.Grid.Ny) '.mat'];
             
            if exist(path,'file')
                S = load(path);
                
                HashInfo.ScatParamToGridAng = S.HashInfo.ScatParamToGridAng;
                HashInfo.RAtScatParamToGridAng = S.HashInfo.RAtScatParamToGridAng;
                HashInfo.nrml_t = S.HashInfo.nrml_t;
                HashInfo.Loaded = true;
                
                %HashInfo.ScatParamToGridAng;
                %HashInfo.RAtScatParamToGridAng;
            end
        end
        
        function Save(obj,HashInfo,ScatParamToGridAng, RAtScatParamToGridAng,nrml_t)            
            HashInfo.ScatParamToGridAng     = ScatParamToGridAng;
            HashInfo.RAtScatParamToGridAng  = RAtScatParamToGridAng;
            HashInfo.nrml_t = nrml_t;
            
            path = [pwd filesep 'SavedData' filesep 'StarShapedBody' filesep obj.TypeOfGrid];
            
            if ~exist(path,'dir')
                mkdir(path);
            end
            
            save([ path  filesep 'h' num2str(HashInfo.HashVal),...
                'x' num2str(obj.Grid.x1) '_' num2str(obj.Grid.xn) 'y' num2str(obj.Grid.y1) '_' num2str(obj.Grid.yn), ...
                 'g' num2str(obj.Grid.Nx) 'x' num2str(obj.Grid.Ny) '.mat'] ,'HashInfo');
        end
        
        function CreateScatParamToGridAng(obj)
            
            options = [];%optimset('TolX',10^-10);
            
            InitialGuess = obj.Th;
            GridPointToTest = obj.Th;
            tmp = zeros(size(obj.Th));
            
            parfor  indx = 1:numel(obj.Th)
                [tmp(indx)...
                    ... ,fval,exitflag,output
                    ]...
                    =   fzero(@(arg) FindMyZeros(obj,arg,InitialGuess(indx)), GridPointToTest(indx),options);
            end
                       
            obj.ScatParamToGridAng = reshape(tmp,obj.Grid.Nx,obj.Grid.Ny);
            obj.RAtScatParamToGridAng = sqrt(obj.XHandle.Derivatives(obj.ScatParamToGridAng).^2 + obj.YHandle.Derivatives(obj.ScatParamToGridAng).^2);

        end
        
        function dDist2dTheta = DerivativeOfDistance(obj,theta,x1,y1)
            % for given (x1,y1) we want to find closest (x0,y0) on the
            % scatterer i.e. the line trough these 2 point will be the normal
            % the distanse^2 is given by (x0-x1)^2+(y0-y1)^2 and we looking for
            % an extremum of this equation, i.e. we have to find the root\zero
            % of the first derivative of the distance equation....we also change
            % unknow (x0,y0) with it's parametrization on the curve....
            
            [x0,dx0dtheta] = obj.XHandle.Derivatives(theta);
            [y0,dy0dtheta] = obj.YHandle.Derivatives(theta);
            
            dDist2dTheta  = (x0 - x1).*dx0dtheta + (y0-y1).*dy0dtheta;
        end
        
        function GridOnScatterer(obj,HashInfo)
            
            obj.r  = obj.R(obj.GridGamma);
            InitialGuess = obj.ScatParamToGridAng(obj.GridGamma); %  obj.th
			obj.th = obj.Th(obj.GridGamma);
            
            x1 = obj.r.*cos(obj.th);
            y1 = obj.r.*sin(obj.th);
            
            options = [];%optimset('TolX',10^-10);
            
            switch obj.FZeroAlg
                case 1 %regular choice
                    if HashInfo.Loaded
                        obj.nrml_t = HashInfo.nrml_t;
                    else
                        if 0
                            tmp = zeros(size(obj.th));
                            parfor indx = 1:numel(obj.th)
                            [tmp(indx),fval,exitflag,output] = fzero(@(arg) DerivativeOfDistance(obj,arg,x1(indx),y1(indx)),InitialGuess(indx),options);
                            end
                            obj.nrml_t = tmp;
                        else
                            [obj.nrml_t,fval,exitflag,output] ...
                                = arrayfun(@(indx) ...
                                fzero(@(arg) obj.DerivativeOfDistance(arg,x1(indx),y1(indx)),InitialGuess(indx),options), ...
                                1:numel(obj.th));%,'UniformOutput',false);
                            obj.nrml_t = obj.nrml_t.';
                        end
                    end
                case 2 %TBD
                case 3
                    if ~strcmpi(obj.Shape,'circle')
                        error('wrong use of debugging parameter');
                    end
                    obj.nrml_t = obj.th;
                    
            end
                         
            if 0
                % a debug test
                x0=r0.*cos(obj.nrml_t);
                y0=r0.*sin(obj.nrml_t);
                %	zero = (y0-y1).*dr0 - (x0-x1);
                
                dx0   = dr0.*cos(obj.nrml_t) - r0.*sin(obj.nrml_t);
                dy0   = dr0.*sin(obj.nrml_t) + r0.*cos(obj.nrml_t);
                
                nrm = sqrt(dx0.^2 + dy0.^2);
                
                mdn = abs((x1+1i*y1)-(x0+1i*y0)).*sign(obj.r - r0);
                
                mx1 = x0 + mdn.*dy0./nrm;
                my1 = y0 - mdn.*dx0./nrm;
                norm(mx1-x1);
                norm(my1-y1);
                
                
            end
            %             obj.dr = obj.r - r0;
            
            %distance between (obj.r, obj.th) and (MyShape(obj.nrml_t),obj.nrml_t)
            
            x0 = obj.XHandle.Derivatives(obj.nrml_t);
            y0 = obj.YHandle.Derivatives(obj.nrml_t);
            
            obj.dn = sqrt((x0-x1).^2+(y0-y1).^2); % distance along normal to grid bdy nodes
            
            % sanity check 1:  check that no normal is longer than the diagonal of the grid
            %  assert(max(obj.dn) <= sqrt(obj.Grid.dx.^2 + obj.Grid.dy.^2));
            
            % sanity check 2: check that no normal is longer than curvature
            %             c = obj.curvature(0:0.00001:2*pi);
            %             if max(obj.dn) > min(abs(c))
            %              warning('max(obj.dn)=%d <= min(c)=%d', max(obj.dn) ,min(abs(c)))
            %             end
            
            %eta0 = acosh(x0./obj.FocalDistance./cos(obj.nrml_t));
            %obj.Eta0 =  acosh(x0./obj.FocalDistance./cos(obj.nrml_t));
            %obj.phi =  angle(x0+1i*y0);
            
            %now we set the direction of dn (delta n)
            %sgn = sign(obj.eta - eta0);
            %sgn = sign(obj.eta - obj.Eta0);
            
            
            sgn = ones(size(obj.GridGamma));
            sgn(obj.Inside(obj.GridGamma)) = -1;
            obj.dn = obj.dn.*sgn;
            
        end
        
        function MyGrid(obj)
            obj.Th = angle(obj.Grid.Z); 
            obj.R  = abs(obj.Grid.Z);
        end
    end
    
    
end

