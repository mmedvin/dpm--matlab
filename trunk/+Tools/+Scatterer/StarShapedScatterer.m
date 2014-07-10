classdef StarShapedScatterer < Tools.Scatterer.SingleScatterer
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
   properties
       
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
            
            obj.SplitGrid();
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
            TS = 'TBD';%struct('R',obj.r0,'Th',obj.th);
            % error('TBD');
        end
    
        function res = Expansion(obj,Xi0,Xi1,~,WaveNumber) % ~ was F
            
            % Xi1.Derivatives() are evaluated at obj.nrml_t rather than obj.th??
            [xi0,xi0t,xi0tt,xi0ttt,xi0tttt] = Xi0.Derivatives();
            [xi1,xi1t,xi1tt,xi1ttt,xi1tttt] = Xi1.Derivatives();            
            
            k = WaveNumber.Derivatives(); %[k,kn,kf,knn,kff] = WaveNumber.Derivatives();

            [h,ht,htt,h3t] = obj.MetricsAtScatterer.metrics(obj.nrml_t);			
            [curv,ds_curv,dsds_curv] = obj.MetricsAtScatterer.curvature(obj.nrml_t);
            
            xi0s    = xi0t./h;
            xi0ss   = (xi0tt./h - ht.*xi0t./h.^2)./h;
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
        
        function HashInfo = hash_func(obj)
            HashInfo.Loaded = 0;
            HashInfo.X = obj.XHandle;
            HashInfo.Y = obj.YHandle;
            HashInfo.Xmeta = metaclass(HashInfo.X);
            HashInfo.Ymeta = metaclass(HashInfo.Y);
            
            
            HashInfo.HashVal = obj.Grid.x1*obj.Grid.xn*obj.Grid.Nx + obj.Grid.y1*obj.Grid.yn*obj.Grid.Ny;
            for i = 1:numel(HashInfo.Xmeta.Properties)
                
                nm1 = HashInfo.Xmeta.Properties{i}.Name;
                val1 = eval(['obj.XHandle.',nm1]);
                if isobject(val1), continue,end
                
                for j = 1:numel(HashInfo.Ymeta.Properties)
                    
                    nm2 = HashInfo.Ymeta.Properties{j}.Name;
                    val2 = eval(['obj.YHandle.',nm2]);
                    
                    if isobject(val2), continue,end
                    
                   HashInfo.HashVal = HashInfo.HashVal + (val1.^i) * (val2.^j);
                end
            end
                        
            HashInfo.HashVal = fix(HashInfo.HashVal *1e4);
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
           
           [obj.ScatParamToGridAng,fval,exitflag,output] = ...
               arrayfun(@(indx) ...
                                fzero(@(arg) obj.FindMyZeros(arg,obj.Th(indx)), obj.Th(indx),options), ... % obj.Th(indx) is an initial guess of root finding algorithm of 'fzero'
                        1:numel(obj.Th));%,'UniformOutput',false);
                                                  
           obj.ScatParamToGridAng = reshape(obj.ScatParamToGridAng,obj.Grid.Nx,obj.Grid.Ny);                    
          
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
                    [obj.nrml_t,fval,exitflag,output] ...
                        = arrayfun(@(indx) ...
                        fzero(@(arg) obj.DerivativeOfDistance(arg,x1(indx),y1(indx)),InitialGuess(indx),options), ... 
                        1:numel(obj.th));%,'UniformOutput',false);
                    
                    obj.nrml_t = obj.nrml_t.';
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

