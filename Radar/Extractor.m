classdef Extractor

    properties(SetAccess = public)
        %M;
        
        
        Basis;
        
        %parameterization of the scatterer boundary
        XHandle;
        YHandle;
        
        %scattererd data on circles defined by r and theta
        r;
        theta;
        scattData;
        
        Scatterer;
        spec;
              
        %ellipse axes
        a;
        b;
    end
    
    methods
        function obj=Extractor(data)
            
            if isa(data,'struct')
                M = data;
            elseif isa(data,'char')
                load(data,'M','a','b');
            end
            
            %obj.M = M;
            
            obj.a = a;
            obj.b = b;
            
            obj.r = M.r;
            obj.theta = M.theta;
            
            obj.spec = M.spec;
            obj.scattData   = M.scattData;
            
            obj.Scatterer = M.Scatterer;
            
            obj.XHandle     = M.Scatterer.Parameterization.XHandle;
            obj.YHandle     = M.Scatterer.Parameterization.YHandle;
            
            obj.Basis       = M.Basis;            
            
            if ~isfield(obj.Basis,'AddParams')
                obj.Basis.AddParams=[];
            end
            
        end
        
        function curve = Curve(obj,theta)
            [x,xt] = obj.XHandle.Derivatives(theta);
            [y,yt] = obj.YHandle.Derivatives(theta);
            curve.z ={x;y} ;
            
            n=abs(yt+1i*xt);
            curve.nz = {yt./n;-xt./n};
        end

        function [scattField,k,phi] = FieldOnCircles(obj,im,il,ind)
            %return the value from an upper circle!!!

            k   = obj.scattData{im,il}.k0;
            phi = obj.scattData{im,il}.phi;
            
            scattField.value        = obj.scattData{im,il}.u{ind};
            scattField.normal_deriv = obj.scattData{im,il}.un{ind};
        end
        
        function curve = Circle(obj,ind)
            curve.z  = obj.r(ind) * [cos(obj.theta); sin(obj.theta)];
            curve.nz =              [cos(obj.theta); sin(obj.theta)];
        end
        
        
        function [scattField,k,phi] = FieldOnScatterer(obj,im,il,theta)
            
            k   = obj.scattData{im,il}.k0;
            phi = obj.scattData{im,il}.phi;
            
            cn0 = obj.scattData{im,il}.cn0;
            cn1 = obj.scattData{im,il}.cn1;
            
            u=0;
            un=0;
            
            for j=1:numel(obj.Basis.Indices0)
                uj = obj.Basis.Handle(theta,obj.Basis.Indices0(j),obj.Basis.AddParams);
                u = u + cn0(j).*uj.xi0;
            end
            
            for j=1:numel(obj.Basis.Indices1)
                uj = obj.Basis.Handle(theta,obj.Basis.Indices1(j),obj.Basis.AddParams);
                un = un + cn1(j).*uj.xi0;
            end
            
            
            uinc = Tools.Common.IncidentWave.PlaneWave(obj.Scatterer,theta,obj.spec.phi_range(im),k);
            
            Duinc = Tools.Common.IncidentWave.dPlaneWave(obj.Scatterer,theta,obj.spec.phi_range(im),k);
            
            scattField.value        = u - uinc;
            scattField.normal_deriv = un - Duinc;
        end

        
    end
end
