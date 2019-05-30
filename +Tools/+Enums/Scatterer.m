classdef Scatterer < Tools.Enums.AbstractEnum
    enumeration
        Circle
        Ellipse
        StarShaped
    end
    
    
    
    methods
        function Print(obj,varargin)
            switch obj
                case Tools.Enums.Scatterer.Circle
                    R0=varargin{1};
                    fprintf('Radius=%-5.2f\n',R0)
                    
                case Tools.Enums.Scatterer.Ellipse
                    a=varargin{1};
                    b=varargin{2};
                    
                    FocalDistance = sqrt(a^2-b^2);
                    Eta0 = acosh(a/FocalDistance);% elliptical radii, e.i. a= FocalDistance*cosh(Eta0)
                    
                    fprintf('Params: FD=%-7.4f, ,Eta0=%-7.4fd, Axes: a=%-5.2f, b=%-5.2f, AR=%-5.2f ,\n',FocalDistance, Eta0 , a,b,a/b);
                    
                case Tools.Enums.Scatterer.StarShaped
                    Parameterization = varargin{5};
                    Parameterization.Print();
            end
        end
        
        function [Handle,Params,Extension,ExParams] = Helper(obj,varargin)
            
            switch(obj)
                case Tools.Enums.Scatterer.Circle
                    
                    R0=varargin{1};
                    %dummy = varargin{2};
                    HankelIndex = varargin{3};
                    HankelType = varargin{4};
                    
                    Handle  = @Tools.Scatterer.PolarScatterer;
                    Params  = struct('r0',R0,'ExpansionType',15,'Stencil',9);
                    Extension        = @Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension;%EBPolarHomoHelmholtz7OrderExtension;
                    ExParams  = struct('ScattererType',obj,'r',R0, 'HankelIndex', HankelIndex,'HankelType',HankelType);
                    
                case Tools.Enums.Scatterer.Ellipse
                    
                    a=varargin{1};
                    b=varargin{2};
                    HankelIndex = varargin{3};
                    HankelType = varargin{4};
                    
                    FocalDistance = sqrt(a^2-b^2);
                    Eta0 = acosh(a/FocalDistance);% elliptical radii, e.i. a= FocalDistance*cosh(Eta0)
                    
                    Handle  = @Tools.Scatterer.EllipticScatterer;
                    Params  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'Stencil',9);
                    Extension = @Tools.Extensions.TwoTupleExtension;
                    ExParams  = struct('ScattererType',obj, ...
                        'eta',Eta0, ...
                        'FocalDistance',FocalDistance, ...
                        'HankelIndex', HankelIndex, ...
                        'HankelType',  HankelType );
                    
                    
                case Tools.Enums.Scatterer.StarShaped
                    
                    Parameterization = varargin{5};
                    
                    Handle  = @Tools.Scatterer.StarShapedScatterer;
                    ExParams = struct('ScattererType',obj,'Parameterization',Parameterization);
                    Params  = ExParams;
                    Params.Stencil = 9;
                    Extension = @Tools.Extensions.TwoTupleExtension;
                    
            end
            
        end
        
        function [UincParams,Ang] = UincOnFieldHelper(obj,Prb,ExParams)
            UincParams = ExParams;
            
            switch(obj)
                case Tools.Enums.Scatterer.Circle
                    UincParams.r = Prb.Scatterer.R;
                    Ang = Prb.Scatterer.Th;
                    
                case Tools.Enums.Scatterer.Ellipse
                    UincParams.eta = Prb.Scatterer.Eta;
                    Ang = Prb.Scatterer.Phi;
                    
                case Tools.Enums.Scatterer.StarShaped 
                    UincParams = rmfield(UincParams,'Parameterization');
                    UincParams.r = Prb.Scatterer.R;
                    Ang = Prb.Scatterer.Th;

            end
            
        end
    end
end

