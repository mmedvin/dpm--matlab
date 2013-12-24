classdef WaveNumberStarShaped < Tools.WaveNumber.WaveNumberPolarR
	%should be a varK wavenumber for starshaped body
		% unfinished for future use

    properties
        kn;kf;knn;kff;
    end
    
    methods(Static = true)
        function [k,kn] = kkn(ellipse,tower,theta,k0,r0)
            
         %   [r,rn] = Tools.WaveNumber.WaveNumberElliptical.chngcoord(ellipse,tower,theta);
            SmB=SubmarineBody(ellipse,tower);
            [r,rn] = SmB.Derivatives(theta);
            
            [k,kr] = Tools.WaveNumber.WaveNumberPolarR.kkr(r,r0,k0);
            kn = kr.*rn ;
        end
        
        function [r,rn,rf,rnn,rff] = chngcoord(theta)
            
            [r,dr] = obj.Submarine.Derivatives(theta); %obj.MyShape(theta);
            
            [dtds,dttds,d3tds,d4tds] = obj.DtDs(r0,dr,d2r,d3r,d4r);
            [curv,ds_curv,dsds_curv] = obj.curvature(obj.nrml_th);
                
            x = r.*cos(theta);
            y = r.*sin(theta);
            xt = dr_.*cos(theta) - r_.*sin(theta);
            yt = dr_.*sin(theta) + r_.*cos(theta);
                
            xs    = xt./dtds;
            ys    = yt./dtds;
            
            
            
            
                rn =(xn.*x + yn.*y)./r;
                rs = (xs.*x + ys.*y)./r;
                
                rnn = (cosh(2*eta).*FocDist.^2 - rn.^2)./r;
                rff = (-cos(2*phi).*FocDist.^2  - rs.^2)./r;
            end
    end
    
    methods
        %function [k,kn,kf,knn,kff, k3n,k3f,k4n,k4f,knf,knff,knnf,knnff] = Derivatives(obj)
        function [k,kn,kf,knn,kff] = Derivatives(obj)
            k       = obj.k;
            kn      = obj.kn;
            kf      = obj.kf;
            knn     = obj.knn;
            kff     = obj.kff;
        end
        
        function obj = WaveNumberStarShaped(Scatterer,AddParams)            
            k0 = AddParams.k0;
            r0 = AddParams.r0;
                                   
            [r,rn,rf,rnn,rff] = ...
                Tools.WaveNumber.WaveNumberElliptical.chngcoord(Scatterer.FocalDistance,Scatterer.Eta,Scatterer.Phi);
            
            PolarScatterer.r  = r;
            PolarScatterer.r0 = r0;
            
            obj=obj@Tools.WaveNumber.WaveNumberPolarR(PolarScatterer,k0);%r,r0);
            
            obj.kn = obj.kr.*rn ;
            obj.kf = obj.kr.*rf ;
            
            obj.knn = obj.krr.*rn.^2 + obj.kr.*rnn ;
            obj.kff = obj.krr.*rf.^2 + obj.kr.*rff ;                
        end        
    end   
end
