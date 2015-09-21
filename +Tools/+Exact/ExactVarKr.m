classdef ExactVarKr < Tools.Exact.SuperExact

    properties%(Access = private)
        u       = 0;
        dudr    = 0; 
%         d2udr2    = 0;
    end
    
    
    properties(Access = private)
        R;
        Th;
    end
    
    methods(Static = true)
        
        function e = Exact(r,th,k)
            x = r.*cos(th);
            e = exp(1i.*k.*x);
        end
        
        function dr = drExact(r,th,k,kr)
            e = Tools.Exact.ExactVarKr.Exact(r,th,k);
            dr = 1i.*cos(th).*(k + kr.*r).*e;
        end
    end    
    
    methods
        
        function obj = ExactVarKr(Scatterer, WaveNumber)
            obj = obj@Tools.Exact.SuperExact(Scatterer, WaveNumber);
            
            try
                obj.R         = Scatterer.R;
                obj.Th        = Scatterer.Th;
            catch
                obj.R         = Scatterer.r;
                obj.Th        = Scatterer.th;
            end
            obj.u  = obj.Exact(obj.R,obj.Th,obj.Coeffs.k);
            obj.dudr = obj.drExact(obj.R,obj.Th,obj.Coeffs.k,obj.Coeffs.kr);                       
        end
        
        function [u,ur,urr,ut,utt,urt] = Derivatives(obj)%[u,ur,urr,u3r,ut,utt,urt,urtt] = Derivatives(obj)
            u      = obj.u;
            
            ur     = obj.dudr;
            
            [k,kr,krr] = obj.Coeffs.Derivatives('r');
            r=obj.R;
            th=obj.Th;
            
            urr = 1i.*cos(th).*(k + kr.*r).*ur + 1i.*cos(th).*(2*kr + krr.*r).*u;
           % u3r = ((k3r.*r+3.*krr).*1i.*cos(th) - 2.*(kr.*r+k).*(krr.*r + 2.*kr).*(cos(th).^2)).*u +...
           %     ((krr.*r + 2.*kr).*1i.*cos (th) - ((kr.*r + k).*cos (th)).^2).*ur;
            ut = -1i.*k.*r.*sin(th).*u;
            utt = -1i.*k.*r.*cos(th).*u - 1i.*k.*r.*sin(th).*ut;
            urt = -1i.*sin(th).*(u.*(k+r.*kr) + ur.*r.*k);
           % urtt = 1i.*(kr.*r + k).*(cos(th).*utt - 2.*sin(th).*ut - cos(th).*u);
        end
        
        function [s,sr,srr,s3r,st,stt,srt,srtt] = calc_s(obj)
            
            [k,kr,krr,k3r,k4r,k5r] = obj.Coeffs.Derivatives('r');
             r=obj.R;
            th=obj.Th;
            
            s   = (krr.*r + 3*kr).*1i.*cos (th) -((kr.* r).^2 + 2*k.*kr.*r).*(cos(th).^2);
            if nargout > 1
                sr  = (k3r.*r + 4*krr).*1i.*cos (th)-2*((kr.*(r.^2) + k.*r).*krr + 2*(kr.^2).*r + k.*kr).*(cos (th).^2);
            end
            if nargout > 2
                srr = (k4r.*r + 5*k3r).*1i.*cos (th)-2*((7*kr.*r + krr.*(r.^2) + 2*k).*krr + (kr.*(r.^2) + k.*r).*k3r +3*(kr.^2)).*(cos (th).^2);
            end
            if nargout > 3
                s3r = (k5r.*r + 6*k4r).*1i.*cos(th) -...
                    2*((9*krr.*r+15*kr+k3r.*(r.^2)).*krr + (10*kr.*r + 2*krr.*(r.^2)+3.*k).*k3r + (kr.*(r.^2) + k.*r).*k4r).*(cos(th).^2);
            end
            if nargout > 4
                st  = -(krr.*r + 3*kr).*1i.*sin (th) +   ((kr.*r).^2 + 2*k.*kr.*r).*sin (2*th);
            end
            if nargout > 5
                stt = -(krr.*r + 3*kr).*1i.*cos (th) + 2*((kr.*r).^2 + 2*k.*kr.*r).*cos (2*th);
            end
            if nargout > 6
                srt = -(k3r.*r + 4*krr).*1i.*sin (th)+2*((kr.*(r.^2) + k.*r).*krr + 2*(kr.^2).*r + k.*kr).*sin(2*th);
            end
            if nargout > 7
                srtt = -(k3r.*r + 4*krr).*1i.*cos(th)+4*((kr.*(r.^2) + k.*r).*krr + 2*(kr.^2).*r + k.*kr).*cos(2*th);
            end
        end
    end
    
end
