classdef SuperNavierStokesExact < handle
    
    methods(Static, Abstract)
        toString();
       
        Psi(theta,Params);
        
        DrDPsi(theta,Params);
        
        DPsiDThetaTheta(theta,Params);
        
        DPsiD4Theta(theta,Params);
        
        DPsiDrThetaTheta(theta,Params);
        
        O = Omega(theta,Params);
        
        DrDOmega(theta,Params);
        
        LO = LaplacianOmega(theta,Params);
        
    end
    
    methods
        function F = NSTimeSource(obj,theta,Params, RN,Ft,t,UseConvTerm)
            % implements the Source F = Omega_t - Laplacian Omega/RN + D(Omega,Psi)
            % of the undergoing Exact
            
            O  = obj.Omega(theta,Params);
            LO = obj.LaplacianOmega(theta,Params);
            
            [f,df] = Ft(t);
            F = O*df - f*LO/RN;
            
            if UseConvTerm
                warning('redo me!');
                % D= -64*xy.*(x.^2-y.^2).*((r.^2 - 3*(Params.r0^2)/4).^2 - (Params.r0^4)/16);
                % g=g+D.*ft.*ft;
            end
        end

        function F = NSSource(obj,theta,Params, Sigma)
            O  = obj.Omega(theta,Params);
            LO = obj.LaplacianOmega(theta,Params);
            F = LO - Sigma*O;
        end
        
        function L =  L_np1(obj,theta,Params, k,Fn,n)
            % returns Laplacian Omega - k Omega at time n
            O  = obj.Omega(theta,Params);
            LO = obj.LaplacianOmega(theta,Params);
            
            L = (LO - k*O)*Fn(n);
        end
        
        function L = L_n(obj,theta,Params, RN,ht,Ft,n,UseConvTerm)
            if nargin<8,UseConvTerm=false; end
            
            O  = obj.Omega(theta,Params);
            LO = obj.LaplacianOmega(theta,Params);
            
            Src = @(n) obj.NSTimeSource(theta,Params,RN,Ft,n*ht,UseConvTerm);
            
            k=2*RN/ht;
            Fn  = @(n) Ft(n*ht);
            
            L = -((LO + k*O)*Fn(n) + RN*(Src(n+1) + Src(n)) );
        end
      
        function L = DL_n(obj,gn,DO,theta,Params, RN,ht,Ft,n,UseConvTerm)
            %Discrete L_n
            if nargin<10,UseConvTerm=false; end
             
            Src = @(n) obj.NSTimeSource(theta,Params,RN,Ft,n*ht,UseConvTerm);
            
            k=2*RN/ht;
            
            L = -((gn + 2*k*DO) + RN*(Src(n+1) + Src(n)) );
        end
        
        
        function Res = Test(obj,theta,Params, RN,ht,Ft,n,UseConvTerm)
            if nargin<8,UseConvTerm=false; end

            k=2*RN/ht;
            
            Fn  = @(n) Ft(n*ht);
                       
            LNp1 = obj.L_np1(theta,Params, k,Fn,n+1); 
            LN   = obj.L_n(theta,Params, RN,ht,Ft,n,UseConvTerm);
            
            Res = LNp1 - LN;
            
        end
        
        function Res = Test2(obj,theta,Params, RN,ht,Ft,n,UseConvTerm)
            if nargin<8,UseConvTerm=false; end

            k=2*RN/ht;
            
            Fn  = @(n) Ft(n*ht);
                       
            LNp1 = obj.L_np1(theta,Params, k,Fn,n+1); 
            
            DO  = obj.Omega(theta,Params)*Fn(n);
            gn = obj.L_np1(theta,Params, k,Fn,n);
            DLN   = obj.DL_n(gn,DO,theta,Params, RN,ht,Ft,n,UseConvTerm);
            
            Res = LNp1 - DLN;
            
        end
    end
    
end
