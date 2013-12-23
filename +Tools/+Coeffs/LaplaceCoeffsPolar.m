classdef LaplaceCoeffsPolar < Tools.Coeffs.AbstractCoeffs
    % Implemntation of Laplacian coefficient for DPM
    % in this case sigma=0; a=b;
    properties
        B;
        
        a; ar; arr; a3r; a4r; a5r;
        b; br; brr; b3r; b4r; b5r;
        sigma; sigma_r; sigma_rr; sigma_3r; sigma_4r; sigma_5r;
    end
    
    methods
        
        function [d,dr,drr,d3r,d4r,d5r] = Derivatives(obj,kind)
            switch kind
                case 'a'
                    d   = obj.a;
                    dr  = obj.ar;
                    drr = obj.arr;
                    d3r = obj.a3r;
                    d4r = obj.a4r;
                    d5r = obj.a5r;
                case 'b'
                    d   = obj.b;
                    dr  = obj.br;
                    drr = obj.brr;
                    d3r = obj.b3r;
                    d4r = obj.b4r;
                    d5r = obj.b5r;
                case 'sigma'
                    d   = obj.sigma;
                    dr  = obj.sigma_r;
                    drr = obj.sigma_rr;
                    d3r = obj.sigma_3r;
                    d4r = obj.sigma_4r;
                    d5r = obj.sigma_5r;
                    
            end
        end
        
        function obj=LaplaceCoeffsPolar(Scatterer,Params)
            r  = Scatterer.r;
            r0 = Scatterer.r0;
            
            obj.B = Params.B;
            
            
            p   = r.^2 + 1;
            pr  = 2*r;
            prr = 2;
            
            if length(r)==1 && r>r0
                p   = obj.B;
                pr  = 0;
                prr = 0;
            else
                p(r>r0)     =obj.B;
                pr(r>r0)    =0;
                prr(r>r0)   =0;
            end
            
            obj.a   = p;
            obj.ar  = pr;
            obj.arr = prr;
            obj.a3r = 0;  obj.a4r = 0;  obj.a5r = 0;
            
            obj.b = obj.a; obj.br = obj.ar; obj.brr = obj.arr; obj.b3r = obj.a3r; obj.b4r = obj.a4r; obj.b5r = obj.a5r;
            obj.sigma=0; obj.sigma_r=0; obj.sigma_rr=0; obj.sigma_3r=0; obj.sigma_4r=0; obj.sigma_5r=0;
        end
    end
end
