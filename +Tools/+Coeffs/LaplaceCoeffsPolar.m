classdef LaplaceCoeffsPolar < Tools.Coeffs.AbstractCoeffs
	% Implemntation of Laplacian coefficient for DPM
	% in this case sigma=0; a=b;
	properties
		B;
		
		a; ar; arr; a3r; a4r; a5r;
		b; bt; btt; b3t; b4t; b5t;
		br;btr;
		sigma; sigma_r; sigma_rr; sigma_3r; sigma_4r; sigma_5r;
	end
	
	methods
		
		function [d,dr,drr,d3r,d4r,d5r] = Derivatives(obj,WhichOne)
			switch WhichOne
				case 'a'
					d   = obj.a;
					dr  = obj.ar;
					drr = obj.arr;
					d3r = obj.a3r;
					d4r = obj.a4r;
					d5r = obj.a5r;
				case 'b'
					d   = obj.b;
					dr  = obj.bt;
					drr = obj.btt;
					d3r = obj.b3t;
					d4r = obj.b4t;
					d5r = obj.b5t;
				case 'sigma'
					d   = obj.sigma;
					dr  = obj.sigma_r;
					drr = obj.sigma_rr;
					d3r = obj.sigma_3r;
					d4r = obj.sigma_4r;
					d5r = obj.sigma_5r;
				otherwise
					error('unknown coeff')
			end
		end
		
		function obj=LaplaceCoeffsPolar(Scatterer,Params)
			r  = Scatterer.r;
			
			if size(r)==1
				r = ones(size(Scatterer.th)).*r;
			end
			
			p   = r.^2 + 1;
			pr  = 2*r;
			prr = 2;
			
			obj.a   = p;
			obj.ar  = pr;
			obj.arr = prr;
			obj.a3r = 0;  obj.a4r = 0;  obj.a5r = 0;
			
			obj.b	= p; 
			obj.bt	= 0; obj.btt = 0; obj.b3t = 0; obj.b4t = 0; obj.b5t = 0;
			obj.br = pr;
			obj.btr = 0;
			
			obj.sigma=0; obj.sigma_r=0; obj.sigma_rr=0; obj.sigma_3r=0; obj.sigma_4r=0; obj.sigma_5r=0;
		end
	end
end
