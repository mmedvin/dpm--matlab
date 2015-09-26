classdef FourierBasis < Tools.Basis.BasisFunctionWD
    
    methods(Static)
		function FBasis = BasisHelper(f,dfdr,err,~)%range)

            if ~exist('err','var'),err = 10^(-10);end
			
            if err < 1
				[cn0,cn1,M] = Tools.Basis.FourierBasis.FourierCoeffNew(f,dfdr,err);
			else
				M=err;
				
				M1 = 2*M+1;
				th=linspace(0,2*pi,M1+1);
				th=th(1:M1);
				fth   = f(th);
				dfth  = dfdr(th);
				
				%calc coeff
				cn0 = Tools.Basis.FourierBasis.FftCoefs(fth,M1).';				
				cn1 = Tools.Basis.FourierBasis.FftCoefs(dfth,M1).';
				
			end
			me = metaclass(Tools.Basis.FourierBasis);
			FBasis = struct('type','Fourier','Handle',str2func(me.Name),...
							'Indices',  -M:M, ...
							'cn0',cn0,'cn1',cn1,'M',M,'MoreParams',[],'NBss', 2*M+1);
		end
        
		function Coefs = FftCoefs(vec, siz)
			if ~exist('siz','var'),siz = numel(vec);end
			
			Coefs = fft(vec);
            Coefs = fftshift(Coefs);
            Coefs = Coefs/siz;  %normalization        
		end
		
        function [cn0,cn1,M] = FourierCoeffNew(f,dfdr,err)
            M0=1020;
            M1 = 2*M0+1;
            %th=linspace(-pi,pi,M1+1);
	     th=linspace(0,2*pi,M1+1);
            th=th(1:M1);
            
            fth   = f(th);
            dfth  = dfdr(th);
            
            %calc coeff
			cn0 = Tools.Basis.FourierBasis.FftCoefs(fth,M1);
            
            cn1 = Tools.Basis.FourierBasis.FftCoefs(dfth,M1);           
            
            Ml = find(abs(cn0) > err,1,'last');
            Mf = find(abs(cn0) > err,1,'first');
            Mm1 = min(Mf,M1-Ml+1);
            if isempty(Mm1),Mm1 = M1;end

            Ml = find(abs(cn1) > err,1,'last');
            Mf = find(abs(cn1) > err,1,'first');
            Mm2 = min(Mf,M1-Ml+1);
            if isempty(Mm2),Mm2 = M1;end
            
            Mm = min(Mm1,Mm2);
            cn0=cn0(Mm:M1-Mm+1).';
            cn1=cn1(Mm:M1-Mm+1).';
            
            M=(numel(cn0)-1)/2;                        
        end
        
        function [cn0,cn1,M] = FourierCoeff(f,dfdr)
            cn0 =   quadgk(f,-pi,pi);
            cn1  =   quadgk(dfdr,-pi,pi);
            j=0;
            notstop=true;
            err=2*pi*(1e-10);
            while notstop
                j=j+1;
                g = @(x) f(x).*exp(-1i*j*x);
                tcn0p =   quadgk(g,-pi,pi);
                
                g = @(x) f(x).*exp(-1i*(-j)*x);
                tcn0m =   quadgk(g,-pi,pi);
                
                dg = @(x) dfdr(x).*exp(-1i*j*x);
                tcn1p  =   quadgk(dg,-pi,pi);
                
                dg = @(x) dfdr(x).*exp(-1i*(-j)*x);
                tcn1m  =   quadgk(dg,-pi,pi);
                
                
                cn0 = [tcn0m; cn0 ;tcn0p];
                cn1 = [tcn1m; cn1 ;tcn1p];
                
                % notstop  = abs(cn0(1))>err;
                notstop  = abs(cn0(1))>err | abs(cn0(2))>err | ...
                    abs(cn0(end))>err | abs(cn0(end-1))>err;
                
            end
            cn0=cn0/(2*pi);
            cn1=cn1/(2*pi);
            M=j;
            
        end
    end
    
    
    methods
        
        function obj    = FourierBasis(x,n,NoMetrics,NoMoreParams)
            if nargin == 0
                obj.xi0     = 0;
                obj.xi0t    = 0;
                obj.xi0tt   = 0;
                obj.xi0ttt  = 0;
                obj.xi0tttt = 0;
                obj.xi0tttttt = 0;
            else
                obj.xi0         = exp(1i.*n.*x);
                obj.xi0t        = 1i.*n.*obj.xi0;
                obj.xi0tt       =-(n^2).*obj.xi0;
                obj.xi0ttt      =-1i*(n^3).*obj.xi0;
                obj.xi0tttt     = (n^4).*obj.xi0;
                obj.xi0tttttt   =-(n^6).*obj.xi0;
            end
        end
        
    end
end
