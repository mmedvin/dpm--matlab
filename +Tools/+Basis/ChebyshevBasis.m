classdef ChebyshevBasis < Tools.Basis.BasisFunctionWD
    properties(Constant)
        epsilon = 2*0.00390625; %magic number, it is 2pi/512 translated to the [-1,1]
        use_epsilon = 1; % 0 or 1
    end
    methods(Static)
        function CBasis = BasisHelper(f,dfdr,range)
            [cn0,cn1,N] = Tools.Basis.ChebyshevBasis.ChebychevCoeff(f,dfdr,range);
            me = metaclass(Tools.Basis.ChebyshevBasis);
            CBasis = struct('type','Chebyshev','Handle',str2func(me.Name),...
                            'Indices',  0:N-1, ...
                            'cn0',cn0.','cn1',cn1.','M',N,'AddParams',range,'NBss', N);
            
        end
        function [cn0,cn1,N] = ChebychevCoeff(f,dfdr,range)
            %N=1024;
            % x = cos(pi / N * (N-0.5 : -1 : 0.5));
            
            cn0 = Tools.Basis.ChebyshevBasis.chebyshevExpansion(f,range);
            cn1 = Tools.Basis.ChebyshevBasis.chebyshevExpansion(dfdr,range);
            
            err=(1e-10);
            cn1(abs(cn1) < err) = 0;
            N = find(cn1,1,'last');
            cn0 = cn0(1:N).';
            cn1 = cn1(1:N).';
            
        end
        
        function coeffs=chebyshevExpansion(func,range)
            
            % func is a function of theta
            % range is the range of theta for func taken as a 1x2 vector [a,b]
            % N is the number of coefficients returned
            
            N=2^10;
            x= Tools.Basis.ChebyshevBasis.chebroot(N); % roots of a large chebyshev polynomial
            
            %x = cos(pi / N * (N-0.5 : -1 : 0.5));
            
            % Linear change of variables: x = (2*theta - a - b) / (b - a);
            % Inverse change of variables:
            eps_= Tools.Basis.ChebyshevBasis.epsilon* Tools.Basis.ChebyshevBasis.use_epsilon;
            
            theta = ((range.b - range.a)*x/(1-eps_) + range.a + range.b)/2;
%             assert(all(theta >= range.a));
%             assert(all(theta <= range.b));
            
            vector=func(theta); % sample the function at those roots
            
            coeffs=dct(vector(:)); % Matlab's built-in "discrete cosine transform"
            coeffs=[coeffs(1)/sqrt(N); coeffs(2:N)*sqrt(2/N)]; % Normalizing
            %coeffs=coeffs(1:N); % Takes first N coefficients
            
            % coeffs=coeffs*sqrt(2/N);
        end
        
        function x = chebroot(n)
            %CHEBROOT Roots of Chebychev polynomial of the first kind.
            %
            %   X = CHEBROOT(N) returns the roots of the Chebychev polynomial of the first
            %   kind of degree N.
            %
            %   Because the extreme values of Chebychev polynomials of the first kind are
            %   either -1 or 1, their roots are often used as starting values for the nodes
            %   in mimimax approximations.
            %
            %   See also CHEBPOLY, CHEBEXTR.
            
            %   Author:      Peter John Acklam
            %   Time-stamp:  2004-09-22 09:21:37 +0200
            %   E-mail:      pjacklam@online.no
            %   URL:         http://home.online.no/~pjacklam
            
            %error(nargchk(1, 1, nargin));
            narginchk(1, 1);
            
            if ~isa(n, 'double') || any(size(n) ~= 1) || ~isreal(n)
                error('N must be double, scalar, and real.');
            end
            if (n < 1) | (n ~= round(n))
                error('N must be a positive integer');
            end
            
            x = cos(pi / n * (n-0.5 : -1 : 0.5));
            
            % adjust for the fact that cos(pi / 2) is not exactly zero
            if rem(n, 2)
                x((n+1) / 2) = 0;
            end
            
            % this is faster if `n' is very large
            % x = zeros(1, n);
            % if rem(n, 2)
            %    t = cos(pi / n * (0.5 : n/2-1));
            %    x( 1 :      (n-1)/2 ) = -t;
            %    x( n : -1 : (n+3)/2 ) =  t;
            % else
            %    t = cos(pi / n * (0.5 : n/2-0.5));
            %    x( 1 :      n/2   ) = -t;
            %    x( n : -1 : n/2+1 ) =  t;
            % end
        end
        
    end
    
    methods
        
        function obj    = ChebyshevBasis(theta,n,Range)
            % theta should be in [a,b], so that x below is in [-1,1]
            if nargin == 0
                obj.xi0     = 0;
                obj.xi0t    = 0;
                obj.xi0tt   = 0;
                obj.xi0ttt  = 0;
                obj.xi0tttt = 0;
                obj.xi0tttttt = 0;
            else
                eps_ = obj.epsilon * obj.use_epsilon;
                x = (1-eps_)*(theta - Range.a) / (Range.b-Range.a) + (-1+eps_)*(theta - Range.b)/(Range.a-Range.b);
                assert(all(theta >= Range.a));
                assert(all(theta <= Range.b));
                assert(all(abs(x)<=1));
                
                obj.xi0=cos(n*acos(x));
                txi0t =(n.*sin(n.*acos(x)))./sqrt(1 - x.^2);                
                txi0tt = (-(n.^2).* obj.xi0 + x.* txi0t)./(1 - x.^2);
                txi0ttt = ((1 - n.^2).* txi0t + 3* x.*txi0tt)./(1 - x.^2);                               
                txi0tttt = ((4 - n.^2).* txi0tt + 5* x.* txi0ttt)./(1 - x.^2);
                txi0ttttt = ((9 - n.^2).* txi0ttt + 7* x.* txi0tttt)./(1 - x.^2);
                txi0tttttt = ((16 - n.^2).* txi0tttt + 9* x.* txi0ttttt)./(1 - x.^2);
                
                if eps_ == 0                    
                    indx = find ((x~=0) & abs(x-sign(x))<obj.epsilon);
                    
                    txi0t(indx) =         DerOfCheNearEnds(obj,x(indx),n,1);
                    txi0tt(indx) =        DerOfCheNearEnds(obj,x(indx),n,2);
                    txi0ttt(indx)   = DerOfCheNearEnds(obj,x(indx),n,3);
                    txi0tttt(indx) =      DerOfCheNearEnds(obj,x(indx),n,4);
                    txi0tttttt(indx) =    DerOfCheNearEnds(obj,x(indx),n,6);
                end
                
%                 txi0t(x==1) = n^2;
%                 txi0t(x==-1) = -((-1)^n) * (n^2);
%                 txi0tt(x==1) = (1/3) * (n^2) *(-1 + n^2);
%                 txi0tt(x==-1) = (1/3) * ((-1)^n) * (n^2) * (-1 + n^2);                
%                 txi0ttt(x==1) =1/15 * (n^2)* (4 - 5* n^2 + n^4);
%                 txi0ttt(x==-1) = -(1/15)* ((-1)^n) *(n^2)* (4 - 5 *n^2 + n^4);               
%                 txi0tttt(x==1) = (1/105) * (n^2) *(-36 + 49 * n^2 - 14* (n^4) + n^6);
%                 txi0tttt(x==-1) =(1/105) * ((-1)^n) * (n^2) *(-36 + 49 * (n^2)
%                 - 14 *(n^4) + n^6);                
% u_xxxxx(1) =1/945 n^2 (576 - 820 n^2 + 273 n^4 - 30 n^6 + n^8)
% 
% u_xxxxx(-1) =-(1/945) (-1)^n n^2 (576 - 820 n^2 + 273 n^4 - 30 n^6 + n^8)
%                 
                dxdtheta = 2*(1-eps_) / (Range.b-Range.a);               
                
                obj.xi0t = txi0t*dxdtheta;
                obj.xi0tt = txi0tt*(dxdtheta^2);
                obj.xi0ttt = txi0ttt*(dxdtheta^3);
                obj.xi0tttt = txi0tttt*(dxdtheta^4);
                obj.xi0tttttt = txi0tttttt*(dxdtheta^6);
                
            end
        end
        
        function dChe=DerOfCheAtEnds(obj,x,n,p)
            %calculate p'th derivatie of n'th order chebyshev function at 1 or -1
            assert(all(abs(x)==1));
            
            
            
            k=0:(p-1);
            t = (n^2 - k.^2)./(2*k+1);
            dChe = prod(t) .* ((x).^(n+p));            
        end
        
        function dChe=DerOfCheNearEnds(obj,x,n,p)
            
            h=x-sign(x);
            dChe = obj.DerOfCheAtEnds(sign(x),n,p) ...
                 + obj.DerOfCheAtEnds(sign(x),n,p+1).*h ...
                 + obj.DerOfCheAtEnds(sign(x),n,p+2).*(h.^2) ...
                 + obj.DerOfCheAtEnds(sign(x),n,p+3).*(h.^3) ...
                 + obj.DerOfCheAtEnds(sign(x),n,p+4).*(h.^4);
        end
        
        function [xi0,xi0f,xi0ff,xi0fff,xi0ffff,xi0ffffff] = Derivatives(obj)
            xi0     = obj.xi0;
            xi0f    = obj.xi0t;
            xi0ff   = obj.xi0tt;
            xi0fff  = obj.xi0ttt;
            xi0ffff = obj.xi0tttt;
            xi0ffffff = obj.xi0tttttt;
        end
        
    end
end
