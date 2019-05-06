classdef Basis < Tools.Enums.AbstractEnum
    enumeration
        Chebyshev
        Fourier
    end
    
    methods
        function B = Helper(obj,varargin)%f,fn,ChebyshevRange)
            
            switch(obj)
                case Tools.Enums.Basis.Chebyshev
                    B = Tools.Basis.ChebyshevBasis.BasisHelper(varargin{:});%(f,fn,ChebyshevRange);
                case Tools.Enums.Basis.Fourier
                    B = Tools.Basis.FourierBasis.BasisHelper(varargin{:});%(f,fn);%,[10,10]);%,[15,15]);%[1e-14,1e-14]);%20);%
            end
            
        end
    end
end

