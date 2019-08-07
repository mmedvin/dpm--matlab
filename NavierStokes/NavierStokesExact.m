classdef NavierStokesExact < Tools.Enums.AbstractEnum
    enumeration
        Exact1
        Exact2
        Exact1DiscreteSrc
        Exact2DiscreteSrc
        Exact1Time
        Exact2Time
    end
    
    methods
        
        function str = toString(obj,varargin)
            Exact = obj.Helper();
            ST = obj.Subtype();
            if ST == Tools.Enums.Category.Time
                str = ['Exact Psi: F(x,y) = ' Exact.toString() ',  Ft=' func2str(varargin{1})];
            else
                str = ['Exact: F(x,y) = ' Exact.toString() ', using ' ST.toString() ' source'];
            end
        end
        
        function ST = Subtype(obj)
            switch(obj)
                case NavierStokesExact.Exact1
                    ST = Tools.Enums.Category.Exact;
                case NavierStokesExact.Exact2
                    ST = Tools.Enums.Category.Exact;
                case NavierStokesExact.Exact1DiscreteSrc
                    ST = Tools.Enums.Category.Discrete;
                case NavierStokesExact.Exact2DiscreteSrc
                    ST = Tools.Enums.Category.Discrete;
                case NavierStokesExact.Exact1Time
                    ST = Tools.Enums.Category.Time;
                case NavierStokesExact.Exact2Time
                    ST = Tools.Enums.Category.Time;
            end
        end
        
        function [Exact,SourceHandle] = Helper(obj,varargin)%f,fn,ChebyshevRange)
            
             switch(obj)
                case NavierStokesExact.Exact1
                    Exact = Tools.Exact.NSExact1;
                    SourceHandle = @Tools.Source.NavierStokesSourceRTh;
                case NavierStokesExact.Exact2
                    Exact = Tools.Exact.NSExact2;
                    SourceHandle = @Tools.Source.NavierStokesSourceRTh2;
                case {NavierStokesExact.Exact1DiscreteSrc, NavierStokesExact.Exact1Time}
                    Exact = Tools.Exact.NSExact1;
                    SourceHandle = @Tools.Source.NavierStokesDescreteSrc;
                case {NavierStokesExact.Exact2DiscreteSrc, NavierStokesExact.Exact2Time}
                    Exact = Tools.Exact.NSExact2;
                    SourceHandle = @Tools.Source.NavierStokesDescreteSrc;
            end
            
        end
    end
end

 