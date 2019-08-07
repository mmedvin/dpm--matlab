classdef PolarScatterer < Tools.Scatterer.SingleScatterer
    properties(SetAccess=protected)
        R;
        Th;
        r;
        th;
        
        dr;
        
        r0;
        
        Inside;
        Outside;
%         ScattererForSource;
        BasisArg;
        TheScatterer;
            
        MetricsAtScatterer = Tools.Metrics.PolarMetrics;
        
        GridLinesIntersection;
    end
    
    methods
        
        function L = get.Inside(obj)
            L=obj.R <= obj.r0;
        end
        function L = get.Outside(obj)
            L = obj.R>obj.r0;
        end
        
        function BA = get.BasisArg(obj)
            BA = obj.th;
        end
        
        function TS = get.TheScatterer(obj)
            TS = struct('r',obj.r0,'th',obj.th);%*ones(size(obj.th))
        end
                     
        function obj = PolarScatterer(Grid,Params)
            obj = obj@Tools.Scatterer.SingleScatterer(Grid);
            obj.r0 = Params.r0;
            
            %obj.ExpansionType = Params.ExpansionType;            
            
            obj.MyGrid();
			
            obj.SplitGrid(Params.Stencil);
            	
            obj.GridGamma = intersect(obj.Np,obj.Nm)';
            
            obj.GridOnScatterer();
            
            obj.IntersectionWithGridLines();
            
        end     
                 
          function res = Expansion(obj,Xi0,Xi1,F,Coeffs)
            %   dummy, need to implement abstract inheritage, leftovers from older version....
          end
    end
    
    
    methods(Access = protected)
        
        function GridOnScatterer(obj)
            obj.r  = obj.R(obj.GridGamma);
            obj.dr = obj.r - obj.r0;
            obj.th = obj.Th(obj.GridGamma);
        end
        
        function MyGrid(obj)
            Z = obj.Grid.Z();
            obj.R=abs(Z);
            obj.Th=angle(Z);
        end	
        
        function IntersectionWithGridLines(obj)

            y = obj.r0^2 - obj.Grid.x.^2;
            msky = y>=0;
            
            tmp = sqrt(y(msky));
            obj.GridLinesIntersection.Y.y = [-tmp,tmp];
            obj.GridLinesIntersection.Y.x = [obj.Grid.x(msky),obj.Grid.x(msky)];
            obj.GridLinesIntersection.Y.msk = msky;
            tmp = find(msky);
            obj.GridLinesIntersection.Y.indx = [tmp,tmp];

            obj.GridLinesIntersection.Y.th = atan2(obj.GridLinesIntersection.Y.y, obj.GridLinesIntersection.Y.x);
            
            
            x = obj.r0^2 - obj.Grid.y.^2;
            mskx = x>=0;
            tmp = sqrt(x(mskx));
            obj.GridLinesIntersection.X.x = [-tmp,tmp];
            obj.GridLinesIntersection.X.y = [obj.Grid.y(mskx),obj.Grid.y(mskx)];
            obj.GridLinesIntersection.X.msk = mskx;
            tmp = find(mskx);
            obj.GridLinesIntersection.X.indx = [tmp,tmp];
            obj.GridLinesIntersection.X.th = atan2(obj.GridLinesIntersection.X.y,obj.GridLinesIntersection.X.x);
        end
        
    end
end
