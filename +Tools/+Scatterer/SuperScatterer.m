classdef SuperScatterer < handle
            
    properties(Abstract = true,SetAccess=protected)%,Access = public)
        Grid;
        
        Mp;
        Mm;
        Np;
        Nm;
        
        GridGamma;               
        
        In; %can't find better name yet...
        Out;
        
        %%%%%%%%%%%%%
        BasisArg;
        TheScatterer;  
        MetricsAtScatterer;
        %ScattererForSource;
        
        Inside;
        Outside;
    end
   
    properties
        Size;
    end
    
    methods(Abstract = true, Access = public)
        %  Grid = AnotherGrid(obj,Grid);
        %res =
        Expansion(obj,Xi0,Xi1,F,WaveNumber);
        
    end   
    
    methods
        function sz = get.Size(obj)%What is this????
            
            sz = obj.Grid.Size;
            
        end
    end
end