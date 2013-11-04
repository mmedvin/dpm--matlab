classdef SupperScatterer < handle
            
    properties(Abstract = true,Access = public)
        Grid;
        
        Mp;
        Mm;
        Np;
        Nm;
        
        GridGamma;
        
        Size;
        
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
   
    
    
    methods(Abstract = true, Access = public)
        %  Grid = AnotherGrid(obj,Grid);
        %res =
        Expansion(obj,Xi0,Xi1,F,WaveNumber);
        
    end     
end