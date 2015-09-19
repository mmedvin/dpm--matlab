classdef SuperExtension < handle
    
    properties
        W;
        Wf;
    end
    
    properties(Access = protected)        
        NoSource = Tools.Source.SuperSource();%SuperHelmholtzSource();  %construct empty source
    end
    
    methods(Access = public,Abstract=true)
        Expand(obj);
        ExpandSource(obj,Source);
        ExpandedBasis(obj,n);
    end
end