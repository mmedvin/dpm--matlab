classdef PolarMetrics < Tools.Metrics.AbstractMetrics
    methods(Static)
        function [h, hn, hnn,h3n,h4n,hf,hff,h3f,h4f] = metrics()%nothing else is used atm....
            h=1;
            hn=0; 
            hnn=0;
            h3n=0;
            h4n=0;
            hf=0; 
            hff=0;
            h3f=0;
            h4f=0;
        end
    end
end
