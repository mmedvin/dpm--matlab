classdef ExactCoeffs < Tools.Coeffs.AbstractCoeffs
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Params;
    end
    
    methods 
        function [d,dx,dxx,d3x] = Derivatives(obj,WhichOne,x,y)
            d=0;dx=0;dxx=0;d3x=0;
            
            switch WhichOne
                case 'c'
                    d   = (y + obj.Params.a).^obj.Params.b;
                    dx  = obj.Params.b*(y + obj.Params.a).^(obj.Params.b-1);
                    dxx = (obj.Params.b-1)*obj.Params.b*(y + obj.Params.a).^(obj.Params.b-2);
                    d3x = (obj.Params.b-2)*(obj.Params.b-1)*obj.Params.b*(y + obj.Params.a).^(obj.Params.b-3);
                case 'd'
                    d   = (x + obj.Params.c).^obj.Params.d;
                    dx  = obj.Params.d*(x + obj.Params.c).^(obj.Params.d-1);
                    dxx = (obj.Params.d-1)*obj.Params.d*(x + obj.Params.c).^(obj.Params.d-2);
                    d3x = (obj.Params.d-2)*(obj.Params.d-1)*obj.Params.d*(x + obj.Params.c).^(obj.Params.d-3);
            end
        end
        

        function obj = ExactCoeffs(Params)
            obj.Params = Params;
        end
    end
    
end

