classdef NavierStokesDescreteSrc < Tools.Source.SuperSource
    properties %(Dependent = true)
        Source;
        
    end
    
    properties%(Access = protected)       
        Scatterer; 
        Params;
        OnNp;
        OnMp;
        OnGG;
        Src;
        fn;
        n;
    end
    
    methods
        function [F,Fr,Frr,Ft,Ftt,Frt] = Derivatives(obj,Ex)
            
            if obj.IsDummy
                F=0;
                %                 Fr=0;
                %                 Frr=0;
                %                 Ft=0;
                %                 Ftt=0;
                %
                %                 Ftt=0;
            else
                
%                 try
%                     r = obj.Scatterer.R;
%                     F = obj.OnNp;
%                 catch exception
%                     if strcmp(exception.identifier,'MATLAB:nonExistentField')

                        %F = obj.OnGG;
                       
                        if 1 %interpolation
                            x = obj.Scatterer.r.*cos(obj.Scatterer.th);
                            y = obj.Scatterer.r.*sin(obj.Scatterer.th);
                            
                            F = griddata(obj.Params.Grid.X(obj.Params.GridGamma),obj.Params.Grid.Y(obj.Params.GridGamma),obj.OnGG,x,y,'cubic');
                        else
                            r  = obj.Scatterer.r;
                            th = obj.Scatterer.th;
                            r0 = obj.Params.r0;
                            ht = obj.Params.ht;
                            
                            [Fn,dFn] = obj.Params.Fn(obj.Params.n);
                            [Fnp1,dFnp1] = obj.Params.Fn(obj.Params.n+1);
                            
                            F = (r.^2).*sin(2*th).*(192*Fnp1 - 4*(4*r.^2-3*r0^2).*(2*Fn/ht + dFn + dFnp1));
                        end
                        
                        %                     else
%                         rethrow(exception);
%                     end
%                end
                
                %F=F - obj.Src.Derivatives()*(obj.fn(obj.n)+obj.fn(obj.n-1))  ;%???
            end
            
            
        end
            
        
            
%         function obj = NavierStokesDescreteSrc(Scatterer,CoeffsClsrHndl,CoeffsParams,Params)
%             CoeffsParams.sigma=1;
%             obj.Src = Params.Handle(Scatterer, CoeffsClsrHndl,CoeffsParams,Params.SourceParams);
         function obj = NavierStokesDescreteSrc(Scatterer,~,~,Params)            
            obj.Scatterer = Scatterer;
            obj.OnNp = Params.Src(Params.Np);
            obj.OnMp = Params.Src(Params.Mp);
            obj.OnGG = Params.Src(Params.GridGamma);
            obj.Src = Params.Src;
            %obj.fn   = Params.fn;
            %obj.n    = Params.n;
            obj.Params = Params;
            obj.IsDummy = false;
        end
        
        
        function S = get.Source(obj)
            S=obj.Src;
            %S=obj.OnNp;
            %S=obj.OnMp;
            
           % S = -obj.Src.Source()*(obj.fn(obj.n)+obj.fn(obj.n-1)) ;
           % S(obj.Scatterer.Np) = S(obj.Scatterer.Np) + obj.OnNp;
            %S(obj.Scatterer.Mp) = S(obj.Scatterer.Mp) + obj.OnMp;
        end
           
    end
    

end

