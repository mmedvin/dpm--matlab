classdef NavierStokesDescreteSrc < Tools.Source.SuperSource
    properties %(Dependent = true)
        Source;
        
    end
    
    properties(SetAccess = protected)       
        Scatterer; 
        Params;
        OnNp;
        OnMp;
        OnGG;
        OnIntGG;
        IntGG
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
                                
                
                if 1 %extrapolation
                    x = obj.Scatterer.r.*cos(obj.Scatterer.th);%good only for polar scatterrer!!!
                    y = obj.Scatterer.r.*sin(obj.Scatterer.th);
                    
%                     F = griddata(obj.Params.Grid.X(obj.Params.GridGamma),obj.Params.Grid.Y(obj.Params.GridGamma),obj.OnGG,x,y,'cubic');
%                     
%                     x = (obj.Scatterer.r-0.1).*cos(obj.Scatterer.th);
%                     y = (obj.Scatterer.r-0.1).*sin(obj.Scatterer.th);
%                     

                    %F = griddata(obj.Params.Grid.X(obj.Params.Mp),obj.Params.Grid.Y(obj.Params.Mp),obj.OnMp,x,y,'v4');
                    
                    F = griddata(obj.Params.Grid.X(obj.IntGG),obj.Params.Grid.Y(obj.IntGG),obj.OnIntGG ,x,y,'v4');
                    
                    if 0 
                        %doesn't work, not sure why, 
                        % to use uncomment few things in PolarScatterer and one in TwoTupleExtension
                        % GridLinesIntersection gives all intersection with circle\scatterer
                        % here it first extrapolates from y, x grid lines to circle\scatterer (GridLinesIntersection points), 
                        % next it interpolates from these points to points on Scatterer.th, i.e. normals from GridGamma, used for Extension operator
                        sz=[obj.Params.Grid.Nx,obj.Params.Grid.Ny];
                        [I,J]=ind2sub(sz,obj.Params.Mp);
                        
                        for i=1:numel(obj.Params.GridLinesIntersection.Y.y)
                            %msk = obj.Params.GridLinesIntersection.Y.msk(i);
                            indx = obj.Params.GridLinesIntersection.Y.indx(i);
                            x = obj.Params.Grid.Y(I(J==indx),indx);
                            y = obj.Src(I(J==indx),indx);
                            newx=obj.Params.GridLinesIntersection.Y.y(i);
                            %sx(i) = spline(x,y,newx);
                            sx(i) = interp1(x,y,newx,'pchip','extrap');
                            %ys = fnxtr(s,3);
                        end
                        
                        for i=1:numel(obj.Params.GridLinesIntersection.X.x)
                            %msk = obj.Params.GridLinesIntersection.X.msk(i);
                            indx = obj.Params.GridLinesIntersection.X.indx(i);
                            x = obj.Params.Grid.X(indx,J(I==indx));
                            y = obj.Src(indx,J(I==indx));
                            newx=obj.Params.GridLinesIntersection.X.x(i);
                            %sy(i) = spline(x,y,newx);
                            sy(i) = interp1(x,y,newx,'pchip','extrap');
                            %ys = fnxtr(s,3);
                        end
                        
                        [x,I] = sort([obj.Params.GridLinesIntersection.X.th, obj.Params.GridLinesIntersection.Y.th]);
                        y = [ sx, sy ];
                        %
                        %                     F = spline(x,y(I),obj.Scatterer.th);
                        F = interp1(x,y(I),obj.Scatterer.th,'pchip','extrap');
                    end
                else
                    r  = obj.Scatterer.r;
                    th = obj.Scatterer.th;
                    r0 = obj.Params.r0;
                    ht = obj.Params.ht;
                    
                    [Fn,dFn] = obj.Params.Fn(obj.Params.n);
                    [Fnp1,dFnp1] = obj.Params.Fn(obj.Params.n+1);
                    
                    F = (r.^2).*sin(2*th).*(192*Fnp1 - 4*(4*r.^2-3*r0^2).*(2*Fn/ht + dFn + dFnp1));
                end
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
            obj.IntGG = intersect(Params.GridGamma,Params.Mp);
            obj.OnIntGG = Params.Src(obj.IntGG);
            obj.Src = Params.Src;
            %obj.fn   = Params.fn;
            %obj.n    = Params.n;
            obj.Params = Params;
            obj.IsDummy = false;
        end
        
        
        function S = get.Source(obj)
            S=obj.Src;
            %S=obj.OnNp;
            %=obj.OnMp;
            
           % S = -obj.Src.Source()*(obj.fn(obj.n)+obj.fn(obj.n-1)) ;
           % S(obj.Scatterer.Np) = S(obj.Scatterer.Np) + obj.OnNp;
            %S(obj.Scatterer.Mp) = S(obj.Scatterer.Mp) + obj.OnMp;
        end
           
    end
    

end

