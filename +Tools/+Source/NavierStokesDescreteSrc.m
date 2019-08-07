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
        %fn;
        %n;
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
                    if 0
                        x = obj.Scatterer.r.*cos(obj.Scatterer.th);%good only for polar scatterrer!!!
                        y = obj.Scatterer.r.*sin(obj.Scatterer.th);
                        
                        %                     F = griddata(obj.Params.Grid.X(obj.Params.GridGamma),obj.Params.Grid.Y(obj.Params.GridGamma),obj.OnGG,x,y,'cubic');
                        %
                        %                     x = (obj.Scatterer.r-0.1).*cos(obj.Scatterer.th);
                        %                     y = (obj.Scatterer.r-0.1).*sin(obj.Scatterer.th);
                        %
                        F=zeros(size(x));
                        X = obj.Params.Grid.X(obj.Params.Mp);
                        Y = obj.Params.Grid.Y(obj.Params.Mp);
                        Z = obj.OnMp;
                        %                     parfor i=1:numel(x)
                        %                         F(i) = griddata(X,Y,Z,x(i),y(i),'v4');
                        %                     end
                        
                        Foo = scatteredInterpolant(X,Y,Z,'linear','nearest');
                        F2=Foo(x,y);
                        %F = griddata(X,Y,Z,x,y,'v4');
                        
                        
                        
                        %F = griddata(obj.Params.Grid.X(obj.IntGG),obj.Params.Grid.Y(obj.IntGG),obj.OnIntGG ,x,y,'v4');
                        
                    else
                        
                        %                         x = [obj.Params.GridLinesIntersection.Y.x , obj.Params.GridLinesIntersection.X.x];
                        %                         y = [obj.Params.GridLinesIntersection.Y.y , obj.Params.GridLinesIntersection.X.y];
                        %
                        %                         %F0 = griddata(obj.Params.Grid.X(obj.Params.Mp),obj.Params.Grid.Y(obj.Params.Mp),obj.OnMp,x,y,'v4');
                        %                         F0 = griddata(obj.Params.Grid.X(obj.IntGG),obj.Params.Grid.Y(obj.IntGG),obj.OnIntGG ,x,y,'v4');
                        %
                        %                         th = [obj.Params.GridLinesIntersection.Y.th , obj.Params.GridLinesIntersection.X.th];
                        %                         %F = interp1(th,F0,obj.Scatterer.th,'pchip');
                        %                         F = spline(th,F0,obj.Scatterer.th);
                        %
                        %                     else
                        % also need few things in PolarScatterer and one in TwoTupleExtension
                        % GridLinesIntersection gives all intersection with circle\scatterer
                        % here it first extrapolates from y, x grid lines to circle\scatterer (GridLinesIntersection points),
                        % next it interpolates from these points to points on Scatterer.th, i.e. normals from GridGamma, used for Extension operator
                        
                        interpType = 'pcchip';%'makima';
                        
                        sz=[obj.Params.Grid.Nx,obj.Params.Grid.Ny];
                        [I,J]=ind2sub(sz,obj.Params.Mp);
                        N = numel(obj.Params.GridLinesIntersection.Y.y);
                        sx=zeros(1,N);
                        parfor i=1:N
                            %msk = obj.Params.GridLinesIntersection.Y.msk(i);
                            indx = obj.Params.GridLinesIntersection.Y.indx(i);
                            x = obj.Params.Grid.Y(I(J==indx),indx);
                            y = obj.Src(I(J==indx),indx);
                            newx=obj.Params.GridLinesIntersection.Y.y(i);
                            %sx(i) = spline(x,y,newx);
                            sx(i) = interp1(x,y,newx,interpType,'extrap');
                            %sx(i) = ppval(fnxtr(csaps(x,y)),newx);
                        end
                        
                        N = numel(obj.Params.GridLinesIntersection.X.x);
                        sy=zeros(1,N);
                        
                        parfor i=1:N
                            %msk = obj.Params.GridLinesIntersection.X.msk(i);
                            indx = obj.Params.GridLinesIntersection.X.indx(i);
                            x = obj.Params.Grid.X(indx,J(I==indx));
                            y = obj.Src(indx,J(I==indx));
                            newx=obj.Params.GridLinesIntersection.X.x(i);
                            %sy(i) = spline(x,y,newx);
                            sy(i) = interp1(x,y,newx,interpType,'extrap');
                            %sy(i) = ppval(fnxtr(csaps(x,y)),newx);
                        end
                        
                        x=[obj.Params.GridLinesIntersection.Y.th, obj.Params.GridLinesIntersection.X.th];
                        %[x,I] = sort(x);
                        y = [ sx, sy ];
                        
                        F = spline(x,y,obj.Scatterer.th);
                        %F = interp1(x,y,obj.Scatterer.th,interpType);
                        %F =  ppval(csaps(x,y(I)),obj.Scatterer.th);
                        % F =  ppval(csaps(x,y),obj.Scatterer.th);
                    end
                else % exact (no Convention Term?)
                    r  = obj.Scatterer.r;
                    th = obj.Scatterer.th;
                    r0 = obj.Params.r0;
                    ht = obj.Params.ht;
                    
                    n = obj.Params.tn;
                    
                    RN = obj.Params.RN;
                    ExParams = struct('r0',r0,'r',r);
                    
                    [Fn,dFn] = obj.Params.Fn(n);
                    
                    Omega = obj.Params.Exact.Omega(th,ExParams).*Fn;
                    gn = obj.Params.Exact.L_np1(th,ExParams, obj.Params.k,obj.Params.Fn,n);
                    
                    D=0;%?
                    
                    F =gn;% obj.Params.TStpSrc(th,r,gn,Omega,n,D);
                    
                    %                    [Fn,dFn] = obj.Params.Fn(obj.Params.n-1);
                    %                     [Fnp1,dFnp1] = obj.Params.Fn(obj.Params.n);
                    %
                    %
                    %                     Omega =  4*(r.^2).*(4*r.^2-3*r0^2).*sin(2*th).*Fn;
                    %                     LaplcianOmega = 192.*sin(2*th).*(r.^2).*Fn; %=384xy
                    %
                    %                     %ftilde = obj.Params.FTILDE(obj.Params.n,G) + obj.Params.FTILDE(obj.Params.n-1,G) ;
                    %                     ftilde = Omega*dFn - LaplcianOmega*Fn/RN ...
                    %                            + Omega*dFnp1 - LaplcianOmega*Fnp1/RN;
                    %
                    %                     F = -(LaplcianOmega +(2*RN/ht)*Omega + RN*ftilde);
                    %
                    %                     %F = (r.^2).*sin(2*th).*(192*Fnp1 - 4*(4*r.^2-3*r0^2).*(2*Fn/ht + dFn + dFnp1));
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
            
            obj.Src = zeros(size(Params.Src));
            obj.Src(Params.Mp) = Params.Src(Params.Mp);
            
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
