classdef SingleScatterer < Tools.Scatterer.SupperScatterer
      properties
        Grid;
        
        Mp;
        Mm;
        Np;
        Nm;
        
        GridGamma;
        
        Size;
        
        In; %can't find better name yet...
        Out;
        
        
              
    end
    
    methods(Abstract = true, Access = public)
        %  Grid = AnotherGrid(obj,Grid);
        %res =
        Expansion(obj,Xi0,Xi1,F,WaveNumber);
        
    end
    
    methods
        function obj = SingleScatterer(Grid)
            if exist('Grid','var')
                obj.Grid = Grid;                
            else
                error('Costructor called without arg Grid');
            end                                 
        end 
        
        function sz = get.Size(obj)
            
            sz = obj.Grid.Size;
            
        end
        
        function in = get.In(obj)
            
            in = setdiff(obj.Mp,obj.GridGamma);            
            
        end
        
        function out = get.Out(obj)
            
            out = setdiff(obj.Mm,obj.GridGamma);
            
        end
        
        
    end 
    
    methods(Access =protected)
        
        function SplitGrid(obj, Stencil)          
           switch Stencil
                case 5
                    obj.SplitGrid5();
                case 9
                    obj.SplitGrid9();
                case 13
                    obj.SplitGrid7x7();
            end
        end
        function SplitGrid9(obj)%,R,r0)
            
            L= obj.Inside();%R <= r0;
            obj.Mp = find(L);
            
            obj.Mp = find(obj.Inside);

            N=obj.Size;
            
            obj.Np=zeros(N);
            obj.Nm=zeros(N);
            
            %%%%%%%%
            [i,j] = find(L);
            
            Nn = [  msub2ind(N,i-1,j-1), msub2ind(N,i-1,j), msub2ind(N,i-1,j+1),...
                msub2ind(N,i,j-1)  , msub2ind(N,i,j)  ,	msub2ind(N,i,j+1),...
                msub2ind(N,i+1,j-1), msub2ind(N,i+1,j),	msub2ind(N,i+1,j+1)];
            
            obj.Np(Nn(Nn>0))=1;
            obj.Np=find(obj.Np)';
            
            %%%%%%%%%%%%%%
            L = obj.Outside;%R>r0;
            obj.Mm = find(L);
            
            [i,j] = find(L);
            
            Nn = [  msub2ind(N,i-1,j-1), msub2ind(N,i-1,j),	msub2ind(N,i-1,j+1),...
                msub2ind(N,i,j-1)  , msub2ind(N,i,j)  ,	msub2ind(N,i,j+1),...
                msub2ind(N,i+1,j-1), msub2ind(N,i+1,j),	msub2ind(N,i+1,j+1)];
            
            obj.Nm(Nn(Nn>0))=1;
            obj.Nm=find(obj.Nm)';
            %%%%%%%%%%%%%%%%%
            
            %             Split.GridGamma = intersect(Np,Nm)';
            %             Split.Inside = setdiff(Mp,Split.GridGamma);
            % obj.Outside= setdiff(Mm,Split.GridGamma);
            %             Split.Mp=Mp;
            %             Split.Np=Np;
            %             Split.Nm=Nm;
            % obj.Mm=Mm;
            
            %[Split.cols,Split.rows] = size(R);
            
            function ind=msub2ind(N,i,j)
                
                %     outside = i<1 | i > N(1) | j<1 | j > N(2);
                              
                j(j==0)=N(2);
                j(j==N(2)+1)=1;

                inside = i>=1 & i<=N(1) & j>=1 & j<=N(2);
                
                % ind = (i+(j-1)*N(2));
                ind=zeros(size(i));
                ind(inside) = sub2ind(N,i(inside),j(inside));
                % ind(outside)=0;
                
                % if ind ~= sub2ind(N,i,j)
                %     error(['msub2ind error, N=' num2str(N) ',i=' num2str(i) ',j=' num2str(j) ',ind=' num2str(ind)])
                % end
            end
		end
		
	function SplitGrid5(obj)%,R,r0)
            
            L= obj.Inside();%R <= r0;
            obj.Mp = find(L);
            
            obj.Mp = find(obj.Inside);

            N=obj.Size;
            
            obj.Np=zeros(N);
            obj.Nm=zeros(N);
            
            %%%%%%%%
            [i,j] = find(L);
            
            Nn = [					 msub2ind(N,i-1,j),...
                msub2ind(N,i,j-1)  , msub2ind(N,i,j)  ,	msub2ind(N,i,j+1),...
									 msub2ind(N,i+1,j)];
            
            obj.Np(Nn(Nn>0))=1;
            obj.Np=find(obj.Np)';
            
            %%%%%%%%%%%%%%
            L = obj.Outside;%R>r0;
            obj.Mm = find(L);
            
            [i,j] = find(L);
            
            Nn = [					 msub2ind(N,i-1,j),	...
                msub2ind(N,i,j-1)  , msub2ind(N,i,j)  ,	msub2ind(N,i,j+1),...
									 msub2ind(N,i+1,j)];
            
            obj.Nm(Nn(Nn>0))=1;
            obj.Nm=find(obj.Nm)';
            %%%%%%%%%%%%%%%%%
            
            %             Split.GridGamma = intersect(Np,Nm)';
            %             Split.Inside = setdiff(Mp,Split.GridGamma);
            % obj.Outside= setdiff(Mm,Split.GridGamma);
            %             Split.Mp=Mp;
            %             Split.Np=Np;
            %             Split.Nm=Nm;
            % obj.Mm=Mm;
            
            %[Split.cols,Split.rows] = size(R);
            
            function ind=msub2ind(N,i,j)
                
                %     outside = i<1 | i > N(1) | j<1 | j > N(2);
                              
                j(j==0)=N(2);
                j(j==N(2)+1)=1;

                inside = i>=1 & i<=N(1) & j>=1 & j<=N(2);
                
                % ind = (i+(j-1)*N(2));
                ind=zeros(size(i));
                ind(inside) = sub2ind(N,i(inside),j(inside));
                % ind(outside)=0;
                
                % if ind ~= sub2ind(N,i,j)
                %     error(['msub2ind error, N=' num2str(N) ',i=' num2str(i) ',j=' num2str(j) ',ind=' num2str(ind)])
                % end
            end
        end	
        
    function SplitGrid7x7(obj)%,R,r0)
            
            L= obj.Inside();%R <= r0;
            obj.Mp = find(L);
            
            obj.Mp = find(obj.Inside);

            N=obj.Size;
            
            obj.Np=zeros(N);
            obj.Nm=zeros(N);
            
            %%%%%%%%
            [i,j] = find(L);
            
            Nn = [  
                                                                           msub2ind(N,i-3,j),                                                          ...
                                                                           msub2ind(N,i-2,j),                                                          ...
                                                                           msub2ind(N,i-1,j),                                                          ...
                msub2ind(N,i,j-3),   msub2ind(N,i,j-2), msub2ind(N,i,j-1), msub2ind(N,i,j)  , msub2ind(N,i,j+1), msub2ind(N,i,j+2), msub2ind(N,i,j+3), ...
                                                                           msub2ind(N,i+1,j),                                                          ...
                                                                           msub2ind(N,i+2,j),                                                          ...
                                                                           msub2ind(N,i+3,j)                                                          ...                                                                           
                  ];
            
            obj.Np(Nn(Nn>0))=1;
            obj.Np=find(obj.Np)';
            
            %%%%%%%%%%%%%%
            L = obj.Outside;%R>r0;
            obj.Mm = find(L);
            
            [i,j] = find(L);
            
            Nn = [  
                                                                           msub2ind(N,i-3,j),                                                          ...
                                                                           msub2ind(N,i-2,j),                                                          ...
                                                                           msub2ind(N,i-1,j),                                                          ...
                msub2ind(N,i,j-3),   msub2ind(N,i,j-2), msub2ind(N,i,j-1), msub2ind(N,i,j)  , msub2ind(N,i,j+1), msub2ind(N,i,j+2), msub2ind(N,i,j+3), ...
                                                                           msub2ind(N,i+1,j),                                                          ...
                                                                           msub2ind(N,i+2,j),                                                          ...
                                                                           msub2ind(N,i+3,j)                                                          ...                                                                           
                  ];

            
            obj.Nm(Nn(Nn>0))=1;
            obj.Nm=find(obj.Nm)';
            %%%%%%%%%%%%%%%%%
            
            %             Split.GridGamma = intersect(Np,Nm)';
            %             Split.Inside = setdiff(Mp,Split.GridGamma);
            % obj.Outside= setdiff(Mm,Split.GridGamma);
            %             Split.Mp=Mp;
            %             Split.Np=Np;
            %             Split.Nm=Nm;
            % obj.Mm=Mm;
            
            %[Split.cols,Split.rows] = size(R);
            
            function ind=msub2ind(N,i,j)
                
                %     outside = i<1 | i > N(1) | j<1 | j > N(2);
                              
                j(j==0)=N(2);
                j(j==N(2)+1)=1;

                inside = i>=1 & i<=N(1) & j>=1 & j<=N(2);
                
                % ind = (i+(j-1)*N(2));
                ind=zeros(size(i));
                ind(inside) = sub2ind(N,i(inside),j(inside));
                % ind(outside)=0;
                
                % if ind ~= sub2ind(N,i,j)
                %     error(['msub2ind error, N=' num2str(N) ',i=' num2str(i) ',j=' num2str(j) ',ind=' num2str(ind)])
                % end
            end
		end
		
    
	end
end
