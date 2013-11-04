classdef SuperHomoSolver < handle
    
    properties(Access = public)
        Q;
        Q0; 
        Q1;
        W;
        W0;
        W1;
        GridGamma;
        Np;
        Nm;
        
        Scatterer; %cosider to make it protected
    end
    
    
    properties(Access = protected)
        
        IsReadyQnW = false;
        myQ0;
        myQ1;
        
        Grid;
        
        Basis;
        
%         BasisIndices;
%         BasisHandle;
%         BasisAddParams;
%         NBss=0;
                
        
        ScattererClsHandle;
        ScattererAddParams;
        
        WaveNumber;
        WaveNumberClsHandle;
        WaveNumberAddParams;
        
        NoSource = Tools.Source.SuperHelmholtzSource();  %construct empty source           
        
        rhs0;
        rhs1;
    end
    
    
    methods(Access = public,Abstract=true)
        %u = 
            P_Omega(xi_gamma);
    end
    
    methods(Access = protected,Abstract=true)
%        f = 
            Lu(obj,u,msk);
%        u = 
            Gf(obj,f);
%        Qj = 
            Qcol(obj,GLW,w);
              
    end

    
    methods      
        
        function q = get.Q(obj)           
            
            if obj.IsReadyQnW == false
               obj.calc_QnW();
            end            
            q=[obj.myQ0,obj.myQ1];
        end 
        
        function q0 = get.Q0(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            q0=obj.myQ0;
            
        end
        
        function q1 = get.Q1(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            q1=obj.myQ1;
            
        end
        
        function w = get.W(obj)           
            
            if obj.IsReadyQnW == false
               obj.calc_QnW();
            end            
            w=[obj.W0,obj.W1];
        end 
        
        
        function GG = get.GridGamma(obj)
            GG = obj.Scatterer.GridGamma;
        end
        
         function N = get.Nm(obj)
            N=obj.Scatterer.Nm;
         end
        
          function N = get.Np(obj)
            N=obj.Scatterer.Np;
        end
        
        
        function obj = SuperHomoSolver( ...
                Basis,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams)
            if nargin == 0
                error('Costructor called without args');
            end
            
            obj.Grid = Grid;
                    
            if isstruct(Basis)
                obj.Basis=Basis;
%                 obj.BasisIndices    = Basis.BasisIndices;
%                 obj.BasisHandle     = Basis.BasisHandle;
%                 obj.BasisAddParams  = Basis.AddParams;
            else % support older versions
                warning('you are using an old version of this class constructor') %#ok<WNTAG>
                obj.Basis.Indices = Basis;
                obj.Basis.Handle  = @FourierBasis;
                obj.Basis.NBss = numel(obj.Basis.Indices);
                obj.Basis.AddParams = [];
%                 obj.BasisIndices = Basis;
%                 obj.BasisHandle  = @FourierBasis;
            end
%             obj.NBss = numel(obj.BasisIndices);
            
            
            obj.ScattererClsHandle = ScattererClsHandle;
            obj.ScattererAddParams = ScattererAddParams;            
            obj.Scatterer = ScattererClsHandle(Grid,obj.ScattererAddParams);
            
            obj.WaveNumberClsHandle = WaveNumberClsHandle;
            obj.WaveNumberAddParams = WaveNumberAddParams;            
            obj.WaveNumber = WaveNumberClsHandle(obj.Scatterer.TheScatterer,obj.WaveNumberAddParams);
            
            
            tmp=spalloc( obj.Grid.Nx*obj.Grid.Ny,obj.Basis.NBss,numel(obj.GridGamma)*obj.Basis.NBss);
            obj.W0= tmp; 
            obj.W1=tmp; 
            obj.rhs0=tmp;
            obj.rhs1=tmp;            
        end        
    end
    
    methods(Access = protected)
    
         function Rhs(obj)
             obj.Expand();
             
%               
%             tmp = obj.Lu(obj.W0);% obj.A*obj.W0;
%             obj.rhs0(obj.Scatterer.Mp,:) = tmp(obj.Scatterer.Mp,:);
%             
%             tmp = obj.Lu(obj.W1) %A*W1;
%             obj.rhs1(obj.Scatterer.Mp,:) = tmp(obj.Scatterer.Mp,:);

             obj.rhs0(obj.Scatterer.Mp,:) = obj.Lu(obj.W0,obj.Scatterer.Mp);
             obj.rhs1(obj.Scatterer.Mp,:) = obj.Lu(obj.W1,obj.Scatterer.Mp);
         end
         
         function calc_QnW(obj)
             obj.Rhs();
             
             GLW = obj.Gf([obj.rhs0,obj.rhs1]);
             
             obj.myQ0 = obj.Qcol( GLW(:,           1:obj.Basis.NBss   ) );
             obj.myQ1 = obj.Qcol( GLW(:,(obj.Basis.NBss+1):2*obj.Basis.NBss ) );
             obj.IsReadyQnW = true;
         end
        
        function Expand(obj)
%             M=obj.BasisIndices(end);
%             for j = -M:M
%                 [obj.W0(obj.Scatterer.GridGamma,j+M+1), ...
%                  obj.W1(obj.Scatterer.GridGamma,j+M+1)] = obj.ExpandedBasis(j) ;
%             end
            
%              for j = 1:obj.Basis.NBss
%                 [obj.W0(obj.GridGamma,j), obj.W1(obj.GridGamma,j)] = obj.ExpandedBasis(obj.BasisIndices(j)) ;
%              end
           
          [tmp1,tmp2] = arrayfun(@(n) obj.ExpandedBasis(n),obj.Basis.Indices,'UniformOutput',false);
          obj.W0(obj.GridGamma,:) =  cell2mat(tmp1);
          obj.W1(obj.GridGamma,:) =  cell2mat(tmp2);
             
        end
           
        function [xi0j,xi1j] = ExpandedBasis(obj,n)
            
            Xi   = obj.Basis.Handle(obj.Scatterer.BasisArg,n,obj.Basis.AddParams);
            NoXi = obj.Basis.Handle();            
                        
            xi0j = obj.Scatterer.Expansion(Xi,NoXi,obj.NoSource,obj.WaveNumber);
            xi1j = obj.Scatterer.Expansion(NoXi,Xi,obj.NoSource,obj.WaveNumber);
            
        end
    end
end
