classdef SuperHomoSolver < handle
% SuperHomoSolver is an abstract class for Difference Potentials Method (to be found http://www.math.utah.edu/~mmedvin/MedvinskyCV.html#publications)
% This class implements main functionality\algorithm and deligate the
% concrete functionality outside
% The terminology used here relate to the terminology in: 
% M. Medvinsky, S. Tsynkov, E. Turkel, The Method of Difference Potentials for the Helmholtz Equation Using Compact High Order Schemes, Journal of Scientific Computing, 53, No. 1 (2012) pp. 150-193. This paper also has an erratum .
% and also in 
% M. Medvinsky, Numerical solution of Maxwell's equations, Ph.D. Dissertation Tel Aviv University, 2013.
% M. Medvinsky, S. Tsynkov, E. Turkel,High Order Numerical Simulation of the Transmission and Scattering of Waves Using the Method of Difference Potentials, Journal of Computational Physics, 243 (2013) pp. 305-322.
    
    properties(Access = public)
        %
        % % Q is a discrete operator, a matrix with columns (each per basis function) computed by Q_gamma = P_gamma-I, where
        % P_gamma = Tr P_omega is a discrete Calderon Projection (Tr - stands for vector trace and P_omega is Calderon Potential )
        % Q = [Q0 | Q1], i.e. Q is matrix concatenation of matricex Q0 and Q1,
        % no data saved in Q,Q0,Q1 they are implemented via get method, the actual data is hidden in myQ0 and myQ1
        %
        Q; Q0; Q1;
        
        % W is a choice of arbitrary function that have to satisfy Tr W = W|_gamma = xi_gamma, see ection 2.1.3. 
        % W = 0 outside of gamma
        % W = [W0 | W1], i.e. W is matrix concatenation of matricex W0 and W1,
        W; W0; W1; 
        
        GridGamma; % the numerical counterpart of the interface shape
        Np; % is the N-plus set
        Nm; % is the N-minus set
        
        Scatterer; %cosider to make it protected
    end
    
    
    properties(Access = protected)
        
        IsReadyQnW = false; % used to prevent multiple calculation of Q and W
        
        myQ0; myQ1; % see Q, Q0, Q1 above
		myW0; myW1; % see W, W0, W1 above
        
        
        Grid;   %Link to Grid Class instance
        Basis;  %Link to Basis Class instance
                
        
        % class handle for Scatterer class (information about the interface shape and more) + the parameters to be sent to its constructor
        ScattererHandle;
        ScattererParams;
        
        % an instance  + class handle for Coeffs class (coefficients of the equation) + the parameters to be sent to its constructor
        Coeffs;
        CoeffsHandle;
        CoeffsParams;
        
        %%%%%%%%%%%%%%%%%%%
        NoSource = Tools.Source.SuperHelmholtzSource();  %construct empty source           
        
        %%%%%%%%%%%%%%%%%%
        % artificial rhs, not the rhs of the original problem
        rhs0; rhs1;
		CollectRhs = 1; % for old versions compatibility
        
        f; % temporary variable used as 'rhs' in solve
    end
    
    
    methods(Access = public,Abstract=true)
        % one need to implement P_omega, i.e. calderon potential
        %u = 
            P_Omega(xi_gamma);
    end
    
    methods(Access = protected,Abstract=true)
        
        % one need to implement the linear operator L
%        f = 
            Lu(obj,u,msk);
        
        % one need to implement inverse operator G = L^-1
%        u = 
            Gf(obj,f);
            
        % one computes columns of Q (=P-I) using the following function, see Q,Q0,Q1 above
%        Qj = 
            Qcol(obj,GLW,w);
              
    end

    % constructor + the get methods for the properties
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
            
            w=[obj.myW0,obj.myW1];
		end 
        
       function w0 = get.W0(obj)           
            if obj.IsReadyQnW == false
               obj.calc_QnW();
            end            
            
            w0=obj.myW0;
	   end 

	   function w1 = get.W1(obj)
		   if obj.IsReadyQnW == false
			   obj.calc_QnW();
		   end
		   
		   w1=obj.myW1;
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
        
        
        function obj = SuperHomoSolver(Basis,Grid,CoeffsHandle,CoeffsParams,ScattererHandle,ScattererParams,CollectRhs)
            if nargin == 0, error('Costructor called without args'); end
			if exist('CollectRhs','var'), obj.CollectRhs = CollectRhs;	end
            
            obj.Grid = Grid;
                    
            if isstruct(Basis)
                obj.Basis=Basis;
            else % support older versions
                warning('you are using an old version of this class constructor') %#ok<WNTAG>
                obj.Basis.Indices = Basis;
                obj.Basis.Handle  = @FourierBasis;
                obj.Basis.NBss = numel(obj.Basis.Indices);
                obj.Basis.AddParams = [];
            end
            
            obj.ScattererHandle = ScattererHandle;
            obj.ScattererParams = ScattererParams;            
            obj.Scatterer = ScattererHandle(Grid,obj.ScattererParams);
            
            obj.CoeffsHandle = CoeffsHandle;
            obj.CoeffsParams = CoeffsParams;            
            obj.Coeffs = CoeffsHandle(obj.Scatterer.TheScatterer,obj.CoeffsParams);
            
            tmp=spalloc( obj.Grid.Nx*obj.Grid.Ny,obj.Basis.NBss,numel(obj.GridGamma)*obj.Basis.NBss);
            obj.myW0= tmp; 
            obj.myW1=tmp; 
            obj.f=zeros(Grid.Nx,Grid.Ny);
        end        
    end
    
    methods(Access = protected)
    
         function Rhs(obj)
             obj.Expand();
            
            obj.rhs0=spalloc( obj.Grid.Nx*obj.Grid.Ny,obj.Basis.NBss,numel(obj.GridGamma)*obj.Basis.NBss);
            obj.rhs1=spalloc( obj.Grid.Nx*obj.Grid.Ny,obj.Basis.NBss,numel(obj.GridGamma)*obj.Basis.NBss);
            
             obj.rhs0(obj.Scatterer.Mp,:) = obj.Lu(obj.myW0,obj.Scatterer.Mp);
             obj.rhs1(obj.Scatterer.Mp,:) = obj.Lu(obj.myW1,obj.Scatterer.Mp);
         end
         
		 function calc_QnW(obj)
			 if obj.CollectRhs
				 obj.Rhs();
				 
				 GLW = obj.Gf([obj.rhs0,obj.rhs1]);
				 
				 obj.myQ0 = obj.Qcol( GLW(:,           1:obj.Basis.NBss   )			,obj.myW0 );
				 obj.myQ1 = obj.Qcol( GLW(:,(obj.Basis.NBss+1):2*obj.Basis.NBss )	,obj.myW1 );
			 else
				 Ngg = numel(obj.GridGamma);
				 w0 = spalloc(obj.Grid.Nx,obj.Grid.Ny,Ngg);
				 w1 = spalloc(obj.Grid.Nx,obj.Grid.Ny,Ngg);
				 
				 for j = 1:obj.Basis.NBss
					 
					 [w0(obj.GridGamma),w1(obj.GridGamma)] = obj.ExpandedBasis(obj.Basis.Indices(j)) ;
					 
					 obj.myW0(obj.GridGamma,j) = w0(obj.GridGamma);
					 obj.myW1(obj.GridGamma,j) = w1(obj.GridGamma);
					 
					 %obj.myQ0(:,j) = obj.Qcol2(w0);
					 %obj.myQ1(:,j) = obj.Qcol2(w1);
                     
                     GLW = obj.Solve(w0);
                     obj.myQ0(:,j) = obj.Qcol(GLW,w0(:));
                     
                     GLW = obj.Solve(w1);
                     obj.myQ1(:,j) = obj.Qcol(GLW,w1(:));
                     
                     
				 end
				 
			 end
			 obj.IsReadyQnW = true;
		 end
        
         function u = Solve(obj,x)
             obj.f(obj.Scatterer.Mp) = obj.Lu(x(:),obj.Scatterer.Mp);
             %obj.Truncate(f);
             u = obj.Gf(obj.f);%(:));
         end
         
         
        function Expand(obj)
          [tmp1,tmp2] = arrayfun(@(n) obj.ExpandedBasis(n),obj.Basis.Indices,'UniformOutput',false);
          obj.myW0(obj.GridGamma,:) =  cell2mat(tmp1);
          obj.myW1(obj.GridGamma,:) =  cell2mat(tmp2);
        end
           
        function [xi0j,xi1j] = ExpandedBasis(obj,n)
            
            Xi   = obj.Basis.Handle(obj.Scatterer.BasisArg,n,obj.Basis.AddParams);
            NoXi = obj.Basis.Handle();            
                        
            xi0j = obj.Scatterer.Expansion(Xi,NoXi,obj.NoSource,obj.Coeffs);
            xi1j = obj.Scatterer.Expansion(NoXi,Xi,obj.NoSource,obj.Coeffs);
            
        end
    end
end
