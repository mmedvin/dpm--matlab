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
        NewQ;
        
        % W is a choice of arbitrary function that have to satisfy Tr W = W|_gamma = xi_gamma, see ection 2.1.3. 
        % W = 0 outside of gamma
        % W = [W0 | W1], i.e. W is matrix concatenation of matricex W0 and W1,
        W; W0; W1; 
        %NewW; %experiment
        
        GridGamma; % the numerical counterpart of the interface shape
        Np; % is the N-plus set
        Nm; % is the N-minus set
        
        Scatterer; %cosider to make it protected
    end
    
    
    properties(Access = protected)
        
        IsReadyQnW = false; % used to prevent multiple calculation of Q and W
        
        %myQ0; myQ1; % see Q, Q0, Q1 above
		%myW0; myW1; % see W, W0, W1 above
        
        
        Grid;   %Link to Grid Class instance
        %Basis;  %Link to Basis Class instance
                
        
        % class handle for Scatterer class (information about the interface shape and more) + the parameters to be sent to its constructor
        ScattererHandle;
        ScattererParams;
        
        % an instance  + class handle for Coeffs class (coefficients of the equation) + the parameters to be sent to its constructor
        %Coeffs;
        CoeffsHandle;
        CoeffsParams;
        
        %%%%%%%%%%%%%%%%%%%
        %NoSource = Tools.Source.SuperHelmholtzSource();  %construct empty source           
        
        %%%%%%%%%%%%%%%%%%
        % artificial rhs, not the rhs of the original problem
        %rhs0; rhs1;
        rhs;
		CollectRhs = 1; % for old versions compatibility
        
        f; % temporary variable used as 'rhs' in solve
        
        Extension;
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
        function obj = SuperHomoSolver(Arguments)		
            if nargin == 0, error('Costructor called without args'); end
			if isfield(Arguments,'CollectRhs'), obj.CollectRhs = Arguments.CollectRhs;	end
            
            obj.Grid            = Arguments.Grid;
            %obj.Basis           = Arguments.Basis;
                    
            obj.ScattererHandle = Arguments.ScattererHandle;
            obj.ScattererParams = Arguments.ScattererParams;            
            obj.Scatterer 		= obj.ScattererHandle(obj.Grid,obj.ScattererParams);
            
            obj.CoeffsHandle 	= Arguments.CoeffsHandle;
            obj.CoeffsParams 	= Arguments.CoeffsParams;            
            %obj.Coeffs 			= obj.CoeffsHandle(obj.Scatterer.TheScatterer,obj.CoeffsParams);%TODO: Consider to remove TheScatterer
            
            Arguments.ExtensionParams.Grid          = obj.Grid;
            Arguments.ExtensionParams.Scatterer     = obj.Scatterer;
            Arguments.ExtensionParams.Basis         = Arguments.Basis;
            Arguments.ExtensionParams.CoeffsHandle 	= Arguments.CoeffsHandle;
             Arguments.ExtensionParams.CoeffsParams	= Arguments.CoeffsParams;            

            %Arguments.ExtensionParams.Coeffs    = obj.Coeffs;
            obj.Extension                       = Arguments.Extension(Arguments.ExtensionParams);

            obj.f=zeros(obj.Grid.Nx,obj.Grid.Ny);
        end        


        function q = get.Q(obj)           
            if obj.IsReadyQnW == false
               obj.calc_QnW();
            end            
            %q=[obj.myQ0,obj.myQ1];
            q=obj.NewQ;
        end 
      
        function qj = Qj(obj,j)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            qj=obj.NewQ{j};%obj.myQ0;
            
        end
        
        function q0 = get.Q0(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            q0=obj.NewQ{1};%obj.myQ0;
            
        end
        
        function q1 = get.Q1(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            q1=obj.NewQ{2};%obj.myQ1;
            
        end
        
        function w = get.W(obj)           
            if obj.IsReadyQnW == false
               obj.calc_QnW();
            end            
            
            %wold=[obj.myW0,obj.myW1];
            w=obj.Extension.W;            
        end 
        
        function wj = Wj(obj,j)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            wj=obj.Extension.W{j}; %obj.myW0;
        end
 
        
       function w0 = get.W0(obj)           
            if obj.IsReadyQnW == false
               obj.calc_QnW();
            end            
            
            w0=obj.Extension.W{1}; %obj.myW0;
	   end 

	   function w1 = get.W1(obj)
		   if obj.IsReadyQnW == false
			   obj.calc_QnW();
		   end
		   
		   w1=obj.Extension.W{2};%myW1;
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
        
    end
    
    methods(Access = protected)
    
        function Rhs(obj)
            obj.Expand();
            tmp=cellfun(@(arg) obj.Lu(arg,obj.Scatterer.Mp),obj.Extension.W,'UniformOutput',false);
            
            obj.rhs = cell(size(tmp));
            for indx=1:numel(tmp)
                [n,m]=size(obj.Extension.W{indx});
                NNZ = nnz(obj.Extension.W{indx});
                obj.rhs{indx} = spalloc( n,m,NNZ);
                obj.rhs{indx}(obj.Scatterer.Mp,:) = tmp{indx};
            end
            
        end
         
		 function calc_QnW(obj)
			 if obj.CollectRhs
				 obj.Rhs();
				                                   
                 GLW        = cellfun(@(arg) obj.Gf(arg),obj.rhs,'UniformOutput',false);
                 obj.NewQ   = cellfun(@(arg1,arg2) obj.Qcol(arg1,arg2),GLW, obj.Extension.W,'UniformOutput',false);
 
             else
                 obj.Expand();                 
                 for indx=1:numel(obj.Extension.W)
                     for j = 1:size(obj.Extension.W{indx},2)
                         GLW                    = obj.Solve(obj.Extension.W{indx}(:,j));
                         obj.NewQ{indx}(:,j)    = obj.Qcol(GLW,obj.Extension.W{indx}(:,j));
                     end
                 end
                 
			 end
			 obj.IsReadyQnW = true;
		 end
        
         function u = Solve(obj,x)
             obj.f(obj.Scatterer.Mp) = obj.Lu(x(:),obj.Scatterer.Mp);
             u = obj.Gf(obj.f);%(:));
         end
         
         function Expand(obj)
             obj.Extension.Expand();
         end
    end
end
