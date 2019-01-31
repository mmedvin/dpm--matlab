classdef SuperHomoNavierStokesSolver < Solvers.SuperHomoSolver
    %SUPERHOMONAVIERSTOCKSSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        P;
    end
    
    properties(Access = protected)
        mP;
    end
    
    methods(Access = protected,Abstract=true)
        Pcol(obj,GLW,w);
    end
    
    methods
        
        function q = get.P(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            q=obj.mP;
        end
        
        function qj = Pj(obj,j)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            qj=obj.mP{j};
            
        end
        
        function obj = SuperHomoNavierStokesSolver(Arguments)
            obj = obj@Solvers.SuperHomoSolver(Arguments);
        end
        
  
    end
    
    methods(Access = protected)
        function calc_QnW(obj)
            if obj.CollectRhs
                obj.Rhs();
                
                GLW        = cellfun(@(arg) obj.Gf(arg),obj.rhs,'UniformOutput',false);
                obj.NewQ   = cellfun(@(arg1,arg2) obj.Qcol(arg1,arg2),GLW, obj.Extension.W,'UniformOutput',false);
                obj.mP   = cellfun(@(arg1,arg2) obj.Pcol(arg1,arg2),GLW, obj.Extension.W,'UniformOutput',false);
            else
                obj.Expand();
                for indx=1:numel(obj.Extension.W)
                    for j = 1:size(obj.Extension.W{indx},2)
                        GLW                    = obj.Solve(obj.Extension.W{indx}(:,j));
                        obj.NewQ{indx}(:,j)    = obj.Qcol(GLW,obj.Extension.W{indx}(:,j));
                        obj.mP{indx}(:,j)    = obj.Pcol(GLW,obj.Extension.W{indx}(:,j));
                    end
                end
                
            end
            obj.IsReadyQnW = true;
        end
    end
    
end

