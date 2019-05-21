classdef SuperHomoNavierStokesSolver < Solvers.SuperHomoSolver
    %SUPERHOMONAVIERSTOCKSSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)       
        Qpsi;
        QpsiOmega;
        
        WpsiOmega;
        Wpsi;
        GPO;
    end
    
    properties(Access = protected)
        rhsPsi;
        ExtensionPsi;
        Op;
        OpPsi;
        mQpsi;
        mQpsiOmega;

    end    
    
    methods(Access = protected,Abstract=true)
        % one need to implement the linear operator L
        %        f =
        %       LuPsi(obj,u,msk);
        
        % one need to implement inverse operator G = L^-1
        %        u =
        %       GfPsi(obj,f);
        
        % one computes columns of Q (=P-I) using the following function, see Q,Q0,Q1 above
        %        Pj =
        %      Pcol(obj,GLW,w
        
        TrGpsiPOmega(obj,GLW,w);
    end
    
    methods
        
        function q = get.Qpsi(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            q=obj.mQpsi;
        end
        
        function q = get.QpsiOmega(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            q=obj.mQpsiOmega;
        end

        function w = get.Wpsi(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            w=obj.ExtensionPsi.Wpsi; %obj.myW0;
        end
        
        function w = get.WpsiOmega(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            w=obj.ExtensionPsi.W; %obj.myW0;
        end
        
        function obj = SuperHomoNavierStokesSolver(Arguments)
            obj = obj@Solvers.SuperHomoSolver(Arguments);
        
            Arguments.ExtensionParamsPsi.Grid           = obj.Grid;
            Arguments.ExtensionParamsPsi.Scatterer      = obj.Scatterer;
            Arguments.ExtensionParamsPsi.Basis          = Arguments.Basis;
            Arguments.ExtensionParamsPsi.CoeffsHandle 	= Arguments.CoeffsHandle;
            Arguments.ExtensionParamsPsi.CoeffsParams	= Arguments.CoeffsParams;
            Arguments.ExtensionParamsPsi.CoeffsParams.sigma=0;
            
            %Arguments.ExtensionParams.Coeffs    = obj.Coeffs;
            obj.ExtensionPsi                        = Arguments.ExtensionPsi(Arguments.ExtensionParamsPsi);

            Arguments.DiffOpParams.Grid         = Arguments.Grid;
            Arguments.DiffOpParams.CoeffsHandle = Arguments.CoeffsHandle;
            Arguments.DiffOpParams.CoeffsParams = Arguments.CoeffsParams;
            
            if isfield(Arguments.ScattererParams,'FocalDistance')
                Arguments.DiffOpParams.CoeffsParams.FocalDistance = Arguments.ScattererParams.FocalDistance;
            end
            
            obj.Op = Arguments.DiffOp(Arguments.DiffOpParams);
            
            ArgPsi = Arguments.DiffOpParams;
            ArgPsi.CoeffsParams.sigma=0;
            %obj.OpPsi = Tools.DifferentialOps.LaplacianOpBCinRhs(ArgPsi);
            %Tools.DifferentialOps.LaplacianOp4(ArgPsi);
            %obj.OpPsi = Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs(ArgPsi);
            obj.OpPsi = Arguments.DiffOp(ArgPsi);

            
            
        end
        
  
    end
    
    methods(Access = protected)
        
        function Rhs(obj)
            Rhs@Solvers.SuperHomoSolver(obj);
            obj.CreateRhsPsi();
        end
        
        function UpdatePsi(obj,PsiBC)
            obj.ExtensionPsi.Update(PsiBC);
            tmp = obj.Lu(obj.ExtensionPsi.Wpsi{1},obj.Scatterer.Mp,obj.OpPsi);
            obj.rhsPsi{end} = obj.doitPsi(obj.ExtensionPsi.Wpsi{1},tmp);
            
            obj.calc_QnWpsi();
        end
        
        function CreateRhsPsi(obj)
                       
            tmp=cellfun(@(arg) obj.Lu(arg,obj.Scatterer.Mp,obj.OpPsi),[obj.ExtensionPsi.W, obj.ExtensionPsi.Wpsi],'UniformOutput',false);
            
            obj.rhsPsi = cell(size(tmp));
            for indx=1:(numel(tmp)-1)
                obj.rhsPsi{indx} = obj.doitPsi(obj.ExtensionPsi.W{indx},tmp{indx});
            end
            
            obj.rhsPsi{end} = obj.doitPsi(obj.ExtensionPsi.Wpsi{1},tmp{end});
            
        end
        
        function r = doitPsi(obj,W,t)
            [n,m]=size(W);
            NNZ = nnz(W);
            r = spalloc( n,m,NNZ);
            r(obj.Scatterer.Mp,:) =t;
        end
        
        function calc_QnW(obj)
            if obj.CollectRhs
                obj.Rhs();
                
                GLW        = cellfun(@(arg) obj.Gf(arg),obj.rhs,'UniformOutput',false);
                obj.NewQ   = cellfun(@(arg1,arg2) obj.Qcol(arg1,arg2),GLW, obj.Extension.W,'UniformOutput',false);
                
                obj.calc_QnWpsi();
                
                b = cellfun(@(arg1,arg2) obj.TrGpsiPOmega(arg1,arg2),GLW,obj.Extension.W,'UniformOutput',false);
                obj.GPO = b;
                
            else
                warning('not yet supposed to work')
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
        
        function calc_QnWpsi(obj)
            GLWPsi  = cellfun(@(arg) obj.Gf(arg,obj.OpPsi),obj.rhsPsi,'UniformOutput',false);
            
            
            a = cellfun(@(arg1,arg2) obj.Qcol(arg1,arg2),GLWPsi, [obj.ExtensionPsi.W, obj.ExtensionPsi.Wpsi],'UniformOutput',false);
            
            obj.mQpsiOmega   = a(1:end-1);%cellfun(@(arg1,arg2) plus(arg1,arg2),a(1:end-1),b,'UniformOutput',false);
            obj.mQpsi=a(end);
            
        end
        
        function Expand(obj)
            Expand@Solvers.SuperHomoSolver(obj);
            
            obj.ExtensionPsi.Expand(); 
        end
        
    end
    
end

