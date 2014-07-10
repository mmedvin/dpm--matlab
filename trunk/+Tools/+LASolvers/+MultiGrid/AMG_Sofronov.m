classdef AMG_Sofronov < Tools.LASolvers.AbstractLA_Solver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        D;
    end
    
    methods
        function obj = AMG_Sofronov(A,k,dec,meth,smothop,a,b,lbase,mbase,problem,variant)
            % m2d_preprop - pre-processor
            
            obj.D = cell(5,k+1); % cell of matrices P,A,R of each grid, they will be stored in the memory
            obj.D{1,j}=A;
            
            for j=k+1:-1:1
                
                % number of grid points on grid j-1
                
                L = lbase*(dec^(j-1));
                M = mbase*(dec^(j-1));
                N = L*M;
                
                %if ( (meth == 'G') || ( (meth == 'A') && (j == k+1) ) )
                %    C{1,j} = matr_A(j,lbase,mbase,a,b,dec,problem); % matrix of the system
                %    %         C{1,j} = A_geom2(j,lbase,mbase,a,b,dec,problem);% another variant
                %end
                
                if (j > 1)
                    D{3,j} = restr(j,lbase,mbase,dec,variant); % restriction
                    D{2,j} = prolong(j,lbase,mbase,dec,variant); % interpolation
                end
                
                if( (meth == 'A') && (j ~= k+1) )
                    D{1,j} = A_algeb(D{3,j+1},D{1,j+1},D{2,j+1}); % algebraic coarse operator
                end
                
                D{4,j} = N;
                
                if(strcmp(smothop,'SI') == 1)
                    D{5,j} = eigs(D{1,j}); % eigenvalues for fixed point iteration method
                end
                
            end
            
        end
    end
    
end

