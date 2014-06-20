classdef FlexibleGmres < handle
	
	properties(Access = protected)
        A; 
    end
		
	properties
	end

	methods
	end
	
	methods(Access = protected)
		function [solX,Residual,iter]=my_flexible_gmres(Rhs,InitialGuess,matvec,prec,m,tol,maxit,print)			
			iter=0;
			iouter=0;
			%n=numel(Rhs);
			Residual0=norm(Rhs);
			%b=Rhs;
			solX = InitialGuess;
			Rhs=Rhs-obj.A*solX;
			while  iter < maxit
				%while iter < maxit
				iouter=iouter+1;
				%clear h v w rmat rs c s
				rs(1,1)=1;
				c(1,1)=0;
				s(1,1)=0;
				dnb=norm(Rhs);
				w(:,1)=Rhs/dnb;
				%iinner(iouter) = 0;
				for j=1:m
					%iinner(iouter) = iinner(iouter)+1;
					v(:,j)=feval(prec,w(:,j));
					w(:,j+1)=obj.A*v(:,j); %feval(matvec,v(:,j));
					iter=iter+1;
					for i=1:j
						h(i,j)=w(:,i)'*w(:,j+1);
						w(:,j+1)=w(:,j+1)-h(i,j)*w(:,i);
					end
					h(j+1,j)=norm(w(:,j+1));
					w(:,j+1)=w(:,j+1)/h(j+1,j);
					rmat(1:j+1,j)=h(1:j+1,j);
					[rmat(:,j),rs,c,s]=zgivens(rmat(:,j),j,rs,c,s);
					Residual=abs(rs(j+1))*dnb/Residual0;
					resrt(iter)=Residual;
					if print==1
						fprintf(' iteration and residual %2.0f%8.5f\n',j,Residual/Residual0);
						[j Residual ]
					end
					if Residual < tol
						break
					end
				end
				yy=dnb*ztrian(rmat,rs,j);
				solX=solX+v(:,1:j)*yy;
				Rhs=Rhs-w*(h*yy);
				if Residual < tol
					return
				end
			end
			% fprintf('outer number of GMRES iterations %6.0f\n',iouter)
			% fprintf('total number of GMRES iterations %6.0f\n',iter);
			% resrt  = log10(resrt/res0);
			% xax   = 1:iter;
			% xao   = 1:iouter;
			%
			% figure(10);
			% plot(xax,resrt,'k-'); hold on;
			% title(' The iterations using flexible GMRES + Right Pre')
			% xlabel(' Iterations')
			% ylabel(' ^{10}log(res_k)')
			%
			% figure(11);
			% plot(xao,iinner,'r'); hold on;
			% title(' The number of inner iterations using flexible GMRES')
			% xlabel(' Outer Iterations')
			% ylabel(' Inner Iterations)')
		end
	end
end

