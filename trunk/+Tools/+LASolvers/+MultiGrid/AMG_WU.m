classdef AMG_WU < handle
	%UNTITLED Summary of this class goes here
	%   Detailed explanation goes here
	
	properties		
		A; % matrix
				
		level;
		max_level;
		relax_it;				% number of smoothing steps
		relaxation_parameter;	% over-relaxation parameter for SOR, 1: usual gauss-seidel and jacobi
		post_smoothing;			% 1: enable, 0: disable postsmoothing step
		max_iter;				% maximal AMG iteration steps
		tollerance;				% stopping tolerance on resdiual
		pc_type;				% 1: jacobi, 2: gauss-seidel, 3: ilu, 4: line jacobi, 5: line gauss-seidel
		connection_threshold;	% threshold value for building strong connection. Set to 0.25 typically
		
		A_stack;
		P_stack;
		MPC;		
		m_cindx;
		l_cindx;
		lpindx;
		
		DefaultGuess;
						
	end
	
	methods
		function obj = AMG_WU(A,Params)
			obj.A = A;
			
			obj.level					= Params.level;
			obj.max_level				= Params.level-1;
			obj.relax_it				= Params.relax_it;
			obj.relaxation_parameter	= Params.relaxation_parameter;
			obj.post_smoothing			= Params.post_smoothing;
			obj.max_iter				= Params.max_iter;
			obj.tollerance				= Params.tollerance;
			obj.pc_type					= Params.pc_type;
			obj.connection_threshold	= Params.connection_threshold;
			
			obj.DefaultGuess=zeros(size(A,1),1);
			
			obj.amg_corsening();
			obj.multi_pcond();
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%    output parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% w : solution
		% error : vector of the ratio of residual and the initial righthand side along iterations           %
		% iter: number of iteration steps for AMG to converge to a given tol
		% flag: 1: converge, 0: not converge
		%%%%%%%%%%%%%%%%%%%%%%%    output parameters     %%%%%%%% %%%%%%%%%%%%%%%%%%%%%%		
		function [Result,resd, nIters] = Solve(obj,Rhs, InitialGuess) 
			
			if nargin>2
				[Result, resd, nIters, flag,clevel]=obj.MG_solver(InitialGuess, Rhs, 1);
			else
				[Result, resd, nIters, flag,clevel]=obj.MG_solver(obj.DefaultGuess, Rhs, 1);
			end
		end
	end
	
	methods(Access=protected)
		
		function amg_corsening(obj)
			
			%%% lpindx(i): number of points in level(i)=lpindx(i+1)-lpindx(i)
			%%% m_cindx(i,:) : global indedx of points in corser grids
			%%% l_cindx(i,:) : pre-level fine grid index of points in current corser grid
			
			n=size(obj.A,2);
			
			obj.P_stack=sparse(2*n,n);
			obj.A_stack=sparse(n,2*n);
			
			obj.m_cindx=sparse(obj.level,n);
			obj.l_cindx=sparse(obj.level,n);
						
			obj.A_stack(1:n,1:n)=obj.A;
			obj.P_stack(1:n,1:n)=speye(n);
			
			obj.m_cindx(1,:)=1:n;
			obj.l_cindx(1,:)=1:n;
			
			lstart=1;
			obj.lpindx(1)=1;
			obj.lpindx(2)=obj.lpindx(1)+n;
			
			fn=n;
			cn=fn;
			
			for i=2:obj.level
				%strong_connect_time=cputime;
				[ni,s,st,lamda]=obj.full_strong_connection(obj.A_stack(1:cn,lstart:lstart+cn-1));

				%[c,f,u] = new_pre_cpoint(A_stack(1:cn,lstart:lstart+cn-1),s,st,lamda);
				[c,f] = obj.pre_cpoint(obj.A_stack(1:cn,lstart:lstart+cn-1),s,st,lamda);
				
				oldc=c;
				oldf=f;
				[c,f,c_to_f]=obj.post_cpoint(obj.A_stack(1:cn,lstart:lstart+cn-1),s,oldc,oldf,ni);

				%cpindx=find(c(:)>0);
				cpindx=(c(:)>0);
				fpindx= f(:)>0;
				lstart=lstart+cn;
				
				%[cn,cm] = size(cpindx);
				cn = nnz(cpindx);
				%[fn,fm]=size(fpindx);
				prolong=sparse(n,cn);
				
				%for j= f(fpindx)%fpindx %1:fn,
					%prolong(f(1,fpindx(j)),:)=c_to_f(f(fpindx(j)),c(cpindx));
					%prolong(f(j),:)=c_to_f(f(j),c(cpindx));
				%	prolong(j,:)=c_to_f(j,c(cpindx));
				%end

				
				%for j=1:cn,
				%	prolong(c(1,cpindx(j)),j)=1.0;
				%	mmm(c(1,cpindx(j)),j)=1.0;
				%end
				
				%[m1,m2] = meshgrid(f(fpindx),1:cn);			
				%prolong = sparse(m1,m2,c_to_f(f(fpindx),c(cpindx)).',n,cn) + sparse(c(cpindx),1:cn,1,n,cn);
				
				prolong(f(fpindx),:)=c_to_f(f(fpindx),c(cpindx));
				
				prolong =prolong + sparse(c(cpindx),1:cn,1,n,cn);
								
				
				%[t1,t2]=size(prolong);
				obj.P_stack(1:n,lstart:lstart+cn-1)=prolong;
				%corser_A=sparse(cn,cn);
				corser_A=prolong'*obj.A_stack(1:n,lstart-n:lstart-1)*prolong;
				
				obj.m_cindx(i,1:cn)=obj.m_cindx(i-1,c(cpindx));
				obj.l_cindx(i,1:cn)=c(cpindx);
				obj.lpindx(i+1)=obj.lpindx(i)+cn;
				
				%obj.A_stack(1:cn,lstart:lstart+cn-1)=sparse(cn,cn);
				obj.A_stack(1:cn,lstart:lstart+cn-1)=corser_A;
				n=cn;	
			end
			
		end
		
		function multi_pcond(obj)
			
			%%% pc_type=1 => Jacobi, pc_type=2 => Gauss-Seidel, pc_type=3 => LU %%%
			%%% w: over-relaxation parameter for SOR
			
			prnum= 1;
			%MPC=sparse(obj.lpindx(obj.level+1)-1,prnum*(obj.lpindx(2)-1));
			
			if(obj.pc_type==3)
				prnum=1;
			end
			iind=[];jind=[];kvalue=[];
			for i=1:obj.level,
				%actual_block_num=0;
				knum=obj.lpindx(i+1)-obj.lpindx(i);
				%PK=sparse(knum,knum);
				PK=obj.A_stack(1:knum,obj.lpindx(i):obj.lpindx(i+1)-1);
				%pindx=obj.m_cindx(i,1:knum);
				%pindx1=1:knum;
				
% 				if(i<=obj.level-1)
% 					%cknum=obj.lpindx(i+2)-obj.lpindx(i+1);
% 					%cpindx=obj.l_cindx(i+1,1:cknum);
% 					%fpindx=full(setdiff(pindx1,cpindx));
% 					%fcount=0;
% 					% fpindx=[];
% 					% for j=1:knum,
% 					% 	%mask=[];
% 					% 	%mask=find(cpindx==pindx1(j));
% 					% 	%if(isempty(mask))
% 					% 	assert(pindx1(j)==j);
% 					% 	if(~any(cpindx==pindx1(j)))
% 					% 		fcount=fcount+1;
% 					% 		fpindx(fcount)=pindx1(j);
% 					% 	end
% 					% end
% 				else
% 					%cknum=knum;
% 					%cpindx=[1:knum];
% 				end
				
				%clear aa bb cc
				for k=1:prnum,
					
					if(obj.pc_type<=3)
						if(obj.pc_type==0)
							%MPC(obj.lpindx(i):obj.lpindx(i+1)-1,(k-1)*knum+1:k*knum)=speye(knum);
							II=speye(knum);
							[aa,bb,cc]=find(II);
						end
						if(obj.pc_type==1)
							JPC=1/obj.relaxation_parameter.*spdiags(diag(PK,0),0,knum,knum);
							[aa,bb,cc]=find(JPC);
							%MPC(obj.lpindx(i):obj.lpindx(i+1)-1,(k-1)*knum+1:k*knum)=JPC;%(iper(k,:),iper(k,:));
							clear JPC;
						end
						if(obj.pc_type==2)
							GSPC=sparse(tril(PK,-1))+1/obj.relaxation_parameter.*spdiags(spdiags(PK,0),0,knum,knum);
							[aa,bb,cc]=find(GSPC);
							%iind=[iind;aa];jind=[jind;bb];kvalue=[kvalue;cc];
							%MPC(obj.lpindx(i):obj.lpindx(i+1)-1,(k-1)*knum+1:k*knum)=GSPC;%(iper(k,:),iper(k,:));
							%clear GSPC;
						end
						if(obj.pc_type==3)
							%[L,U,P]=luinc(PK,1e-04);
							[L,U]=ilu(PK,1e-04);
							LUPC=L*U;
							[aa,bb,cc]=find(LUPC);
							%MPC(obj.lpindx(i):obj.lpindx(i+1)-1,(k-1)*knum+1:k*knum)=LUPC;%(iper(k,:),iper(k,:));
						end
						
					else
						if (obj.pc_type>3 && obj.pc_type<=5)
							
							if(obj.pc_type==4)
								JPC=1/obj.relaxation_parameter.*(diag(diag(PK,0),0)+diag(diag(PK,1),1)+diag(diag(PK,-1),-1));
								[aa,bb,cc]=find(JPC);
								%MPC(obj.lpindx(i):obj.lpindx(i+1)-1,(k-1)*knum+1:k*knum)=JPC;%(iper(k,:),iper(k,:));
							end
							if(obj.pc_type==5)
								GSPC=sparse(tril(PK,-2))+1/obj.relaxation_parameter.*(spdiags(spdiags(PK,-1),-1,knum,knum)+spdiags(spdiags(PK,0),0,knum,knum)+...
									spdiags(spdiags(PK,1),1,knum,knum));
								[aa,bb,cc]=find(GSPC);
								%MPC(obj.lpindx(i):obj.lpindx(i+1)-1,(k-1)*knum+1:k*knum)=GSPC;%(iper(k,:),iper(k,:));
								clear GSPC;
							end
						else
							disp('error: undefined smoother')
							return
						end
					end
					iind=[iind;obj.lpindx(i)+aa-1];
					jind=[jind;(k-1)*knum+bb];
					kvalue=[kvalue;cc];
					%clear per iper secs this_K this_point PK
				end
			end
			obj.MPC=sparse(iind,jind,kvalue);
		end

		function [n,s,st,lamda]=full_strong_connection(obj,Matrix)
			
			[row,col]=size(Matrix);
			%s=sparse(row,col);
			%st=sparse(col,row);
			%n=sparse(row,col);
			s	= cell(row,1);
			n	= cell(row,1);
			st	= cell(col,1); 
			
			%lamda=sparse(1,col);
			lamda=zeros(1,col);
			
			% A_plus=abs(Matrix);
			% for i=1:row,
			% 	A_plus(i,i)=0;
			% 	[max_value(i),max_ind(i)]=max(A_plus(i,:));
			% end
			
			A_plus=abs(Matrix - diag(diag(Matrix)));

			[NzMat_i,NzMat_j] = find(Matrix~=0);
			
			
				max_value=max(A_plus,[],2);
			trhsld = obj.connection_threshold*max_value;		
			
			for indx=1:numel(NzMat_i)
				if A_plus(NzMat_i(indx),NzMat_j(indx))~=0
					n{NzMat_i(indx)}= [n{NzMat_i(indx)} , NzMat_j(indx)];
				end
				if (NzMat_i(indx)~=NzMat_j(indx)) && (A_plus(NzMat_i(indx),NzMat_j(indx)) > trhsld(NzMat_i(indx)))
					s{NzMat_i(indx)}= [s{NzMat_i(indx)} , NzMat_j(indx)];
				end
			end
			

% 			for i=1:row,
% 				col=find(Matrix(i,:));
% 				%col = find(MatrixT(:,i)).';
% 				%col=NzMat_j(NzMat_i==i).';
% 				
%  				indx = A_plus(i,col)~=0;
% 				n{i} = col(indx);
% 				
% 				indx = find(A_plus(i,col) > trhsld(i));
% 				isnoteq_i = col(indx)~=i;
% 				indx = indx(isnoteq_i);
% 				s{i} = col(indx);
% 
% 			end

		

			
			col=row;
			for i=1:col,
				row=find(Matrix(:,i));
				%rsize=size(row);
				%rnum=rsize(1);
				
				%row=NzMat_i(NzMat_j==i);
				
				%assert(all(row==find(Matrix(:,i))));
				
				indx = find(A_plus(row,i) >= trhsld(row));
				isnoteq_i = row(indx)~=i;
				indx = indx(isnoteq_i);
				%st1(i,1:numel(indx)) = row(indx);
				st{i} = row(indx);
				%lamda(1,i)=lamda(1,i)+numel(indx);
				assert(lamda(i)==0);
				lamda(i)=lamda(i)+numel(indx);
				
% 				k=1;l=1;
% 				for r=1:rnum,
% 					j=row(r);
% 					if(j~=i && A_plus(j,i) >= trhsld(j))
% 						st(i,l)=j;
% 						l=l+1;
% 						lamda(1,i)=lamda(1,i)+1;
% 					end
% 				end
				%row=[];
			end
		end
		
		function [InitialGuess, resd, iter, flag,clevel]=MG_solver( obj,InitialGuess, Rhs, cycle)
			
			%%%
			% tol_type: 1: discrete l2-norm, 2:contineous l2-norm
			
			%%% mw = initial guess %%%%
			iter=0;
			pnum=numel(InitialGuess);
			%itsize=size(obj.tollerance,2);
			err(1)=0;
			if(obj.max_iter>1)
				resd=Rhs-obj.A_stack(1:pnum,1:pnum)*InitialGuess;
				err(1)=norm(resd)/norm(Rhs);
			else
				err(1)=1;   %%% for 1 step v-cycle multigrid preconditioner in GMRES
			end
			
			ini_iteration=1;
			clevel=obj.max_level;
			%pre_err=0;
			flag=0;
			%true_x=sparse(pnum,obj.max_level+1);
			%while (ini_iteration==1 || (err(iter+1) > obj.tollerance && iter<obj.max_iter)),
			while (ini_iteration==1 || (err > obj.tollerance && iter<obj.max_iter)),
				ini_iteration=0;
				%pre_err=err; %pre_err=err(iter+1);
				
				[InitialGuess, resd, iter, flag, clevel] = obj.MG(InitialGuess, Rhs, iter,  clevel, cycle);
				
				err=norm(resd)/norm(Rhs);%err(iter+1)=norm(resd)/norm(Rhs);
				%if(err(iter+1)<obj.tollerance)
				%	flag=1;
				%end
				flag = err<obj.tollerance;
				%conv_rate=err(iter+1)/pre_err;
			end
		end
		
		function [InitialGuess, resd, iter, flag, clevel] = MG(obj,InitialGuess, Rhs, iter, clevel,cycle)
			
			% obj.m_cindx point to global index for nodes in each level %
			% obj.lpindx record the start indx and end indx of each level for extracting
			% matrix A, Prolongation operator and Restriction operator from obj.A_stack,P_stack,R_stack
			% obj.m_cindx and Mrs are no use at this stage %%%
			% obj.relax_it: smoothing steps
			% cycle: multigrid cycle types
			% iter: multigrid iteration counts
			% clevel: current grid level on calling MG
			
			pernum= 1;
			head=obj.lpindx(obj.max_level-clevel+1);
			tail=obj.lpindx(obj.max_level-clevel+2)-1;
			delta_head_tail=tail-head+1;
			mA=obj.A_stack(1:delta_head_tail,head:tail);
			
			%[L,U,perm] = lu(mA,'vector');
			%[L,U] = lu(mA);
			n = size(mA,2);
			%delta_x=zeros(n,1);
			
			if(obj.max_level==0)
				
				Ml=obj.MPC(head:tail,1:pernum*delta_head_tail);
				%prolong=obj.P_stack(head:tail,head:tail);
				%cx=zeros(n,1);
				prelax=0;
				if(size(obj.relax_it,2)==1)
					[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it*(2^prelax),n);
				else
					[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it(1),n);
				end
				iter=iter+1;
				clevel=obj.max_level;
				return
			end
			
			if(clevel>0)
				
				Ml=obj.MPC(head:tail,1:pernum*delta_head_tail);
				prelax=0;
				clevel=clevel-1;
				chead=obj.lpindx(obj.max_level-clevel+1);
				ctail=obj.lpindx(obj.max_level-clevel+2)-1;
				cdelta_head_tail=obj.lpindx(obj.max_level-clevel+2)-obj.lpindx(obj.max_level-clevel+1);
				prolong=obj.P_stack(1:delta_head_tail,chead:ctail);
				restrict=prolong';
				%resd=Rhs-mA*InitialGuess;
				
				if(size(obj.relax_it,2)==1)
					[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it-(obj.max_level-clevel)+1,n);
				else
					[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it(obj.max_level-clevel),n);
				end
				
				corser_b=restrict*resd;
				cini=zeros(cdelta_head_tail,1);
				[cx, resd, iter, flag, clevel]=obj.MG( cini, corser_b, iter,clevel,cycle);
				
			end
			
			if(clevel==0&&obj.max_level>0)
				InitialGuess=mA\Rhs;
				
				%InitialGuess = U\(L\(Rhs(perm,:)));
				%InitialGuess = U\(L\(Rhs));
				
				clevel=clevel-1;
				resd=zeros(n,1);
				%iter=iter;
				flag=0;
				return
			end
			
			if(clevel<0 && clevel + obj.max_level>0)
				InitialGuess=InitialGuess+prolong*cx;
				if(obj.post_smoothing==1)
					if(size(obj.relax_it,2)==1)
						[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it*(2^prelax),n);
					else
						[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it(obj.max_level+clevel+1),n);
					end
				else
					resd=Rhs-mA*InitialGuess;
				end
				%iter=iter;
				flag=0;
				clevel=clevel-1;
				return
			end
			
			if((clevel<0 && clevel+obj.max_level==0) || obj.max_level==0)
				InitialGuess=InitialGuess+prolong*cx;
				%resd=Rhs-mA*InitialGuess;
				if(obj.post_smoothing==1)
					if(size(obj.relax_it,2)==1)
						[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it*(2^prelax),n);
					else
						[InitialGuess,resd]=precond_inv(InitialGuess,pernum,obj.relax_it(obj.max_level+clevel+1),n);
					end
				else
					resd=Rhs-mA*InitialGuess;
				end
				iter=iter+1;
				clevel=obj.max_level;
			end
			
			function [InitialGuess,r]=precond_inv(InitialGuess,pernum,relax_iter,n)
				
				for i=1:relax_iter,
					for piter=1:pernum,
						M=Ml(1:n,(piter-1)*n+1:piter*n);
						c = M\(Rhs-(mA*InitialGuess));
						InitialGuess=InitialGuess+c;
					end
				end
				r=Rhs-mA*InitialGuess;
			end

			
		end
		
		function [c,f,u]=pre_cpoint(obj,A,s,st,lamda)
			
			[row,col]=size(A);
			c=sparse(1,col);
			f=sparse(1,col);
			u=[1:row];
			cnum=1;
			l=0;
			%m=0;
			%start=1;
			while (any(u)),
				[max_value,max_lamda]=max(lamda);
				c(cnum)=max_lamda; %%% C = C union {vi}
				cnum=cnum+1;
				u(max_lamda)=0;    %%% U= U - {vi}
				lamda(max_lamda)=-9999999;
				%temp=find(st(max_lamda,:)~=0);
				%sit_col=size(temp,2);
				%clear temp
				
				%sit_col=nnz(st1(max_lamda,:));
				sit_col = numel(st{max_lamda});
				%temp=find(s(max_lamda,:)~=0);
				%si_col=size(temp,2);
				%clear temp
				%sit_col=count_nonzero(st(max_lamda,:));
				%si_col=count_nonzero(s(max_lamda,:));
				
				%si_col=nnz(s1(max_lamda,:));
				si_col=numel(s{max_lamda});
				%add_to_f=0;
				for j=1:sit_col,
					%mask=[];
					%mask=find(u==st(max_lamda,j));
					%if(~isempty(mask)) %%% (S_i)T and U ~= empty
					%if any(u==st1(max_lamda,j))
					if any(u==st{max_lamda}(j))
						%%%step 4
						l=l+1;
						%f(l)=st1(max_lamda,j);  %%% F=F union {vj=st(max_lamda,j)}
						f(l)=st{max_lamda}(j);  %%% F=F union {vj=st(max_lamda,j)}
						u(f(l))=0;             %%% U=U-{vj}
						lamda(f(l))=-9999999;
						%%% step 5
						%sj_col=nnz(s1(f(l),:)); %count_nonzero(s(f(l),:));
						sj_col=numel(s{f(1)});
						for k=1:sj_col,
							%mask=[];
							%mask=find(u==s(f(l),k));
							%if(~isempty(mask))  %%%(S_j) and U ~= empty
							%if any(u==s1(f(l),k))
							if any(u==s{f(l)}(k))
								%lamda(s1(f(l),k))=lamda(s1(f(l),k))+1; %%% zk=zk+1
								lamda(s{f(l)}(k))=lamda(s{f(l)}(k))+1; %%% zk=zk+1
							end
						end
					end
				end
				%%% step 6
				for k=1:si_col,
					%mask=[];
					%mask=find(u==s(max_lamda,k));
					%if(~isempty(mask))
					%if any(u==s1(max_lamda,k))
					if any(u==s{max_lamda}(k))
						%lamda(s1(max_lamda,k))=lamda(s1(max_lamda,k))-1; %%%zj=zj-1
						lamda(s{max_lamda}(k))=lamda(s{max_lamda}(k))-1; %%%zj=zj-1
					end
				end
			end %%% end while
		end
		
		function [C1,F1,W1]=post_cpoint(obj,A,S,C1,F1,N)
			
			[row,col]=size(A);
			old_F=F1;
			W1=sparse(row,col);
			T1=zeros(1,col);
			%temp=find(F1~=0);
			%fnum=size(temp,2);
			%clear temp
			%temp=find(C~=0);
			%cnum=size(temp,2);
			cnum=nnz(C1);
			C1=sort(C1);
			l=0;
			i=0;
			while (any(F1))
				CT=sparse(1,col);
				%CI=sparse(1,col);
				d=sparse(1,col);
				%mask=[];
				in_T=1;
				i=i+1;
				while(in_T==1 && i<=col)
					%mask=find(T==F(i));
					%if(~isempty(mask)|F(i)==0)
					if(F1(i)==0 || any(T1==F1(i)) )
					%if( F(i)==0 || ismember(F(i),T))
						i=i+1;
					else
						in_T=0;
						T1(i)=F1(i);
						break
					end
				end
				if(i>=col)
					break
				end
				%sindx=find(S(F(i),:));
				%nindx=find(N(F(i),:));
				sindx=find(S{F1(i)});
				nindx=find(N{F1(i)});
				snum=numel(sindx);
				nnum=numel(nindx);
				if(snum==0)
					disp('There is no strong connection for this fine point')
					cnum=cnum+1;
					C1(cnum)=F1(i);
					F1(i)=0;
				end
				
				%m=0;
				dwm=0;
				%dsm=0;
				%%% setup Ci,Di_strong==DS %%%
				%DS=sparse(1,snum);
				%DW=sparse(1,nnum);
% 				for j=1:snum,
% 					if(any(C==S{F1(i)}(sindx(j))))
% 						m=m+1;   %%% m = number of elements in CI %%%%
% 						CI(m)=S{F1(i)}(sindx(j));
% 						d(CI(m))=A(F1(i),CI(m));
% 					end
% 				end
 				
				CI = [];
				DS = [];
				for j=S{F1(i)}
					if ismember(j,C1)
						CI = [CI j];
					else
						DS = [DS j];
					end
				end

				%CI=intersect(C(C~=0),S{F1(i)});
				%CI=intersect(S{F1(i)},C);%???
				d(CI)=A(F1(i),CI);
				m=numel(CI);
				
% 				for j=1:snum,
% 					if(~any(CI==S{F1(i)}(sindx(j))))
% 						dsm=dsm+1;
% 						DS(dsm)=S{F1(i)}(sindx(j));
% 					end
% 				end
				
				%DS = setdiff(CI,S{F1(i)});
				dsm=numel(DS);
				
				%%% setup Di_weak==DW
				if (F1(i)~=0)
					d(F1(i))=A(F1(i),F1(i));
% 					for j=1:nnum,
% 						%mask=[];
% 						%mask=find(S1(F(i),:)==N(F(i),nindx(j)),1);
% 						mask=find(S{F1(i)}==N{F1(i)}(nindx(j)),1);
% 						if(isempty(mask))
% 							dwm=dwm+1; %%% dwm = number of elements in DW %%%%
% 							%DW(dwm)=N1(F(i),nindx(j));
% 							DW(dwm)=N{F1(i)}(nindx(j));
% 							d(F1(i))=d(F1(i))+A(F1(i),DW(dwm));
% 						end
% 					end
					
					DW = setdiff(S{F1(i)}, N{F1(i)} );
					for j=1:numel(DW)
						d(F1(i))=d(F1(i))+A(F1(i),DW(dwm));
					end
				end
				%%%% start step 5 %%%%
				cm=0;
				done=0;
				while(~done&snum~=0)
					nextf=0;
					if dsm==0, done=1; end %temporary
					for j=1:dsm %snum,
						if(DS(j)~=0)
							for k=1:m,
								if(any(S{DS(j)}==CI(k)))            %%% if S_j ^ C_i not empty %%%
									nextf=nextf+1;
									sum_ajl=0;
									for k=1:m,
										sum_ajl=sum_ajl+A(DS(j),CI(k));
									end
									for k=1:m,
										if(sum_ajl==0)
											disp('error in post_cpoint: divided by 0 when computing interpolation weights')
										end
										d(CI(k))=d(CI(k))+A(F1(i),DS(j))*A(DS(j),CI(k))/sum_ajl;
									end
									break;
								end
							end
							%if(isempty(mask))               %%% if S_j ^ C_i is empty %%%
							%if(any(S1(DS(j),:)==CI(k)))               %%% if S_j ^ C_i is empty %%%
							if(any(S{DS(j)}==CI(k)))               %%% if S_j ^ C_i is empty %%%
								if(any(CT)==0)
									cm=cm+1;
									m=m+1;
									dsm=dsm-1;
									CT(cm)=DS(j);
									CI(m)=DS(j);
									DS(j)=0;
									d(CI(m))=A(F1(i),CI(m));
									done=0;
									break;
								else
									cnum=cnum+1;
									C1(cnum)=F1(i);
									F1(i)=0;
									done=1;
									break;
								end
							end %if(isempty(mask))
						elseif(j==snum)
							done=1;
							break;
						else
							done=0;
						end  %if(DS(j)~=0)
						
					end %for j=1:snum
				end %while(~done)
				
				if(F1(i)~=0)
					if(any(CT))
						cnum=cnum+1;
						C1(cnum)=CT(cm);
						remove_f= F1==CT(cm);
						F1(remove_f)=0;
					end
					
					if(d(F1(i))==0)
						disp('error in final c-point selection: interpolation weight are not well defined')
					end
					
					W1(F1(i),CI(1:m))=-d(CI(1:m))/d(F1(i));
					%for k=1:m,						
					%	W1(F1(i),CI(k))=-d(CI(k))/d(F1(i));
					%end
				end
				if(T1==old_F)
					break
				end
			end %%% end while loop %%%
		end

		
	end
	
end

